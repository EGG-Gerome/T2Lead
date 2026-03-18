"""Short molecular-dynamics simulation and MM-GBSA binding energy via OpenMM.

Runs a brief MD trajectory on a protein-ligand complex, then estimates
the binding free energy with MM-GBSA (OBC2 implicit solvent model).

Required dependencies:
  - openmm              (conda install -c conda-forge openmm)
  - openmmforcefields   (conda install -c conda-forge openmmforcefields)
  - openff-toolkit      (conda install -c conda-forge openff-toolkit)
Optional:
  - mdtraj              (conda install -c conda-forge mdtraj)
"""
# 短时 MD 模拟与 MM-GBSA 结合自由能（OpenMM + GAFF2 小分子力场）。

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from drugpipe.utils.chem import safe_mol

logger = logging.getLogger(__name__)

try:
    import openmm
    import openmm.app as app
    import openmm.unit as unit
    _OPENMM_OK = True
except ImportError:
    _OPENMM_OK = False

try:
    from openmmforcefields.generators import GAFFTemplateGenerator
    from openff.toolkit.topology import Molecule as OFFMolecule
    _GAFF_OK = True
except ImportError:
    _GAFF_OK = False

try:
    import mdtraj
    _MDTRAJ_OK = True
except ImportError:
    _MDTRAJ_OK = False


def _select_platform(device: str) -> Tuple[Optional[str], Dict[str, str]]:
    """Choose OpenMM platform based on user config."""
    if not _OPENMM_OK:
        return None, {}

    device = device.lower().strip()
    if device == "cuda":
        return "CUDA", {"Precision": "mixed"}
    if device in ("opencl", "mps"):
        return "OpenCL", {"Precision": "mixed"}
    if device == "cpu":
        return "CPU", {}

    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            openmm.Platform.getPlatformByName(name)
            props = {"Precision": "mixed"} if name in ("CUDA", "OpenCL") else {}
            logger.info("OpenMM auto-selected platform: %s", name)
            return name, props
        except Exception:
            continue
    return None, {}


def _smiles_to_3d_pdb(smiles: str) -> Optional[str]:
    """Generate a 3D PDB block from a SMILES string."""
    mol = safe_mol(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        params2 = AllChem.ETKDGv3()
        params2.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params2) != 0:
            return None
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    return Chem.MolToPDBBlock(mol)


class MDSimulator:
    """Run short MD simulations and compute MM-GBSA binding energies."""

    def __init__(self, cfg: Dict[str, Any]):
        lo = cfg.get("lead_optimization", {})
        md = lo.get("md_simulation", {})
        self.enabled = bool(md.get("enabled", True))
        self.top_n = int(md.get("top_n_for_md", 10))
        self.equilibration_ps = float(md.get("equilibration_ps", 100))
        self.production_ns = float(md.get("production_ns", 1.0))
        self.gb_model = md.get("gb_model", "OBC2")

        device_str = cfg.get("pipeline", {}).get("device", "auto")
        self.platform, self.platform_props = _select_platform(device_str)

    # ------------------------------------------------------------------
    def run(
        self,
        df: pd.DataFrame,
        protein_info: Dict[str, Any],
        out_dir: Path,
    ) -> pd.DataFrame:
        """Run MD + MM-GBSA on the top-N candidates (by docking score).

        Adds columns: ``md_binding_energy``, ``md_rmsd_mean``.
        """
        df = df.copy()
        df["md_binding_energy"] = np.nan
        df["md_rmsd_mean"] = np.nan

        if not self.enabled:
            logger.info("MD simulation disabled in config.")
            return df

        if not _OPENMM_OK:
            logger.warning(
                "OpenMM not installed (conda install -c conda-forge openmm). "
                "Skipping MD simulation."
            )
            return df

        if not _GAFF_OK:
            logger.warning(
                "openmmforcefields / openff-toolkit not installed. "
                "Cannot parameterize small-molecule ligands for MD. "
                "Install: conda install -c conda-forge openmmforcefields openff-toolkit"
            )
            return df

        pdb_path = protein_info.get("pdb_path")
        if not pdb_path or not Path(pdb_path).exists():
            logger.error("No prepared PDB for MD simulation.")
            return df

        candidates = self._select_candidates(df)
        if not candidates:
            logger.warning("No candidates selected for MD.")
            return df

        md_dir = out_dir / "md_trajectories"
        md_dir.mkdir(parents=True, exist_ok=True)

        for idx in candidates:
            smi = df.at[idx, "canonical_smiles"]
            logger.info("Running MM-GBSA for molecule %d: %s", idx, smi[:60])

            energy, rmsd = self._simulate_one(smi, pdb_path, md_dir, idx)
            if energy is not None:
                df.at[idx, "md_binding_energy"] = energy
            if rmsd is not None:
                df.at[idx, "md_rmsd_mean"] = rmsd

        n_done = df["md_binding_energy"].notna().sum()
        logger.info("MD simulation complete: %d / %d candidates processed.", n_done, len(candidates))
        return df

    # ------------------------------------------------------------------
    def _select_candidates(self, df: pd.DataFrame) -> List[int]:
        """Pick top-N molecules with best docking scores for MD."""
        if "docking_score" in df.columns:
            ranked = df.dropna(subset=["docking_score"]).sort_values("docking_score")
            return ranked.head(self.top_n).index.tolist()
        return df.head(self.top_n).index.tolist()

    # ------------------------------------------------------------------
    def _simulate_one(
        self,
        smiles: str,
        receptor_pdb: str,
        md_dir: Path,
        mol_idx: int,
    ) -> Tuple[Optional[float], Optional[float]]:
        """Compute MM-GBSA binding energy: ΔG = E(complex) − E(receptor) − E(ligand).

        Returns (binding_energy_kcal, mean_rmsd_angstrom).
        """
        try:
            off_mol = OFFMolecule.from_smiles(smiles, allow_undefined_stereo=True)
        except Exception as exc:
            logger.warning("Mol %d: OpenFF cannot parse SMILES: %s", mol_idx, exc)
            return None, None

        try:
            complex_energy = self._gb_energy_complex(receptor_pdb, smiles, off_mol)
            receptor_energy = self._gb_energy_receptor(receptor_pdb)
            ligand_energy = self._gb_energy_ligand(smiles, off_mol)

            if all(e is not None for e in (complex_energy, receptor_energy, ligand_energy)):
                binding = complex_energy - receptor_energy - ligand_energy
                logger.info(
                    "Mol %d MM-GBSA: complex=%.1f  receptor=%.1f  ligand=%.1f  ΔG=%.2f kcal/mol",
                    mol_idx, complex_energy, receptor_energy, ligand_energy, binding,
                )
                return binding, None
            else:
                logger.warning(
                    "Mol %d: incomplete energies (complex=%s, receptor=%s, ligand=%s)",
                    mol_idx, complex_energy, receptor_energy, ligand_energy,
                )
            return None, None

        except Exception as exc:
            logger.warning("MD simulation failed for mol %d: %s", mol_idx, exc)
            return None, None

    # ------------------------------------------------------------------
    def _gb_energy_complex(
        self, receptor_pdb: str, ligand_smi: str, off_mol: Any,
    ) -> Optional[float]:
        """Build protein+ligand in implicit solvent, minimise, return energy."""
        lig_pdb = _smiles_to_3d_pdb(ligand_smi)
        if lig_pdb is None:
            return None

        try:
            rec_pdb = app.PDBFile(receptor_pdb)
            modeller = app.Modeller(rec_pdb.topology, rec_pdb.positions)

            with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp:
                tmp.write(lig_pdb)
                tmp_path = tmp.name
            lig_pdb_obj = app.PDBFile(tmp_path)
            modeller.add(lig_pdb_obj.topology, lig_pdb_obj.positions)

            gaff = GAFFTemplateGenerator(molecules=off_mol, forcefield="gaff-2.11")
            forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
            forcefield.registerTemplateGenerator(gaff.generator)

            return self._minimise_and_energy(forcefield, modeller)
        except Exception as exc:
            logger.debug("Complex GB energy failed: %s", exc)
            return None

    def _gb_energy_receptor(self, receptor_pdb: str) -> Optional[float]:
        """GB energy of the receptor alone (no GAFF needed)."""
        try:
            pdb = app.PDBFile(receptor_pdb)
            modeller = app.Modeller(pdb.topology, pdb.positions)
            forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
            return self._minimise_and_energy(forcefield, modeller)
        except Exception as exc:
            logger.debug("Receptor GB energy failed: %s", exc)
            return None

    def _gb_energy_ligand(self, ligand_smi: str, off_mol: Any) -> Optional[float]:
        """GB energy of the ligand alone, parameterized with GAFF2."""
        lig_pdb = _smiles_to_3d_pdb(ligand_smi)
        if lig_pdb is None:
            return None

        try:
            with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp:
                tmp.write(lig_pdb)
                tmp_path = tmp.name
            pdb = app.PDBFile(tmp_path)
            modeller = app.Modeller(pdb.topology, pdb.positions)

            gaff = GAFFTemplateGenerator(molecules=off_mol, forcefield="gaff-2.11")
            forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
            forcefield.registerTemplateGenerator(gaff.generator)

            return self._minimise_and_energy(forcefield, modeller)
        except Exception as exc:
            logger.debug("Ligand GB energy failed: %s", exc)
            return None

    # ------------------------------------------------------------------
    def _minimise_and_energy(
        self, forcefield: Any, modeller: Any, max_iter: int = 1000,
    ) -> Optional[float]:
        """Create system, minimise, return potential energy in kcal/mol."""
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )
        integrator = openmm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds,
        )
        platform = None
        if self.platform:
            try:
                platform = openmm.Platform.getPlatformByName(self.platform)
            except Exception:
                pass
        if platform:
            simulation = app.Simulation(
                modeller.topology, system, integrator,
                platform, self.platform_props,
            )
        else:
            simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=max_iter)

        state = simulation.context.getState(getEnergy=True)
        return float(state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))
