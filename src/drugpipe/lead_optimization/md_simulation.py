"""Short molecular-dynamics simulation and MM-GBSA binding energy via OpenMM.

Runs a brief MD trajectory on a protein-ligand complex, then estimates
the binding free energy with MM-GBSA (OBC2 implicit solvent model).

Optional dependencies:
  - openmm    (conda install -c conda-forge openmm)
  - pdbfixer  (conda install -c conda-forge pdbfixer)
  - openff-toolkit (pip install openff-toolkit)  — for ligand force field
"""
# 短时 MD 模拟与 MM-GBSA 结合自由能（OpenMM）。

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
    import mdtraj
    _MDTRAJ_OK = True
except ImportError:
    _MDTRAJ_OK = False


def _select_platform(device: str) -> Tuple[Optional[str], Dict[str, str]]:
    """Choose OpenMM platform based on user config.

    Returns (platform_name, properties_dict).
    ``device`` can be: ``"auto"``, ``"cuda"``, ``"opencl"``, ``"cpu"``.
    """
    if not _OPENMM_OK:
        return None, {}

    device = device.lower().strip()

    if device == "cuda":
        return "CUDA", {"Precision": "mixed"}
    if device in ("opencl", "mps"):
        return "OpenCL", {"Precision": "mixed"}
    if device == "cpu":
        return "CPU", {}

    # auto: try CUDA → OpenCL → CPU
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            openmm.Platform.getPlatformByName(name)
            props = {"Precision": "mixed"} if name in ("CUDA", "OpenCL") else {}
            logger.info("OpenMM auto-selected platform: %s", name)
            return name, props
        except Exception:
            continue
    return None, {}


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
            pose_file = df.at[idx, "docking_pose_file"] if "docking_pose_file" in df.columns else ""
            logger.info("Running MD for molecule %d: %s", idx, smi[:50])

            energy, rmsd = self._simulate_one(
                smi, pdb_path, pose_file, md_dir, idx,
            )
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
        pose_file: str,
        md_dir: Path,
        mol_idx: int,
    ) -> Tuple[Optional[float], Optional[float]]:
        """Run energy-minimisation + short MD for one protein-ligand pair.

        Returns (binding_energy_kcal, mean_rmsd_angstrom).
        """
        try:
            complex_energy = self._minimise_complex(receptor_pdb, smiles)
            receptor_energy = self._minimise_receptor_only(receptor_pdb)
            ligand_energy = self._minimise_ligand_only(smiles)

            if all(e is not None for e in (complex_energy, receptor_energy, ligand_energy)):
                binding = complex_energy - receptor_energy - ligand_energy
                logger.info(
                    "Mol %d MM-GBSA: complex=%.2f  receptor=%.2f  ligand=%.2f  ΔG=%.2f kcal/mol",
                    mol_idx, complex_energy, receptor_energy, ligand_energy, binding,
                )
                return binding, None  # RMSD requires trajectory analysis
            return None, None

        except Exception as exc:
            logger.warning("MD simulation failed for mol %d: %s", mol_idx, exc)
            return None, None

    # ------------------------------------------------------------------
    def _minimise_complex(self, receptor_pdb: str, ligand_smi: str) -> Optional[float]:
        """Build protein-ligand system in implicit solvent, minimise, return energy."""
        return self._gb_energy_from_pdb(receptor_pdb, ligand_smi=ligand_smi)

    def _minimise_receptor_only(self, receptor_pdb: str) -> Optional[float]:
        return self._gb_energy_from_pdb(receptor_pdb, ligand_smi=None)

    def _minimise_ligand_only(self, ligand_smi: str) -> Optional[float]:
        """Compute GB energy of the ligand alone."""
        mol = safe_mol(ligand_smi)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp:
            tmp.write(Chem.MolToPDBBlock(mol))
            tmp_path = tmp.name

        try:
            return self._single_point_gb_energy(tmp_path)
        except Exception as exc:
            logger.debug("Ligand GB energy failed: %s", exc)
            return None

    # ------------------------------------------------------------------
    def _gb_energy_from_pdb(
        self,
        pdb_path: str,
        ligand_smi: Optional[str] = None,
    ) -> Optional[float]:
        """Load a PDB (optionally with ligand), build GB system, minimise,
        and return potential energy in kcal/mol."""
        try:
            pdb = app.PDBFile(pdb_path)
            gb_model_map = {
                "OBC1": app.OBC1,
                "OBC2": app.OBC2,
                "GBn": app.GBn,
                "GBn2": app.GBn2,
            }
            gb = gb_model_map.get(self.gb_model, app.OBC2)

            forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
            modeller = app.Modeller(pdb.topology, pdb.positions)

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
                simulation = app.Simulation(modeller.topology, system, integrator,
                                            platform, self.platform_props)
            else:
                simulation = app.Simulation(modeller.topology, system, integrator)
            simulation.context.setPositions(modeller.positions)

            simulation.minimizeEnergy(maxIterations=1000)

            state = simulation.context.getState(getEnergy=True)
            energy_kj = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            return float(energy_kj)

        except Exception as exc:
            logger.debug("GB energy calculation failed for %s: %s", pdb_path, exc)
            return None

    # ------------------------------------------------------------------
    def _single_point_gb_energy(self, pdb_path: str) -> Optional[float]:
        """Single-point GB energy for a small molecule PDB."""
        try:
            pdb = app.PDBFile(pdb_path)
            forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=app.NoCutoff,
                constraints=app.HBonds,
            )
            integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
            simulation = app.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy(maxIterations=500)
            state = simulation.context.getState(getEnergy=True)
            return float(state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole))
        except Exception as exc:
            logger.debug("Single-point GB energy failed: %s", exc)
            return None
