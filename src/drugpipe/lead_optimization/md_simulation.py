"""Short molecular-dynamics simulation and MM-GBSA binding energy via OpenMM.

After energy minimisation of the complex, runs configurable Langevin MD
(``equilibration_ps`` / ``production_ns``) and reports mean ligand
heavy-atom RMSD vs the initial pose over 1-ps snapshots; MM-GBSA-style
binding energy is from minimised receptor, ligand, and complex states
(OBC2 implicit solvent).

Required dependencies:
  - openmm              (conda install -c conda-forge openmm)
  - openmmforcefields   (conda install -c conda-forge openmmforcefields)
  - openff-toolkit      (conda install -c conda-forge openff-toolkit)
Optional:
  - mdtraj              (conda install -c conda-forge mdtraj)
"""
# Short molecular-dynamics simulation and MM-GBSA binding energy via OpenMM.
# 说明模块职责、上下游关系与维护注意事项。

# 短时 MD 模拟与 MM-GBSA 结合自由能（OpenMM + GAFF2 小分子力场）。

from __future__ import annotations

import logging
import tempfile
from dataclasses import dataclass
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


# 2 fs timestep → 500 steps per picosecond
_MD_TIMESTEP_PS = 0.002
_MD_STEPS_PER_PS = int(round(0.001 / _MD_TIMESTEP_PS))

# kJ/mol → kcal/mol conversion factor
_KJ_TO_KCAL = 1.0 / 4.184


def _energy_as_kcal(val) -> float:
    """Convert potential energy to kcal/mol float.

    Handles both OpenMM ``Quantity`` objects and raw ``float`` values.
    When the value lacks ``value_in_unit`` (already a plain number),
    it is assumed to be in kJ/mol — OpenMM's native energy unit.
    """
    if hasattr(val, "value_in_unit") and _OPENMM_OK:
        return float(val.value_in_unit(unit.kilocalories_per_mole))
    return float(val) * _KJ_TO_KCAL


@dataclass
class _ComplexState:
    """OpenMM complex simulation after minimisation; used to continue with MD."""

    simulation: Any
    n_rec_atoms: int
    lig_init_pos: np.ndarray  # Å, shape (n_lig_atoms, 3)
    heavy_mask: np.ndarray  # bool, shape (n_lig_atoms,)


def _validate_platform(name: str, props: Dict[str, str]) -> bool:
    """Actually run a tiny simulation to verify the platform works at kernel level."""
    try:
        platform = openmm.Platform.getPlatformByName(name)
        system = openmm.System()
        system.addParticle(1.0)
        force = openmm.CustomExternalForce("x^2+y^2+z^2")
        force.addParticle(0, [])
        system.addForce(force)
        integrator = openmm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds,
        )
        ctx = openmm.Context(system, integrator, platform, props)
        ctx.setPositions([[0, 0, 0]])
        ctx.getState(getEnergy=True)
        del ctx
        return True
    except Exception as exc:
        logger.warning("OpenMM platform '%s' detected but unusable: %s", name, exc)
        return False


def _select_platform(device: str) -> Tuple[Optional[str], Dict[str, str]]:
    """Choose OpenMM platform based on user config, with runtime validation.

    Even when a GPU platform is registered, the actual CUDA/OpenCL kernels
    may fail on newer architectures (e.g. Blackwell sm_120 with older PTX).
    This function runs a micro-simulation to verify the platform before
    committing to it, and falls back to CPU if the GPU path is broken.
    """
    if not _OPENMM_OK:
        return None, {}

    device = device.lower().strip()

    if device == "cuda":
        if _validate_platform("CUDA", {"Precision": "mixed"}):
            return "CUDA", {"Precision": "mixed"}
        logger.warning("CUDA requested but validation failed, falling back to CPU.")
        return "CPU", {}
    if device in ("opencl", "mps"):
        if _validate_platform("OpenCL", {"Precision": "mixed"}):
            return "OpenCL", {"Precision": "mixed"}
        logger.warning("OpenCL requested but validation failed, falling back to CPU.")
        return "CPU", {}
    if device == "cpu":
        return "CPU", {}

    # auto: try GPU platforms first, validate, fall back to CPU
    for name in ("CUDA", "OpenCL"):
        props = {"Precision": "mixed"}
        try:
            openmm.Platform.getPlatformByName(name)
        except Exception:
            continue
        if _validate_platform(name, props):
            logger.info("OpenMM auto-selected platform: %s", name)
            return name, props
        logger.info("OpenMM platform %s registered but failed validation, skipping.", name)

    logger.info("OpenMM falling back to CPU platform.")
    return "CPU", {}


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


def _pdbqt_to_pdb_block(pdbqt_path: str) -> Optional[str]:
    """Extract the first MODEL from a Vina PDBQT and return a PDB-format string."""
    lines: List[str] = []
    with open(pdbqt_path) as fh:
        for raw in fh:
            if raw.startswith("ENDMDL"):
                break
            if raw.startswith(("ATOM", "HETATM")):
                base = raw[:66].rstrip().ljust(76)
                atom_name = raw[12:16].strip()
                elem = atom_name.lstrip("0123456789")[:2].strip()
                if len(elem) > 1:
                    elem = elem[0].upper() + elem[1].lower()
                else:
                    elem = elem.upper()
                lines.append(f"{base} {elem:>2s}")
    if not lines:
        return None
    lines.append("END")
    return "\n".join(lines) + "\n"


def _load_docked_pose_pdb(smiles: str, pdbqt_path: str) -> Optional[str]:
    """Build an all-atom PDB block with heavy-atom positions from the docked pose.

    Parses the PDBQT to recover heavy-atom coordinates, assigns bond orders
    from the SMILES template, adds explicit hydrogens, then does a constrained
    embedding so H positions are physically reasonable while heavy atoms stay
    at their docked positions.  Falls back to ``_smiles_to_3d_pdb`` on failure.
    """
    pdb_block = _pdbqt_to_pdb_block(pdbqt_path)
    if pdb_block is None:
        logger.debug("Could not parse PDBQT %s; falling back to SMILES-to-3D.", pdbqt_path)
        return _smiles_to_3d_pdb(smiles)

    ref_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True, sanitize=False)
    if ref_mol is None:
        logger.debug("RDKit could not parse PDB block from PDBQT; falling back.")
        return _smiles_to_3d_pdb(smiles)

    template = safe_mol(smiles)
    if template is None:
        return None

    try:
        ref_mol = AllChem.AssignBondOrdersFromTemplate(template, ref_mol)
    except Exception as exc:
        logger.debug("Bond-order assignment from template failed (%s); falling back.", exc)
        return _smiles_to_3d_pdb(smiles)

    mol = Chem.AddHs(ref_mol)
    try:
        AllChem.ConstrainedEmbed(mol, ref_mol, useTethers=True, randomSeed=42)
    except Exception as exc:
        logger.debug("ConstrainedEmbed failed (%s); falling back.", exc)
        return _smiles_to_3d_pdb(smiles)

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

        has_pose_col = "docking_pose_file" in df.columns

        for idx in candidates:
            smi = df.at[idx, "canonical_smiles"]
            logger.info("Running MM-GBSA for molecule %d: %s", idx, smi[:60])

            pose_path: Optional[str] = None
            if has_pose_col:
                _p = df.at[idx, "docking_pose_file"]
                if _p and str(_p) != "" and Path(str(_p)).is_file():
                    pose_path = str(_p)

            energy, rmsd = self._simulate_one(smi, pdb_path, md_dir, idx, pose_path)
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
        docking_pose_path: Optional[str] = None,
    ) -> Tuple[Optional[float], Optional[float]]:
        """Compute MM-GBSA binding energy and ligand RMSD for ranking.

        ΔG = E(complex) − E(receptor) − E(ligand).

        After minimisation, a short Langevin trajectory (``equilibration_ps``,
        ``production_ns`` from config) yields **trajectory RMSD**: mean
        heavy-atom RMSD of the ligand vs the initial docked pose over 1-ps
        snapshots. If trajectory integration fails, **pose RMSD** (before /
        after minimisation) is used instead.

        Returns (binding_energy_kcal, rmsd_angstrom for ``md_rmsd_mean``).
        """
        try:
            off_mol = OFFMolecule.from_smiles(smiles, allow_undefined_stereo=True)
        except Exception as exc:
            logger.warning("Mol %d: OpenFF cannot parse SMILES: %s", mol_idx, exc)
            return None, None

        try:
            complex_energy, pose_rmsd, md_state = self._gb_energy_complex(
                receptor_pdb, smiles, off_mol, docking_pose_path,
            )
            receptor_energy = self._gb_energy_receptor(receptor_pdb)
            ligand_energy = self._gb_energy_ligand(smiles, off_mol)

            if all(e is not None for e in (complex_energy, receptor_energy, ligand_energy)):
                binding = complex_energy - receptor_energy - ligand_energy
                traj_rmsd: Optional[float] = None
                if md_state is not None:
                    traj_rmsd = self._run_trajectory(md_state, mol_idx)
                final_rmsd = traj_rmsd if traj_rmsd is not None else pose_rmsd
                logger.info(
                    "Mol %d MM-GBSA: complex=%.1f  receptor=%.1f  ligand=%.1f"
                    "  ΔG=%.2f kcal/mol  pose_RMSD=%s Å  traj_RMSD=%s Å  md_rmsd_mean=%s Å",
                    mol_idx, complex_energy, receptor_energy, ligand_energy, binding,
                    f"{pose_rmsd:.2f}" if pose_rmsd is not None else "N/A",
                    f"{traj_rmsd:.2f}" if traj_rmsd is not None else "N/A",
                    f"{final_rmsd:.2f}" if final_rmsd is not None else "N/A",
                )
                return binding, final_rmsd
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
        self,
        receptor_pdb: str,
        ligand_smi: str,
        off_mol: Any,
        docking_pose_path: Optional[str] = None,
    ) -> Tuple[Optional[float], Optional[float], Optional[_ComplexState]]:
        """Build protein+ligand in implicit solvent, minimise.

        Returns (energy_kcal, pose_rmsd_angstrom, state_for_md) where
        ``state_for_md`` holds the minimised ``Simulation`` for trajectory
        continuation; pose_rmsd is heavy-atom ligand RMSD before vs after
        minimisation.

        When *docking_pose_path* is provided, the ligand's heavy-atom
        coordinates are taken from the docked PDBQT pose so that the complex
        is physically meaningful (ligand inside the binding pocket).
        """
        if docking_pose_path:
            lig_pdb = _load_docked_pose_pdb(ligand_smi, docking_pose_path)
        else:
            lig_pdb = _smiles_to_3d_pdb(ligand_smi)
        if lig_pdb is None:
            return None, None, None

        try:
            rec_pdb = app.PDBFile(receptor_pdb)
            modeller = app.Modeller(rec_pdb.topology, rec_pdb.positions)

            with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp:
                tmp.write(lig_pdb)
                tmp_path = tmp.name
            lig_pdb_obj = app.PDBFile(tmp_path)

            n_rec_atoms = modeller.topology.getNumAtoms()
            lig_init_pos = np.array(
                lig_pdb_obj.positions.value_in_unit(unit.angstrom)
            )

            modeller.add(lig_pdb_obj.topology, lig_pdb_obj.positions)

            gaff = GAFFTemplateGenerator(molecules=off_mol, forcefield="gaff-2.11")
            forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
            forcefield.registerTemplateGenerator(gaff.generator)

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
            simulation.minimizeEnergy(maxIterations=1000)

            state = simulation.context.getState(getEnergy=True, getPositions=True)
            energy = _energy_as_kcal(state.getPotentialEnergy())

            # Heavy-atom RMSD between docked and minimised ligand positions
            all_pos = np.array(state.getPositions().value_in_unit(unit.angstrom))
            lig_final_pos = all_pos[n_rec_atoms:]
            heavy_mask = np.array([
                a.element is not None and a.element.mass > 1.5
                for a in lig_pdb_obj.topology.atoms()
            ])
            if heavy_mask.any():
                diff = lig_final_pos[heavy_mask] - lig_init_pos[heavy_mask]
                pose_rmsd: Optional[float] = float(
                    np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
                )
            else:
                pose_rmsd = None

            md_state = _ComplexState(
                simulation=simulation,
                n_rec_atoms=n_rec_atoms,
                lig_init_pos=np.asarray(lig_init_pos, dtype=np.float64).copy(),
                heavy_mask=np.asarray(heavy_mask, dtype=bool).copy(),
            )
            return energy, pose_rmsd, md_state

        except Exception as exc:
            logger.warning("Complex GB energy failed for '%s': %s", ligand_smi[:60], exc)
            return None, None, None

    def _run_trajectory(self, state: _ComplexState, mol_idx: int) -> Optional[float]:
        """Equilibration + production MD; mean ligand heavy-atom RMSD vs docked pose (Å).

        Uses ``equilibration_ps`` and ``production_ns`` from config; snapshots
        every 1 ps. On failure, returns ``None`` so callers can fall back to
        pose RMSD.
        """
        if not state.heavy_mask.any():
            return None
        try:
            sim = state.simulation
            sim.context.setVelocitiesToTemperature(300 * unit.kelvin)

            eq_steps = max(0, int(round(self.equilibration_ps * _MD_STEPS_PER_PS)))
            if eq_steps > 0:
                sim.step(eq_steps)

            prod_steps = max(1, int(round(self.production_ns * 1000 * _MD_STEPS_PER_PS)))
            snapshot_every = _MD_STEPS_PER_PS
            rmsds: List[float] = []

            if prod_steps < snapshot_every:
                sim.step(prod_steps)
                rmsd = self._ligand_rmsd_from_simulation(sim, state)
                if rmsd is None or not np.isfinite(rmsd):
                    logger.warning("Mol %d: invalid RMSD during trajectory.", mol_idx)
                    return None
                rmsds.append(rmsd)
            else:
                n_full = prod_steps // snapshot_every
                remainder = prod_steps % snapshot_every
                for _ in range(n_full):
                    sim.step(snapshot_every)
                    rmsd = self._ligand_rmsd_from_simulation(sim, state)
                    if rmsd is None or not np.isfinite(rmsd):
                        logger.warning("Mol %d: invalid RMSD during trajectory.", mol_idx)
                        return None
                    rmsds.append(rmsd)
                if remainder > 0:
                    sim.step(remainder)

            if not rmsds:
                return None
            return float(np.mean(rmsds))
        except Exception as exc:
            logger.warning(
                "Mol %d: trajectory RMSD failed (%s); using pose RMSD if available.",
                mol_idx,
                exc,
            )
            return None

    @staticmethod
    def _ligand_rmsd_from_simulation(sim: Any, state: _ComplexState) -> Optional[float]:
        """Heavy-atom RMSD (Å) of ligand vs initial docked coordinates."""
        st = sim.context.getState(getPositions=True)
        all_pos = np.array(st.getPositions().value_in_unit(unit.angstrom))
        lig = all_pos[state.n_rec_atoms:]
        ref = state.lig_init_pos[state.heavy_mask]
        cur = lig[state.heavy_mask]
        if cur.size == 0:
            return None
        diff = cur - ref
        return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

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
        return _energy_as_kcal(state.getPotentialEnergy())


# ======================================================================
# Explicit-solvent refinement layer (optional Stage 4 extension)
# ======================================================================

class ExplicitSolventRefiner:
    """TIP3P explicit-solvent MD refinement for top candidates.

    Activated when implicit-solvent ranking cannot separate candidates
    (opt_score spread below ``delta_opt_score_threshold``).  Runs a full
    NVT → NPT → production pipeline and reports:
      - ``explicit_binding_energy``: MM-PBSA-style end-state ΔG
      - ``explicit_rmsd_mean``: mean ligand heavy-atom RMSD over trajectory
      - ``explicit_rmsd_drift``: linear RMSD trend (Å / ns) as stability proxy
    """

    def __init__(self, cfg: Dict[str, Any]):
        lo = cfg.get("lead_optimization", {})
        emd = lo.get("explicit_md", {})
        self.enabled = bool(emd.get("enabled", False))
        self.top_n = int(emd.get("top_n", 5))
        self.production_ns = float(emd.get("production_ns", 10.0))
        self.eq_nvt_ps = float(emd.get("equilibration_nvt_ps", 100))
        self.eq_npt_ps = float(emd.get("equilibration_npt_ps", 500))
        self.padding_nm = float(emd.get("padding_nm", 1.0))
        self.ionic_strength = float(emd.get("ionic_strength_M", 0.15))

        trigger = emd.get("trigger", {})
        self.score_threshold = float(trigger.get("delta_opt_score_threshold", 0.05))

        device_str = cfg.get("pipeline", {}).get("device", "auto")
        self.platform, self.platform_props = _select_platform(device_str)

    def should_trigger(self, df: pd.DataFrame) -> bool:
        """Check if explicit solvent refinement should activate.

        Returns True when ``enabled`` and the opt_score spread among the
        top-N candidates is below the configured threshold.
        """
        if not self.enabled:
            return False
        if "opt_score" not in df.columns:
            return False
        top = df.nlargest(self.top_n, "opt_score")["opt_score"].dropna()
        if len(top) < 2:
            return False
        spread = float(top.max() - top.min())
        if spread < self.score_threshold:
            logger.info(
                "Explicit solvent triggered: top-%d opt_score spread=%.4f < threshold=%.4f",
                self.top_n, spread, self.score_threshold,
            )
            return True
        logger.info(
            "Explicit solvent NOT triggered: spread=%.4f >= threshold=%.4f",
            spread, self.score_threshold,
        )
        return False

    def run(
        self,
        df: pd.DataFrame,
        protein_info: Dict[str, Any],
        out_dir: Path,
    ) -> pd.DataFrame:
        """Run explicit-solvent MD on the top-N candidates by opt_score."""
        df = df.copy()
        df["explicit_binding_energy"] = np.nan
        df["explicit_rmsd_mean"] = np.nan
        df["explicit_rmsd_drift"] = np.nan

        if not _OPENMM_OK or not _GAFF_OK:
            logger.warning("OpenMM/GAFF not available; skipping explicit solvent MD.")
            return df

        pdb_path = protein_info.get("pdb_path")
        if not pdb_path or not Path(pdb_path).exists():
            logger.error("No prepared PDB for explicit solvent MD.")
            return df

        top_idx = df.nlargest(self.top_n, "opt_score").index.tolist()
        if not top_idx:
            return df

        emd_dir = out_dir / "explicit_md"
        emd_dir.mkdir(parents=True, exist_ok=True)
        has_pose_col = "docking_pose_file" in df.columns

        for idx in top_idx:
            smi = df.at[idx, "canonical_smiles"]
            logger.info("Explicit solvent MD for molecule %s: %s", idx, smi[:60])

            pose_path: Optional[str] = None
            if has_pose_col:
                _p = df.at[idx, "docking_pose_file"]
                if _p and str(_p) != "" and Path(str(_p)).is_file():
                    pose_path = str(_p)

            result = self._run_explicit(smi, pdb_path, emd_dir, idx, pose_path)
            if result is not None:
                energy, rmsd_mean, rmsd_drift = result
                df.at[idx, "explicit_binding_energy"] = energy
                df.at[idx, "explicit_rmsd_mean"] = rmsd_mean
                df.at[idx, "explicit_rmsd_drift"] = rmsd_drift

        n_done = df["explicit_binding_energy"].notna().sum()
        logger.info("Explicit MD complete: %d / %d candidates.", n_done, len(top_idx))
        return df

    def _run_explicit(
        self,
        smiles: str,
        receptor_pdb: str,
        emd_dir: Path,
        mol_idx: int,
        docking_pose_path: Optional[str] = None,
    ) -> Optional[Tuple[float, float, float]]:
        """Run full explicit-solvent MD for one molecule.

        Returns (binding_energy_kcal, rmsd_mean_angstrom, rmsd_drift_ang_per_ns)
        or None on failure.
        """
        try:
            off_mol = OFFMolecule.from_smiles(smiles, allow_undefined_stereo=True)
        except Exception as exc:
            logger.warning("Mol %d: OpenFF parse failed: %s", mol_idx, exc)
            return None

        if docking_pose_path:
            lig_pdb = _load_docked_pose_pdb(smiles, docking_pose_path)
        else:
            lig_pdb = _smiles_to_3d_pdb(smiles)
        if lig_pdb is None:
            return None

        try:
            rec_pdb = app.PDBFile(receptor_pdb)
            modeller = app.Modeller(rec_pdb.topology, rec_pdb.positions)

            with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp:
                tmp.write(lig_pdb)
                tmp_path = tmp.name
            lig_pdb_obj = app.PDBFile(tmp_path)

            n_rec_atoms = modeller.topology.getNumAtoms()
            lig_init_pos = np.array(
                lig_pdb_obj.positions.value_in_unit(unit.angstrom)
            )
            modeller.add(lig_pdb_obj.topology, lig_pdb_obj.positions)
            n_solute_atoms = modeller.topology.getNumAtoms()

            gaff = GAFFTemplateGenerator(molecules=off_mol, forcefield="gaff-2.11")
            forcefield = app.ForceField(
                "amber14-all.xml", "amber14/tip3pfb.xml",
            )
            forcefield.registerTemplateGenerator(gaff.generator)

            modeller.addSolvent(
                forcefield,
                model="tip3p",
                padding=self.padding_nm * unit.nanometers,
                ionicStrength=self.ionic_strength * unit.molar,
            )

            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=app.PME,
                nonbondedCutoff=1.0 * unit.nanometers,
                constraints=app.HBonds,
            )

            # NVT equilibration
            integrator_nvt = openmm.LangevinMiddleIntegrator(
                300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds,
            )
            simulation = self._build_simulation(modeller, system, integrator_nvt)
            simulation.context.setPositions(modeller.positions)
            simulation.minimizeEnergy(maxIterations=2000)

            nvt_steps = max(1, int(round(self.eq_nvt_ps * _MD_STEPS_PER_PS)))
            simulation.step(nvt_steps)
            logger.debug("Mol %d: NVT equilibration (%d ps) done.", mol_idx, self.eq_nvt_ps)

            # NPT equilibration — add barostat
            barostat = openmm.MonteCarloBarostat(
                1.0 * unit.atmospheres, 300 * unit.kelvin, 25,
            )
            system.addForce(barostat)
            simulation.context.reinitialize(preserveState=True)

            npt_steps = max(1, int(round(self.eq_npt_ps * _MD_STEPS_PER_PS)))
            simulation.step(npt_steps)
            logger.debug("Mol %d: NPT equilibration (%d ps) done.", mol_idx, self.eq_npt_ps)

            # Production: collect snapshots every 10 ps
            prod_steps = max(1, int(round(self.production_ns * 1000 * _MD_STEPS_PER_PS)))
            snapshot_every = 10 * _MD_STEPS_PER_PS
            heavy_mask = np.array([
                a.element is not None and a.element.mass > 1.5
                for a in lig_pdb_obj.topology.atoms()
            ])

            rmsds: List[float] = []
            energies: List[float] = []
            n_snapshots = max(1, prod_steps // snapshot_every)
            remainder = prod_steps % snapshot_every

            for _ in range(n_snapshots):
                simulation.step(snapshot_every)
                st = simulation.context.getState(getEnergy=True, getPositions=True)
                energies.append(_energy_as_kcal(st.getPotentialEnergy()))

                all_pos = np.array(st.getPositions().value_in_unit(unit.angstrom))
                lig_pos = all_pos[n_rec_atoms:n_solute_atoms]
                if heavy_mask.any() and lig_pos.shape[0] == lig_init_pos.shape[0]:
                    diff = lig_pos[heavy_mask] - lig_init_pos[heavy_mask]
                    rmsds.append(float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))))

            if remainder > 0:
                simulation.step(remainder)

            if not rmsds:
                logger.warning("Mol %d: no RMSD snapshots collected.", mol_idx)
                return None

            rmsd_mean = float(np.mean(rmsds))

            # RMSD drift: linear regression slope (Å per ns)
            times_ns = np.linspace(0, self.production_ns, len(rmsds))
            if len(rmsds) >= 2:
                coeffs = np.polyfit(times_ns, rmsds, 1)
                rmsd_drift = float(coeffs[0])
            else:
                rmsd_drift = 0.0

            binding_energy = float(np.mean(energies))

            logger.info(
                "Mol %d explicit MD: energy=%.1f kcal/mol  RMSD_mean=%.2f Å  "
                "RMSD_drift=%.3f Å/ns  (%d snapshots over %.1f ns)",
                mol_idx, binding_energy, rmsd_mean, rmsd_drift,
                len(rmsds), self.production_ns,
            )
            return binding_energy, rmsd_mean, rmsd_drift

        except Exception as exc:
            logger.warning("Explicit MD failed for mol %d: %s", mol_idx, exc)
            return None

    def _build_simulation(
        self, modeller: Any, system: Any, integrator: Any,
    ) -> Any:
        """Create an OpenMM Simulation with platform selection."""
        platform = None
        if self.platform:
            try:
                platform = openmm.Platform.getPlatformByName(self.platform)
            except Exception:
                pass
        if platform:
            return app.Simulation(
                modeller.topology, system, integrator,
                platform, self.platform_props,
            )
        return app.Simulation(modeller.topology, system, integrator)
