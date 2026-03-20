"""Molecular docking module using AutoDock Vina.

Converts SMILES to 3-D conformers, docks them against a prepared receptor,
and returns docking scores (kcal/mol).  Supports molecule-level parallelism
to fully utilize multi-core machines.

Optional dependencies:
  - vina   (pip install vina)
  - meeko  (pip install meeko)
"""
# EN: Module overview and key intent for maintainers.
# 中文：模块总览与关键设计意图，便于后续维护。

# 分子对接模块（AutoDock Vina）：SMILES → 3D 构象 → 对接 → 打分，支持分子级并行。

from __future__ import annotations

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

from drugpipe.utils.chem import safe_mol

logger = logging.getLogger(__name__)

try:
    from vina import Vina
    _VINA_OK = True
except ImportError:
    _VINA_OK = False

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    _MEEKO_OK = True
except ImportError:
    _MEEKO_OK = False


# EN: VinaDocking core behavior and intent.
# 中文：VinaDocking 的核心行为与设计意图。
class VinaDocking:
    """Batch molecular docking with AutoDock Vina.

    When the machine has more CPU cores than ``exhaustiveness``, multiple
    molecules are docked **in parallel** so that total CPU utilization
    stays near 100%.
    """

    # EN: __init__ core behavior and intent.
    # 中文：__init__ 的核心行为与设计意图。
    def __init__(self, cfg: Dict[str, Any]):
        lo = cfg.get("lead_optimization", {})
        dk = lo.get("docking", {})
        self.enabled = bool(dk.get("enabled", True))
        self.exhaustiveness = int(dk.get("exhaustiveness", 32))
        self.n_poses = int(dk.get("n_poses", 5))
        self.score_threshold = float(dk.get("score_threshold", -6.0))
        self.n_cpus = int(dk.get("n_cpus", 0))  # 0 = auto-detect

    # ------------------------------------------------------------------
    # EN: dock_leads core behavior and intent.
    # 中文：dock_leads 的核心行为与设计意图。
    def dock_leads(
        self,
        df: pd.DataFrame,
        protein_info: Dict[str, Any],
        out_dir: Path,
    ) -> pd.DataFrame:
        """Dock every lead and add ``docking_score`` column."""
        df = df.copy()
        df["docking_score"] = np.nan
        df["docking_pose_file"] = ""

        if not self.enabled:
            logger.info("Docking disabled in config.")
            return df

        if not _VINA_OK:
            logger.warning("AutoDock Vina not installed (pip install vina). Skipping docking.")
            return df

        receptor_pdbqt = protein_info.get("pdbqt_path")
        if not receptor_pdbqt or not Path(receptor_pdbqt).exists():
            logger.error("No receptor PDBQT available for docking. "
                         "Ensure protein preparation produced a .pdbqt file.")
            return df

        center = protein_info["center"]
        box_size = protein_info["box_size"]

        total_cpus = self.n_cpus or os.cpu_count() or 4
        cpus_per_mol = min(total_cpus, self.exhaustiveness)
        n_parallel = max(1, total_cpus // cpus_per_mol)

        logger.info(
            "Docking config: exhaustiveness=%d, cpus_per_mol=%d, "
            "parallel_workers=%d (total_cpus=%d)",
            self.exhaustiveness, cpus_per_mol, n_parallel, total_cpus,
        )

        poses_dir = out_dir / "docking_poses"
        poses_dir.mkdir(parents=True, exist_ok=True)

        jobs = [
            (idx, df.at[idx, "canonical_smiles"])
            for idx in df.index
        ]
        n_total = len(jobs)
        n_docked = 0
        n_failed = 0

        if n_parallel <= 1:
            for i, (idx, smi) in enumerate(jobs):
                score, pose_path = _dock_one_standalone(
                    smi, receptor_pdbqt, center, box_size,
                    str(poses_dir), idx, cpus_per_mol,
                    self.exhaustiveness, self.n_poses,
                )
                self._store_result(df, idx, score, pose_path)
                n_docked, n_failed = self._update_counts(
                    score, n_docked, n_failed, i, n_total,
                )
        else:
            with ProcessPoolExecutor(max_workers=n_parallel) as pool:
                futures = {
                    pool.submit(
                        _dock_one_standalone,
                        smi, receptor_pdbqt, center, box_size,
                        str(poses_dir), idx, cpus_per_mol,
                        self.exhaustiveness, self.n_poses,
                    ): idx
                    for idx, smi in jobs
                }
                for i, future in enumerate(as_completed(futures)):
                    idx = futures[future]
                    try:
                        score, pose_path = future.result()
                    except Exception as exc:
                        logger.debug("Docking worker error for idx %d: %s", idx, exc)
                        score, pose_path = None, None
                    self._store_result(df, idx, score, pose_path)
                    n_docked, n_failed = self._update_counts(
                        score, n_docked, n_failed, i, n_total,
                    )

        passed = (df["docking_score"] <= self.score_threshold).sum()
        logger.info(
            "Docking complete: %d/%d docked, %d passed threshold (%.1f kcal/mol).",
            n_docked, n_total, passed, self.score_threshold,
        )
        return df

    # ------------------------------------------------------------------
    @staticmethod
    # EN: _store_result core behavior and intent.
    # 中文：_store_result 的核心行为与设计意图。
    def _store_result(
        df: pd.DataFrame, idx: int,
        score: Optional[float], pose_path: Optional[str],
    ) -> None:
        if score is not None:
            df.at[idx, "docking_score"] = score
            df.at[idx, "docking_pose_file"] = pose_path or ""

    @staticmethod
    # EN: _update_counts core behavior and intent.
    # 中文：_update_counts 的核心行为与设计意图。
    def _update_counts(
        score: Optional[float],
        n_docked: int, n_failed: int, i: int, n_total: int,
    ) -> Tuple[int, int]:
        if score is not None:
            n_docked += 1
        else:
            n_failed += 1
        if (i + 1) % 5 == 0 or (i + 1) == n_total:
            logger.info("Docking progress: %d docked, %d failed / %d total",
                        n_docked, n_failed, n_total)
        return n_docked, n_failed


# ======================================================================
# Top-level function for multiprocessing (must be picklable)
# ======================================================================

# EN: _dock_one_standalone core behavior and intent.
# 中文：_dock_one_standalone 的核心行为与设计意图。
def _dock_one_standalone(
    smiles: str,
    receptor_pdbqt: str,
    center: List[float],
    box_size: List[float],
    poses_dir: str,
    mol_idx: int,
    cpus: int,
    exhaustiveness: int,
    n_poses: int,
) -> Tuple[Optional[float], Optional[str]]:
    """Dock a single SMILES. Returns (best_score, pose_path_str)."""
    mol = safe_mol(smiles)
    if mol is None:
        return None, None

    mol3d = _embed_3d(mol)
    if mol3d is None:
        return None, None

    lig_pdbqt = _mol_to_pdbqt(mol3d)
    if lig_pdbqt is None:
        return None, None

    try:
        v = Vina(sf_name="vina", cpu=cpus, verbosity=0)
        v.set_receptor(rigid_pdbqt_filename=receptor_pdbqt)
        v.set_ligand_from_string(lig_pdbqt)
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        energies = v.energies()
        best_score = float(energies[0][0]) if len(energies) > 0 else None

        pose_path = str(Path(poses_dir) / f"pose_{mol_idx}.pdbqt")
        v.write_poses(pdbqt_filename=pose_path, n_poses=1, overwrite=True)

        return best_score, pose_path

    except Exception:
        return None, None


# ======================================================================
# Standalone helpers (module-level for pickling by ProcessPoolExecutor)
# ======================================================================

# EN: _embed_3d core behavior and intent.
# 中文：_embed_3d 的核心行为与设计意图。
def _embed_3d(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Generate a 3-D conformer via ETKDG."""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        params2 = AllChem.ETKDGv3()
        params2.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params2) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    return mol


# EN: _mol_to_pdbqt core behavior and intent.
# 中文：_mol_to_pdbqt 的核心行为与设计意图。
def _mol_to_pdbqt(mol: Chem.Mol) -> Optional[str]:
    """Convert an RDKit Mol with 3-D coords to a PDBQT string for Vina."""
    if _MEEKO_OK:
        try:
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            if mol_setups:
                pdbqt_str, is_ok, _ = PDBQTWriterLegacy.write_string(
                    mol_setups[0]
                )
                if is_ok:
                    return pdbqt_str
        except Exception:
            pass

    try:
        pdb_block = Chem.MolToPDBBlock(mol)
        if pdb_block:
            return _pdb_to_pdbqt_naive(pdb_block)
    except Exception:
        pass

    return None


# EN: _pdb_to_pdbqt_naive core behavior and intent.
# 中文：_pdb_to_pdbqt_naive 的核心行为与设计意图。
def _pdb_to_pdbqt_naive(pdb_block: str) -> str:
    """Minimal PDB → PDBQT with correct fixed-width column layout."""
    ad4 = {"C": "C", "N": "NA", "O": "OA", "S": "SA", "H": "HD",
           "F": "F", "P": "P", "I": "I", "CL": "Cl", "BR": "Br"}
    lines = []
    for line in pdb_block.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            elem = line[76:78].strip() if len(line) >= 78 else ""
            if not elem:
                elem = line[12:16].strip().lstrip("0123456789")[:1]
            ad_type = ad4.get(elem.upper(), elem.upper()[:2] if len(elem) > 1 else "C")
            base = line[:54].ljust(54)
            occ_bfact = line[54:66].ljust(12) if len(line) >= 66 else "  1.00  0.00"
            new_line = f"{base}{occ_bfact}    {0.000:+.3f} {ad_type:<2s}"
            lines.append(new_line)
        elif line.startswith("END"):
            lines.append(line)
    return "\n".join(lines)
