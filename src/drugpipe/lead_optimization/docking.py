"""Molecular docking module using AutoDock Vina.

Converts SMILES to 3-D conformers, docks them against a prepared receptor,
and returns docking scores (kcal/mol).

Optional dependencies:
  - vina   (pip install vina)
  - meeko  (pip install meeko)
"""
# 分子对接模块（AutoDock Vina）：SMILES → 3D 构象 → 对接 → 打分。

from __future__ import annotations

import logging
import os
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


class VinaDocking:
    """Batch molecular docking with AutoDock Vina."""

    def __init__(self, cfg: Dict[str, Any]):
        lo = cfg.get("lead_optimization", {})
        dk = lo.get("docking", {})
        self.enabled = bool(dk.get("enabled", True))
        self.exhaustiveness = int(dk.get("exhaustiveness", 32))
        self.n_poses = int(dk.get("n_poses", 5))
        self.score_threshold = float(dk.get("score_threshold", -6.0))
        self.n_cpus = int(dk.get("n_cpus", 0))  # 0 = auto-detect

    # ------------------------------------------------------------------
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

        cpus = self.n_cpus or os.cpu_count() or 4
        logger.info("Docking config: exhaustiveness=%d, n_poses=%d, cpus=%d",
                     self.exhaustiveness, self.n_poses, cpus)

        poses_dir = out_dir / "docking_poses"
        poses_dir.mkdir(parents=True, exist_ok=True)

        n_total = len(df)
        n_docked = 0
        n_failed = 0
        for i, idx in enumerate(df.index):
            smi = df.at[idx, "canonical_smiles"]
            score, pose_path = self._dock_one(
                smi, receptor_pdbqt, center, box_size, poses_dir, idx, cpus,
            )
            if score is not None:
                df.at[idx, "docking_score"] = score
                df.at[idx, "docking_pose_file"] = str(pose_path) if pose_path else ""
                n_docked += 1
            else:
                n_failed += 1

            if (i + 1) % 5 == 0 or (i + 1) == n_total:
                logger.info("Docking progress: %d/%d docked, %d failed",
                            n_docked, n_failed, n_total)

        passed = (df["docking_score"] <= self.score_threshold).sum()
        logger.info(
            "Docking complete: %d/%d docked, %d passed threshold (%.1f kcal/mol).",
            n_docked, n_total, passed, self.score_threshold,
        )
        return df

    # ------------------------------------------------------------------
    def _dock_one(
        self,
        smiles: str,
        receptor_pdbqt: str,
        center: List[float],
        box_size: List[float],
        poses_dir: Path,
        mol_idx: int,
        cpus: int,
    ) -> Tuple[Optional[float], Optional[Path]]:
        """Dock a single SMILES. Returns (best_score, pose_path)."""
        mol = safe_mol(smiles)
        if mol is None:
            logger.debug("Mol %d: invalid SMILES", mol_idx)
            return None, None

        mol3d = self._embed_3d(mol)
        if mol3d is None:
            logger.debug("Mol %d: 3D embedding failed", mol_idx)
            return None, None

        lig_pdbqt = self._mol_to_pdbqt(mol3d)
        if lig_pdbqt is None:
            logger.debug("Mol %d: PDBQT conversion failed", mol_idx)
            return None, None

        try:
            v = Vina(sf_name="vina", cpu=cpus, verbosity=0)
            v.set_receptor(rigid_pdbqt_filename=receptor_pdbqt)
            v.set_ligand_from_string(lig_pdbqt)
            v.compute_vina_maps(center=center, box_size=box_size)
            v.dock(exhaustiveness=self.exhaustiveness, n_poses=self.n_poses)

            energies = v.energies()
            best_score = float(energies[0][0]) if len(energies) > 0 else None

            pose_path = poses_dir / f"pose_{mol_idx}.pdbqt"
            v.write_poses(pdbqt_filename=str(pose_path), n_poses=1, overwrite=True)

            logger.debug("Mol %d: score=%.2f kcal/mol", mol_idx, best_score)
            return best_score, pose_path

        except Exception as exc:
            logger.debug("Docking failed for mol %d (%s…): %s",
                         mol_idx, smiles[:30], exc)
            return None, None

    # ------------------------------------------------------------------
    @staticmethod
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

    @staticmethod
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
            except Exception as exc:
                logger.debug("Meeko ligand prep failed: %s", exc)

        # Fallback: manual PDB→PDBQT
        try:
            pdb_block = Chem.MolToPDBBlock(mol)
            if pdb_block:
                return _pdb_to_pdbqt_naive(pdb_block)
        except Exception:
            pass

        return None


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
