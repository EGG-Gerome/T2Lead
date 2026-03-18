"""Protein structure retrieval and preparation for molecular docking / MD.

Fetches a PDB structure from RCSB, cleans it with PDBFixer (optional),
identifies the binding site from a co-crystallized ligand, and converts
to PDBQT format for AutoDock Vina.

Optional dependencies:
  - pdbfixer  (conda install -c conda-forge pdbfixer)
  - meeko     (pip install meeko)
"""
# 蛋白结构获取与准备：从 RCSB 下载 PDB → PDBFixer 清洗 → 识别结合位点 → 转换为 PDBQT。

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import requests

logger = logging.getLogger(__name__)

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    _PDBFIXER_OK = True
except ImportError:
    _PDBFIXER_OK = False


# AutoDock 4 atom-type mapping (element → AD4 type)
_AD4_TYPE = {
    "C": "C", "N": "N", "O": "OA", "S": "SA", "H": "HD",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
    "ZN": "Zn", "FE": "Fe", "MG": "Mg", "CA": "Ca", "MN": "Mn",
}


class ProteinPreparator:
    """Fetch, fix, and prepare a protein structure for docking / MD."""

    def __init__(self, cfg: Dict[str, Any]):
        lo = cfg.get("lead_optimization", {})
        self.pdb_id = (lo.get("pdb_id", "") or "").strip().upper()
        bs = lo.get("binding_site", {})
        self.auto_detect = bool(bs.get("auto_detect", True))
        self._cfg_center: List[float] = bs.get("center", [0.0, 0.0, 0.0])
        self._cfg_box: List[float] = bs.get("box_size", [25.0, 25.0, 25.0])

    # ------------------------------------------------------------------
    def prepare(self, out_dir: Path) -> Optional[Dict[str, Any]]:
        """Download PDB, fix, detect binding site, write PDBQT receptor.

        Returns a dict with paths and binding-site geometry, or None on
        failure / missing PDB ID.
        """
        if not self.pdb_id:
            logger.warning("No pdb_id configured — skipping protein preparation.")
            return None

        pdb_path = self._fetch_pdb(out_dir)
        if pdb_path is None:
            return None

        center, box_size = self._detect_binding_site(pdb_path)

        fixed_pdb = self._fix_structure(pdb_path, out_dir)

        pdbqt_path = self._pdb_to_pdbqt(fixed_pdb, out_dir)

        return {
            "pdb_path": str(fixed_pdb),
            "pdbqt_path": str(pdbqt_path) if pdbqt_path else None,
            "center": center,
            "box_size": box_size,
        }

    # ------------------------------------------------------------------
    def _fetch_pdb(self, out_dir: Path) -> Optional[Path]:
        """Download structure from RCSB."""
        url = f"https://files.rcsb.org/download/{self.pdb_id}.pdb"
        dest = out_dir / f"{self.pdb_id}.pdb"
        if dest.exists():
            logger.info("PDB file already cached: %s", dest)
            return dest

        logger.info("Downloading PDB %s from RCSB ...", self.pdb_id)
        try:
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
        except Exception as exc:
            logger.error("Failed to download PDB %s: %s", self.pdb_id, exc)
            return None

        dest.write_text(resp.text, encoding="utf-8")
        logger.info("PDB saved: %s (%d bytes)", dest, len(resp.text))
        return dest

    # ------------------------------------------------------------------
    def _fix_structure(self, pdb_path: Path, out_dir: Path) -> Path:
        """Use PDBFixer to add missing *side-chain* atoms and hydrogens.

        Skips missing-residue (loop) rebuilding to avoid the very slow
        modelling step that can take 10+ minutes for large structures.
        """
        fixed = out_dir / f"{self.pdb_id}_fixed.pdb"
        if fixed.exists():
            logger.info("Fixed PDB already cached: %s", fixed)
            return fixed

        if not _PDBFIXER_OK:
            logger.warning("PDBFixer not installed — using raw PDB without repair.")
            return pdb_path

        logger.info("Fixing PDB with PDBFixer (skip loop building) ...")
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.removeHeterogens(keepWater=False)
        # Only fix missing side-chain atoms, NOT missing residues (loops)
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        with open(fixed, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh)
        logger.info("Fixed PDB saved: %s", fixed)
        return fixed

    # ------------------------------------------------------------------
    def _detect_binding_site(
        self, pdb_path: Path
    ) -> Tuple[List[float], List[float]]:
        """Return (center_xyz, box_size_xyz).

        If ``auto_detect`` is True, parse HETATM records of the original
        PDB to find co-crystallized ligands and compute a bounding box.
        """
        if not self.auto_detect:
            return list(self._cfg_center), list(self._cfg_box)

        coords = self._hetatm_coords(pdb_path)
        if len(coords) == 0:
            logger.warning(
                "No HETATM ligand atoms found for auto-detect; "
                "using config center/box."
            )
            return list(self._cfg_center), list(self._cfg_box)

        arr = np.array(coords)
        cmin = arr.min(axis=0)
        cmax = arr.max(axis=0)
        center = ((cmin + cmax) / 2.0).tolist()
        padding = 10.0
        box = (cmax - cmin + padding).tolist()
        box = [max(min(b, 30.0), 20.0) for b in box]  # clamp to [20, 30] Å

        logger.info("Auto-detected binding site center=%.1f,%.1f,%.1f  box=%.1f,%.1f,%.1f",
                     *center, *box)
        return center, box

    @staticmethod
    def _hetatm_coords(pdb_path: Path) -> List[List[float]]:
        """Extract XYZ of HETATM atoms that are not water (HOH/WAT)."""
        skip = {"HOH", "WAT", "DOD"}
        coords: List[List[float]] = []
        for line in pdb_path.read_text().splitlines():
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if resname in skip:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue
        return coords

    # ------------------------------------------------------------------
    def _pdb_to_pdbqt(self, pdb_path: Path, out_dir: Path) -> Optional[Path]:
        """Convert fixed PDB receptor to PDBQT for AutoDock Vina.

        PDBQT format (fixed-width):
          cols  1-6  : record type
          cols  7-66 : same as PDB (atom name, residue, coords, occ, bfact)
          cols 67-76 : blank (unused in PDBQT) or padding
          cols 71-76 : partial charge (Gasteiger, here 0.000)
          cols 77-78 : blank
          cols 79-80 : AutoDock4 atom type (right-justified)
        """
        pdbqt = out_dir / f"{self.pdb_id}_receptor.pdbqt"
        if pdbqt.exists():
            logger.info("Receptor PDBQT already cached: %s", pdbqt)
            return pdbqt

        logger.info("Converting receptor PDB → PDBQT ...")
        lines_out: List[str] = []

        for line in pdb_path.read_text().splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                if line.startswith(("END", "TER")):
                    lines_out.append(line)
                continue

            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                atom_name = line[12:16].strip()
                element = atom_name.lstrip("0123456789")[:1].upper()

            ad_type = _ad4_atom_type(element, line[12:16].strip())

            # PDBQT: first 66 chars from PDB, then q and type at fixed positions
            base = line[:54].ljust(54)               # coords end at col 54
            occ_bfact = line[54:66] if len(line) >= 66 else "  1.00  0.00"
            occ_bfact = occ_bfact.ljust(12)

            # cols 67-76: partial charge (right-justified, 6.3f)
            # cols 77-78: space + AD type (left-justified, 2 chars)
            pdbqt_line = f"{base}{occ_bfact}    {0.000:+.3f} {ad_type:<2s}"
            lines_out.append(pdbqt_line)

        if not lines_out:
            logger.error("No ATOM records found in %s", pdb_path)
            return None

        pdbqt.write_text("\n".join(lines_out) + "\n", encoding="utf-8")
        logger.info("Receptor PDBQT written: %s (%d atom lines)", pdbqt, len(lines_out))
        return pdbqt


def _ad4_atom_type(element: str, atom_name: str) -> str:
    """Map element + atom name to a valid AutoDock 4 atom type."""
    el = element.upper()
    if el == "H":
        # Polar hydrogen: bonded to N or O (name starts with H + N/O letter)
        if any(atom_name.startswith(p) for p in ("HN", "HO", "HE", "HH", "HG1")):
            return "HD"
        return "H"
    if el == "C":
        return "C"
    if el == "N":
        return "NA"
    if el == "O":
        return "OA"
    if el == "S":
        return "SA"
    return _AD4_TYPE.get(el, el[:2] if len(el) >= 2 else "C")
