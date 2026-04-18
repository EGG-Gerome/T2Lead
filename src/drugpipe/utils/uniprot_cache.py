"""Local UniProt Swiss-Prot FASTA + ID-mapping cache.

Eliminates redundant network requests when the same UniProt accession or
gene symbol is queried across multiple variants.  The cache lazily loads
from on-disk files (download once, reuse forever).

Expected files under *db_dir*:
  - uniprot_sprot.fasta(.gz)         Swiss-Prot FASTA (~89 MB gzipped)
  - HUMAN_9606_idmapping.dat(.gz)    Human gene→UniProt→PDB mapping
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


class UniProtLocalCache:
    """Index-free, memory-mapped local lookup for Swiss-Prot sequences."""

    def __init__(self, db_dir: Optional[str] = None):
        self._seq_index: Dict[str, str] = {}
        self._gene_to_uniprot: Dict[str, List[str]] = {}
        self._uniprot_to_gene: Dict[str, str] = {}
        self._loaded = False
        self._db_dir = Path(db_dir) if db_dir else None

    @property
    def available(self) -> bool:
        if self._db_dir is None:
            return False
        return (
            self._find_file("uniprot_sprot.fasta") is not None
            or self._find_file("HUMAN_9606_idmapping.dat") is not None
        )

    def get_sequence(self, uniprot_id: str) -> Optional[str]:
        self._ensure_loaded()
        return self._seq_index.get(uniprot_id.strip().upper())

    def gene_to_uniprot_ids(self, gene_symbol: str) -> List[str]:
        self._ensure_loaded()
        return self._gene_to_uniprot.get(gene_symbol.strip().upper(), [])

    def uniprot_to_gene(self, uniprot_id: str) -> Optional[str]:
        self._ensure_loaded()
        return self._uniprot_to_gene.get(uniprot_id.strip().upper())

    def _ensure_loaded(self) -> None:
        if self._loaded or self._db_dir is None:
            return
        self._loaded = True
        self._load_fasta()
        self._load_idmapping()

    def _load_fasta(self) -> None:
        path = self._find_file("uniprot_sprot.fasta")
        if path is None:
            return
        logger.info("Loading Swiss-Prot FASTA from %s ...", path)
        opener = gzip.open if str(path).endswith(".gz") else open
        current_id = ""
        chunks: List[str] = []
        count = 0
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if current_id and chunks:
                        self._seq_index[current_id] = "".join(chunks)
                        count += 1
                    parts = line[1:].split("|")
                    current_id = parts[1].strip().upper() if len(parts) >= 2 else ""
                    chunks = []
                else:
                    chunks.append(line.strip())
            if current_id and chunks:
                self._seq_index[current_id] = "".join(chunks)
                count += 1
        logger.info("Indexed %d Swiss-Prot sequences.", count)

    def _load_idmapping(self) -> None:
        path = self._find_file("HUMAN_9606_idmapping.dat")
        if path is None:
            return
        logger.info("Loading human ID mapping from %s ...", path)
        opener = gzip.open if str(path).endswith(".gz") else open
        gene_rows = 0
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                uid, db, value = parts[0].strip().upper(), parts[1].strip(), parts[2].strip()
                if db == "Gene_Name":
                    gene = value.upper()
                    self._gene_to_uniprot.setdefault(gene, [])
                    if uid not in self._gene_to_uniprot[gene]:
                        self._gene_to_uniprot[gene].append(uid)
                    if uid not in self._uniprot_to_gene:
                        self._uniprot_to_gene[uid] = gene
                    gene_rows += 1
        logger.info("Loaded %d gene→UniProt mappings.", gene_rows)

    def _find_file(self, basename: str) -> Optional[Path]:
        if self._db_dir is None:
            return None
        for suffix in ("", ".gz"):
            p = self._db_dir / f"{basename}{suffix}"
            if p.is_file():
                return p
        return None
