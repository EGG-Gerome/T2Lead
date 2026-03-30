"""Generate mutant protein sequences from somatic variant annotations.

Takes ``SomaticVariant`` entries (from ``vcf_parser``) and fetches the
wildtype protein sequence from UniProt or Ensembl, then applies the
amino-acid substitution to produce a mutant FASTA ready for structure
prediction (ESMFold / ColabFold).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from drugpipe.variant_analysis.vcf_parser import SomaticVariant

logger = logging.getLogger(__name__)

_UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
_ENSEMBL_API = "https://rest.ensembl.org"

# Standard 3-letter → 1-letter amino-acid mapping
_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}


@dataclass
class MutantProtein:
    """A mutant protein sequence derived from a somatic variant."""

    variant: SomaticVariant
    gene_symbol: str
    uniprot_id: Optional[str]
    wildtype_sequence: str
    mutant_sequence: str
    mutation_label: str

    @property
    def fasta_header(self) -> str:
        return f">{self.gene_symbol}_{self.mutation_label} uniprot={self.uniprot_id or 'unknown'}"

    def to_fasta(self) -> str:
        lines = [self.fasta_header]
        seq = self.mutant_sequence
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
        return "\n".join(lines) + "\n"


class MutantSequenceBuilder:
    """Build mutant protein sequences from somatic variants."""

    def __init__(self, cfg: Optional[Dict[str, Any]] = None):
        va = (cfg or {}).get("variant_analysis", {})
        self._timeout = int(va.get("api_timeout_s", 30))
        self._gene_to_uniprot: Dict[str, str] = {}

    def build(self, variants: List[SomaticVariant]) -> List[MutantProtein]:
        """Generate mutant sequences for a list of somatic variants."""
        results: List[MutantProtein] = []
        for var in variants:
            mp = self._build_one(var)
            if mp is not None:
                results.append(mp)
        logger.info("Built %d mutant sequences from %d variants.", len(results), len(variants))
        return results

    def write_fastas(
        self, mutants: List[MutantProtein], out_dir: Path,
    ) -> List[Path]:
        """Write each mutant sequence to an individual FASTA file."""
        out_dir.mkdir(parents=True, exist_ok=True)
        paths: List[Path] = []
        for mp in mutants:
            fname = f"{mp.gene_symbol}_{mp.mutation_label}.fasta"
            path = out_dir / fname
            path.write_text(mp.to_fasta())
            paths.append(path)
        return paths

    def _build_one(self, var: SomaticVariant) -> Optional[MutantProtein]:
        """Build a single mutant protein from a variant."""
        ref_aa_1 = self._to_single_letter(var.ref_aa)
        alt_aa_1 = self._to_single_letter(var.alt_aa)
        pos = var.protein_position

        if not all([ref_aa_1, alt_aa_1, pos]):
            logger.debug(
                "Skipping %s: incomplete AA change (ref=%s alt=%s pos=%s)",
                var.short_label, var.ref_aa, var.alt_aa, pos,
            )
            return None

        uniprot_id = self._resolve_uniprot(var.gene_symbol)
        wt_seq = self._fetch_sequence(var.gene_symbol, uniprot_id, var.transcript_id)
        if not wt_seq:
            logger.warning("Could not fetch wildtype sequence for %s", var.gene_symbol)
            return None

        if pos > len(wt_seq):
            logger.warning(
                "%s: position %d exceeds sequence length %d",
                var.short_label, pos, len(wt_seq),
            )
            return None

        actual_wt_aa = wt_seq[pos - 1]
        if actual_wt_aa != ref_aa_1:
            logger.warning(
                "%s: expected ref AA '%s' at pos %d but found '%s' in sequence",
                var.short_label, ref_aa_1, pos, actual_wt_aa,
            )

        mut_seq = wt_seq[: pos - 1] + alt_aa_1 + wt_seq[pos:]
        label = f"{ref_aa_1}{pos}{alt_aa_1}"

        return MutantProtein(
            variant=var,
            gene_symbol=var.gene_symbol,
            uniprot_id=uniprot_id,
            wildtype_sequence=wt_seq,
            mutant_sequence=mut_seq,
            mutation_label=label,
        )

    def _resolve_uniprot(self, gene_symbol: str) -> Optional[str]:
        """Map gene symbol to UniProt accession (human, reviewed)."""
        if gene_symbol in self._gene_to_uniprot:
            return self._gene_to_uniprot[gene_symbol]
        try:
            url = (
                f"{_UNIPROT_API}/search?"
                f"query=gene_exact:{gene_symbol}+AND+organism_id:9606+AND+reviewed:true"
                f"&fields=accession&size=1&format=json"
            )
            resp = requests.get(url, timeout=self._timeout)
            resp.raise_for_status()
            data = resp.json()
            results = data.get("results", [])
            if results:
                acc = results[0]["primaryAccession"]
                self._gene_to_uniprot[gene_symbol] = acc
                return acc
        except Exception as exc:
            logger.debug("UniProt lookup failed for %s: %s", gene_symbol, exc)
        return None

    def _fetch_sequence(
        self,
        gene_symbol: str,
        uniprot_id: Optional[str],
        transcript_id: str,
    ) -> Optional[str]:
        """Fetch protein sequence from UniProt (preferred) or Ensembl."""
        if uniprot_id:
            seq = self._fetch_uniprot_sequence(uniprot_id)
            if seq:
                return seq

        if transcript_id:
            seq = self._fetch_ensembl_sequence(transcript_id)
            if seq:
                return seq

        logger.debug("No sequence source available for %s", gene_symbol)
        return None

    def _fetch_uniprot_sequence(self, accession: str) -> Optional[str]:
        """Fetch protein sequence from UniProt REST API."""
        try:
            url = f"{_UNIPROT_API}/{accession}.fasta"
            resp = requests.get(url, timeout=self._timeout)
            resp.raise_for_status()
            lines = resp.text.strip().split("\n")
            return "".join(l for l in lines if not l.startswith(">"))
        except Exception as exc:
            logger.debug("UniProt sequence fetch failed for %s: %s", accession, exc)
            return None

    def _fetch_ensembl_sequence(self, transcript_id: str) -> Optional[str]:
        """Fetch protein sequence from Ensembl REST API."""
        clean_id = transcript_id.split(".")[0]
        try:
            url = f"{_ENSEMBL_API}/sequence/id/{clean_id}?type=protein"
            resp = requests.get(
                url, headers={"Content-Type": "application/json"},
                timeout=self._timeout,
            )
            resp.raise_for_status()
            data = resp.json()
            return data.get("seq")
        except Exception as exc:
            logger.debug("Ensembl sequence fetch failed for %s: %s", transcript_id, exc)
            return None

    @staticmethod
    def _to_single_letter(aa: Optional[str]) -> Optional[str]:
        """Convert amino acid to single-letter code (handles 1-letter or 3-letter)."""
        if not aa:
            return None
        if len(aa) == 1 and aa.isalpha():
            return aa.upper()
        if len(aa) == 3:
            return _AA3TO1.get(aa.capitalize())
        return _AA3TO1.get(aa.capitalize(), aa[0].upper() if aa else None)
