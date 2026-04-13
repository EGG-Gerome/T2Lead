"""Route mutant protein variants to the best structure-prediction strategy.

Decision tree:
  1. If a wildtype crystal structure exists in RCSB PDB → use it as base
     and apply point mutation (placeholder for FoldX integration).
  2. Otherwise → use ESMFold (or existing T2Lead protein_prep path)
     with the full mutant sequence.

This module produces PDB files that feed into T2Lead Stage 4 (docking + MD).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from drugpipe.utils.http import HTTPClient
from drugpipe.variant_analysis.mutant_sequence import MutantProtein

logger = logging.getLogger(__name__)

_RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core"
_ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
_ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"


class StructureBridge:
    """Select and execute the best structure-prediction route per mutation.

    Returns a list of dicts with keys: gene, mutation, pdb_path, method.
    """

    def __init__(self, cfg: Optional[Dict[str, Any]] = None):
        va = (cfg or {}).get("variant_analysis", {})
        self._timeout = int(va.get("api_timeout_s", 30))
        self._prefer_experimental = bool(va.get("prefer_experimental_structure", True))
        self._pdb_cache: Dict[str, List[str]] = {}
        retries = max(1, int(va.get("structure_http_retries", 5)))
        polite = float(va.get("structure_polite_sleep", 0.05))
        self._http = HTTPClient(
            timeout=max(self._timeout, 30),
            retries=retries,
            polite_sleep=polite,
        )
        self._http_esm = HTTPClient(
            timeout=max(self._timeout, 120),
            retries=retries,
            polite_sleep=polite,
        )

    def resolve_structures(
        self,
        mutants: List[MutantProtein],
        out_dir: Path,
    ) -> List[Dict[str, Any]]:
        """For each mutant, find or predict a 3D structure.

        Returns a list of structure records ready for Stage 4.
        """
        out_dir.mkdir(parents=True, exist_ok=True)
        results: List[Dict[str, Any]] = []

        for mp in mutants:
            record = self._resolve_one(mp, out_dir)
            if record is not None:
                results.append(record)

        logger.info(
            "Resolved %d / %d mutant structures.", len(results), len(mutants),
        )
        return results

    def _resolve_one(
        self, mp: MutantProtein, out_dir: Path,
    ) -> Optional[Dict[str, Any]]:
        """Resolve structure for a single mutant protein."""
        gene = mp.gene_symbol
        label = mp.mutation_label

        if self._prefer_experimental:
            pdb_ids = self._find_experimental_pdbs(gene, mp.uniprot_id)
            for pdb_id in pdb_ids:
                logger.info(
                    "%s %s: trying experimental PDB %s (will apply point mutation).",
                    gene, label, pdb_id,
                )
                pdb_path = self._download_pdb(pdb_id, out_dir)
                if pdb_path:
                    return {
                        "gene": gene,
                        "mutation": label,
                        "pdb_path": str(pdb_path),
                        "method": f"experimental+mutation({pdb_id})",
                        "base_pdb_id": pdb_id,
                        "mutant_protein": mp,
                    }
            if pdb_ids:
                logger.warning(
                    "%s %s: all candidate experimental PDB downloads failed (%s).",
                    gene, label, ", ".join(pdb_ids[:5]),
                )

        # ESMFold public API rejects many long proteins (413); use AlphaFold DB
        # wildtype model as a robust fallback for long sequences.
        if len(mp.mutant_sequence) > 400 and mp.uniprot_id:
            af_path = self._download_alphafold(mp.uniprot_id, out_dir)
            if af_path:
                return {
                    "gene": gene,
                    "mutation": label,
                    "pdb_path": str(af_path),
                    "method": f"alphafold_db_wildtype({mp.uniprot_id})",
                    "mutant_protein": mp,
                }

        pdb_path = self._predict_esmfold(mp, out_dir)
        if pdb_path:
            return {
                "gene": gene,
                "mutation": label,
                "pdb_path": str(pdb_path),
                "method": "esmfold_mutant",
                "mutant_protein": mp,
            }

        if mp.uniprot_id:
            af_path = self._download_alphafold(mp.uniprot_id, out_dir)
            if af_path:
                return {
                    "gene": gene,
                    "mutation": label,
                    "pdb_path": str(af_path),
                    "method": f"alphafold_db_wildtype({mp.uniprot_id})",
                    "mutant_protein": mp,
                }

        logger.warning("No structure resolved for %s %s", gene, label)
        return None

    def _find_experimental_pdbs(
        self, gene_symbol: str, uniprot_id: Optional[str],
    ) -> List[str]:
        """Query RCSB for candidate experimental structures."""
        cache_key = uniprot_id or gene_symbol
        if cache_key in self._pdb_cache:
            return self._pdb_cache[cache_key]

        pdb_ids: List[str] = []

        if uniprot_id:
            pdb_ids.extend(self._rcsb_search_by_uniprot(uniprot_id))
        if not pdb_ids:
            pdb_ids.extend(self._rcsb_search_by_gene(gene_symbol))

        # De-duplicate while preserving order.
        seen = set()
        uniq: List[str] = []
        for pid in pdb_ids:
            if pid and pid not in seen:
                seen.add(pid)
                uniq.append(pid)

        self._pdb_cache[cache_key] = uniq
        return uniq

    def _rcsb_search_by_uniprot(self, uniprot_id: str) -> List[str]:
        """Search RCSB by UniProt accession, prefer highest resolution X-ray."""
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers."
                                 "reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id,
                },
            },
            "return_type": "entry",
            "request_options": {
                "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
                "paginate": {"start": 0, "rows": 10},
            },
        }
        return self._run_rcsb_search(query)

    def _rcsb_search_by_gene(self, gene_symbol: str) -> List[str]:
        """Search RCSB by gene name, human organism, best resolution."""
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                            "operator": "exact_match",
                            "value": gene_symbol,
                        },
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.ncbi_taxonomy_id",
                            "operator": "exact_match",
                            "value": 9606,
                        },
                    },
                ],
            },
            "return_type": "entry",
            "request_options": {
                "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
                "paginate": {"start": 0, "rows": 10},
            },
        }
        return self._run_rcsb_search(query)

    def _run_rcsb_search(self, query: dict) -> List[str]:
        """Execute an RCSB search query and return candidate PDB IDs."""
        try:
            data = self._http.post_json(_RCSB_SEARCH_API, query)
            if not data:
                return []
            results = data.get("result_set", [])
            return [r.get("identifier") for r in results if r.get("identifier")]
        except Exception as exc:
            logger.debug("RCSB search failed: %s", exc)
        return []

    def _download_pdb(self, pdb_id: str, out_dir: Path) -> Optional[Path]:
        """Download a PDB file from RCSB."""
        pdb_path = out_dir / f"{pdb_id.lower()}.pdb"
        if pdb_path.exists():
            return pdb_path
        try:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            text = self._http.get_text(url)
            pdb_path.write_text(text)
            logger.info("Downloaded PDB %s → %s", pdb_id, pdb_path)
            return pdb_path
        except Exception as exc:
            logger.warning("PDB download failed for %s: %s", pdb_id, exc)
            return None

    def _download_alphafold(self, uniprot_id: str, out_dir: Path) -> Optional[Path]:
        """Download AlphaFold DB wildtype model by UniProt ID."""
        pdb_path = out_dir / f"af_{uniprot_id.lower()}.pdb"
        if pdb_path.exists():
            return pdb_path
        try:
            url = _ALPHAFOLD_PDB_URL.format(uniprot=uniprot_id)
            text = self._http.get_text(url)
            if "ATOM" not in text:
                logger.warning("AlphaFold DB returned invalid PDB for %s", uniprot_id)
                return None
            pdb_path.write_text(text)
            logger.info("Downloaded AlphaFold DB model %s → %s", uniprot_id, pdb_path)
            return pdb_path
        except Exception as exc:
            logger.warning("AlphaFold DB download failed for %s: %s", uniprot_id, exc)
            return None

    def _predict_esmfold(self, mp: MutantProtein, out_dir: Path) -> Optional[Path]:
        """Predict structure via ESMFold API (or local ESM if available)."""
        pdb_path = out_dir / f"{mp.gene_symbol}_{mp.mutation_label}_esmfold.pdb"
        if pdb_path.exists():
            return pdb_path

        seq = mp.mutant_sequence
        if len(seq) > 400:
            logger.info(
                "%s %s: sequence length %d > 400 — ESMFold API may be slow or fail; "
                "consider local ESM or ColabFold.",
                mp.gene_symbol, mp.mutation_label, len(seq),
            )

        try:
            pdb_text = self._http_esm.post_text(
                _ESMFOLD_API,
                seq,
                headers={"Content-Type": "text/plain"},
            )
            if "ATOM" not in pdb_text:
                logger.warning("ESMFold returned empty structure for %s", mp.mutation_label)
                return None
            pdb_path.write_text(pdb_text)
            logger.info("ESMFold prediction saved → %s", pdb_path)
            return pdb_path
        except Exception as exc:
            logger.warning(
                "ESMFold prediction failed for %s %s: %s",
                mp.gene_symbol, mp.mutation_label, exc,
            )
            return None
