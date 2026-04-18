"""Route variant proteins to robust structure generation for Stage 4.

Decision tree (variant path):
  1. Prefer experimental RCSB structure for target gene/UniProt.
  2. Build mutant model by point-mutation on the experimental template.
  3. Run local relaxation (OpenMM restrained minimisation).
  4. If any of the above fails, fall back to ESMFold (local first, then remote).
  5. If ESMFold fails or seq too long, use AlphaFold DB WT as last resort.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, ClassVar, Dict, List, Optional, Tuple

from drugpipe.utils.http import HTTPClient
from drugpipe.variant_analysis.mutant_sequence import MutantProtein

logger = logging.getLogger(__name__)

_RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
_ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
_ESMFOLD_LOCAL_DEFAULT = "http://localhost:8000/predict"
_ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"

_AA3_TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}
_AA1_TO3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    "U": "SEC", "O": "PYL",
}

_MUT_LABEL = re.compile(r"^(?P<ref>[A-Z\*])(?P<pos>\d+)(?P<alt>[A-Z\*])$")


class StructureBridge:
    """Select and execute the best structure route per mutation."""

    _global_structure_cache: ClassVar[Dict[str, Dict[str, Any]]] = {}

    def __init__(self, cfg: Optional[Dict[str, Any]] = None):
        self.cfg = cfg or {}
        va = self.cfg.get("variant_analysis", {})
        self._timeout = int(va.get("api_timeout_s", 30))
        self._prefer_experimental = bool(va.get("prefer_experimental_structure", True))
        self._pdb_cache: Dict[str, List[str]] = {}

        esm_cfg = va.get("esmfold", {}) if isinstance(va.get("esmfold"), dict) else {}
        self._esmfold_skip_threshold = int(esm_cfg.get("skip_above_length", 800))
        self._esmfold_local_url = str(esm_cfg.get("local_url", "")).strip() or ""
        self._esmfold_local_max_len = int(esm_cfg.get("local_max_length", 700))
        self._esmfold_remote_enabled = bool(esm_cfg.get("remote_enabled", False))

        retries = max(1, int(va.get("structure_http_retries", 5)))
        polite = float(va.get("structure_polite_sleep", 0.05))
        self._http = HTTPClient(
            timeout=max(self._timeout, 30),
            retries=retries,
            polite_sleep=polite,
        )
        self._http_esm = HTTPClient(
            timeout=max(self._timeout, 120),
            retries=max(1, int(esm_cfg.get("remote_retries", 2))),
            polite_sleep=polite,
        )

        mm = va.get("mutation_modeling", {}) if isinstance(va.get("mutation_modeling", {}), dict) else {}
        self._mutation_enabled = bool(mm.get("enabled", True))
        self._local_opt_enabled = bool(mm.get("local_optimize", True))
        self._local_opt_window = max(0, int(mm.get("local_optimize_window", 2)))
        self._local_opt_k = float(mm.get("local_optimize_k_kcal_per_a2", 10.0))
        self._local_opt_max_iter = max(10, int(mm.get("local_optimize_max_iterations", 500)))

    def resolve_structures(
        self,
        mutants: List[MutantProtein],
        out_dir: Path,
    ) -> List[Dict[str, Any]]:
        """For each mutant, resolve an actionable PDB path."""
        out_dir.mkdir(parents=True, exist_ok=True)
        results: List[Dict[str, Any]] = []
        cache_hits = 0

        for mp in mutants:
            cache_key = f"{mp.uniprot_id or mp.gene_symbol}_{mp.mutation_label}"
            cached = self._global_structure_cache.get(cache_key)
            if cached and Path(cached["pdb_path"]).exists():
                logger.info(
                    "Global cache hit: %s → %s (%s)",
                    cache_key, cached["pdb_path"], cached.get("method", "?"),
                )
                result = dict(cached)
                result["mutant_protein"] = mp
                results.append(result)
                cache_hits += 1
                continue

            record = self._resolve_one_uncached(mp, out_dir)
            if record is not None:
                self._global_structure_cache[cache_key] = {
                    k: v for k, v in record.items() if k != "mutant_protein"
                }
                results.append(record)

        logger.info(
            "Resolved %d / %d mutant structures (%d from global cache).",
            len(results), len(mutants), cache_hits,
        )
        return results

    def _resolve_one_uncached(self, mp: MutantProtein, out_dir: Path) -> Optional[Dict[str, Any]]:
        """Resolve structure for a single mutant protein (no cache)."""
        gene = mp.gene_symbol
        label = mp.mutation_label

        if self._prefer_experimental:
            pdb_ids = self._find_experimental_pdbs(gene, mp.uniprot_id)
            for pdb_id in pdb_ids:
                logger.info(
                    "%s %s: trying experimental PDB %s (point mutation + local optimisation).",
                    gene, label, pdb_id,
                )
                base_pdb = self._download_pdb(pdb_id, out_dir)
                if not base_pdb:
                    continue
                if self._mutation_enabled:
                    mut_pdb = self._build_mutant_from_experimental(mp, base_pdb, out_dir, pdb_id)
                    if mut_pdb:
                        return {
                            "gene": gene,
                            "mutation": label,
                            "pdb_path": str(mut_pdb),
                            "method": f"experimental_mutated_localopt({pdb_id})",
                            "base_pdb_id": pdb_id,
                            "mutant_protein": mp,
                        }
                else:
                    return {
                        "gene": gene,
                        "mutation": label,
                        "pdb_path": str(base_pdb),
                        "method": f"experimental_no_mutation_model({pdb_id})",
                        "base_pdb_id": pdb_id,
                        "mutant_protein": mp,
                    }
            if self._prefer_experimental and self._find_experimental_pdbs(gene, mp.uniprot_id):
                logger.warning("%s %s: experimental mutation modelling failed, falling back to ESMFold.", gene, label)

        # Fallback 1: ESMFold with mutant sequence (patient-specific sequence route).
        pdb_path = self._predict_esmfold(mp, out_dir)
        if pdb_path:
            return {
                "gene": gene,
                "mutation": label,
                "pdb_path": str(pdb_path),
                "method": "esmfold_mutant",
                "mutant_protein": mp,
            }

        # Fallback 2: AlphaFold WT model (last resort).
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

    # ---------------------------------------------------------------------
    # Experimental PDB search and download
    # ---------------------------------------------------------------------
    def _find_experimental_pdbs(self, gene_symbol: str, uniprot_id: Optional[str]) -> List[str]:
        """Query RCSB for candidate experimental structures."""
        cache_key = uniprot_id or gene_symbol
        if cache_key in self._pdb_cache:
            return self._pdb_cache[cache_key]

        pdb_ids: List[str] = []
        if uniprot_id:
            pdb_ids.extend(self._rcsb_search_by_uniprot(uniprot_id))
        if not pdb_ids:
            pdb_ids.extend(self._rcsb_search_by_gene(gene_symbol))

        seen = set()
        uniq: List[str] = []
        for pid in pdb_ids:
            if pid and pid not in seen:
                seen.add(pid)
                uniq.append(pid)

        self._pdb_cache[cache_key] = uniq
        return uniq

    def _rcsb_search_by_uniprot(self, uniprot_id: str) -> List[str]:
        """Search RCSB by UniProt accession; higher quality first."""
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
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
        """Search RCSB by gene symbol in human proteins."""
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
        """Execute RCSB search query and return candidate PDB IDs."""
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
            logger.info("Downloaded PDB %s -> %s", pdb_id, pdb_path)
            return pdb_path
        except Exception as exc:
            logger.warning("PDB download failed for %s: %s", pdb_id, exc)
            return None

    # ---------------------------------------------------------------------
    # Mutation modelling on experimental template
    # ---------------------------------------------------------------------
    def _build_mutant_from_experimental(
        self,
        mp: MutantProtein,
        base_pdb_path: Path,
        out_dir: Path,
        pdb_id: str,
    ) -> Optional[Path]:
        """Create mutant structure from experimental template and relax locally."""
        try:
            from pdbfixer import PDBFixer
            import openmm
            import openmm.app as app
            import openmm.unit as unit
        except Exception as exc:
            logger.warning("Mutation modelling unavailable (pdbfixer/openmm missing): %s", exc)
            return None

        parsed = self._extract_variant_mutation(mp)
        if not parsed:
            logger.warning("Cannot parse mutation for %s %s; skip template mutation.", mp.gene_symbol, mp.mutation_label)
            return None
        ref_aa, pos, alt_aa = parsed
        alt3 = _AA1_TO3.get(alt_aa)
        if alt3 is None:
            logger.warning("Unsupported target amino acid %s for %s", alt_aa, mp.mutation_label)
            return None

        try:
            fixer = PDBFixer(filename=str(base_pdb_path))
            fixer.removeHeterogens(keepWater=False)
        except Exception as exc:
            logger.warning("Failed loading template PDB %s: %s", base_pdb_path, exc)
            return None

        selected = self._select_chain_residue_for_mutation(fixer, mp, ref_aa, pos)
        if selected is None:
            logger.warning(
                "%s %s: cannot map mutation site to template %s; fallback required.",
                mp.gene_symbol,
                mp.mutation_label,
                pdb_id,
            )
            return None

        chain_id, residue_id, orig_resname, chain_residues, mut_chain_idx = selected
        mut_str = f"{orig_resname}-{residue_id}-{alt3}"

        try:
            fixer.applyMutations([mut_str], chain_id)
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.0)
            positions = fixer.positions
            if self._local_opt_enabled:
                positions = self._local_relax_mutation_region(
                    topology=fixer.topology,
                    positions=positions,
                    chain_id=chain_id,
                    chain_residues=chain_residues,
                    mut_chain_idx=mut_chain_idx,
                    openmm_mod=openmm,
                    app_mod=app,
                    unit_mod=unit,
                )

            out_path = out_dir / f"{base_pdb_path.stem}_{mp.mutation_label}_mut_localopt.pdb"
            with open(out_path, "w") as fh:
                app.PDBFile.writeFile(fixer.topology, positions, fh, keepIds=True)
            logger.info(
                "%s %s: template mutation model built (%s, %s -> %s) -> %s",
                mp.gene_symbol,
                mp.mutation_label,
                pdb_id,
                orig_resname,
                alt3,
                out_path,
            )
            return out_path
        except Exception as exc:
            logger.warning(
                "%s %s: template mutation modelling failed on %s (%s)",
                mp.gene_symbol,
                mp.mutation_label,
                pdb_id,
                exc,
            )
            return None

    def _extract_variant_mutation(self, mp: MutantProtein) -> Optional[Tuple[str, int, str]]:
        """Extract (refAA, position, altAA) from variant annotations."""
        var = mp.variant
        pos = int(var.protein_position) if var.protein_position else None

        ref = str(var.ref_aa or "").strip().upper()[:1]
        alt = str(var.alt_aa or "").strip().upper()[:1]

        if (not ref or not alt or not pos) and mp.mutation_label:
            m = _MUT_LABEL.match(str(mp.mutation_label).strip().upper())
            if m:
                ref = ref or m.group("ref")
                alt = alt or m.group("alt")
                pos = pos or int(m.group("pos"))

        if not ref or not alt or not pos:
            return None
        if ref == "*" or alt == "*":
            # Truncation is not a simple side-chain substitution.
            return None
        return ref, int(pos), alt

    def _select_chain_residue_for_mutation(
        self,
        fixer: Any,
        mp: MutantProtein,
        ref_aa: str,
        protein_pos: int,
    ) -> Optional[Tuple[str, str, str, List[Any], int]]:
        """Map canonical mutation position to best template chain residue."""
        wt_seq = str(mp.wildtype_sequence or "").strip().upper()
        if not wt_seq or protein_pos < 1 or protein_pos > len(wt_seq):
            return None

        best: Optional[Tuple[float, str, str, str, List[Any], int]] = None

        for ci, chain in enumerate(fixer.topology.chains()):
            chain_id = chain.id if str(chain.id).strip() else str(ci)
            chain_seq, chain_residues = self._chain_sequence_and_residues(chain)
            if not chain_seq:
                continue
            aln_score, q_to_t = self._needleman_wunsch_map(wt_seq, chain_seq)
            q_idx = protein_pos - 1
            if q_idx not in q_to_t:
                continue
            t_idx = q_to_t[q_idx]
            if t_idx < 0 or t_idx >= len(chain_residues):
                continue
            residue = chain_residues[t_idx]
            aa1 = _AA3_TO1.get(str(residue.name).upper(), "")
            if aa1 != ref_aa:
                continue
            cand = (aln_score, chain_id, str(residue.id), str(residue.name).upper(), chain_residues, t_idx)
            if best is None or cand[0] > best[0]:
                best = cand

        if best is not None:
            _, chain_id, residue_id, resname, chain_residues, mut_idx = best
            return chain_id, residue_id, resname, chain_residues, mut_idx

        # Fallback: direct residue id match if numbering is canonical and AA matches.
        rid = str(protein_pos)
        for ci, chain in enumerate(fixer.topology.chains()):
            chain_id = chain.id if str(chain.id).strip() else str(ci)
            chain_seq, chain_residues = self._chain_sequence_and_residues(chain)
            if not chain_seq:
                continue
            for idx, residue in enumerate(chain_residues):
                aa1 = _AA3_TO1.get(str(residue.name).upper(), "")
                if str(residue.id) == rid and aa1 == ref_aa:
                    return chain_id, rid, str(residue.name).upper(), chain_residues, idx
        return None

    @staticmethod
    def _chain_sequence_and_residues(chain: Any) -> Tuple[str, List[Any]]:
        residues: List[Any] = []
        chars: List[str] = []
        for residue in chain.residues():
            aa = _AA3_TO1.get(str(residue.name).upper())
            if aa:
                residues.append(residue)
                chars.append(aa)
        return "".join(chars), residues

    @staticmethod
    def _needleman_wunsch_map(query: str, target: str) -> Tuple[float, Dict[int, int]]:
        """Global alignment map: query index -> target index for matched columns."""
        m, n = len(query), len(target)
        if m == 0 or n == 0:
            return float("-inf"), {}

        match = 2
        mismatch = -1
        gap = -2

        score = [[0] * (n + 1) for _ in range(m + 1)]
        trace = [[""] * (n + 1) for _ in range(m + 1)]

        for i in range(1, m + 1):
            score[i][0] = i * gap
            trace[i][0] = "U"
        for j in range(1, n + 1):
            score[0][j] = j * gap
            trace[0][j] = "L"

        for i in range(1, m + 1):
            qi = query[i - 1]
            for j in range(1, n + 1):
                tj = target[j - 1]
                s_diag = score[i - 1][j - 1] + (match if qi == tj else mismatch)
                s_up = score[i - 1][j] + gap
                s_left = score[i][j - 1] + gap
                if s_diag >= s_up and s_diag >= s_left:
                    score[i][j] = s_diag
                    trace[i][j] = "D"
                elif s_up >= s_left:
                    score[i][j] = s_up
                    trace[i][j] = "U"
                else:
                    score[i][j] = s_left
                    trace[i][j] = "L"

        q_to_t: Dict[int, int] = {}
        i, j = m, n
        while i > 0 or j > 0:
            move = trace[i][j]
            if move == "D":
                i -= 1
                j -= 1
                q_to_t[i] = j
            elif move == "U":
                i -= 1
            else:
                j -= 1

        return float(score[m][n]), q_to_t

    def _local_relax_mutation_region(
        self,
        topology: Any,
        positions: Any,
        chain_id: str,
        chain_residues: List[Any],
        mut_chain_idx: int,
        openmm_mod: Any,
        app_mod: Any,
        unit_mod: Any,
    ) -> Any:
        """Restrained minimisation keeping movement local to mutation window."""
        try:
            forcefield = app_mod.ForceField("amber14-all.xml", "implicit/obc2.xml")
            system = forcefield.createSystem(
                topology,
                nonbondedMethod=app_mod.NoCutoff,
                constraints=app_mod.HBonds,
            )

            lo = max(0, mut_chain_idx - self._local_opt_window)
            hi = min(len(chain_residues) - 1, mut_chain_idx + self._local_opt_window)
            mobile_ids = {
                str(chain_residues[i].id)
                for i in range(lo, hi + 1)
            }

            k_kj_nm2 = self._local_opt_k * 418.4
            restraint = openmm_mod.CustomExternalForce(
                "0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
            )
            restraint.addGlobalParameter("k", float(k_kj_nm2))
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")

            added = 0
            for atom in topology.atoms():
                elem = atom.element
                if elem is None or elem.symbol == "H":
                    continue
                rid = str(atom.residue.id)
                cid = str(atom.residue.chain.id).strip()
                if cid == chain_id and rid in mobile_ids:
                    continue
                p = positions[atom.index].value_in_unit(unit_mod.nanometer)
                restraint.addParticle(atom.index, [float(p.x), float(p.y), float(p.z)])
                added += 1

            if added > 0:
                system.addForce(restraint)

            integrator = openmm_mod.LangevinMiddleIntegrator(
                300 * unit_mod.kelvin,
                1.0 / unit_mod.picosecond,
                0.002 * unit_mod.picoseconds,
            )
            simulation = app_mod.Simulation(topology, system, integrator)
            simulation.context.setPositions(positions)
            simulation.minimizeEnergy(maxIterations=self._local_opt_max_iter)
            st = simulation.context.getState(getPositions=True)
            logger.info(
                "Local mutation relaxation done: chain=%s residue_idx=%d window=%d restrained_atoms=%d",
                chain_id,
                mut_chain_idx,
                self._local_opt_window,
                added,
            )
            return st.getPositions()
        except Exception as exc:
            logger.warning("Local mutation relaxation failed; using pre-relax positions: %s", exc)
            return positions

    # ---------------------------------------------------------------------
    # Sequence-based structure fallback
    # ---------------------------------------------------------------------
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
            logger.info("Downloaded AlphaFold DB model %s -> %s", uniprot_id, pdb_path)
            return pdb_path
        except Exception as exc:
            logger.warning("AlphaFold DB download failed for %s: %s", uniprot_id, exc)
            return None

    def _predict_esmfold(self, mp: MutantProtein, out_dir: Path) -> Optional[Path]:
        """Predict structure via ESMFold — local-first, skip-long, remote-optional."""
        pdb_path = out_dir / f"{mp.gene_symbol}_{mp.mutation_label}_esmfold.pdb"
        if pdb_path.exists():
            return pdb_path

        seq = mp.mutant_sequence
        seq_len = len(seq)

        if seq_len > self._esmfold_skip_threshold:
            logger.info(
                "%s %s: seq length %d > skip threshold %d — "
                "skipping ESMFold entirely, proceeding to next fallback.",
                mp.gene_symbol, mp.mutation_label, seq_len,
                self._esmfold_skip_threshold,
            )
            return None

        if self._esmfold_local_url and seq_len <= self._esmfold_local_max_len:
            result = self._call_esmfold_local(mp, seq, pdb_path)
            if result:
                return result
            logger.info(
                "%s %s: local ESMFold failed, trying next option.",
                mp.gene_symbol, mp.mutation_label,
            )

        if not self._esmfold_remote_enabled:
            logger.info(
                "%s %s: remote ESMFold disabled (seq=%d) — next fallback.",
                mp.gene_symbol, mp.mutation_label, seq_len,
            )
            return None

        return self._call_esmfold_remote(mp, seq, pdb_path)

    def _call_esmfold_local(
        self, mp: MutantProtein, seq: str, pdb_path: Path,
    ) -> Optional[Path]:
        """Call locally deployed ESMFold (e.g. LiteFold)."""
        import requests as _req

        try:
            logger.info(
                "%s %s: calling local ESMFold (%d aa) ...",
                mp.gene_symbol, mp.mutation_label, len(seq),
            )
            resp = _req.post(
                self._esmfold_local_url,
                data=seq,
                headers={"Content-Type": "text/plain"},
                timeout=300,
            )
            resp.raise_for_status()
            pdb_text = resp.text
            if "ATOM" not in pdb_text:
                return None
            pdb_path.write_text(pdb_text)
            logger.info("Local ESMFold prediction saved -> %s", pdb_path)
            return pdb_path
        except Exception as exc:
            logger.warning("Local ESMFold failed for %s: %s", mp.mutation_label, exc)
            return None

    def _call_esmfold_remote(
        self, mp: MutantProtein, seq: str, pdb_path: Path,
    ) -> Optional[Path]:
        """Call remote ESMFold API (esmatlas.com) — only if enabled."""
        try:
            pdb_text = self._http_esm.post_text(
                _ESMFOLD_API,
                seq,
                headers={"Content-Type": "text/plain"},
            )
            if "ATOM" not in pdb_text:
                logger.warning("Remote ESMFold returned empty structure for %s", mp.mutation_label)
                return None
            pdb_path.write_text(pdb_text)
            logger.info("Remote ESMFold prediction saved -> %s", pdb_path)
            return pdb_path
        except Exception as exc:
            logger.warning(
                "Remote ESMFold prediction failed for %s %s: %s",
                mp.gene_symbol, mp.mutation_label, exc,
            )
            return None
