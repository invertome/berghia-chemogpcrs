#!/usr/bin/env python3
# build_structural_channel.py
# Purpose: PRODUCER for the Foldseek structural-evidence channel (Glue task
#   G1, docs/plans/2026-07-01-ml-plm-chemoreceptor-ranking.md). Bridges the
#   raw Unity output (scripts/unity/run_foldseek_candidates.sh: three
#   per-DB Foldseek easy-search tabs) to the TSV
#   scripts/rank_aggregation.py's merge_evidence_channels() consumes.
"""Merge Foldseek hits across DBs, classify, and write the struct channel TSV.

scripts/structural_evidence.py is the SCORER: it parses one Foldseek
easy-search tab file and classifies a single candidate's best hit. This
module is the PRODUCER one level up the pipeline:

    1. merge_best_hits    -- across the PDB/AFDB50/GPCRdb hit files, keep
                              the highest-alntmscore hit per candidate.
    2. build_family_map   -- target id -> family label, from the anchor set
                              (+ optional GPCRdb-specific metadata).
    3. build_structural_channel -- classify every candidate's merged best
                              hit via structural_evidence.classify_hit.
    4. write_channel_tsv  -- emit the TSV rank_aggregation.
                              merge_evidence_channels() left-joins by id.

Same council rule as structural_evidence.py (bead berghia-chemogpcrs-875):
cross-species structural resemblance is used ONLY for RECALL (novelty --
no confident hit to anything characterized) and EXCLUSION (corroborating a
known NON-chemoreceptor GPCR family). A confident hit to a family that is
not a recognized non-chemoreceptor family ("known_other") is never surfaced
as positive evidence for chemoreceptor identity.

This module is a pure parser/merger/writer: no Foldseek binary is invoked
here, and there are no import-time side effects.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from typing import Dict, List, Optional, Set

from structural_evidence import (
    _strip_structure_suffixes, classify_hit, parse_foldseek,
)

# The exact columns scripts/rank_aggregation.py's merge_evidence_channels()
# reads from a struct TSV: the id_col join key ("id") plus its struct
# value_cols (["struct_novelty", "struct_nonchemo_corrob", "struct_state"]).
# has_struct_data is written too for human-readability / parity with
# structural_channel()'s own per-candidate schema, even though
# merge_evidence_channels recomputes that flag itself from row presence.
CHANNEL_TSV_COLUMNS = [
    "id", "struct_state", "struct_novelty", "struct_nonchemo_corrob",
    "has_struct_data",
]


# --------------------------------------------------------------------------- #
# Foldseek QUERY-identifier normalisation
#
# The target side already normalises through structural_evidence.target_keys() /
# _strip_structure_suffixes(). The query side used the raw foldseek field
# VERBATIM as the channel join key, so any decoration foldseek adds to the query
# made the channel row join to nothing -- silently, because an all-miss left
# join just reads as has_struct_data=0 for every candidate.
#
# MEASURED, not assumed (foldseek 10.941cd33, this project's berghia-gpcr env on
# Unity, srun 61999969; flat query dir staged <candidate_id>.<ext> exactly as
# scripts/unity/run_foldseek_candidates.sh stages it):
#
#     BersteEVm000001t1.cif  (1 chain)   -> query "BersteEVm000001t1"
#     BersteEVm000002t1.pdb  (1 chain)   -> query "BersteEVm000002t1"
#     BersteEVm000003t1.cif  (2 chains)  -> queries "BersteEVm000003t1_A"
#                                                   "BersteEVm000003t1_B"
#
# So foldseek STRIPS the extension (structcreatedb.cpp: Util::remove_extension,
# applied twice for .gz/.zst) -- "<id>.cif" and hence "<id>.cif_A" are forms it
# never emits -- and APPENDS "_<chain>" once a file holds more than one chain,
# emitting one row PER CHAIN. The chain suffix is the real break: a multi-chain
# AF3 model (receptor + peptide/ligand/G-protein, or any complex) yields
# <cand_id>_A / <cand_id>_B, neither of which joins to the bare candidate id,
# and one candidate silently becomes two orphan channel rows.
#
# target_keys() alone does not repair that: its _PDB_RE wants a 4-character id
# starting with a digit, so "BersteEVm000003t1_A" normalises to itself. And the
# chain strip cannot be applied blindly, because a candidate id may legitimately
# end in "_<alnum>" ("Berghia_scaffold_12"). It is therefore applied ONLY when
# the stripped form is corroborated against the real candidate id universe --
# resolution by identity, never by pattern-guess.
# --------------------------------------------------------------------------- #

# Trailing "_<chain>" as foldseek appends it (mmCIF/PDB chain ids are short
# alphanumerics). Only ever consulted against a known candidate id set.
_CHAIN_SUFFIX_RE = re.compile(r"^(.+)_([A-Za-z0-9]{1,4})$")


def query_keys(query: Optional[str]) -> List[str]:
    """Ordered candidate join keys for one Foldseek query id.

    Mirrors structural_evidence.target_keys() for the query side, reusing the
    same _strip_structure_suffixes() helper so the two sides cannot drift:
    raw -> path basename -> extension-stripped basename -> chain-stripped.

    Deduplicated and order-stable. An empty/None query yields [].
    """
    raw = (query or "").strip()
    if not raw:
        return []

    keys: List[str] = []

    def add(key: str) -> None:
        if key and key not in keys:
            keys.append(key)

    add(raw)
    basename = raw.rsplit("/", 1)[-1]
    add(basename)
    stem = _strip_structure_suffixes(basename)
    add(stem)

    chain = _CHAIN_SUFFIX_RE.match(stem)
    if chain:
        add(chain.group(1))
    return keys


def canonical_query_id(query: Optional[str],
                        candidate_ids: Optional[Set[str]] = None) -> Optional[str]:
    """Resolve a Foldseek query id into the candidate id namespace.

    With `candidate_ids`, returns the first of query_keys() that is a REAL
    candidate id, or None when none of them is -- an unresolvable query is
    reported, never rewritten into something that merely looks joinable.

    Without `candidate_ids` there is nothing to corroborate against, so only
    the unambiguous extension strip is applied and the chain suffix is left
    alone (stripping it on a guess could truncate a legitimate id).
    """
    keys = query_keys(query)
    if not keys:
        return None
    if candidate_ids is None:
        return _strip_structure_suffixes(keys[0].rsplit("/", 1)[-1]) or None
    for key in keys:
        if key in candidate_ids:
            return key
    return None


def read_candidate_ids(fasta_path: str) -> Set[str]:
    """Candidate id universe from a FASTA: the first whitespace token of each
    header. Missing/empty path yields an empty set (the caller decides)."""
    ids: Set[str] = set()
    if not fasta_path or not os.path.exists(fasta_path):
        return ids
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                token = line[1:].split()
                if token:
                    ids.add(token[0])
    return ids


def merge_best_hits(foldseek_tsv_paths: List[str]) -> Dict[str, dict]:
    """One best hit per query across all given Foldseek hit files.

    Reuses structural_evidence.parse_foldseek per file (so missing files,
    blank lines, and malformed rows are all already handled gracefully
    there) and keeps the hit with the highest alntmscore per query across
    every file -- mirroring parse_foldseek's own within-file tie-break.

    Args:
        foldseek_tsv_paths: paths to Foldseek easy-search tab files
            (typically 3: PDB, AFDB50, GPCRdb). A path that doesn't exist
            simply contributes no hits. An empty/None list returns {}.

    Returns:
        {query_id: {"target": str, "fident": float, "alntmscore": float,
                    "evalue": float}}
    """
    merged: Dict[str, dict] = {}
    for path in foldseek_tsv_paths or []:
        for query, hit in parse_foldseek(path).items():
            current = merged.get(query)
            if current is None or hit["alntmscore"] > current["alntmscore"]:
                merged[query] = hit
    return merged


def build_family_map(anchor_set_tsv: str, gpcrdb_meta: Optional[str] = None) -> Dict[str, str]:
    """Best-effort target id -> family label map for classify_hit().

    Primary source is the anchor set (references/anchors/anchor_set.tsv,
    columns accession,tier,taxid,species,family,class,evidence -- written
    by scripts/build_anchor_set.py): accession -> family. If `gpcrdb_meta`
    is given (a TSV with "target"/"family" columns -- e.g. an
    accession-less GPCRdb-DB-specific mapping), it is merged in on top.

    Both inputs are optional-in-effect: a missing/nonexistent path (or
    gpcrdb_meta=None) simply contributes nothing rather than raising --
    targets absent from the returned map make classify_hit() fall back to
    "known_other", the graceful/conservative default. Rows missing either
    field are skipped.

    Returns:
        {target_id: family_label}
    """
    family_map: Dict[str, str] = {}

    def _merge_from(path: str, id_field_candidates: tuple) -> None:
        if not path or not os.path.exists(path):
            return
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                target = ""
                for field in id_field_candidates:
                    target = (row.get(field) or "").strip()
                    if target:
                        break
                family = (row.get("family") or "").strip()
                if target and family:
                    family_map[target] = family

    _merge_from(anchor_set_tsv, ("accession",))
    _merge_from(gpcrdb_meta, ("target", "accession"))
    return family_map


def build_structural_channel(foldseek_tsv_paths: List[str], family_map: Dict[str, str],
                              tm_threshold: float = 0.5,
                              candidate_ids: Optional[Set[str]] = None) -> Dict[str, dict]:
    """Per-candidate Foldseek structural-evidence channel across all DBs.

    Merges the given Foldseek hit files (merge_best_hits), normalises each
    Foldseek query into the candidate id namespace (canonical_query_id --
    symmetric with the target side), and classifies each candidate's merged
    best hit (structural_evidence.classify_hit), exactly mirroring
    structural_channel()'s per-candidate schema.

    `candidate_ids` is the real candidate id universe (e.g. the class-A
    candidate FASTA's headers). Passing it is strongly preferred: it is what
    lets a multi-chain query's "_<chain>" suffix be resolved rather than
    guessed, collapses a model's chains into ONE candidate row (best
    alntmscore wins), and enables the join assertion below.

    Raises:
        ValueError: if `candidate_ids` is given but empty, or if there are
            hits and NOT ONE of them resolves to a candidate. Zero key
            overlap means the channel would left-join to nothing and the
            structural voter would go dormant while everything exits 0 --
            the exact silent failure this guard exists to make loud.
    """
    if candidate_ids is not None and not candidate_ids:
        raise ValueError(
            "build_structural_channel: candidate_ids is empty -- the candidate "
            "universe could not be read, so every Foldseek query would fail to "
            "resolve and the structural channel would silently join to nothing. "
            "Check the --candidate-fasta path."
        )

    hits = merge_best_hits(foldseek_tsv_paths)

    resolved: Dict[str, dict] = {}
    unresolved: List[str] = []
    for query, best in hits.items():
        canonical = canonical_query_id(query, candidate_ids)
        if canonical is None:
            unresolved.append(query)
            continue
        # Chains of one model collapse here; keep the strongest hit.
        current = resolved.get(canonical)
        if current is None or best["alntmscore"] > current["alntmscore"]:
            resolved[canonical] = best

    if hits and candidate_ids is not None and not resolved:
        raise ValueError(
            "build_structural_channel: ZERO of "
            f"{len(hits)} Foldseek queries resolved to a known candidate id. "
            "The structural channel would left-join to nothing and go silently "
            "dormant.\n"
            f"  saw (up to 5):      {sorted(unresolved)[:5]}\n"
            f"  expected ids like:  {sorted(candidate_ids)[:5]}\n"
            "  Likely cause: the Foldseek query dir was staged under the wrong "
            "id scheme (e.g. AF3 job names instead of candidate ids -- bead "
            "5ubd), or --candidate-fasta points at a different candidate set."
        )

    if unresolved:
        print(
            f"[build_structural_channel] WARNING: {len(unresolved)} of "
            f"{len(hits)} Foldseek queries did not resolve to a candidate id "
            f"and were dropped: {sorted(unresolved)[:5]}"
            f"{' ...' if len(unresolved) > 5 else ''}",
            file=sys.stderr,
        )

    channel: Dict[str, dict] = {}
    for candidate_id, best in resolved.items():
        state = classify_hit(best, family_map, tm_threshold=tm_threshold)
        channel[candidate_id] = {
            "struct_state": state,
            "struct_novelty": 1 if state == "novel" else 0,
            "struct_nonchemo_corrob": 1 if state == "known_non_chemoreceptor" else 0,
            "has_struct_data": True,
        }
    return channel


def write_channel_tsv(channel: Dict[str, dict], path: str) -> None:
    """Write the struct channel TSV that merge_evidence_channels() reads.

    Columns are exactly CHANNEL_TSV_COLUMNS (id + struct_state +
    struct_novelty + struct_nonchemo_corrob + has_struct_data). Rows are
    sorted by candidate id for deterministic output. An empty channel still
    writes a valid header-only TSV (pandas reads it as a zero-row frame
    with the right columns, which merge_evidence_channels handles fine).
    """
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(CHANNEL_TSV_COLUMNS)
        for candidate_id in sorted(channel):
            row = channel[candidate_id]
            writer.writerow([
                candidate_id,
                row["struct_state"],
                row["struct_novelty"],
                row["struct_nonchemo_corrob"],
                row["has_struct_data"],
            ])


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def build_args_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Merge PDB/AFDB50/GPCRdb Foldseek hits, classify vs the "
                     "anchor-set family map, and write the structural-evidence "
                     "channel TSV that rank_aggregation.merge_evidence_channels "
                     "consumes.")
    p.add_argument("--foldseek-tsvs", nargs="+", required=True,
                    help="Foldseek easy-search hit TSVs (typically 3: PDB, "
                         "AFDB50, GPCRdb; from "
                         "scripts/unity/run_foldseek_candidates.sh)")
    p.add_argument("--anchor-set", default="references/anchors/anchor_set.tsv",
                    help="anchor_set.tsv (accession,tier,taxid,species,family,"
                         "class,evidence) for the target->family map "
                         "(default: references/anchors/anchor_set.tsv)")
    p.add_argument("--gpcrdb-meta", default=None,
                    help="Optional target->family TSV for GPCRdb-DB-specific "
                         "targets not covered by --anchor-set")
    p.add_argument("--candidate-fasta", default=None,
                    help="Candidate FASTA whose headers give the candidate id "
                         "universe. Used to resolve Foldseek query ids into the "
                         "ranking table's namespace (a multi-chain model emits "
                         "'<cand_id>_<chain>' per chain) and to assert the join "
                         "is non-empty. Strongly recommended: without it a "
                         "chain-suffixed query cannot be resolved and its row "
                         "joins to nothing.")
    p.add_argument("--out", required=True,
                    help="Output structural-channel TSV path")
    return p


def main(argv=None) -> int:
    args = build_args_parser().parse_args(argv)
    family_map = build_family_map(args.anchor_set, args.gpcrdb_meta)

    candidate_ids = None
    if args.candidate_fasta:
        candidate_ids = read_candidate_ids(args.candidate_fasta)
        if not candidate_ids:
            print(f"[build_structural_channel] ERROR: no candidate ids read "
                  f"from {args.candidate_fasta}", file=sys.stderr)
            return 2
        print(f"[build_structural_channel] candidate universe: "
              f"{len(candidate_ids)} ids from {args.candidate_fasta}",
              file=sys.stderr)
    else:
        print("[build_structural_channel] WARNING: no --candidate-fasta given; "
              "Foldseek query ids cannot be resolved against the candidate "
              "namespace and the join is unverified.", file=sys.stderr)

    channel = build_structural_channel(args.foldseek_tsvs, family_map,
                                        candidate_ids=candidate_ids)
    write_channel_tsv(channel, args.out)
    n_novel = sum(1 for r in channel.values() if r["struct_state"] == "novel")
    n_excl = sum(1 for r in channel.values()
                 if r["struct_state"] == "known_non_chemoreceptor")
    print(f"[build_structural_channel] {len(channel)} candidates -> {args.out} "
          f"({n_novel} novel, {n_excl} known_non_chemoreceptor)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
