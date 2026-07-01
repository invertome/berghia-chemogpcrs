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
import sys
from typing import Dict, List, Optional

from structural_evidence import classify_hit, parse_foldseek

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
                              tm_threshold: float = 0.5) -> Dict[str, dict]:
    """Per-candidate Foldseek structural-evidence channel across all DBs.

    Merges the given Foldseek hit files (merge_best_hits) and classifies
    each candidate's merged best hit (structural_evidence.classify_hit),
    exactly mirroring structural_channel()'s per-candidate schema.

    Returns:
        {candidate_id: {"struct_state": "novel"|"known_non_chemoreceptor"|
                        "known_other",
                        "struct_novelty": 1 if novel else 0,
                        "struct_nonchemo_corrob": 1 if known_non_chemoreceptor
                                                  else 0,
                        "has_struct_data": True}}
    """
    hits = merge_best_hits(foldseek_tsv_paths)
    channel: Dict[str, dict] = {}
    for query, best in hits.items():
        state = classify_hit(best, family_map, tm_threshold=tm_threshold)
        channel[query] = {
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
    p.add_argument("--out", required=True,
                    help="Output structural-channel TSV path")
    return p


def main(argv=None) -> int:
    args = build_args_parser().parse_args(argv)
    family_map = build_family_map(args.anchor_set, args.gpcrdb_meta)
    channel = build_structural_channel(args.foldseek_tsvs, family_map)
    write_channel_tsv(channel, args.out)
    n_novel = sum(1 for r in channel.values() if r["struct_state"] == "novel")
    n_excl = sum(1 for r in channel.values()
                 if r["struct_state"] == "known_non_chemoreceptor")
    print(f"[build_structural_channel] {len(channel)} candidates -> {args.out} "
          f"({n_novel} novel, {n_excl} known_non_chemoreceptor)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
