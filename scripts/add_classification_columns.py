#!/usr/bin/env python3
"""add_classification_columns.py — Augment ranked CSV with non-chemoreceptor
classification columns from the 06c consensus TSV.

Phase 5 Task 5.1 of the non-chemoreceptor classification feature.

Mirrors the existing add_hcr_columns.py / add_og_coverage_columns.py
pattern — reads the ranked candidates CSV (rank_candidates.py output),
joins on candidate_id with the 06c consensus TSV, and appends:

    classification             chemoreceptor-candidate (default) |
                               likely-non-chemoreceptor (medium) |
                               non-chemoreceptor (high)
    classification_confidence  NA / medium / high
    classification_family      coarse family if classified, else ''
    classification_subfamily   medium-granularity (aminergic/peptide only)
                               or '' if no subfamily consensus
    classification_evidence    semicolon-separated source:label string,
                               e.g. 'hmm:aminergic;og:aminergic;placement:aminergic'

Candidates not in the consensus TSV (e.g. 06c hasn't run, or the
candidate had no orthogroup membership and thus no entry) default to
'chemoreceptor-candidate' / 'NA' / '' / '' / ''.

Also applies classification-based rank suppression (bead f2e): 06c-classified
non-chemoreceptors are down-weighted in rank_score and the CSV re-sorted
(configurable via NONCHEMO_RANK_FACTOR / LIKELY_NONCHEMO_RANK_FACTOR; suppress,
don't drop — a mislabelled divergent candidate stays recoverable).

Usage:
    python3 add_classification_columns.py \\
        --ranked-csv results/ranking/ranked_candidates_sorted.csv \\
        --consensus-tsv results/classification/candidate_classifications.tsv \\
        --out results/ranking/ranked_candidates_sorted.csv
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

NEW_COLUMNS = [
    "classification",
    "classification_confidence",
    "classification_family",
    "classification_subfamily",
    "classification_evidence",
]

# cl2: a source (hmm/og) is "dark" when it made no real family call.
_DARK_SOURCE_LABELS = {"", "na", "unclassified-hmm", "unclassified-og",
                       "unclassified-placement"}


def _source_is_dark(evidence: str, source: str) -> bool:
    """True if `source` (e.g. 'hmm'/'og') made no real family call in the
    evidence string 'hmm:X;og:Y;placement:Z'. A source absent from the string
    counts as dark (no signal)."""
    for part in str(evidence).split(";"):
        src, sep, label = part.partition(":")
        if sep and src.strip() == source:
            return label.strip().lower() in _DARK_SOURCE_LABELS
    return True


def needs_manual_review(classification: str, evidence: str) -> bool:
    """cl2: flag a default 'chemoreceptor-candidate' whose HMM and OG sources
    BOTH went dark — the divergent/unclassifiable case where the two
    annotation-based classifiers give no signal, so a human should eyeball it
    before a divergent real chemoreceptor is silently deprioritized. Placement
    alone does not rescue (least-reliable source)."""
    return (str(classification) == "chemoreceptor-candidate"
            and _source_is_dark(evidence, "hmm")
            and _source_is_dark(evidence, "og"))


def add_classification_columns(ranked_csv_path: str,
                                consensus_tsv_path: str,
                                out_path: str,
                                id_column: str = "id",
                                nonchemo_rank_factor: float = 0.1,
                                likely_nonchemo_rank_factor: float = 0.5) -> int:
    """Read ranked CSV + consensus TSV, join, write annotated CSV. Returns
    number of rows annotated. Also applies classification-based rank suppression
    (bead f2e) — see the suppression block below."""
    df = pd.read_csv(ranked_csv_path, keep_default_na=False, dtype=str)

    if (consensus_tsv_path and os.path.exists(consensus_tsv_path)
            and os.path.getsize(consensus_tsv_path) > 0):
        try:
            consensus = pd.read_csv(consensus_tsv_path, sep="\t",
                                     keep_default_na=False, dtype=str)
        except pd.errors.EmptyDataError:
            consensus = pd.DataFrame(columns=["candidate_id"] + NEW_COLUMNS)
    else:
        print(f"WARN: consensus TSV not found or empty: {consensus_tsv_path}; "
              f"defaulting all candidates to chemoreceptor-candidate",
              file=sys.stderr)
        consensus = pd.DataFrame(columns=["candidate_id"] + NEW_COLUMNS)

    # Build lookup dict: candidate_id -> dict of new columns
    lookup: dict[str, dict[str, str]] = {}
    if not consensus.empty and "candidate_id" in consensus.columns:
        for _, row in consensus.iterrows():
            cid = str(row.get("candidate_id", ""))
            lookup[cid] = {col: str(row.get(col, "")) for col in NEW_COLUMNS}

    # Defaults for candidates not in consensus
    defaults = {
        "classification": "chemoreceptor-candidate",
        "classification_confidence": "NA",
        "classification_family": "",
        "classification_subfamily": "",
        "classification_evidence": "",
    }

    for col in NEW_COLUMNS:
        df[col] = df[id_column].astype(str).map(
            lambda cid: lookup.get(cid, defaults).get(col, defaults[col]))

    # cl2: flag candidates the automated classifiers can't decide (HMM + OG both
    # dark) for manual review, so a divergent real chemoreceptor isn't lost.
    df["needs_manual_review"] = [
        "yes" if needs_manual_review(c, e) else ""
        for c, e in zip(df["classification"], df["classification_evidence"])
    ]
    n_flagged = int((df["needs_manual_review"] == "yes").sum())
    if n_flagged:
        print(f"cl2: {n_flagged} candidate(s) flagged needs_manual_review "
              f"(chemoreceptor-candidate with HMM+OG both dark)", file=sys.stderr)

    # f2e: classification-based rank suppression. The 06c classification is a
    # post-hoc annotation, not a hard filter, so without this a non-chemoreceptor
    # with strong non-phylo signals could surface in the shortlist (more likely
    # now that o98 dropped the crude phylo-absence penalty). Multiply rank_score
    # by a per-classification factor (default: non-chemoreceptor x0.1,
    # likely-non-chemoreceptor x0.5, chemoreceptor-candidate x1.0), preserve the
    # unsuppressed score, and re-sort. This SUPPRESSES rather than drops, so a
    # mislabelled divergent candidate stays recoverable from the full list; set
    # both factors to 1.0 to disable.
    #
    # The re-sort keys on `final_rank`, the PRODUCTION order, NOT on
    # `rank_score`. Under RANK_METHOD=rankagg (the production default)
    # `rank_score` stays the DEMOTED weighted composite while row order carries
    # the rank-aggregation result, so sorting by rank_score silently restored
    # the weighted ordering -- even with nothing reclassified and every factor
    # 1.0. That matters downstream because 08_structural_analysis.sh picks the
    # AlphaFold set by ROW POSITION (`head -n N | cut -d, -f1`).
    #
    # Suppression is therefore applied as a demotion of the production order:
    # candidates sort by (suppression factor descending, final_rank ascending),
    # so a classified non-chemoreceptor sinks below the unclassified ones while
    # everything keeps its production order WITHIN its class, under either
    # RANK_METHOD. Both factors at 1.0 is then exactly order-preserving.
    if "rank_score" in df.columns:
        factor = df["classification"].map({
            "non-chemoreceptor": float(nonchemo_rank_factor),
            "likely-non-chemoreceptor": float(likely_nonchemo_rank_factor),
        }).fillna(1.0)
        df["rank_score_prefilter"] = df["rank_score"]
        df["classification_rank_factor"] = factor
        df["rank_score"] = pd.to_numeric(df["rank_score"], errors="coerce") * factor

        if "final_rank" in df.columns:
            df["final_rank_prefilter"] = pd.to_numeric(
                df["final_rank"], errors="coerce")
            df = (df.assign(_factor=factor)
                    .sort_values(["_factor", "final_rank_prefilter"],
                                 ascending=[False, True], kind="mergesort")
                    .drop(columns="_factor")
                    .reset_index(drop=True))
            df["final_rank"] = range(1, len(df) + 1)
        else:
            # Legacy ranked CSV predating final_rank: the weighted composite is
            # the only order available.
            df = df.sort_values("rank_score", ascending=False,
                                kind="mergesort").reset_index(drop=True)
        if "rank" in df.columns:
            df["rank"] = range(1, len(df) + 1)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return len(df)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True)
    ap.add_argument("--consensus-tsv", required=True,
                    help="06c output: candidate_classifications.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--id-column", default="id",
                    help="Column in ranked CSV that holds the candidate ID "
                         "(default 'id')")
    args = ap.parse_args()
    n = add_classification_columns(
        args.ranked_csv, args.consensus_tsv, args.out, args.id_column,
        nonchemo_rank_factor=float(os.getenv("NONCHEMO_RANK_FACTOR", 0.1)),
        likely_nonchemo_rank_factor=float(os.getenv("LIKELY_NONCHEMO_RANK_FACTOR", 0.5)))
    print(f"Annotated {n} rows -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
