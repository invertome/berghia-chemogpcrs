#!/usr/bin/env python3
"""validate_against_aplysia.py — Retrospective validation harness.

Bead -bdu. Given a ranked-candidates CSV produced by running the pipeline
on the *Aplysia californica* genome (NCBI accession GCF_000002075.2),
compute recall@N and precision@N against the labeled chemoreceptor set
in references/cummins2009_aplysia_chemoreceptors.csv.

Output: a TSV recall-curve plus a small JSON summary. NOT a gate — this
is a method-validation diagnostic for the manuscript / paper review.

Usage:
    python validate_against_aplysia.py \\
        --ranked-csv results/aplysia_run/ranking/ranked_candidates_sorted.csv \\
        --labeled-csv references/cummins2009_aplysia_chemoreceptors.csv \\
        --out results/validation/aplysia_recall.tsv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _validation_lib import average_rank, recall_curve  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True)
    ap.add_argument("--labeled-csv", required=True)
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--ns", default="10,25,50,100,250,500",
                    help="Comma-separated N values for the recall curve")
    ap.add_argument("--id-column", default="aplysia_protein_id",
                    help="Which column of the labeled CSV holds IDs that "
                         "match the ranked CSV's 'id' column")
    args = ap.parse_args()

    ranked = pd.read_csv(args.ranked_csv)
    if "id" not in ranked.columns:
        print("ERROR: ranked CSV missing 'id' column", file=sys.stderr)
        return 1
    if "rank_score" in ranked.columns:
        ranked = ranked.sort_values("rank_score", ascending=False)
    ranked_ids = ranked["id"].astype(str).tolist()

    labeled = pd.read_csv(args.labeled_csv)
    if args.id_column not in labeled.columns:
        print(f"ERROR: labeled CSV missing column {args.id_column!r}", file=sys.stderr)
        return 1
    truth_ids = (
        labeled[args.id_column]
        .dropna()
        .astype(str)
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .tolist()
    )
    if not truth_ids:
        print(f"WARN: no usable IDs in column {args.id_column!r} of "
              f"{args.labeled_csv} — nothing to validate against. "
              f"Edit the CSV to add concrete IDs first.",
              file=sys.stderr)

    ns = [int(x) for x in args.ns.split(",") if x.strip()]
    rc_df = recall_curve(ranked_ids, truth_ids, ns)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    rc_df.to_csv(args.out, sep="\t", index=False)

    summary = {
        "n_ranked": len(ranked_ids),
        "n_labeled": len(truth_ids),
        "average_rank_of_truth": average_rank(ranked_ids, truth_ids),
        "recall_curve": rc_df.to_dict(orient="records"),
    }
    summary_path = Path(args.out).with_suffix(".summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote {args.out} and {summary_path}", file=sys.stderr)
    if truth_ids:
        for r in rc_df.to_dict(orient="records"):
            print(f"  N={r['N']}: recall={r['recall']:.3f}, "
                  f"precision={r['precision']:.3f} ({r['n_truth_found_in_topN']}/{r['n_truth_total']})",
                  file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
