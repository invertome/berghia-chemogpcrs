#!/usr/bin/env python3
"""check_positive_controls.py — sanity-check known chemoreceptor markers.

Bead -edx. Reads ``references/hcr_positive_controls.csv`` and looks each
control gene up in the ranked-candidates CSV. Reports each control's
rank + percentile + composite score. Raises a non-fatal alert if any
control falls below the 50th percentile (pipeline drift detector).

The control set starts with one entry (Gα_olf, the user's example —
present by HCR in fed slugs but absent from this dataset's RNA-seq).
It grows as wet-lab data accumulates.

Output: ``results/ranking/positive_controls_check.tsv`` with one row
per control: gene_name, found, candidate_id, rank, percentile, score, alert
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def find_control_in_ranked(controls_row: pd.Series,
                           ranked_df: pd.DataFrame) -> dict:
    """Look up a single positive-control row against the ranked CSV.

    Match strategies (priority order):
      1. exact match on `id` column to ranked candidate id
      2. case-insensitive match against any of `gene_name`, `aliases`,
         `protein_id`, `ncbi_gene_id`, `refseq_protein` (if columns exist)
      3. substring match (fallback) — only on columns that exist
    """
    name = str(controls_row.get("gene_name", "")).strip()
    aliases_raw = str(controls_row.get("aliases", "") or "")
    aliases = [a.strip() for a in aliases_raw.split(",") if a.strip()]
    needle_set = {n.lower() for n in [name] + aliases if n}
    refseq = str(controls_row.get("refseq_protein", "") or "").strip()
    if refseq:
        needle_set.add(refseq.lower())

    # Strategy 1: exact id match
    if name in ranked_df["id"].astype(str).values:
        idx = ranked_df.index[ranked_df["id"].astype(str) == name][0]
        return {"matched_strategy": "exact_id", "row_index": int(idx)}

    # Strategy 2 + 3: search through name-bearing columns for full-token
    # case-insensitive match; substring match as last resort.
    name_cols = [c for c in ("id", "gene_name", "description",
                             "protein_id", "name") if c in ranked_df.columns]
    if not name_cols:
        return {"matched_strategy": "none", "row_index": None}

    # Build a single concatenated search string per row
    haystack = ranked_df[name_cols].astype(str).agg(" ".join, axis=1).str.lower()

    for needle in needle_set:
        m = haystack == needle
        if m.any():
            return {"matched_strategy": "exact_alias", "row_index": int(m.idxmax())}
    for needle in needle_set:
        m = haystack.str.contains(rf"\b{needle}\b", regex=True, na=False)
        if m.any():
            return {"matched_strategy": "word_boundary", "row_index": int(m.idxmax())}

    return {"matched_strategy": "none", "row_index": None}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True,
                    help="Path to ranked_candidates_sorted.csv from rank_candidates.py")
    ap.add_argument("--controls-csv", required=True,
                    help="Path to references/hcr_positive_controls.csv")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--alert-percentile", type=float, default=50.0,
                    help="Alert if a control falls below this percentile (default 50)")
    args = ap.parse_args()

    ranked = pd.read_csv(args.ranked_csv)
    controls = pd.read_csv(args.controls_csv)

    # Ensure ranked is sorted by rank_score descending; assign rank + percentile
    if "rank_score" in ranked.columns:
        ranked = ranked.sort_values("rank_score", ascending=False).reset_index(drop=True)
    ranked["rank"] = ranked.index + 1
    n = len(ranked)
    ranked["percentile"] = (1.0 - (ranked["rank"] - 1) / max(n, 1)) * 100.0

    rows = []
    n_alerts = 0
    for _, ctrl in controls.iterrows():
        m = find_control_in_ranked(ctrl, ranked)
        if m["row_index"] is None:
            rows.append({
                "gene_name": ctrl.get("gene_name", ""),
                "found": False,
                "candidate_id": "",
                "rank": "",
                "percentile": "",
                "rank_score": "",
                "matched_strategy": m["matched_strategy"],
                "alert": True,
                "notes": "Not found in ranked CSV",
            })
            n_alerts += 1
            continue
        r = ranked.iloc[m["row_index"]]
        alert = float(r["percentile"]) < args.alert_percentile
        if alert:
            n_alerts += 1
        rows.append({
            "gene_name": ctrl.get("gene_name", ""),
            "found": True,
            "candidate_id": r["id"],
            "rank": int(r["rank"]),
            "percentile": float(r["percentile"]),
            "rank_score": float(r.get("rank_score", float("nan"))),
            "matched_strategy": m["matched_strategy"],
            "alert": alert,
            "notes": ctrl.get("notes", ""),
        })

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    if n_alerts:
        print(f"WARN: {n_alerts}/{len(rows)} positive controls below "
              f"the {args.alert_percentile}th percentile or missing.",
              file=sys.stderr)
    else:
        print(f"OK: all {len(rows)} positive controls present and above "
              f"the {args.alert_percentile}th percentile.", file=sys.stderr)
    return 0  # non-fatal: this is a sanity check, not a gate


if __name__ == "__main__":
    sys.exit(main())
