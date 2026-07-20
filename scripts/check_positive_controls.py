#!/usr/bin/env python3
"""check_positive_controls.py — sanity-check known chemoreceptor markers.

Bead -edx. Reads ``references/hcr_positive_controls.csv`` and looks each
control gene up in the ranked-candidates CSV. Reports each control's
rank + percentile + composite score.

Each control lands in exactly one of three ``status`` values:

``found_healthy``
    Present in the ranking at or above ``--alert-percentile``. ``alert=False``.
``found_below_percentile``
    Present but ranked below the threshold. This is the genuine pipeline-drift
    signal — a known marker sinking in the ranking means the weights moved.
    ``alert=True``.
``not_found``
    Absent from the ranked CSV. INFORMATIONAL, ``alert=False``.

Absence is deliberately not an alert (bead 444). The only control shipped
today is Gα_olf, a G-protein alpha subunit rather than a class-A GPCR, so it
can never appear in a class-A chemoreceptor ranked CSV; flagging that as an
alert fired red on every run and trained the reader to ignore the check.
Absent controls are still reported in full — their own TSV row, their own
``not_found`` status, and their own count in the stderr summary — because the
goal is to stop crying wolf, not to hide information.

The control set grows as wet-lab data accumulates.

Output: ``results/ranking/positive_controls_check.tsv`` with one row per
control: gene_name, found, candidate_id, rank, percentile, rank_score,
matched_strategy, status, alert, notes
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Per-control outcome. Only STATUS_BELOW sets `alert` — see module docstring.
STATUS_HEALTHY = "found_healthy"
STATUS_BELOW = "found_below_percentile"
STATUS_NOT_FOUND = "not_found"

# Pinned so the TSV schema is stable for stage 09, which reads gene_name /
# found / rank / percentile / alert by name. `status` was appended (bead 444);
# nothing was renamed or removed.
OUTPUT_COLUMNS = [
    "gene_name", "found", "candidate_id", "rank", "percentile", "rank_score",
    "matched_strategy", "status", "alert", "notes",
]


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
    n_below = 0      # found but below the alert percentile -> the drift signal
    n_absent = 0     # not in the ranked CSV at all -> informational
    n_healthy = 0    # found and at/above the threshold
    for _, ctrl in controls.iterrows():
        m = find_control_in_ranked(ctrl, ranked)
        if m["row_index"] is None:
            # Absence is NOT an alert: a control that isn't a class-A GPCR
            # (e.g. Ga_olf) can never appear in a class-A ranked CSV.
            rows.append({
                "gene_name": ctrl.get("gene_name", ""),
                "found": False,
                "candidate_id": "",
                "rank": "",
                "percentile": "",
                "rank_score": "",
                "matched_strategy": m["matched_strategy"],
                "status": STATUS_NOT_FOUND,
                "alert": False,
                "notes": "Not found in ranked CSV (informational, not an alert)",
            })
            n_absent += 1
            continue
        r = ranked.iloc[m["row_index"]]
        below = float(r["percentile"]) < args.alert_percentile
        if below:
            n_below += 1
        else:
            n_healthy += 1
        rows.append({
            "gene_name": ctrl.get("gene_name", ""),
            "found": True,
            "candidate_id": r["id"],
            "rank": int(r["rank"]),
            "percentile": float(r["percentile"]),
            "rank_score": float(r.get("rank_score", float("nan"))),
            "matched_strategy": m["matched_strategy"],
            "status": STATUS_BELOW if below else STATUS_HEALTHY,
            "alert": below,
            "notes": ctrl.get("notes", ""),
        })

    out = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    # Machine-readable breakdown first, so absent controls stay visible even
    # when the verdict line is the quiet OK branch.
    print(f"positive controls: total={len(rows)} healthy={n_healthy} "
          f"below_percentile={n_below} not_found={n_absent}", file=sys.stderr)
    if n_below:
        print(f"WARN: {n_below}/{len(rows)} positive controls found but below "
              f"the {args.alert_percentile}th percentile (possible ranking "
              f"drift); {n_absent} not found.", file=sys.stderr)
    elif n_absent:
        print(f"OK: {n_healthy}/{len(rows)} positive controls found and above "
              f"the {args.alert_percentile}th percentile; {n_absent} not found "
              f"in the ranked CSV (informational, not an alert).",
              file=sys.stderr)
    else:
        print(f"OK: all {len(rows)} positive controls present and above "
              f"the {args.alert_percentile}th percentile.", file=sys.stderr)
    return 0  # non-fatal: this is a sanity check, not a gate


if __name__ == "__main__":
    sys.exit(main())
