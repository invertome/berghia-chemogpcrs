#!/usr/bin/env python3
"""rank_confidence.py — signal-bootstrap rank CIs + P(in top-k) + tiers.

Ranks are estimated with heterogeneous precision (Hall & Miller 2009): top/
bottom are stable, the middle is not. Under the production RRA aggregator
(weight-free, so weight-perturbation is uninformative), we resample the SIGNAL
SET with replacement, re-aggregate, and read each candidate's rank distribution
-> a rank CI + P(candidate in top-k). The top-k-with-confidence set is the wet-
lab shortlist; mid-list exact ranks are reported only as tiers. (Mid-list RRA
scores are densely packed and, where candidates share a rho, exactly tied and
tie-broken by id -> their rank CIs are legitimately wide, which is the honest
'mid-ranks not separable' statement. Before bead 8k8e the packing was far worse:
the Bonferroni bound clipped 46.5% of the cohort to a single score of 1.0.)

CLI (stage-07 augmenter): reads a ranked CSV, builds the SIGNAL_SPEC ranklists
via rank_aggregation.build_ranklists_from_df, computes intervals + tiers, and
LEFT-JOINS rank_ci_lo / rank_ci_hi / p_top_k / rank_tier into the ranked CSV in
place (same read-as-str, map-by-id, write idiom as add_embedding_columns.py).
"""
from __future__ import annotations

import argparse
import sys
from typing import Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from rank_aggregation import build_ranklists_from_df, rra_score


def bootstrap_rank_intervals(per_signal: Mapping[str, Mapping[str, float]],
                             k: int, n_boot: int = 1000, seed: int = 0,
                             alpha: float = 0.05) -> Dict[str, Dict[str, float]]:
    names = list(per_signal)
    ids = sorted({c for s in per_signal.values() for c in s})
    rng = np.random.default_rng(seed)
    ranks: Dict[str, list] = {i: [] for i in ids}
    in_top: Dict[str, int] = {i: 0 for i in ids}
    for _ in range(n_boot):
        pick = rng.choice(len(names), len(names), replace=True)
        sub = {f"{names[j]}__{idx}": dict(per_signal[names[j]])
               for idx, j in enumerate(pick)}
        rho = rra_score(sub)
        order = sorted(rho, key=lambda i: (rho[i], i))
        pos = {i: r + 1 for r, i in enumerate(order)}
        for i in ids:
            r = pos.get(i, len(ids))
            ranks[i].append(r)
            if r <= k:
                in_top[i] += 1
    out: Dict[str, Dict[str, float]] = {}
    for i in ids:
        arr = np.array(ranks[i])
        out[i] = {
            "rank_ci_lo": int(np.quantile(arr, alpha / 2)),
            "rank_ci_hi": int(np.quantile(arr, 1 - alpha / 2)),
            "rank_median": float(np.median(arr)),
            "p_top_k": in_top[i] / n_boot,
        }
    return out


def assign_tiers(intervals: Mapping[str, Mapping[str, float]],
                 high: float = 0.8, plausible: float = 0.2) -> Dict[str, str]:
    out = {}
    for i, d in intervals.items():
        p = d["p_top_k"]
        out[i] = "high" if p >= high else ("plausible" if p >= plausible else "tail")
    return out


# Columns this augmenter owns; dropped-then-rejoined so re-runs stay idempotent.
_JOIN_COLUMNS = ("rank_ci_lo", "rank_ci_hi", "p_top_k", "rank_tier")


def annotate_ranked_csv(ranked_csv: str, k: int, n_boot: int = 1000,
                        seed: int = 0, id_column: str = "id") -> int:
    """Join rank_ci_lo/rank_ci_hi/p_top_k/rank_tier into ``ranked_csv`` in place.

    Ranklists are built from a numeric read (matching the production rankagg
    path); the in-place write reads the CSV as strings so every existing column
    is preserved byte-for-byte. Candidates covered by no signal are left blank.
    Returns the number of rows written.
    """
    per_signal = build_ranklists_from_df(pd.read_csv(ranked_csv), id_col=id_column)
    intervals = bootstrap_rank_intervals(per_signal, k=k, n_boot=n_boot, seed=seed)
    tiers = assign_tiers(intervals)

    lo = {i: d["rank_ci_lo"] for i, d in intervals.items()}
    hi = {i: d["rank_ci_hi"] for i, d in intervals.items()}
    ptk = {i: d["p_top_k"] for i, d in intervals.items()}

    df = pd.read_csv(ranked_csv, dtype=str, keep_default_na=False)
    df = df.drop(columns=[c for c in _JOIN_COLUMNS if c in df.columns])
    df["rank_ci_lo"] = df[id_column].map(lo)
    df["rank_ci_hi"] = df[id_column].map(hi)
    df["p_top_k"] = df[id_column].map(ptk)
    df["rank_tier"] = df[id_column].map(tiers)
    for c in _JOIN_COLUMNS:
        df[c] = df[c].fillna("")
    df.to_csv(ranked_csv, index=False)
    return len(df)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True,
                    help="Ranked CSV to annotate in place")
    ap.add_argument("--k", type=int, default=20,
                    help="top-k membership threshold (default 20)")
    ap.add_argument("--n-boot", type=int, default=1000,
                    help="signal-bootstrap draws (default 1000)")
    args = ap.parse_args(argv)

    n_rows = annotate_ranked_csv(args.ranked_csv, k=args.k, n_boot=args.n_boot)
    print(f"Annotated {args.ranked_csv} ({n_rows} rows) with rank CIs + tiers",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
