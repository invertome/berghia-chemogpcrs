#!/usr/bin/env python3
"""compare_rankings.py — Pairwise comparison of ranked-candidate CSVs.

Bead -ab1 (option-G ablation framework). Given two or more ranked CSVs
(produced by rank_candidates.py under different parameter / filter
configurations), compute robustness metrics:

    Spearman rank correlation     — overall rank-order agreement
    Kendall tau                   — pairwise concordance
    Jaccard top-N overlap         — set-overlap of top-N candidates
    Top-N rank shift              — median |Δrank| within top-N

Usage:
    # All-pairs comparison across configs
    python compare_rankings.py \\
        --csv default=results/ranking/ranked_candidates_sorted.csv \\
        --csv no_dnds=results/ablations/no_dnds/ranked.csv \\
        --csv no_phylo=results/ablations/no_phylo/ranked.csv \\
        --top-ns 10,25,50,100 \\
        --out results/ablations/comparison.tsv

The output TSV has one row per (config_a, config_b, N) triple with
all four metrics. JSON summary alongside.
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.stats import kendalltau, spearmanr


def load_ranked(path: str, id_col: str = "id",
                score_col: str = "rank_score") -> list[str]:
    """Read a ranked CSV, sort by score (descending), return list of IDs in order."""
    df = pd.read_csv(path)
    if id_col not in df.columns:
        raise ValueError(f"{path}: missing {id_col!r} column")
    if score_col in df.columns:
        df = df.sort_values(score_col, ascending=False)
    return df[id_col].astype(str).tolist()


def jaccard_top_n(a: list[str], b: list[str], n: int) -> float:
    """Jaccard index of the top-N sets of two ranked lists."""
    sa = set(a[:n])
    sb = set(b[:n])
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


def median_rank_shift_top_n(a: list[str], b: list[str], n: int) -> float:
    """For IDs in the top-N of `a`, what's the median |rank_in_a − rank_in_b|?

    IDs absent from `b` get rank=len(b)+1 (capped). Returns -1 if a has fewer
    than n members."""
    if len(a) < n:
        n = len(a)
    if n == 0:
        return -1.0
    rank_a = {x: i for i, x in enumerate(a)}
    rank_b = {x: i for i, x in enumerate(b)}
    cap = len(b)
    shifts: list[int] = []
    for x in a[:n]:
        ra = rank_a[x]
        rb = rank_b.get(x, cap)
        shifts.append(abs(ra - rb))
    shifts.sort()
    mid = len(shifts) // 2
    if len(shifts) % 2 == 1:
        return float(shifts[mid])
    return (shifts[mid - 1] + shifts[mid]) / 2


def compare_two(a_name: str, a_ranked: list[str],
                b_name: str, b_ranked: list[str],
                ns: list[int]) -> list[dict]:
    """Compute Spearman, Kendall tau, Jaccard top-N, and median rank shift."""
    common = sorted(set(a_ranked) & set(b_ranked))
    if not common:
        return [{
            "config_a": a_name, "config_b": b_name, "N": n,
            "spearman_rho": None, "kendall_tau": None,
            "jaccard_topN": 0.0, "median_rank_shift_topN": -1.0,
            "n_common": 0,
        } for n in ns]

    rank_a = {x: i for i, x in enumerate(a_ranked)}
    rank_b = {x: i for i, x in enumerate(b_ranked)}
    a_ranks = [rank_a[x] for x in common]
    b_ranks = [rank_b[x] for x in common]
    rho = spearmanr(a_ranks, b_ranks).correlation if len(common) > 2 else None
    tau = kendalltau(a_ranks, b_ranks).correlation if len(common) > 2 else None

    rows: list[dict] = []
    for n in ns:
        rows.append({
            "config_a": a_name, "config_b": b_name, "N": n,
            "spearman_rho": float(rho) if rho is not None else None,
            "kendall_tau": float(tau) if tau is not None else None,
            "jaccard_topN": jaccard_top_n(a_ranked, b_ranked, n),
            "median_rank_shift_topN": median_rank_shift_top_n(a_ranked, b_ranked, n),
            "n_common": len(common),
        })
    return rows


def parse_csv_args(args: list[str]) -> dict[str, str]:
    """Parse --csv name=path entries into a dict (preserves order)."""
    out: dict[str, str] = {}
    for entry in args:
        if "=" not in entry:
            raise ValueError(f"--csv expects name=path, got {entry!r}")
        name, path = entry.split("=", 1)
        out[name.strip()] = path.strip()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--csv", action="append", required=True,
                    help="name=path entry (repeat for each config)")
    ap.add_argument("--top-ns", default="10,25,50,100",
                    help="Comma-separated top-N values for set overlap")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--id-col", default="id")
    ap.add_argument("--score-col", default="rank_score")
    args = ap.parse_args()

    configs = parse_csv_args(args.csv)
    if len(configs) < 2:
        print("ERROR: at least 2 --csv entries required for pairwise comparison",
              file=sys.stderr)
        return 1
    ns = [int(x) for x in args.top_ns.split(",") if x.strip()]

    rankings = {name: load_ranked(path, args.id_col, args.score_col)
                for name, path in configs.items()}

    rows: list[dict] = []
    for a_name, b_name in combinations(configs.keys(), 2):
        rows.extend(compare_two(a_name, rankings[a_name],
                                b_name, rankings[b_name], ns))

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    summary_path = Path(args.out).with_suffix(".summary.json")
    summary = {
        "configs": list(configs.keys()),
        "config_paths": configs,
        "ns": ns,
        "n_pairs": len(list(combinations(configs.keys(), 2))),
        "rows": rows,
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote {args.out} and {summary_path}", file=sys.stderr)
    print(f"  {len(configs)} configs, {len(list(combinations(configs.keys(), 2)))} pairs, "
          f"top-N: {ns}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
