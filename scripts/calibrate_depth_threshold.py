#!/usr/bin/env python3
"""calibrate_depth_threshold.py — Null-calibrated branch-distance threshold.

Bead -m6k. The previous pipeline used within-dataset percentiles (75th for
LSE depth, 90th for ASR deep nodes) as thresholds. By construction such
percentiles always flag the same proportion of nodes regardless of biological
signal — they are descriptive cuts, not statistical tests.

This script generates a null distribution of internal-node branch distances
under a Yule birth-death model with no LSEs (constant-rate birth, low or
zero death), then writes a JSON file with:
  - p95_depth: 95th percentile across all simulated trees
  - p99_depth: 99th percentile (more conservative)
  - p_max_per_tree: distribution of maximum depth per tree (for empirical p)

Uses dendropy.simulate.treesim (https://dendropy.org/), the canonical Python
phylogeny simulator. Do not roll our own.

Usage:
    python calibrate_depth_threshold.py \\
        --n-taxa 2400 --n-sims 1000 \\
        --birth-rate 1.0 --death-rate 0.0 \\
        --out results/calibration/depth_thresholds.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def simulate_null_depths(n_taxa: int, n_sims: int,
                         birth_rate: float = 1.0,
                         death_rate: float = 0.0,
                         seed: int = 12345) -> dict:
    """Simulate Yule trees and collect internal-node depth distributions.

    Returns:
        Dict with keys:
          all_depths: flattened list of all internal-node depths across sims
          per_tree_max: list of max depth per simulated tree
    """
    import dendropy
    import random
    from dendropy.simulate import treesim

    # dendropy expects stdlib random.Random (not numpy.random.Generator)
    rng = random.Random(seed)
    all_depths: list[float] = []
    per_tree_max: list[float] = []

    for sim_i in range(n_sims):
        t = treesim.birth_death_tree(
            birth_rate=birth_rate,
            death_rate=death_rate,
            num_extant_tips=n_taxa,
            repeat_until_success=True,
            rng=rng,
        )
        depths = [float(node.distance_from_root())
                  for node in t.preorder_internal_node_iter()
                  if not node.is_leaf()]
        if not depths:
            continue
        all_depths.extend(depths)
        per_tree_max.append(max(depths))

    return {"all_depths": all_depths, "per_tree_max": per_tree_max}


def summarize(d: dict) -> dict:
    arr = np.array(d["all_depths"], dtype=float)
    pt = np.array(d["per_tree_max"], dtype=float)
    return {
        "n_simulations": int(len(pt)),
        "n_internal_nodes": int(len(arr)),
        "depth_mean": float(arr.mean()) if len(arr) else float("nan"),
        "depth_median": float(np.median(arr)) if len(arr) else float("nan"),
        "depth_p75": float(np.percentile(arr, 75)) if len(arr) else float("nan"),
        "depth_p90": float(np.percentile(arr, 90)) if len(arr) else float("nan"),
        "depth_p95": float(np.percentile(arr, 95)) if len(arr) else float("nan"),
        "depth_p99": float(np.percentile(arr, 99)) if len(arr) else float("nan"),
        "per_tree_max_mean": float(pt.mean()) if len(pt) else float("nan"),
        "per_tree_max_p95": float(np.percentile(pt, 95)) if len(pt) else float("nan"),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--n-taxa", type=int, default=2400,
                    help="Number of extant tips to simulate (default: 2400)")
    ap.add_argument("--n-sims", type=int, default=1000,
                    help="Number of trees to simulate (default: 1000)")
    ap.add_argument("--birth-rate", type=float, default=1.0)
    ap.add_argument("--death-rate", type=float, default=0.0,
                    help="Yule (default 0) or birth-death (>0)")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", required=True,
                    help="Output JSON path")
    args = ap.parse_args()

    print(f"Simulating {args.n_sims} Yule trees with {args.n_taxa} tips "
          f"(birth={args.birth_rate}, death={args.death_rate}, seed={args.seed})...",
          file=sys.stderr)
    d = simulate_null_depths(args.n_taxa, args.n_sims,
                             args.birth_rate, args.death_rate, args.seed)
    summary = summarize(d)
    summary["params"] = {
        "n_taxa": args.n_taxa, "n_sims": args.n_sims,
        "birth_rate": args.birth_rate, "death_rate": args.death_rate,
        "seed": args.seed,
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote thresholds to {args.out}", file=sys.stderr)
    print(f"  p95 internal-node depth: {summary['depth_p95']:.4f}", file=sys.stderr)
    print(f"  p99 internal-node depth: {summary['depth_p99']:.4f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
