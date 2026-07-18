#!/usr/bin/env python3
"""a1_compare_filter_trees.py — compare A1-tree filter conditions.

Given the FastTree trees built under each PREQUAL x TAPER filter condition (all
sharing the canonical MAFFT-DASH E-INS-i alignment lineage), quantify what the
filters do to the two things A1 actually cares about:

  (1) the candidate -> nearest-reference patristic-distance RANK (A1 residualizes
      on rank(tree_distance), so rank stability vs the canonical/unfiltered tree
      is the decisive metric — a filter that barely changes the ranks is neutral
      for A1), and
  (2) branch support (mean SH-like local support) — does the filter buy a
      better-supported tree?

Reuses tree_distance_to_refs.nearest_ref_distance (unit-tested). Pure analysis;
emits a JSON verdict table. Decision rule is reported, not hard-coded: prefer
the MINIMAL filtering that does not degrade support, because over-masking removes
the ECL/ICL divergence that IS the novelty signal (same rationale the pipeline
already applies to 04b).
"""
from __future__ import annotations

import argparse
import json
from typing import Dict, List

import numpy as np
from scipy.stats import spearmanr

from tree_distance_to_refs import nearest_ref_distance


def _mean_support(newick: str) -> float:
    """Mean internal-node support (FastTree writes SH-like support as the
    internal node label). Uses ete3; leaves and the root are excluded."""
    from ete3 import Tree

    t = Tree(newick, format=0)
    sup = [n.support for n in t.traverse() if not n.is_leaf() and n.support is not None]
    return float(np.mean(sup)) if sup else float("nan")


def compare(tree_paths: Dict[str, str], ref_ids: List[str],
            baseline: str = "canonical") -> Dict[str, object]:
    """Per-condition distance + support, and rank-correlation vs the baseline."""
    dists: Dict[str, Dict[str, float]] = {}
    support: Dict[str, float] = {}
    for name, path in tree_paths.items():
        with open(path) as fh:
            nwk = fh.read()
        dists[name] = nearest_ref_distance(nwk, ref_ids)
        support[name] = _mean_support(nwk)

    base = dists.get(baseline, {})
    rows = {}
    for name, d in dists.items():
        shared = sorted(set(d) & set(base))
        if name != baseline and len(shared) >= 3:
            rho, _ = spearmanr([d[i] for i in shared], [base[i] for i in shared])
        else:
            rho = 1.0 if name == baseline else float("nan")
        rows[name] = {
            "n_candidates_scored": len(d),
            "mean_support": round(support[name], 4),
            "dist_rank_spearman_vs_baseline": round(float(rho), 4),
            "median_nearest_ref_dist": round(float(np.median(list(d.values()))), 4) if d else None,
        }
    return {"baseline": baseline, "conditions": rows}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--tree", action="append", default=[], metavar="NAME=PATH",
                    help="a condition tree as NAME=PATH (repeatable)")
    ap.add_argument("--ref-ids", required=True, help="reference (anchor) id file, one per line")
    ap.add_argument("--baseline", default="canonical")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    tree_paths = {}
    for spec in args.tree:
        name, path = spec.split("=", 1)
        tree_paths[name] = path
    with open(args.ref_ids) as fh:
        ref_ids = [ln.strip() for ln in fh if ln.strip()]

    result = compare(tree_paths, ref_ids, baseline=args.baseline)
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
