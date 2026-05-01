#!/usr/bin/env python3
# select_deep_nodes.py
# Purpose: Select deep internal nodes containing a focal taxon for ancestral
# sequence reconstruction (ASR), based on a null-calibrated or in-dataset
# percentile threshold of internal-node depths.
#
# Bead -mqt fixes:
#   1. Taxid extraction: previous version used `leaf.name.split('_')[0]` which
#      always returned 'ref' for reference leaves named like 'ref_TAXID_N',
#      so NCBI-tax-id-based filtering only worked for non-reference leaves.
#      Now uses a structured extraction matching update_headers.py /
#      lse_refine.py conventions.
#   2. Percentile threshold: previous version computed the 90th percentile
#      over ALL internal-node depths, which biased the threshold downward
#      (reference-only deep nodes inflated the distribution). Now computed
#      over only nodes containing the focal taxon.
#
# Bead -m6k: optional null-calibrated threshold from
# scripts/calibrate_depth_threshold.py (use --null-threshold-file).
#
# Inputs: tree_file, focal taxon prefix, min distance (positional, legacy)
# Outputs: space-separated list of node IDs
# Author: Jorge L. Perez-Moreno, Ph.D.

import argparse
import json
import os
import sys

import numpy as np
from ete3 import Tree


def extract_focal_token(leaf_name: str) -> str:
    """Extract the focal-taxon token from a leaf name.

    Header conventions in this pipeline (see update_headers.py):
      - Reference: ``ref_TAXID_N`` (or ``ref_outgroup_TAXID_N``) where
        TAXID is an NCBI taxon id (numeric).
      - Berghia / focal: ``BERGHIA_FILE_PREFIX``-prefixed (e.g.
        ``1287507_berghia_stephanieae_TRINITY_DN1``), or ``taxid_berghia_*``
        legacy form. Header check needs to compare the FULL prefix (taxid +
        genus + species), not just split on first underscore.

    For ASR-deep-node filtering we want to know whether a leaf "belongs to"
    the focal taxon. The simplest reliable check is prefix-match: leaf name
    starts with the focal token (the user-supplied --taxon arg, e.g.
    ``1287507_berghia_stephanieae`` or just ``1287507``). We support both.
    """
    return leaf_name


def leaf_belongs_to_focal(leaf_name: str, focal_prefix: str) -> bool:
    """Return True iff this leaf is from the focal taxon.

    Matches on:
      - exact equality
      - prefix match (focal followed by an underscore separator)
      - substring of the form '_<focal>_' (refs like ref_TAXID_N where
        focal is the numeric TAXID)
    """
    if leaf_name == focal_prefix:
        return True
    if leaf_name.startswith(focal_prefix + "_"):
        return True
    if f"_{focal_prefix}_" in leaf_name:
        return True
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("tree_file", help="IQ-TREE / FastTree treefile")
    ap.add_argument("focal_taxon",
                    help="Focal taxon token to look for in leaf names "
                         "(e.g. '1287507_berghia_stephanieae' or '1287507')")
    ap.add_argument("min_distance", type=float,
                    help="Fallback minimum branch distance from root")
    ap.add_argument("--percentile", type=float, default=90.0,
                    help="Percentile threshold over focal-containing-node "
                         "depths (default 90)")
    ap.add_argument("--null-threshold-file", default=None,
                    help="JSON from calibrate_depth_threshold.py with "
                         "depth_p95/depth_p99; use null-calibrated cutoff "
                         "instead of within-dataset percentile when provided")
    ap.add_argument("--null-percentile-key", default="depth_p95",
                    choices=["depth_p75", "depth_p90", "depth_p95", "depth_p99"],
                    help="Which key to read from the null-threshold JSON "
                         "(default depth_p95)")
    args = ap.parse_args()

    t = Tree(args.tree_file, format=1)  # IQ-TREE output with support values

    nodes_with_taxon = []
    for node in t.traverse():
        if node.is_leaf():
            continue
        dist_from_root = t.get_distance(node)
        leaves = node.get_leaves()
        has_focal = any(leaf_belongs_to_focal(l.name, args.focal_taxon)
                        for l in leaves)
        if has_focal:
            nodes_with_taxon.append({
                "node": node,
                "dist": dist_from_root,
                "name": node.name or f"node_{id(node)}",
            })

    if not nodes_with_taxon:
        print("", end="")
        return 0

    focal_dists = [n["dist"] for n in nodes_with_taxon]

    # Decide threshold
    if args.null_threshold_file and os.path.exists(args.null_threshold_file):
        with open(args.null_threshold_file) as f:
            null_summary = json.load(f)
        null_threshold = float(null_summary.get(args.null_percentile_key, 0.0))
        # Use null threshold; fall back to within-dataset percentile only if
        # the null cutoff would select nothing.
        in_dataset_p = (np.percentile(focal_dists, args.percentile)
                        if len(focal_dists) > 1 else focal_dists[0])
        # Use the LESS conservative of the two (whichever cutoff selects more
        # candidates) but never below configured min_distance.
        effective = max(min(null_threshold, in_dataset_p),
                        args.min_distance * 0.5)
        threshold_kind = "null+within-dataset"
    else:
        if len(focal_dists) > 1:
            within = np.percentile(focal_dists, args.percentile)
            effective = max(within, args.min_distance * 0.5)
        else:
            effective = args.min_distance * 0.5
        threshold_kind = "within-dataset percentile"

    deep_nodes = [n["name"] for n in nodes_with_taxon if n["dist"] >= effective]
    if not deep_nodes:
        # Fallback: take the single deepest focal-containing node
        deepest = max(nodes_with_taxon, key=lambda x: x["dist"])
        deep_nodes.append(deepest["name"])

    print(" ".join(deep_nodes))
    print(f"# threshold={effective:.4f} ({threshold_kind}); "
          f"selected {len(deep_nodes)}/{len(nodes_with_taxon)} focal-containing nodes",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
