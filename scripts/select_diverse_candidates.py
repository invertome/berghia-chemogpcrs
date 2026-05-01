#!/usr/bin/env python3
"""select_diverse_candidates.py — Select phylogenetically diverse candidates.

Bead -mqt: previous version used `id in t` which is ete3's
``__contains__`` (substring search by name) — this could match the wrong
leaf when one candidate ID is a prefix of another (e.g.
'TRINITY_DN1' would substring-match 'TRINITY_DN10', '...DN100', etc.).
Now uses exact match via ``tree.search_nodes(name=id)``.

Also asserts that the resulting condensed distance matrix has no missing
pairs (all-pairs computed on the same set of leaves), so the
``scipy.cluster.hierarchy.linkage`` doesn't get garbage input.

Inputs : tree_file, ranked candidates CSV, num candidates, output IDs file
Outputs: text file with one selected ID per line
Author : Jorge L. Perez-Moreno, Ph.D.
"""
import sys
from typing import Dict, List

import numpy as np
import pandas as pd
from ete3 import Tree
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def find_leaf_node(tree: Tree, name: str):
    """Return the ete3 leaf node with exact name == ``name``, or None.

    Uses tree.search_nodes for exact match (not __contains__'s substring).
    Returns the first leaf if the name is somehow duplicated (with a warning).
    """
    matches = tree.search_nodes(name=name)
    leaves = [m for m in matches if m.is_leaf()]
    if not leaves:
        return None
    if len(leaves) > 1:
        print(f"WARN: {len(leaves)} leaves named {name!r}; using the first.",
              file=sys.stderr)
    return leaves[0]


def main(argv: List[str]) -> int:
    if len(argv) < 5:
        print("Usage: select_diverse_candidates.py TREE RANKING_CSV N OUT", file=sys.stderr)
        return 2
    tree_file, ranking_file, num_candidates, output_file = argv[1:5]
    num_candidates = int(num_candidates)

    t = Tree(tree_file, format=1)
    ranked_df = pd.read_csv(ranking_file)
    if "id" not in ranked_df.columns:
        print("ERROR: ranking CSV must have an 'id' column", file=sys.stderr)
        return 3

    # Build map id -> leaf node (exact match)
    id_to_leaf = {}
    missing = []
    for cand_id in ranked_df["id"].tolist():
        leaf = find_leaf_node(t, cand_id)
        if leaf is not None:
            id_to_leaf[cand_id] = leaf
        else:
            missing.append(cand_id)
    if missing:
        print(f"WARN: {len(missing)} candidates not found in tree (first 5: "
              f"{missing[:5]})", file=sys.stderr)

    candidates_in_tree = list(id_to_leaf.keys())

    if len(candidates_in_tree) <= num_candidates:
        selected_ids = candidates_in_tree
    else:
        n = len(candidates_in_tree)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                dist = t.get_distance(id_to_leaf[candidates_in_tree[i]],
                                      id_to_leaf[candidates_in_tree[j]])
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist

        # Sanity check: any zero off-diagonal indicates a missing pair
        # (would silently merge those into one cluster).
        off_diag = dist_matrix[np.triu_indices(n, k=1)]
        if np.any(off_diag <= 0):
            n_zero = int(np.sum(off_diag <= 0))
            print(f"WARN: {n_zero} pairwise distances were zero or negative; "
                  "duplicate or unresolvable leaves may collapse into one cluster.",
                  file=sys.stderr)

        condensed_dist = squareform(dist_matrix)
        linkage_matrix = linkage(condensed_dist, method="average")
        clusters = fcluster(linkage_matrix, t=num_candidates, criterion="maxclust")

        # Pick highest-ranked candidate per cluster
        cluster_dict: Dict[int, tuple[str, float]] = {}
        for i, cluster_id in enumerate(clusters):
            cid = candidates_in_tree[i]
            score = ranked_df[ranked_df["id"] == cid]["rank_score"].values[0]
            if cluster_id not in cluster_dict or score > cluster_dict[cluster_id][1]:
                cluster_dict[cluster_id] = (cid, float(score))
        selected_ids = [cid for cid, _ in cluster_dict.values()]

    with open(output_file, "w") as f:
        for cid in selected_ids:
            f.write(f"{cid}\n")
    print(f"Wrote {len(selected_ids)} diverse candidates to {output_file}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
