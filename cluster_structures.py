#!/usr/bin/env python3
# cluster_structures.py
# Purpose: Cluster protein structures using UPGMA based on TM-align scores.
# Inputs: TM-align scores CSV ($1), output tree file ($2)
# Outputs: UPGMA tree in Newick format
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import numpy as np
import sys
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from ete3 import Tree

# Command-line arguments
tmalign_scores_file = sys.argv[1]  # Input TM-align scores CSV
output_tree_file = sys.argv[2]     # Output tree file (Newick)

# Load TM-align scores and extract IDs
scores = pd.read_csv(tmalign_scores_file, names=['file', 'tm_score'], sep=',')
scores['id1'] = scores['file'].apply(lambda x: x.split('tmalign_')[1].split('_')[0])
scores['id2'] = scores['file'].apply(lambda x: '_'.join(x.split('tmalign_')[1].split('_')[1:]).replace('.txt', ''))

# Get unique IDs and create distance matrix
ids = list(set(scores['id1']).union(scores['id2']))
n = len(ids)

# Guard against insufficient data
if n < 2:
    print(f"Error: Need at least 2 structures for clustering, got {n}", file=sys.stderr)
    # Write single-node tree if only one structure
    if n == 1:
        Tree(f"{ids[0]};").write(outfile=output_tree_file)
    sys.exit(0)

dist_matrix = pd.DataFrame(1.0, index=ids, columns=ids)
for _, row in scores.iterrows():
    dist_matrix.loc[row['id1'], row['id2']] = 1 - row['tm_score']
    dist_matrix.loc[row['id2'], row['id1']] = 1 - row['tm_score']
np.fill_diagonal(dist_matrix.values, 0)  # Diagonal set to 0

# Perform UPGMA clustering
# Convert square matrix to condensed form (1D) for scipy linkage
condensed_dist = squareform(dist_matrix.values, checks=False)
linkage_matrix = linkage(condensed_dist, method='average')
tree = to_tree(linkage_matrix)

# Convert to Newick format and save
def tree_to_newick(node, ids):
    if node.is_leaf():
        return ids[node.id]
    left = tree_to_newick(node.left, ids)
    right = tree_to_newick(node.right, ids)
    return f"({left},{right}):{node.dist/2}"

newick = tree_to_newick(tree, ids)
Tree(newick).write(outfile=output_tree_file)
