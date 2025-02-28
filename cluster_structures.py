#!/usr/bin/env python3
# cluster_structures.py
# Purpose: Cluster protein structures using UPGMA based on TM-align scores.
# Inputs: TM-align scores CSV ($1), output tree file ($2)
# Outputs: UPGMA tree in Newick format
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import sys
from scipy.cluster.hierarchy import linkage, to_tree
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
dist_matrix = pd.DataFrame(1.0, index=ids, columns=ids)
for _, row in scores.iterrows():
    dist_matrix.loc[row['id1'], row['id2']] = 1 - row['tm_score']
    dist_matrix.loc[row['id2'], row['id1']] = 1 - row['tm_score']
dist_matrix.values[[range(n)]*2] = 0  # Diagonal set to 0

# Perform UPGMA clustering
linkage_matrix = linkage(dist_matrix.values, method='average')
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
