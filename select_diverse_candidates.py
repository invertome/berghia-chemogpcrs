#!/usr/bin/env python3
# select_diverse_candidates.py
# Purpose: Select diverse candidates from the phylogeny using hierarchical clustering.
# Inputs: Tree file ($1), ranked candidates CSV ($2), num candidates ($3), output IDs ($4)
# Outputs: File with selected IDs
# Logic: Uses ete3 to compute tree distances, scipy for hierarchical clustering, and selects top-ranked candidate per cluster.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from ete3 import Tree
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster

tree_file = sys.argv[1]
ranking_file = sys.argv[2]
num_candidates = int(sys.argv[3])
output_file = sys.argv[4]

# Load tree and rankings
t = Tree(tree_file)
ranked_df = pd.read_csv(ranking_file)
candidates_in_tree = [id for id in ranked_df['id'].tolist() if id in t]

if len(candidates_in_tree) <= num_candidates:
    selected_ids = candidates_in_tree
else:
    # Compute distance matrix from tree
    dist_matrix = []
    for i in range(len(candidates_in_tree)):
        row = []
        for j in range(len(candidates_in_tree)):
            if i == j:
                row.append(0)
            else:
                row.append(t.get_distance(candidates_in_tree[i], candidates_in_tree[j]))
        dist_matrix.append(row)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(dist_matrix, method='average')
    clusters = fcluster(linkage_matrix, t=num_candidates, criterion='maxclust')
    
    # Select highest-ranked candidate per cluster
    cluster_dict = {}
    for i, cluster_id in enumerate(clusters):
        candidate_id = candidates_in_tree[i]
        rank_score = ranked_df[ranked_df['id'] == candidate_id]['rank_score'].values[0]
        if cluster_id not in cluster_dict or rank_score > cluster_dict[cluster_id][1]:
            cluster_dict[cluster_id] = (candidate_id, rank_score)
    
    selected_ids = [id for id, _ in cluster_dict.values()]

# Write selected IDs
with open(output_file, 'w') as f:
    for id in selected_ids:
        f.write(f"{id}\n")
