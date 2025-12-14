#!/usr/bin/env python3
# select_deep_nodes.py
# Purpose: Select deep nodes in a tree for ASR based on relative depth metrics.
# Inputs: Tree file ($1), taxid ($2), min distance ($3 from config.sh as MIN_ASR_DISTANCE)
# Outputs: Space-separated list of node IDs
# Logic: Traverses tree, selects nodes in the top percentile of depth containing the specified taxid.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from ete3 import Tree
import numpy as np

tree_file = sys.argv[1]
taxid = sys.argv[2]
min_distance = float(sys.argv[3])  # Used as fallback or minimum threshold

t = Tree(tree_file, format=1)  # format=1 for IQ-TREE output with support values

# Calculate distances from root for all internal nodes
node_distances = []
nodes_with_taxid = []

for node in t.traverse():
    if not node.is_leaf():
        # Get distance from root
        dist_from_root = t.get_distance(node)

        # Check if this node has descendants with the specified taxid
        leaves = [leaf.name.split('_')[0] for leaf in node.get_leaves()]
        has_taxid = taxid in leaves

        node_distances.append(dist_from_root)

        if has_taxid:
            nodes_with_taxid.append({
                'node': node,
                'dist': dist_from_root,
                'name': node.name or f"node_{id(node)}"
            })

if not nodes_with_taxid:
    # No nodes with the specified taxid
    print("")
    sys.exit(0)

# Calculate relative threshold: use top 10% deepest nodes OR nodes above min_distance
# whichever selects more nodes (to ensure we get meaningful ASR targets)
all_dists = [n['dist'] for n in nodes_with_taxid]

if len(all_dists) > 1:
    # Use 90th percentile as threshold for "deep" nodes
    percentile_threshold = np.percentile(all_dists, 90)

    # Use the more permissive threshold (lower value selects more nodes)
    # But ensure we use at least the configured minimum
    effective_threshold = max(min(percentile_threshold, min_distance), min_distance * 0.5)
else:
    # Only one node - use it if above half the minimum distance
    effective_threshold = min_distance * 0.5

# Select deep nodes
deep_nodes = []
for n in nodes_with_taxid:
    if n['dist'] >= effective_threshold:
        deep_nodes.append(n['name'])

# If no nodes pass the threshold but we have candidates, take the deepest one
if not deep_nodes and nodes_with_taxid:
    deepest = max(nodes_with_taxid, key=lambda x: x['dist'])
    deep_nodes.append(deepest['name'])

print(" ".join(deep_nodes))
