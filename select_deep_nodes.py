#!/usr/bin/env python3
# select_deep_nodes.py
# Purpose: Select deep nodes in a tree for ASR based on a minimum distance threshold.
# Inputs: Tree file ($1), taxid ($2), min distance ($3 from config.sh as MIN_ASR_DISTANCE)
# Outputs: Space-separated list of node IDs
# Logic: Traverses tree, selects nodes with distance > MIN_ASR_DISTANCE containing the specified taxid in descendants.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from ete3 import Tree

tree_file = sys.argv[1]
taxid = sys.argv[2]
min_distance = float(sys.argv[3])

t = Tree(tree_file)
deep_nodes = []
for node in t.traverse():
    if node.dist > min_distance:
        leaves = [leaf.name.split('_')[0] for leaf in node.get_leaves()]
        if taxid in leaves:
            deep_nodes.append(node.name or f"node_{id(node)}")

print(" ".join(deep_nodes))
