#!/usr/bin/env python3
# select_deep_nodes.py
# Purpose: Select deep nodes for ASR based on subtree size and taxonomic diversity.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree, NCBITaxa

tree_file = sys.argv[1]
berghia_taxid = int(sys.argv[2])

ncbi = NCBITaxa()
t = Tree(tree_file)

deep_nodes = []
for node in t.traverse():
    if not node.is_leaf():
        descendants = node.get_leaves()
        taxids = {leaf.name.split('_')[0] for leaf in descendants}
        berghia_count = sum(1 for leaf in descendants if leaf.name.startswith(str(berghia_taxid)))
        if berghia_count > 2 and len(taxids) > 1 and node.dist > 0.5:
            deep_nodes.append(node.name)

print(" ".join(deep_nodes))
