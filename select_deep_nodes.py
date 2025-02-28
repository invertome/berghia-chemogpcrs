#!/usr/bin/env python3
# select_deep_nodes.py
# Purpose: Select deep nodes in a phylogenetic tree for ancestral sequence reconstruction (ASR) based on subtree size and diversity.
# Inputs: Tree file ($1), Berghia TaxID ($2)
# Outputs: Space-separated list of node names to stdout
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree, NCBITaxa

tree_file = sys.argv[1]  # Input tree file (Newick format)
berghia_taxid = int(sys.argv[2])  # TaxID for Berghia stephanieae

# Load taxonomic data and tree
ncbi = NCBITaxa()
t = Tree(tree_file)

# Identify deep nodes with significant Berghia descendants and taxonomic diversity
deep_nodes = []
for node in t.traverse():
    if not node.is_leaf():
        descendants = node.get_leaves()
        taxids = {leaf.name.split('_')[0] for leaf in descendants}
        berghia_count = sum(1 for leaf in descendants if leaf.name.startswith(str(berghia_taxid)))
        if berghia_count > 2 and len(taxids) > 1 and node.dist > 0.5:  # Criteria: >2 Berghia, diverse taxa, significant branch length
            deep_nodes.append(node.name)

# Output node names
print(" ".join(deep_nodes))
