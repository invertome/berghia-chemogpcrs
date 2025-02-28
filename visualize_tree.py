#!/usr/bin/env python3
# visualize_tree.py
# Purpose: Visualize phylogenetic trees with multiple layouts and receptor counts.
# Inputs: Tree file ($1), output prefix ($2)
# Outputs: Basic (${output_prefix}_basic.png), receptor (${output_prefix}_receptors.png), circular (${output_prefix}_circular.png) plots
# Logic: Plots basic tree, adds receptor counts if metadata available, and provides circular layout.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

tree_file = sys.argv[1]
output_prefix = sys.argv[2]

t = Tree(tree_file)
ts = TreeStyle()
ts.show_leaf_name = True

# Basic tree plot
t.render(f"{output_prefix}_basic.png", w=800, units='px', tree_style=ts)

# Plot with receptor counts (assuming receptor_count in metadata, adjust as needed)
for leaf in t:
    if hasattr(leaf, 'receptor_count'):  # Placeholder for actual metadata
        leaf.add_face(TextFace(f"Receptors: {leaf.receptor_count}"), column=0, position="branch-right")
t.render(f"{output_prefix}_receptors.png", w=800, units='px', tree_style=ts)

# Circular layout
ts_circular = TreeStyle()
ts_circular.mode = "c"
ts_circular.show_leaf_name = True
t.render(f"{output_prefix}_circular.png", w=800, units='px', tree_style=ts_circular)
