#!/usr/bin/env python3
# visualize_tree.py
# Purpose: Generate PNG and SVG visualizations of phylogenetic trees with bootstrap/posterior support.
# Inputs: Tree file ($1), output prefix ($2)
# Outputs: Tree visualizations (${output_prefix}.png, ${output_prefix}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

tree_file = sys.argv[1]  # Input tree file (Newick format)
output_prefix = sys.argv[2]  # Output prefix for visualization files

# Load tree
t = Tree(tree_file)

# Define tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 50  # Scale branch lengths visually

# Add bootstrap/posterior support to nodes
for node in t.traverse():
    if not node.is_leaf() and node.support is not None:
        support_label = f"{node.support:.2f}" if node.support <= 1 else f"{int(node.support)}"
        node.add_face(TextFace(support_label, fsize=8), column=0, position="branch-top")

# Render tree in PNG and SVG formats
t.render(f"{output_prefix}.png", w=800, units='px', tree_style=ts, dpi=300)
t.render(f"{output_prefix}.svg", w=800, units='px', tree_style=ts)
