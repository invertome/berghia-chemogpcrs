#!/usr/bin/env python3
# visualize_tree.py
# Purpose: Generate PNG and HTML visualizations of phylogenetic trees.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree, TreeStyle

tree_file = sys.argv[1]
output_prefix = sys.argv[2]

t = Tree(tree_file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 50
t.render(f"{output_prefix}.png", w=800, units='px', tree_style=ts, dpi=300)
t.render(f"{output_prefix}.html", w=800, units='px', tree_style=ts)
