#!/usr/bin/env python3
# plot_asr.py
# Purpose: Visualize ancestral sequence reconstruction (ASR) results on phylogenetic trees.
# Inputs: Tree file ($1), ASR sequences CSV ($2), output plot file ($3)
# Outputs: Tree with ASR annotations (${output_file}.png, ${output_file}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree, TreeStyle, TextFace
import pandas as pd

tree_file = sys.argv[1]  # Input tree file (Newick format)
asr_file = sys.argv[2]   # CSV file with ASR sequences
output_file = sys.argv[3]  # Output plot file (without extension)

# Load tree and style it
t = Tree(tree_file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 50

# Add ASR sequences to nodes
try:
    asr_seqs = pd.read_csv(asr_file, header=None)
    for node in t.traverse():
        if node.name in asr_seqs[0].values:
            node.add_feature('seq', asr_seqs[asr_seqs[0] == node.name][1].values[0])
            node.add_face(TextFace(f"ASR: {node.seq[:10]}..."), column=0, position="branch-right")
    t.render(f"{output_file}.png", w=800, units='px', tree_style=ts, dpi=300)
    t.render(f"{output_file}.svg", w=800, units='px', tree_style=ts)
except Exception as e:
    print(f"Error plotting ASR for {tree_file}: {e}", file=sys.stderr)
    t.render(f"{output_file}_empty.png", w=800, units='px', tree_style=ts, dpi=300)
    t.render(f"{output_file}_empty.svg", w=800, units='px', tree_style=ts)
