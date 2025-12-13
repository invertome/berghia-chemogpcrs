#!/usr/bin/env python3
# plot_asr.py
# Purpose: Visualize ASR sequences on a tree with sophisticated plotting.
# Inputs: Tree file ($1), ASR FASTA ($2), output prefix ($3)
# Outputs: Basic (${output_prefix}.png), circular (${output_prefix}_circular.png) plots
# Logic: Highlights ASR nodes, adds sequence length annotations, and provides circular layout.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

tree_file = sys.argv[1]
asr_fasta = sys.argv[2]
output_prefix = sys.argv[3]

# Load ASR sequences
asr_ids = [line[1:].strip() for line in open(asr_fasta) if line.startswith('>')]
seq_lengths = {}
with open(asr_fasta, 'r') as f:
    seq_id = None
    seq = ''
    for line in f:
        if line.startswith('>'):
            if seq_id:
                seq_lengths[seq_id] = len(seq)
            seq_id = line[1:].strip()
            seq = ''
        else:
            seq += line.strip()
    if seq_id:
        seq_lengths[seq_id] = len(seq)

# Load and style tree
t = Tree(tree_file)
for node in t.traverse():
    if node.name in asr_ids:
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "red"
        node.set_style(nstyle)
        if node.name in seq_lengths:
            node.add_face(TextFace(f"Len: {seq_lengths[node.name]}"), column=0, position="branch-right")

# Basic plot
t.render(f"{output_prefix}.png", w=800, units='px')

# Circular plot
ts_circular = TreeStyle()
ts_circular.mode = "c"
ts_circular.show_leaf_name = True
t.render(f"{output_prefix}_circular.png", w=800, units='px', tree_style=ts_circular)
