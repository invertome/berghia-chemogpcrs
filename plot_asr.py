#!/usr/bin/env python3
# plot_asr.py
# Purpose: Visualize ASR sequences on a tree with sophisticated plotting.
# Inputs: Tree file ($1), ASR FASTA ($2), output prefix ($3)
# Outputs: Basic (${output_prefix}.png), circular (${output_prefix}_circular.png), PDF plots
# Logic: Highlights ASR nodes, adds sequence length annotations, and provides circular layout.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from Bio import SeqIO
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

tree_file = sys.argv[1]
asr_fasta = sys.argv[2]
output_prefix = sys.argv[3]

# Load ASR sequences using Biopython for robust FASTA parsing
asr_ids = []
seq_lengths = {}
for record in SeqIO.parse(asr_fasta, "fasta"):
    asr_ids.append(record.id)
    seq_lengths[record.id] = len(record.seq)

# Load and style tree
t = Tree(tree_file, format=1)  # format=1 for IQ-TREE output with support values
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
