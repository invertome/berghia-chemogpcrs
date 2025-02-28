#!/usr/bin/env python3
# plot_struct_vs_seq.py
# Purpose: Compare structural and sequence-based trees.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree
import matplotlib.pyplot as plt

struct_tree_file = sys.argv[1]
seq_tree_file = sys.argv[2]
output_file = sys.argv[3]

t_struct = Tree(struct_tree_file)
t_seq = Tree(seq_tree_file)

rf_dist = t_struct.robinson_foulds(t_seq)[0]

plt.figure(figsize=(10, 6))
plt.text(0.5, 0.5, f"Robinson-Foulds Distance: {rf_dist}", ha='center', va='center', fontsize=14)
plt.title("Structural vs. Sequence Tree Comparison")
plt.axis('off')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
