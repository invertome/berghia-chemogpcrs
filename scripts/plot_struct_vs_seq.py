#!/usr/bin/env python3
# plot_struct_vs_seq.py
# Purpose: Compare structural and sequence trees using generalized RF distance, accounting for differing taxon sets.
# Inputs: Structural tree ($1), sequence tree ($2), output prefix ($3)
# Outputs: Comparison plot in ${output_prefix}.png
# Logic: Prunes both trees to common taxa before RF calculation.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from ete3 import Tree
import matplotlib.pyplot as plt

struct_tree_file = sys.argv[1]
seq_tree_file = sys.argv[2]
output_prefix = sys.argv[3]

t_struct = Tree(struct_tree_file)
t_seq = Tree(seq_tree_file)

# Prune to common taxa
common_taxa = set(t_struct.get_leaf_names()).intersection(set(t_seq.get_leaf_names()))
if not common_taxa:
    print("Error: No common taxa between trees")
    sys.exit(1)
t_seq.prune(common_taxa, preserve_branch_length=True)
t_struct.prune(common_taxa, preserve_branch_length=True)

# Calculate generalized RF distance
rf, max_rf, _, _, _, _, _ = t_struct.robinson_foulds(t_seq, unrooted_trees=True)
grf = rf / max_rf if max_rf > 0 else 0

plt.figure(figsize=(10, 6))
plt.text(0.5, 0.5, f"Generalized RF Distance: {grf:.4f}", ha='center', va='center', fontsize=14)
plt.title("Structural vs. Sequence Tree Comparison")
plt.axis('off')
plt.savefig(f"{output_prefix}.png", dpi=300)
plt.close()
