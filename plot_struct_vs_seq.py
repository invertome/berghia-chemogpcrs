#!/usr/bin/env python3
# plot_struct_vs_seq.py
# Purpose: Compare structural and sequence-based phylogenetic trees using Robinson-Foulds distance.
# Inputs: Structural tree ($1), sequence tree ($2), output plot file ($3)
# Outputs: Comparison plot (${output_file}.png, ${output_file}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree
import matplotlib.pyplot as plt

# Command-line arguments
struct_tree_file = sys.argv[1]  # Structural tree file (Newick)
seq_tree_file = sys.argv[2]     # Sequence tree file (Newick)
output_file = sys.argv[3]       # Output plot file (without extension)

# Load trees
t_struct = Tree(struct_tree_file)
t_seq = Tree(seq_tree_file)

# Calculate Robinson-Foulds distance
rf_dist = t_struct.robinson_foulds(t_seq)[0]

# Create comparison plot
plt.figure(figsize=(10, 6))
plt.text(0.5, 0.5, f"Robinson-Foulds Distance: {rf_dist}", ha='center', va='center', fontsize=14)
plt.title("Structural vs. Sequence Tree Comparison")
plt.axis('off')

# Save in both PNG and SVG formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
plt.close()
