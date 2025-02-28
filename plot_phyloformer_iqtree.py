#!/usr/bin/env python3
# plot_phyloformer_iqtree.py
# Purpose: Compare Phyloformer and IQ-TREE trees, optionally including MrBayes, with RF distance.
# Inputs: IQ-TREE file ($1), Phyloformer file ($2), output file ($3)
# Outputs: Comparison plot (${output_file}.png, ${output_file}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree
import matplotlib.pyplot as plt
import os

iqtree_file = sys.argv[1]  # IQ-TREE output file
phyloformer_file = sys.argv[2]  # Phyloformer output file
output_file = sys.argv[3]  # Output plot file (without extension)

# Loadtrees
t_iq = Tree(iqtree_file)
t_ph = Tree(phyloformer_file)

# Calculate RF distance between IQ-TREE and Phyloformer
rf_dist_ph = t_iq.robinson_foulds(t_ph)[0]

# Check for optional MrBayes tree
mrbayes_file = iqtree_file.replace('.treefile', '.nex.con.tre')
text = f"IQ-TREE vs. Phyloformer RF Distance: {rf_dist_ph}"
if os.path.exists(mrbayes_file):
    t_mb = Tree(mrbayes_file)
    rf_dist_mb = t_iq.robinson_foulds(t_mb)[0]
    text += f"\nIQ-TREE vs. MrBayes RF Distance: {rf_dist_mb}"

# Create plot
plt.figure(figsize=(10, 6))
plt.text(0.5, 0.5, text, ha='center', va='center', fontsize=14)
plt.title("Tree Comparison")
plt.axis('off')
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
plt.close()
