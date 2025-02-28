#!/usr/bin/env python3
# plot_phyloformer_iqtree.py
# Purpose: Compare Phyloformer and IQ-TREE trees.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from ete3 import Tree
import matplotlib.pyplot as plt

iqtree_file = sys.argv[1]
phyloformer_file = sys.argv[2]
output_file = sys.argv[3]

t_iq = Tree(iqtree_file)
t_ph = Tree(phyloformer_file)

rf_dist = t_iq.robinson_foulds(t_ph)[0]

plt.figure(figsize=(10, 6))
plt.text(0.5, 0.5, f"Robinson-Foulds Distance: {rf_dist}", ha='center', va='center', fontsize=14)
plt.title("Phyloformer vs. IQ-TREE Comparison")
plt.axis('off')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
