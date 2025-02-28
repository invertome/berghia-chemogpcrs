#!/usr/bin/env python3
# plot_synteny.py
# Purpose: Visualize synteny blocks.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import matplotlib.pyplot as plt
import pandas as pd
import glob
import sys

collinearity_files = glob.glob(sys.argv[1])
output_file = sys.argv[2]

plt.figure(figsize=(12, 8))
for file in collinearity_files:
    data = pd.read_csv(file, sep='\t', comment='#', names=['block', 'gene1', 'gene2', 'score'])
    plt.scatter(data['gene1'], data['gene2'], s=10, label=file.split('/')[-1].replace('_mcscanx.collinearity', ''))
plt.xlabel('Berghia Genes')
plt.ylabel('Reference Genes')
plt.title('Synteny Blocks Across Genomes')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
