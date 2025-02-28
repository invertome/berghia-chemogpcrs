#!/usr/bin/env python3
# plot_synteny.py
# Purpose: Visualize synteny blocks across genomes.
# Inputs: Collinearity files (glob pattern $1), output plot file ($2)
# Outputs: Synteny plot (${output_file}.png, ${output_file}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import matplotlib.pyplot as plt
import pandas as pd
import glob
import sys

# Command-line arguments
collinearity_files = glob.glob(sys.argv[1])  # Glob pattern for collinearity files
output_file = sys.argv[2]  # Output plot file (without extension)

# Initialize plot
plt.figure(figsize=(12, 8))

# Plot each synteny block from collinearity files
for file in collinearity_files:
    data = pd.read_csv(file, sep='\t', comment='#', names=['block', 'gene1', 'gene2', 'score'])
    plt.scatter(data['gene1'], data['gene2'], s=10, label=file.split('/')[-1].replace('_mcscanx.collinearity', ''))

# Customize plot
plt.xlabel('Berghia Genes')
plt.ylabel('Reference Genes')
plt.title('Synteny Blocks Across Genomes')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save in both PNG and SVG formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
plt.close()
