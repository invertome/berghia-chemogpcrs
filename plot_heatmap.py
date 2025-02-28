#!/usr/bin/env python3
# plot_heatmap.py
# Purpose: Generate a heatmap of TM-align scores for structural similarity.
# Inputs: TM-align scores CSV ($1), output plot file ($2)
# Outputs: Heatmap plot (${output_file}.png, ${output_file}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Command-line arguments
tmalign_scores_file = sys.argv[1]  # Input TM-align scores CSV
output_file = sys.argv[2]          # Output plot file (without extension)

# Load and process TM-align scores
scores = pd.read_csv(tmalign_scores_file, names=['file', 'tm_score'], sep=',')
scores['id1'] = scores['file'].apply(lambda x: x.split('tmalign_')[1].split('_')[0])
scores['id2'] = scores['file'].apply(lambda x: '_'.join(x.split('tmalign_')[1].split('_')[1:]).replace('.txt', ''))

# Create similarity matrix
pivot = scores.pivot(index='id1', columns='id2', values='tm_score').fillna(0)

# Plot heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(pivot, cmap='viridis', annot=True, fmt='.2f')
plt.title('TM-align Score Heatmap')
plt.tight_layout()

# Save in both PNG and SVG formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
plt.close()
