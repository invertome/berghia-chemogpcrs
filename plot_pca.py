#!/usr/bin/env python3
# plot_pca.py
# Purpose: Generate a PCA plot of TM-align scores to visualize structural relationships.
# Inputs: TM-align scores CSV ($1), output plot file ($2)
# Outputs: PCA plot (${output_file}.png, ${output_file}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
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

# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(pivot)

# Plot PCA
plt.figure(figsize=(10, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1], c='blue', label=f'Explained variance: {sum(pca.explained_variance_ratio_):.2f}')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA of TM-align Scores')
plt.legend()
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
plt.close()
