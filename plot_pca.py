#!/usr/bin/env python3
# plot_pca.py
# Purpose: Generate PCA plot of TM-align scores.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import sys

tmalign_scores_file = sys.argv[1]
output_file = sys.argv[2]

scores = pd.read_csv(tmalign_scores_file, names=['file', 'tm_score'], sep=',')
scores['id1'] = scores['file'].apply(lambda x: x.split('tmalign_')[1].split('_')[0])
scores['id2'] = scores['file'].apply(lambda x: '_'.join(x.split('tmalign_')[1].split('_')[1:]).replace('.txt', ''))

pivot = scores.pivot(index='id1', columns='id2', values='tm_score').fillna(0)
pca = PCA(n_components=2)
pca_result = pca.fit_transform(pivot)

plt.figure(figsize=(10, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1], c='blue', label=f'Explained variance: {sum(pca.explained_variance_ratio_):.2f}')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA of TM-align Scores')
plt.legend()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
