#!/usr/bin/env python3
# plot_heatmap.py
# Purpose: Generate heatmap of TM-align scores.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

tmalign_scores_file = sys.argv[1]
output_file = sys.argv[2]

scores = pd.read_csv(tmalign_scores_file, names=['file', 'tm_score'], sep=',')
scores['id1'] = scores['file'].apply(lambda x: x.split('tmalign_')[1].split('_')[0])
scores['id2'] = scores['file'].apply(lambda x: '_'.join(x.split('tmalign_')[1].split('_')[1:]).replace('.txt', ''))

pivot = scores.pivot(index='id1', columns='id2', values='tm_score').fillna(0)
plt.figure(figsize=(12, 10))
sns.heatmap(pivot, cmap='viridis', annot=True, fmt='.2f')
plt.title('TM-align Score Heatmap')
plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
