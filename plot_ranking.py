#!/usr/bin/env python3
# plot_ranking.py
# Purpose: Visualize ranked candidates.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

ranking_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(ranking_file)
plt.figure(figsize=(10, 6))
sns.barplot(data=df.head(10), x='id', y='rank_score')
plt.xticks(rotation=45, ha='right')
plt.xlabel('Candidate ID')
plt.ylabel('Rank Score')
plt.title('Top 10 Ranked GPCR Candidates')
plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
