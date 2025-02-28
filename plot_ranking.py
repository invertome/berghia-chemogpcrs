#!/usr/bin/env python3
# plot_ranking.py
# Purpose: Visualize top-ranked GPCR candidates.
# Inputs: Ranked candidates CSV ($1), output plot file ($2)
# Outputs: Bar plot of top candidates (${output_file}.png, ${output_file}.svg)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Command-line arguments
ranking_file = sys.argv[1]  # Input ranked candidates CSV
output_file = sys.argv[2]   # Output plot file (without extension)

# Load ranked data
df = pd.read_csv(ranking_file)

# Plot top 10 candidates
plt.figure(figsize=(10, 6))
sns.barplot(data=df.head(10), x='id', y='rank_score')
plt.xticks(rotation=45, ha='right')
plt.xlabel('Candidate ID')
plt.ylabel('Rank Score')
plt.title('Top 10 Ranked GPCR Candidates')
plt.tight_layout()

# Save in both PNG and SVG formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
plt.close()
