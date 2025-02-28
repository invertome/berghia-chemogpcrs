#!/usr/bin/env python3
# plot_selective_pressure.py
# Purpose: Visualize selective pressure results.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

lrt_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(lrt_file, names=['base', 'chi2_stat', 'p_value'])
if df.empty:
    print(f"Warning: No data in {lrt_file}", file=sys.stderr)
    plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, "No Selective Pressure Data", ha='center', va='center', fontsize=14)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    sys.exit(0)

plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='chi2_stat', y=-df['p_value'].apply(lambda x: pd.np.log10(x) if x > 0 else 0), hue='base', legend=False)
plt.xlabel('Chi-squared Statistic')
plt.ylabel('-log10(p-value)')
plt.title('Selective Pressure Analysis (LRT)')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
