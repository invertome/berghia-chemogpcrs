#!/usr/bin/env python3
# plot_pca.py
# Purpose: Generate PCA visualization of structural similarity from TM-align scores.
# Inputs: TM-align scores CSV ($1), output plot file ($2), [ID mapping file ($3)]
# Outputs: Multi-panel PCA figure (${output_file}.png, ${output_file}.svg, ${output_file}.pdf)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster
import seaborn as sns
import sys
import os
import re

# Set publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

# Color scheme for sequence types
TYPE_COLORS = {
    'Reference': '#1f77b4',
    'Berghia': '#2ca02c',
    'Other': '#ff7f0e'
}

# Command-line arguments
tmalign_scores_file = sys.argv[1]
output_file = sys.argv[2]
id_map_file = sys.argv[3] if len(sys.argv) > 3 else None

# Load and process TM-align scores
try:
    scores = pd.read_csv(tmalign_scores_file, names=['file', 'tm_score'], sep=',')
except Exception as e:
    print(f"Error loading {tmalign_scores_file}: {e}", file=sys.stderr)
    # Create empty plot
    plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, "No TM-align Data", ha='center', va='center', fontsize=14)
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
    plt.close()
    sys.exit(0)

if scores.empty:
    print(f"Warning: No data in {tmalign_scores_file}", file=sys.stderr)
    plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, "No TM-align Data", ha='center', va='center', fontsize=14)
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
    plt.close()
    sys.exit(0)

# Parse IDs from filenames
def parse_ids(filename):
    """Extract pair of IDs from TM-align filename."""
    # Pattern: tmalign_ID1_ID2.txt or similar
    match = re.search(r'tmalign_(.+?)_(.+?)\.txt', filename)
    if match:
        return match.group(1), match.group(2)
    # Fallback pattern
    parts = filename.replace('.txt', '').split('_')
    if len(parts) >= 3:
        return parts[1], '_'.join(parts[2:])
    return None, None

scores['id1'], scores['id2'] = zip(*scores['file'].apply(parse_ids))
scores = scores.dropna(subset=['id1', 'id2'])

if scores.empty:
    print("Warning: Could not parse any IDs from filenames", file=sys.stderr)
    plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, "Could not parse TM-align data", ha='center', va='center', fontsize=14)
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
    plt.close()
    sys.exit(0)

# Create similarity matrix and convert to distance
pivot = scores.pivot(index='id1', columns='id2', values='tm_score').fillna(0)

# Make symmetric if not already
all_ids = list(set(pivot.index) | set(pivot.columns))
full_matrix = pd.DataFrame(0.0, index=all_ids, columns=all_ids)
for idx in pivot.index:
    for col in pivot.columns:
        val = pivot.loc[idx, col]
        full_matrix.loc[idx, col] = val
        full_matrix.loc[col, idx] = val

# Fill diagonal with 1.0 (perfect self-similarity)
np.fill_diagonal(full_matrix.values, 1.0)

# Convert similarity to distance
distance_matrix = 1 - full_matrix

# Determine sequence types
def get_sequence_type(seq_id):
    """Categorize sequence by ID prefix."""
    if seq_id.startswith('ref_'):
        return 'Reference'
    elif 'berghia' in seq_id.lower():
        return 'Berghia'
    else:
        return 'Other'

seq_types = {seq_id: get_sequence_type(seq_id) for seq_id in all_ids}

# Create figure with multiple panels
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.25)

# --- Panel A: PCA Plot ---
ax1 = fig.add_subplot(gs[0, 0])

# Check if we have enough data for 2D PCA
n_samples = len(distance_matrix)
if n_samples < 2:
    ax1.text(0.5, 0.5, f"Insufficient data for PCA\n(n={n_samples})",
             ha='center', va='center', fontsize=12, transform=ax1.transAxes)
    ax1.set_title('A. PCA of Structural Similarity', fontsize=12, fontweight='bold', loc='left')
    pca_result = None
    var_explained = [0, 0]
else:
    # Perform PCA on distance matrix
    n_components = min(2, n_samples)
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(distance_matrix)

    # Handle case where we only have 1 component
    if pca_result.shape[1] == 1:
        # Add a zero column for the second dimension
        pca_result = np.column_stack([pca_result, np.zeros(pca_result.shape[0])])

    # Plot by sequence type
    for seq_type, color in TYPE_COLORS.items():
        mask = [seq_types[seq_id] == seq_type for seq_id in distance_matrix.index]
        if any(mask):
            indices = [i for i, m in enumerate(mask) if m]
            ax1.scatter(pca_result[indices, 0], pca_result[indices, 1],
                       c=color, label=f'{seq_type} (n={sum(mask)})',
                       s=60, alpha=0.7, edgecolors='white', linewidths=0.5)

    # Add explained variance
    var_explained = list(pca.explained_variance_ratio_)
    if len(var_explained) < 2:
        var_explained.append(0)

ax1.set_xlabel(f'PC1 ({var_explained[0]*100:.1f}% variance)', fontsize=11)
ax1.set_ylabel(f'PC2 ({var_explained[1]*100:.1f}% variance)' if var_explained[1] > 0 else 'PC2', fontsize=11)
ax1.set_title('A. PCA of Structural Similarity', fontsize=12, fontweight='bold', loc='left')
ax1.legend(loc='best', fontsize=9, framealpha=0.9)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Add total variance explained
total_var = sum(var_explained) * 100
ax1.text(0.05, 0.95, f'Total variance: {total_var:.1f}%', transform=ax1.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# --- Panel B: PCA with K-means Clustering ---
ax2 = fig.add_subplot(gs[0, 1])

# Determine optimal number of clusters (max 5)
n_clusters = min(5, len(distance_matrix) // 3)
if n_clusters >= 2:
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(distance_matrix)

    # Plot clusters
    cluster_colors = plt.cm.Set1(np.linspace(0, 1, n_clusters))
    for cluster_id in range(n_clusters):
        mask = cluster_labels == cluster_id
        ax2.scatter(pca_result[mask, 0], pca_result[mask, 1],
                   c=[cluster_colors[cluster_id]], label=f'Cluster {cluster_id+1}',
                   s=60, alpha=0.7, edgecolors='white', linewidths=0.5)

    ax2.legend(loc='best', fontsize=9, framealpha=0.9)
else:
    ax2.scatter(pca_result[:, 0], pca_result[:, 1], c='steelblue',
               s=60, alpha=0.7, edgecolors='white', linewidths=0.5)

ax2.set_xlabel(f'PC1 ({var_explained[0]*100:.1f}% variance)', fontsize=11)
ax2.set_ylabel(f'PC2 ({var_explained[1]*100:.1f}% variance)' if len(var_explained) > 1 else 'PC2', fontsize=11)
ax2.set_title('B. K-means Structural Clusters', fontsize=12, fontweight='bold', loc='left')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# --- Panel C: TM-score Distribution ---
ax3 = fig.add_subplot(gs[1, 0])

# Get all TM-scores (excluding diagonal)
tm_vals = full_matrix.values[np.triu_indices_from(full_matrix.values, k=1)]
tm_vals = tm_vals[tm_vals > 0]  # Remove zeros

if len(tm_vals) > 0:
    ax3.hist(tm_vals, bins=30, color='steelblue', edgecolor='white', alpha=0.7)

    # Add interpretation lines
    ax3.axvline(x=0.5, color='orange', linestyle='--', linewidth=1.5,
               label='Same fold (TM > 0.5)')
    ax3.axvline(x=0.17, color='red', linestyle=':', linewidth=1.5,
               label='Random similarity (TM ~ 0.17)')

    # Add statistics
    mean_tm = np.mean(tm_vals)
    median_tm = np.median(tm_vals)
    ax3.axvline(x=mean_tm, color='green', linestyle='-', linewidth=1.5,
               label=f'Mean = {mean_tm:.3f}')

    ax3.legend(loc='upper right', fontsize=9, framealpha=0.9)

ax3.set_xlabel('TM-score', fontsize=11)
ax3.set_ylabel('Count', fontsize=11)
ax3.set_title('C. TM-score Distribution', fontsize=12, fontweight='bold', loc='left')
ax3.set_xlim(0, 1)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

# --- Panel D: Scree Plot ---
ax4 = fig.add_subplot(gs[1, 1])

# Perform full PCA for scree plot
n_components = min(10, len(distance_matrix) - 1)
if n_components > 0:
    pca_full = PCA(n_components=n_components)
    pca_full.fit(distance_matrix)

    var_ratio = pca_full.explained_variance_ratio_
    cum_var = np.cumsum(var_ratio)

    x = np.arange(1, len(var_ratio) + 1)
    ax4.bar(x, var_ratio * 100, color='steelblue', alpha=0.7, label='Individual')
    ax4.plot(x, cum_var * 100, 'ro-', markersize=6, linewidth=2, label='Cumulative')

    # Add 80% threshold line
    ax4.axhline(y=80, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax4.text(len(var_ratio), 81, '80%', fontsize=9, ha='right')

    ax4.legend(loc='center right', fontsize=9, framealpha=0.9)

ax4.set_xlabel('Principal Component', fontsize=11)
ax4.set_ylabel('Variance Explained (%)', fontsize=11)
ax4.set_title('D. Scree Plot', fontsize=12, fontweight='bold', loc='left')
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)

# Add overall title
fig.suptitle('Structural Similarity Analysis (TM-align)', fontsize=14, fontweight='bold', y=0.98)

# Save in multiple formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.pdf", format='pdf', bbox_inches='tight', facecolor='white')
plt.close()

# --- Generate Summary Statistics ---
summary_file = f"{output_file}_summary.txt"
with open(summary_file, 'w') as f:
    f.write("Structural Similarity Analysis Summary\n")
    f.write("=" * 40 + "\n\n")

    f.write(f"Structures analyzed: {len(all_ids)}\n")
    f.write(f"  References: {sum(1 for t in seq_types.values() if t == 'Reference')}\n")
    f.write(f"  Berghia: {sum(1 for t in seq_types.values() if t == 'Berghia')}\n")
    f.write(f"  Other: {sum(1 for t in seq_types.values() if t == 'Other')}\n\n")

    f.write(f"TM-score Statistics:\n")
    if len(tm_vals) > 0:
        f.write(f"  Mean: {np.mean(tm_vals):.4f}\n")
        f.write(f"  Median: {np.median(tm_vals):.4f}\n")
        f.write(f"  Min: {np.min(tm_vals):.4f}\n")
        f.write(f"  Max: {np.max(tm_vals):.4f}\n")
        f.write(f"  Pairs with same fold (TM > 0.5): {sum(tm_vals > 0.5)} ({sum(tm_vals > 0.5)/len(tm_vals)*100:.1f}%)\n")

    f.write(f"\nPCA Results:\n")
    f.write(f"  PC1 variance: {var_explained[0]*100:.2f}%\n")
    if len(var_explained) > 1:
        f.write(f"  PC2 variance: {var_explained[1]*100:.2f}%\n")
    f.write(f"  Total (2 PCs): {sum(var_explained[:2])*100:.2f}%\n")

print(f"Generated PCA visualizations:")
print(f"  PNG: {output_file}.png")
print(f"  SVG: {output_file}.svg")
print(f"  PDF: {output_file}.pdf")
print(f"  Summary: {summary_file}")
