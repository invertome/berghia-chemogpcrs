#!/usr/bin/env python3
# plot_heatmap.py
# Purpose: Generate hierarchically-clustered heatmap of TM-align structural similarity.
# Inputs: TM-align scores CSV ($1), output plot file ($2)
# Outputs: Clustered heatmap (${output_file}.png, ${output_file}.svg, ${output_file}.pdf)
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import sys
import re

# Set publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 1.0

# Color scheme for sequence types
TYPE_COLORS = {
    'Reference': '#1f77b4',
    'Berghia': '#2ca02c',
    'Other': '#ff7f0e'
}

# Command-line arguments
tmalign_scores_file = sys.argv[1]
output_file = sys.argv[2]

# Load and process TM-align scores
try:
    scores = pd.read_csv(tmalign_scores_file, names=['file', 'tm_score'], sep=',')
except Exception as e:
    print(f"Error loading {tmalign_scores_file}: {e}", file=sys.stderr)
    plt.figure(figsize=(10, 8))
    plt.text(0.5, 0.5, "No TM-align Data", ha='center', va='center', fontsize=14)
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
    plt.close()
    sys.exit(0)

if scores.empty:
    print(f"Warning: No data in {tmalign_scores_file}", file=sys.stderr)
    plt.figure(figsize=(10, 8))
    plt.text(0.5, 0.5, "No TM-align Data", ha='center', va='center', fontsize=14)
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight')
    plt.close()
    sys.exit(0)

# Parse IDs from filenames
def parse_ids(filename):
    """Extract pair of IDs from TM-align filename."""
    match = re.search(r'tmalign_(.+?)_(.+?)\.txt', filename)
    if match:
        return match.group(1), match.group(2)
    parts = filename.replace('.txt', '').split('_')
    if len(parts) >= 3:
        return parts[1], '_'.join(parts[2:])
    return None, None

scores['id1'], scores['id2'] = zip(*scores['file'].apply(parse_ids))
scores = scores.dropna(subset=['id1', 'id2'])

if scores.empty:
    print("Warning: Could not parse any IDs from filenames", file=sys.stderr)
    plt.figure(figsize=(10, 8))
    plt.text(0.5, 0.5, "Could not parse TM-align data", ha='center', va='center', fontsize=14)
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight')
    plt.close()
    sys.exit(0)

# Create similarity matrix
pivot = scores.pivot(index='id1', columns='id2', values='tm_score').fillna(0)

# Make symmetric
all_ids = sorted(set(pivot.index) | set(pivot.columns))
full_matrix = pd.DataFrame(0.0, index=all_ids, columns=all_ids)
for idx in pivot.index:
    for col in pivot.columns:
        val = pivot.loc[idx, col]
        full_matrix.loc[idx, col] = val
        full_matrix.loc[col, idx] = val

# Fill diagonal with 1.0
np.fill_diagonal(full_matrix.values, 1.0)

# Determine sequence types for color annotation
def get_sequence_type(seq_id):
    """Categorize sequence by ID prefix."""
    if seq_id.startswith('ref_'):
        return 'Reference'
    elif 'berghia' in seq_id.lower():
        return 'Berghia'
    else:
        return 'Other'

seq_types = pd.Series({seq_id: get_sequence_type(seq_id) for seq_id in all_ids})
row_colors = seq_types.map(TYPE_COLORS)

# Create clustermap with hierarchical clustering
# Use distance (1 - similarity) for clustering
distance_matrix = 1 - full_matrix

# Ensure distance matrix is valid for clustering
distance_matrix = distance_matrix.clip(lower=0)
np.fill_diagonal(distance_matrix.values, 0)

# Check matrix size for appropriate settings
n_samples = len(full_matrix)
show_labels = n_samples <= 50
annot = n_samples <= 20

# Compute linkage
try:
    # Convert to condensed form for linkage
    condensed = squareform(distance_matrix.values, checks=False)
    linkage_matrix = linkage(condensed, method='average')
except Exception as e:
    print(f"Warning: Clustering failed, using simple heatmap: {e}", file=sys.stderr)
    linkage_matrix = None

# Create figure
if linkage_matrix is not None:
    # Clustered heatmap with dendrograms
    g = sns.clustermap(
        full_matrix,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        row_colors=row_colors,
        col_colors=row_colors,
        cmap='YlOrRd',
        vmin=0, vmax=1,
        xticklabels=show_labels,
        yticklabels=show_labels,
        annot=annot,
        fmt='.2f',
        annot_kws={'size': 7} if annot else {},
        figsize=(14, 12),
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        cbar_kws={'label': 'TM-score'}
    )

    # Rotate labels for readability
    if show_labels:
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)

    # Add legend for sequence types
    legend_patches = [mpatches.Patch(color=color, label=label)
                      for label, color in TYPE_COLORS.items()
                      if label in seq_types.values]
    g.ax_heatmap.legend(handles=legend_patches, loc='upper left',
                        bbox_to_anchor=(1.15, 1), framealpha=0.9, fontsize=9)

    # Add title
    g.fig.suptitle('Structural Similarity Heatmap (TM-align)', fontsize=14, fontweight='bold', y=1.02)

    # Add interpretation guide
    guide_text = "TM-score: >0.5 = same fold, >0.17 = similar topology"
    g.fig.text(0.5, -0.02, guide_text, ha='center', fontsize=10, style='italic')

    fig = g.fig

else:
    # Simple heatmap without clustering
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(full_matrix, cmap='YlOrRd', vmin=0, vmax=1,
                xticklabels=show_labels, yticklabels=show_labels,
                annot=annot, fmt='.2f', ax=ax,
                cbar_kws={'label': 'TM-score'})

    if show_labels:
        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=8)

    ax.set_title('Structural Similarity Heatmap (TM-align)', fontsize=14, fontweight='bold')

# Save in multiple formats
plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.svg", format='svg', bbox_inches='tight', facecolor='white')
plt.savefig(f"{output_file}.pdf", format='pdf', bbox_inches='tight', facecolor='white')
plt.close()

# --- Generate separate dendrogram figure ---
if linkage_matrix is not None and n_samples <= 100:
    fig_dend, ax_dend = plt.subplots(figsize=(12, 6))

    # Color function based on sequence type
    def get_leaf_color(label):
        return TYPE_COLORS.get(get_sequence_type(label), '#7f7f7f')

    # Plot dendrogram
    dendro = dendrogram(
        linkage_matrix,
        labels=all_ids,
        ax=ax_dend,
        leaf_rotation=90,
        leaf_font_size=8 if n_samples <= 30 else 6,
        color_threshold=0.5  # Clusters with distance < 0.5 (TM-score > 0.5)
    )

    ax_dend.set_ylabel('Distance (1 - TM-score)', fontsize=11)
    ax_dend.set_title('Hierarchical Clustering of Structural Similarity', fontsize=12, fontweight='bold')

    # Add threshold line
    ax_dend.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.7,
                    label='Same fold threshold (TM=0.5)')
    ax_dend.legend(loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig(f"{output_file}_dendrogram.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{output_file}_dendrogram.pdf", format='pdf', bbox_inches='tight', facecolor='white')
    plt.close()

# --- Generate Summary Statistics ---
summary_file = f"{output_file}_summary.txt"
with open(summary_file, 'w') as f:
    f.write("Structural Similarity Heatmap Summary\n")
    f.write("=" * 40 + "\n\n")

    f.write(f"Matrix dimensions: {n_samples} x {n_samples}\n\n")

    f.write("Sequence types:\n")
    for stype in ['Reference', 'Berghia', 'Other']:
        count = (seq_types == stype).sum()
        if count > 0:
            f.write(f"  {stype}: {count}\n")

    # TM-score statistics (excluding diagonal)
    tm_vals = full_matrix.values[np.triu_indices_from(full_matrix.values, k=1)]

    f.write(f"\nTM-score Statistics:\n")
    f.write(f"  Mean: {np.mean(tm_vals):.4f}\n")
    f.write(f"  Median: {np.median(tm_vals):.4f}\n")
    f.write(f"  Std: {np.std(tm_vals):.4f}\n")
    f.write(f"  Min: {np.min(tm_vals):.4f}\n")
    f.write(f"  Max: {np.max(tm_vals):.4f}\n")

    # Cluster analysis
    same_fold = np.sum(tm_vals > 0.5)
    total_pairs = len(tm_vals)
    f.write(f"\nStructural Relationships:\n")
    f.write(f"  Same fold (TM > 0.5): {same_fold} pairs ({same_fold/total_pairs*100:.1f}%)\n")
    f.write(f"  Similar topology (TM > 0.17): {np.sum(tm_vals > 0.17)} pairs\n")

print(f"Generated heatmap visualizations:")
print(f"  PNG: {output_file}.png")
print(f"  SVG: {output_file}.svg")
print(f"  PDF: {output_file}.pdf")
if linkage_matrix is not None and n_samples <= 100:
    print(f"  Dendrogram: {output_file}_dendrogram.png")
print(f"  Summary: {summary_file}")
