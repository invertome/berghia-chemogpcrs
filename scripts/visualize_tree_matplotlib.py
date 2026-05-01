#!/usr/bin/env python3
"""
visualize_tree_matplotlib.py
Publication-quality GPCR phylogenetic tree visualizations using matplotlib.
Avoids ete3 rendering (problematic on headless systems).

Produces:
  1. Circular/radial tree overview (colored by Berghia vs reference)
  2. Bootstrap support histogram from contree

Author: Jorge L. Perez-Moreno, Ph.D.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # headless backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from Bio import Phylo
from collections import Counter


# ---------------------------------------------------------------------------
# Utility: compute radial coordinates for a phylogenetic tree
# ---------------------------------------------------------------------------

def compute_circular_coords(tree):
    """
    Compute (x, y) coordinates for every clade in a circular/radial layout.

    Returns:
        coords: dict mapping clade -> (x, y)
        leaf_angles: dict mapping terminal clade -> angle (radians)
        max_depth: float, maximum root-to-tip distance
    """
    terminals = tree.get_terminals()
    n_leaves = len(terminals)

    # Assign each leaf an evenly-spaced angle
    leaf_angles = {}
    for i, leaf in enumerate(terminals):
        leaf_angles[leaf] = 2.0 * np.pi * i / n_leaves

    # For internal nodes, angle = mean of descendant leaf angles
    # We also need the depth (distance from root) for the radial distance
    depths = tree.depths(unit_branch_lengths=False)

    # Compute angles bottom-up
    clade_angles = {}
    for leaf in terminals:
        clade_angles[leaf] = leaf_angles[leaf]

    def get_angle(clade):
        if clade in clade_angles:
            return clade_angles[clade]
        child_angles = [get_angle(c) for c in clade.clades]
        # Use circular mean to handle wrap-around
        sin_mean = np.mean([np.sin(a) for a in child_angles])
        cos_mean = np.mean([np.cos(a) for a in child_angles])
        angle = np.arctan2(sin_mean, cos_mean)
        clade_angles[clade] = angle
        return angle

    get_angle(tree.root)

    # Convert (depth, angle) -> (x, y)
    max_depth = max(depths.values()) if depths.values() else 1.0
    coords = {}
    for clade, depth in depths.items():
        angle = clade_angles.get(clade, 0)
        r = depth / max_depth  # normalize to [0, 1]
        coords[clade] = (r * np.cos(angle), r * np.sin(angle))

    return coords, clade_angles, max_depth


def draw_circular_tree(tree, ax, coords, clade_angles, max_depth,
                       berghia_color='#d62728', ref_color='#b0b0b0',
                       line_color='#555555', linewidth=0.15):
    """
    Draw the tree on a matplotlib Axes in circular layout.
    Branches are drawn as straight radial lines + arcs.
    """
    depths = tree.depths(unit_branch_lengths=False)

    # Collect line segments for batch drawing (much faster than individual lines)
    branch_segments = []
    branch_colors = []

    for clade in tree.find_clades(order='level'):
        if clade.clades:
            parent_depth = depths[clade] / max_depth
            parent_angle = clade_angles[clade]

            for child in clade.clades:
                child_depth = depths[child] / max_depth
                child_angle = clade_angles[child]

                # Radial line: from (parent_depth, child_angle) to (child_depth, child_angle)
                x1 = parent_depth * np.cos(child_angle)
                y1 = parent_depth * np.sin(child_angle)
                x2 = child_depth * np.cos(child_angle)
                y2 = child_depth * np.sin(child_angle)
                branch_segments.append([(x1, y1), (x2, y2)])
                branch_colors.append(line_color)

                # Arc at parent_depth from parent_angle toward child_angle
                # Use small arc segments
                angle_diff = child_angle - parent_angle
                # Normalize to [-pi, pi]
                while angle_diff > np.pi:
                    angle_diff -= 2 * np.pi
                while angle_diff < -np.pi:
                    angle_diff += 2 * np.pi

                n_arc_pts = max(2, int(abs(angle_diff) / (np.pi / 180)))  # ~1 degree steps
                arc_angles = np.linspace(parent_angle, parent_angle + angle_diff, n_arc_pts)
                for j in range(len(arc_angles) - 1):
                    ax1 = parent_depth * np.cos(arc_angles[j])
                    ay1 = parent_depth * np.sin(arc_angles[j])
                    ax2 = parent_depth * np.cos(arc_angles[j + 1])
                    ay2 = parent_depth * np.sin(arc_angles[j + 1])
                    branch_segments.append([(ax1, ay1), (ax2, ay2)])
                    branch_colors.append(line_color)

    # Draw all branches at once
    lc = LineCollection(branch_segments, colors=branch_colors,
                        linewidths=linewidth, antialiaseds=True)
    ax.add_collection(lc)

    # Draw leaf tips as colored dots
    berghia_x, berghia_y = [], []
    ref_x, ref_y = [], []

    for leaf in tree.get_terminals():
        x, y = coords[leaf]
        if leaf.name and leaf.name.startswith('Berste'):
            berghia_x.append(x)
            berghia_y.append(y)
        else:
            ref_x.append(x)
            ref_y.append(y)

    # Plot reference tips first (underneath), then Berghia on top
    ax.scatter(ref_x, ref_y, s=2.0, c=ref_color, alpha=0.6, linewidths=0,
               zorder=2, rasterized=True)
    ax.scatter(berghia_x, berghia_y, s=8.0, c=berghia_color, alpha=0.9,
               linewidths=0, zorder=3, rasterized=True)


def plot_circular_tree(treefile, outpath):
    """Generate publication-quality circular tree plot."""
    print(f"Reading tree from {treefile} ...")
    tree = Phylo.read(treefile, 'newick')

    n_leaves = len(tree.get_terminals())
    n_berghia = sum(1 for l in tree.get_terminals()
                    if l.name and l.name.startswith('Berste'))
    n_ref = n_leaves - n_berghia

    print(f"  {n_leaves} leaves ({n_berghia} Berghia, {n_ref} reference)")
    print("Computing circular layout ...")
    coords, clade_angles, max_depth = compute_circular_coords(tree)

    print("Drawing ...")
    fig, ax = plt.subplots(figsize=(20, 20), dpi=300)
    ax.set_aspect('equal')
    ax.axis('off')

    draw_circular_tree(tree, ax, coords, clade_angles, max_depth,
                       berghia_color='#d62728',   # red
                       ref_color='#9ec5e8',        # light blue
                       line_color='#333333',
                       linewidth=0.25)

    # Legend
    legend_patches = [
        mpatches.Patch(color='#d62728',
                       label=f'Berghia stephanieae (n={n_berghia})'),
        mpatches.Patch(color='#9ec5e8',
                       label=f'Reference GPCRs (n={n_ref})'),
    ]
    ax.legend(handles=legend_patches, loc='lower left', fontsize=14,
              frameon=True, fancybox=True, framealpha=0.9,
              edgecolor='#cccccc', borderpad=1.0)

    # Title
    ax.set_title(
        'GPCR Phylogeny of Berghia stephanieae\n'
        f'VT+I+R10 model | {n_leaves:,} sequences | 1,000 UFBoot replicates',
        fontsize=18, fontweight='bold', pad=20
    )

    # Tight layout with padding
    margin = 0.08
    ax.set_xlim(-1 - margin, 1 + margin)
    ax.set_ylim(-1 - margin, 1 + margin)

    fig.tight_layout(pad=2.0)
    fig.savefig(outpath, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_bootstrap_histogram(contree_file, outpath):
    """Generate histogram of UFBoot support values from contree."""
    print(f"Reading contree from {contree_file} ...")
    tree = Phylo.read(contree_file, 'newick')

    # Extract bootstrap support values from internal nodes
    supports = []
    for clade in tree.get_nonterminals():
        if clade.confidence is not None:
            supports.append(float(clade.confidence))

    n_total = len(supports)
    print(f"  {n_total} internal nodes with support values")

    if n_total == 0:
        print("  ERROR: No support values found. Aborting histogram.")
        return

    supports = np.array(supports)

    # Statistics
    n_high = np.sum(supports >= 95)
    n_moderate = np.sum((supports >= 70) & (supports < 95))
    n_low = np.sum(supports < 70)
    median_val = np.median(supports)
    mean_val = np.mean(supports)

    print(f"  High (>=95): {n_high} ({100*n_high/n_total:.1f}%)")
    print(f"  Moderate (70-94): {n_moderate} ({100*n_moderate/n_total:.1f}%)")
    print(f"  Low (<70): {n_low} ({100*n_low/n_total:.1f}%)")
    print(f"  Median: {median_val:.1f}, Mean: {mean_val:.1f}")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # Bins every 5 units from 0 to 100
    bin_edges = np.arange(0, 105, 5)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    counts, _ = np.histogram(supports, bins=bin_edges)

    # Color each bar based on its bin center
    bar_colors = []
    for center in bin_centers:
        if center >= 95:
            bar_colors.append('#2ca02c')   # green
        elif center >= 70:
            bar_colors.append('#f0c929')   # yellow/gold
        else:
            bar_colors.append('#d62728')   # red

    bars = ax.bar(bin_centers, counts, width=4.5, color=bar_colors,
                  edgecolor='white', linewidth=0.5, zorder=3)

    # Threshold lines
    ax.axvline(x=70, color='#555555', linestyle='--', linewidth=1.5,
               zorder=4, label='Moderate threshold (70)')
    ax.axvline(x=95, color='#222222', linestyle='--', linewidth=1.5,
               zorder=4, label='Strong threshold (95)')

    # Labels and title
    ax.set_xlabel('UFBoot Support Value', fontsize=14, fontweight='bold')
    ax.set_ylabel('Number of Internal Nodes', fontsize=14, fontweight='bold')
    ax.set_title(
        'Distribution of Bootstrap Support Values\n'
        f'Consensus tree | {n_total:,} internal nodes | '
        f'Median={median_val:.0f} | Mean={mean_val:.1f}',
        fontsize=15, fontweight='bold'
    )

    # Custom legend
    legend_patches = [
        mpatches.Patch(color='#2ca02c',
                       label=f'Strong (>=95): {n_high} ({100*n_high/n_total:.1f}%)'),
        mpatches.Patch(color='#f0c929',
                       label=f'Moderate (70-94): {n_moderate} ({100*n_moderate/n_total:.1f}%)'),
        mpatches.Patch(color='#d62728',
                       label=f'Weak (<70): {n_low} ({100*n_low/n_total:.1f}%)'),
    ]
    ax.legend(handles=legend_patches, fontsize=11, loc='upper left',
              frameon=True, fancybox=True, framealpha=0.9)

    ax.set_xlim(-2, 102)
    ax.set_xticks(np.arange(0, 105, 10))
    ax.tick_params(axis='both', labelsize=11)

    # Light grid
    ax.yaxis.grid(True, linestyle=':', alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout(pad=1.5)
    fig.savefig(outpath, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    base_dir = '/home/workspace/Desktop/projects/umass/berghia-chemogpcrs/preliminary/results/phylogenies/protein/v2'

    # Allow overriding via CLI args
    if len(sys.argv) >= 2:
        base_dir = sys.argv[1]

    treefile = os.path.join(base_dir, 'gpcrs.treefile')
    contree = os.path.join(base_dir, 'gpcrs.contree')
    circular_out = os.path.join(base_dir, 'gpcrs_tree_circular.png')
    hist_out = os.path.join(base_dir, 'gpcrs_bootstrap_histogram.png')

    # 1. Circular tree
    if os.path.isfile(treefile):
        plot_circular_tree(treefile, circular_out)
    else:
        print(f"WARNING: Tree file not found: {treefile}")

    # 2. Bootstrap histogram
    if os.path.isfile(contree):
        plot_bootstrap_histogram(contree, hist_out)
    else:
        print(f"WARNING: Contree file not found: {contree}")

    print("\nDone.")
