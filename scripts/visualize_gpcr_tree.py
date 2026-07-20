#!/usr/bin/env python3
"""
Visualize GPCR phylogenetic tree with:
  - Midpoint rooting
  - Ladderized descending (largest clades first)
  - Outer ring colored by taxonomic group
  - Candidates highlighted by data-driven tiers (Q90 + natural break)
  - Root marker
  - Circular layout with branch length compression

Usage:
  python3 visualize_gpcr_tree.py \
    --tree results/phylogenies/protein/v2/gpcrs.treefile \
    --ranking results/ranking/ranked_candidates_sorted.csv \
    --output results/phylogenies/protein/v2/gpcrs_tree_taxonomy.png
"""

import argparse
import csv
import math
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from Bio import Phylo

# Taxonomic group -> color
TAXON_COLORS = {
    'Gastropoda':       '#e41a1c',
    'Bivalvia':         '#377eb8',
    'Cephalopoda':      '#4daf4a',
    'Annelida':         '#ff7f00',
    'Platyhelminthes':  '#984ea3',
    'Other Mollusca':   '#a65628',
    'Other Lophotrochozoa': '#f781bf',
    'Berghia':          '#2d2d2d',  # overridden to white in dark theme
    'Tier 1':           '#FFD700',  # gold — Q90
    'Tier 2':           '#C0C0C0',  # silver — above natural break
    'Unknown':          '#cccccc',
}

# Species prefix -> taxonomic group
PREFIX_TO_GROUP = {}
_groups = {
    'Gastropoda': [
        'alvmar', 'amba', 'anvo', 'aplcal', 'arvu', 'baar', 'baare', 'bigl', 'bipf',
        'bist', 'caun', 'chori', 'chsq', 'cobe', 'coco', 'cotr', 'dipe', 'drsu',
        'elch', 'elma', 'giae', 'gima', 'goco', 'gofi', 'goge', 'goku', 'gole',
        'hacra', 'hadi', 'hala', 'haru', 'logi', 'losc', 'lyst', 'metu', 'orid',
        'pade', 'pape', 'pavu', 'phau', 'phli', 'ploc', 'poca', 'raau', 'rave', 'trte',
    ],
    'Bivalvia': [
        'arma', 'arpu', 'atja', 'bapl', 'bofu', 'ceed', 'chfa', 'cobi', 'cofl', 'coku',
        'craan', 'craar', 'crahon', 'crgi', 'crpl', 'crvi', 'drpo', 'drro', 'frfr',
        'frwh', 'gate', 'hihi', 'hycu', 'lael', 'lifo', 'lian2', 'lini', 'maho',
        'mama', 'maqu', 'memer', 'miye', 'moph', 'myar', 'myca', 'myco', 'myed',
        'myga', 'myvi', 'osed', 'oslu', 'page', 'paye', 'pema', 'pevi', 'pifu',
        'pifuma', 'piim', 'post', 'ruph', 'sagl', 'sapu', 'scbr', 'sico', 'spso',
        'tegr', 'trcr', 'trgi', 'unde', 'unpi',
    ],
    'Cephalopoda': [
        'arar', 'ardu', 'eusc', 'hama', 'jadi', 'mule', 'mulo', 'napo', 'ocbi',
        'ocde', 'ocin', 'ocmi', 'ocsi', 'octmim', 'octmy', 'octrub', 'ocvu', 'seph',
        'watsci',
    ],
    'Annelida': [
        'acsq', 'alge', 'alvi', 'amco', 'ampa', 'apocal', 'cate', 'digy', 'eueu',
        'haimp', 'hero', 'himan', 'hime', 'hive', 'hyel', 'lalu', 'lasa', 'lecl',
        'lute', 'mevu', 'owfu', 'paec', 'papa', 'pldu', 'poma', 'ripis', 'ripa',
        'sinu', 'stbe', 'stli', 'tela', 'whpi',
    ],
    'Platyhelminthes': [
        'atwi', 'cafo', 'clsi', 'dican', 'diden', 'dila', 'duja', 'eccap', 'eccan',
        'ecgr', 'ecmu', 'euni', 'fabu', 'fagi', 'fahe', 'gitig', 'gybu', 'gysa',
        'heame', 'hydi', 'hyta', 'macl', 'mahy', 'mali', 'opfe', 'opvi', 'pahe',
        'pask', 'pawe', 'prcro', 'prfa', 'rona', 'scbo', 'sccu', 'scgu', 'scha',
        'scin', 'scja', 'scman', 'scmar', 'scmat', 'scme', 'scro', 'scso', 'sctur',
        'sper', 'sppr', 'tacr', 'tamu', 'tasa', 'trre',
    ],
    'Other Mollusca': [
        'accr', 'acgra', 'epba', 'gato', 'gype', 'laant', 'mosw', 'move', 'wiar',
    ],
    'Other Lophotrochozoa': [
        'bune', 'crmu', 'crpa', 'lian', 'liun', 'lilo', 'meme', 'noge', 'pece',
        # 'phaust' = Phoronis australis. Was 'phau2', which matched no
        # sequence anywhere, so phoronid leaves fell through to the
        # gastropod 'phau' entry and were coloured as molluscs. 'phaust'
        # is the rename target documented in
        # recover_cds_from_assemblies.SPECIES_MAP; it takes effect once the
        # 's/^>phau_/>phaust_/' rename is applied to the Phoronis proteome.
        'phaust', 'phov', 'phpa', 'teme',
    ],
}
for group, prefixes in _groups.items():
    for p in prefixes:
        PREFIX_TO_GROUP[p] = group


def get_taxon_group(leaf_name):
    """Determine taxonomic group from leaf name."""
    if not leaf_name:
        return 'Unknown'
    if leaf_name.startswith('Berste') or leaf_name.startswith('TRINITY_'):
        return 'Berghia'
    name = leaf_name[4:] if leaf_name.startswith('ref_') else leaf_name
    for p in sorted(PREFIX_TO_GROUP.keys(), key=len, reverse=True):
        if name.startswith(p + '_') or name == p:
            return PREFIX_TO_GROUP[p]
    first = name.split('_')[0]
    return PREFIX_TO_GROUP.get(first, 'Unknown')


def compute_candidate_tiers(ranking_csv):
    """Compute data-driven candidate tiers from ranking scores.

    Tier 1 (gold): >= Q90 (top 10th percentile)
    Tier 2 (silver): above largest natural gap in score distribution

    Returns: (tier1_set, tier2_set, q90_threshold, gap_threshold)
    """
    scores = []
    ids_scores = []
    with open(ranking_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            s = float(row['rank_score'])
            scores.append(s)
            ids_scores.append((row['id'], s))

    scores_arr = np.array(scores)
    q90 = float(np.percentile(scores_arr, 90))

    # Find largest gap in descending sorted scores
    sorted_desc = np.sort(scores_arr)[::-1]
    gaps = np.abs(np.diff(sorted_desc))
    # Only look at gaps in top 60% to avoid noise in the tail
    search_range = max(10, int(0.6 * len(gaps)))
    gap_idx = np.argmax(gaps[:search_range])
    gap_threshold = float(sorted_desc[gap_idx + 1])  # score just below the gap

    tier1 = set()
    tier2 = set()
    for cid, s in ids_scores:
        if s >= q90:
            tier1.add(cid)
        elif s >= gap_threshold:
            tier2.add(cid)

    return tier1, tier2, q90, gap_threshold


def draw_tree(tree, tier1, tier2, q90_thresh, gap_thresh, output_path, dpi=300, theme='light'):
    """Draw circular phylogenetic tree with taxonomy ring and tiered candidates."""
    dark = theme == 'dark'
    bg_color = '#000000' if dark else 'white'
    text_color = '#e0e0e0' if dark else 'black'
    branch_alpha = 0.85 if dark else 0.7
    arc_alpha = 0.7 if dark else 0.5
    ring_alpha = 0.95 if dark else 0.8
    root_marker_color = '#ffffff' if dark else 'black'
    root_edge_color = '#000000' if dark else 'white'
    star_edge_color = '#ffffff' if dark else 'black'
    legend_bg = '#111111' if dark else 'white'
    legend_edge = '#555555' if dark else '#cccccc'

    # Override Berghia color for dark theme
    if dark:
        TAXON_COLORS['Berghia'] = '#ffffff'
    else:
        TAXON_COLORS['Berghia'] = '#2d2d2d'

    fig, ax = plt.subplots(1, 1, figsize=(28, 28), subplot_kw={'projection': 'polar'})
    fig.patch.set_facecolor(bg_color)

    terminals = tree.get_terminals()
    n_leaves = len(terminals)

    # Leaf angles (evenly spaced)
    leaf_angles = {}
    for i, leaf in enumerate(terminals):
        leaf_angles[leaf.name or str(i)] = 2 * math.pi * i / n_leaves

    # Depths with clipping
    depths = tree.depths(unit_branch_lengths=False)
    max_depth = max(depths.values()) if depths.values() else 1
    all_depths = sorted(depths.values())
    p95 = all_depths[int(0.95 * len(all_depths))] if all_depths else max_depth
    clip_depth = max(p95 * 1.2, max_depth * 0.5)

    # Node positions (postorder for angles, any order for radii)
    node_angles = {}
    node_radii = {}
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():
            node_angles[id(clade)] = leaf_angles.get(clade.name, 0)
        else:
            child_angles = [node_angles[id(c)] for c in clade.clades if id(c) in node_angles]
            if child_angles:
                mean_sin = sum(math.sin(a) for a in child_angles) / len(child_angles)
                mean_cos = sum(math.cos(a) for a in child_angles) / len(child_angles)
                node_angles[id(clade)] = math.atan2(mean_sin, mean_cos)
            else:
                node_angles[id(clade)] = 0
        d = depths.get(clade, 0)
        node_radii[id(clade)] = min(d, clip_depth) / clip_depth * 0.85

    # Parent map
    parent_map = {}
    for clade in tree.find_clades(order='preorder'):
        for child in clade.clades:
            parent_map[id(child)] = clade

    # --- Draw root marker ---
    root_angle = node_angles.get(id(tree.root), 0)
    root_radius = node_radii.get(id(tree.root), 0)
    ax.scatter([root_angle], [root_radius], c=root_marker_color, s=80, zorder=15,
               marker='D', edgecolors=root_edge_color, linewidths=1.0)

    # --- Draw branches ---
    print(f'  Drawing {n_leaves} leaves...')
    for clade in tree.find_clades(order='preorder'):
        if clade == tree.root:
            continue
        parent = parent_map.get(id(clade))
        if parent is None:
            continue

        angle = node_angles.get(id(clade), 0)
        radius = node_radii.get(id(clade), 0)
        parent_angle = node_angles.get(id(parent), 0)
        parent_radius = node_radii.get(id(parent), 0)

        color = _clade_color(clade)

        ax.plot([angle, angle], [parent_radius, radius],
                color=color, linewidth=0.5 if dark else 0.4, alpha=branch_alpha, solid_capstyle='round')

        if abs(angle - parent_angle) > 0.002:
            a1, a2 = parent_angle, angle
            diff = a2 - a1
            if diff > math.pi:
                a2 -= 2 * math.pi
            elif diff < -math.pi:
                a2 += 2 * math.pi
            arc_angles = np.linspace(a1, a2, max(3, int(abs(a2 - a1) * 30)))
            arc_radii = np.full_like(arc_angles, parent_radius)
            ax.plot(arc_angles, arc_radii,
                    color=color, linewidth=0.5 if dark else 0.4, alpha=arc_alpha, solid_capstyle='round')

    # --- Outer rings ---
    ring_inner = 0.88
    ring_outer = 0.93
    berghia_inner = 0.93
    berghia_outer = 0.96
    half_width = math.pi / n_leaves

    for leaf in terminals:
        name = leaf.name or ''
        group = get_taxon_group(name)
        angle = leaf_angles.get(name, 0)
        color = TAXON_COLORS.get(group, '#cccccc')

        if group == 'Berghia':
            if name in tier1:
                c = TAXON_COLORS['Tier 1']
                ax.bar(angle, berghia_outer - berghia_inner, width=half_width * 2,
                       bottom=berghia_inner, color=c, alpha=1.0, edgecolor='none')
            else:
                ax.bar(angle, berghia_outer - berghia_inner, width=half_width * 2,
                       bottom=berghia_inner, color=color, alpha=0.6 if dark else 0.5, edgecolor='none')
        else:
            ax.bar(angle, ring_outer - ring_inner, width=half_width * 2,
                   bottom=ring_inner, color=color, alpha=ring_alpha, edgecolor='none')

    # Star markers for Tier 1
    for leaf in terminals:
        name = leaf.name or ''
        if name in tier1:
            angle = leaf_angles.get(name, 0)
            ax.scatter([angle], [1.00], c=TAXON_COLORS['Tier 1'],
                       s=60, zorder=10, alpha=1.0, edgecolors=star_edge_color,
                       linewidths=0.6, marker='*')

    # --- Style ---
    ax.set_ylim(0, 1.08)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    ax.set_facecolor(bg_color)

    # --- Legend ---
    group_counts = defaultdict(int)
    for leaf in terminals:
        group_counts[get_taxon_group(leaf.name or '')] += 1

    legend_patches = []
    for group in ['Gastropoda', 'Bivalvia', 'Cephalopoda', 'Annelida',
                   'Platyhelminthes', 'Other Mollusca', 'Other Lophotrochozoa']:
        n = group_counts.get(group, 0)
        if n > 0:
            legend_patches.append(mpatches.Patch(
                color=TAXON_COLORS[group], label=f'{group} (n={n})'))

    n_berghia = group_counts.get('Berghia', 0)
    legend_patches.append(mpatches.Patch(
        color=TAXON_COLORS['Berghia'], label=f'B. stephanieae (n={n_berghia})'))
    legend_patches.append(plt.scatter([], [], c=TAXON_COLORS['Tier 1'],
                                       s=100, marker='*', edgecolors=star_edge_color, linewidths=0.6,
                                       label=f'Tier 1: score \u2265 {q90_thresh:.1f} (P90, n={len(tier1)})'))
    legend_patches.append(mlines.Line2D([], [], color=root_marker_color, marker='D', markersize=8,
                                         markerfacecolor=root_marker_color, markeredgecolor=root_edge_color,
                                         linestyle='None', label='Root'))

    leg = ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(-0.15, 1.10),
                    fontsize=11, ncol=2, handlelength=1.5,
                    facecolor=legend_bg, edgecolor=legend_edge,
                    labelcolor=text_color, framealpha=0.95)

    ax.set_title(f'GPCR Phylogeny \u2014 Berghia stephanieae\n'
                 f'VT+I+R10 | {n_leaves} sequences | 1,000 UFBoot | '
                 f'Midpoint rooted, ladderized',
                 fontsize=18, fontweight='bold', pad=30, color=text_color)

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor=bg_color)
    plt.close()
    print(f'Saved: {output_path}')


def _clade_color(clade):
    """Color for a clade based on majority leaf group (excluding Berghia weight)."""
    if clade.is_terminal():
        return TAXON_COLORS.get(get_taxon_group(clade.name or ''), '#cccccc')
    groups = defaultdict(float)
    for leaf in clade.get_terminals():
        g = get_taxon_group(leaf.name or '')
        groups[g] += 1 if g != 'Berghia' else 0.3
    if groups:
        return TAXON_COLORS.get(max(groups, key=groups.get), '#cccccc')
    return '#cccccc'


def main():
    parser = argparse.ArgumentParser(description='Visualize GPCR tree with taxonomy coloring')
    parser.add_argument('--tree', required=True, help='Newick tree file')
    parser.add_argument('--ranking', required=True, help='Ranked candidates CSV')
    parser.add_argument('--output', required=True, help='Output PNG path')
    parser.add_argument('--theme', choices=['light', 'dark'], default='light')
    parser.add_argument('--dpi', type=int, default=300, help='Output DPI')
    args = parser.parse_args()

    # Compute candidate tiers from ranking data
    print('Computing candidate tiers...')
    tier1, tier2, q90, gap_thresh = compute_candidate_tiers(args.ranking)
    print(f'  Tier 1 (P90, score >= {q90:.2f}): {len(tier1)} candidates')
    print(f'  Tier 2 (natural gap, score >= {gap_thresh:.2f}): {len(tier2)} candidates')

    # Load tree
    print('Loading tree...')
    tree = Phylo.read(args.tree, 'newick')
    print(f'  {len(tree.get_terminals())} tips')

    print('Midpoint rooting...')
    tree.root_at_midpoint()

    print('Ladderizing (descending)...')
    tree.ladderize(reverse=True)

    print(f'Drawing ({args.theme} theme)...')
    draw_tree(tree, tier1, tier2, q90, gap_thresh, args.output, args.dpi, args.theme)


if __name__ == '__main__':
    main()
