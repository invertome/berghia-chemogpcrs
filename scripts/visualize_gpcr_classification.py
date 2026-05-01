#!/usr/bin/env python3
"""
Visualize GPCR phylogenetic tree colored by expansion-based functional classification:
  - Large expansion (>=10 Berghia paralogs): chemoreceptor candidates
  - Moderate expansion (3-9 paralogs): recent duplications
  - Conserved (1-2 Berghia + >=3 refs): signaling GPCRs
  - Reference sequences colored by Nath et al. LSE vs 1:1 ortholog
  - Midpoint rooted (outgroup rooting invalid for multi-subfamily gene tree)
  - Ladderized descending
"""

import argparse
import csv
import math
from collections import defaultdict, Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from Bio import Phylo

# Color scheme
COLORS = {
    'large_expansion':   '#e41a1c',  # red — chemoreceptor candidate
    'moderate_expansion': '#ff7f00',  # orange — recent duplication
    'small_clade':       '#984ea3',  # purple — unresolved
    'conserved':         '#377eb8',  # blue — signaling GPCR
    'lse':               '#e41a1c',  # red — LSE reference
    'one_to_one_ortholog': '#377eb8', # blue — conserved reference
    'unknown':           '#cccccc',
    'tier1':             '#FFD700',  # gold
}


def draw_tree(tree, ref_category, berghia_class, tier1, output_path, dpi=300, theme='light'):
    """Draw circular tree colored by expansion-based classification."""
    dark = theme == 'dark'
    bg_color = '#000000' if dark else 'white'
    text_color = '#e0e0e0' if dark else 'black'
    branch_alpha = 0.85 if dark else 0.7
    arc_alpha = 0.7 if dark else 0.5
    ring_alpha = 0.95 if dark else 0.8
    root_mc = '#ffffff' if dark else 'black'
    root_ec = '#000000' if dark else 'white'
    star_ec = '#ffffff' if dark else 'black'
    legend_bg = '#111111' if dark else 'white'
    legend_ec = '#555555' if dark else '#cccccc'

    fig, ax = plt.subplots(1, 1, figsize=(28, 28), subplot_kw={'projection': 'polar'})
    fig.patch.set_facecolor(bg_color)

    terminals = tree.get_terminals()
    n_leaves = len(terminals)

    leaf_angles = {}
    for i, leaf in enumerate(terminals):
        leaf_angles[leaf.name or str(i)] = 2 * math.pi * i / n_leaves

    depths = tree.depths(unit_branch_lengths=False)
    all_depths = sorted(depths.values())
    p95 = all_depths[int(0.95 * len(all_depths))]
    clip_depth = max(p95 * 1.2, max(depths.values()) * 0.5)

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

    parent_map = {}
    for clade in tree.find_clades(order='preorder'):
        for child in clade.clades:
            parent_map[id(child)] = clade

    # Root marker
    root_angle = node_angles.get(id(tree.root), 0)
    root_radius = node_radii.get(id(tree.root), 0)
    ax.scatter([root_angle], [root_radius], c=root_mc, s=80, zorder=15,
               marker='D', edgecolors=root_ec, linewidths=1.0)

    # Branches
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

        color = _branch_color(clade, ref_category, berghia_class)

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

    # Outer ring
    ring_inner = 0.88
    ring_outer = 0.96
    half_width = math.pi / n_leaves

    for leaf in terminals:
        name = leaf.name or ''
        angle = leaf_angles.get(name, 0)

        if name.startswith('Berste'):
            cat = berghia_class.get(name, 'unknown')
        else:
            cat = ref_category.get(name, 'unknown')

        color = COLORS.get(cat, '#cccccc')
        alpha = 0.9 if name.startswith('Berste') else 0.7

        ax.bar(angle, ring_outer - ring_inner, width=half_width * 2,
               bottom=ring_inner, color=color, alpha=ring_alpha if not name.startswith('Berste') else alpha, edgecolor='none')

    # Tier 1 stars
    for leaf in terminals:
        name = leaf.name or ''
        if name in tier1:
            angle = leaf_angles.get(name, 0)
            ax.scatter([angle], [0.99], c=COLORS['tier1'], s=60, zorder=10,
                       alpha=1.0, edgecolors=star_ec, linewidths=0.6, marker='*')

    # Style
    ax.set_ylim(0, 1.08)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    ax.set_facecolor(bg_color)

    # Counts
    b_cats = Counter(berghia_class.values())
    r_cats = Counter(ref_category.values())

    # Legend
    legend_patches = [
        mpatches.Patch(color=COLORS['large_expansion'],
                       label=f'Large expansion (\u226510 paralogs): '
                             f'Berghia {b_cats.get("large_expansion",0)}, '
                             f'LSE refs {r_cats.get("lse",0)}'),
        mpatches.Patch(color=COLORS['moderate_expansion'],
                       label=f'Moderate expansion (3\u20139 paralogs): '
                             f'Berghia {b_cats.get("moderate_expansion",0)}'),
        mpatches.Patch(color=COLORS['conserved'],
                       label=f'Conserved signaling (\u22642 + \u22653 refs): '
                             f'Berghia {b_cats.get("conserved",0)}, '
                             f'1:1 refs {r_cats.get("one_to_one_ortholog",0)}'),
        mpatches.Patch(color=COLORS['small_clade'],
                       label=f'Small clade (unresolved): '
                             f'Berghia {b_cats.get("small_clade",0)}'),
    ]

    legend_patches.append(plt.scatter([], [], c=COLORS['tier1'], s=100, marker='*',
                                       edgecolors=star_ec, linewidths=0.6,
                                       label=f'Top candidates (P90, n={len(tier1)})'))
    legend_patches.append(mlines.Line2D([], [], color=root_mc, marker='D', markersize=8,
                                         markerfacecolor=root_mc, markeredgecolor=root_ec,
                                         linestyle='None', label='Root (midpoint)'))

    leg = ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(-0.15, 1.10),
                    fontsize=11, ncol=1, handlelength=1.5,
                    facecolor=legend_bg, edgecolor=legend_ec,
                    labelcolor=text_color, framealpha=0.95)

    ax.set_title(f'GPCR Functional Classification \u2014 Berghia stephanieae\n'
                 f'Expansion-based: {b_cats.get("large_expansion",0)} large + '
                 f'{b_cats.get("moderate_expansion",0)} moderate expansions, '
                 f'{b_cats.get("conserved",0)} conserved | Midpoint rooted',
                 fontsize=16, fontweight='bold', pad=30, color=text_color)

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor=bg_color)
    plt.close()
    print(f'Saved: {output_path}')


def _branch_color(clade, ref_category, berghia_class):
    """Color branch by functional category of descendants."""
    if clade.is_terminal():
        name = clade.name or ''
        if name.startswith('Berste'):
            return COLORS.get(berghia_class.get(name, 'unknown'), '#cccccc')
        return COLORS.get(ref_category.get(name, 'unknown'), '#cccccc')

    cats = defaultdict(float)
    for leaf in clade.get_terminals():
        name = leaf.name or ''
        if name.startswith('Berste'):
            cat = berghia_class.get(name, 'unknown')
        else:
            cat = ref_category.get(name, 'unknown')

        if cat in ('large_expansion', 'lse'):
            cats['expansion'] += 1
        elif cat in ('conserved', 'one_to_one_ortholog'):
            cats['conserved'] += 1
        elif cat == 'moderate_expansion':
            cats['moderate'] += 1
        else:
            cats['other'] += 0.3

    if cats:
        top = max(cats, key=cats.get)
        if top == 'expansion':
            return COLORS['large_expansion']
        elif top == 'conserved':
            return COLORS['conserved']
        elif top == 'moderate':
            return COLORS['moderate_expansion']
    return '#cccccc'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', required=True)
    parser.add_argument('--ranking', required=True)
    parser.add_argument('--classification', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--theme', choices=['light', 'dark'], default='light')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()

    # Load Berghia classification
    berghia_class = {}
    with open(args.classification) as f:
        reader = csv.DictReader(f)
        for row in reader:
            berghia_class[row['id']] = row['classification']
    print(f'Berghia: {Counter(berghia_class.values())}')

    # Load reference categories
    short_to_original = {}
    with open('preliminary/results/reference_sequences/subsampled_id_map.csv') as f:
        reader = csv.DictReader(f)
        for row in reader:
            short_to_original[row['short_id']] = row['original_id']

    header_to_category = {}
    with open('preliminary/results/reference_sequences/header_taxid_map.tsv') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                header_to_category[parts[0]] = parts[3]

    ref_category = {}
    tree = Phylo.read(args.tree, 'newick')
    for leaf in tree.get_terminals():
        name = leaf.name or ''
        if not name.startswith('Berste'):
            orig = short_to_original.get(name)
            if orig:
                ref_category[name] = header_to_category.get(orig, 'unknown')
    print(f'References: {Counter(ref_category.values())}')

    # Top candidates
    scores = []
    ids_scores = []
    with open(args.ranking) as f:
        reader = csv.DictReader(f)
        for row in reader:
            s = float(row['rank_score'])
            scores.append(s)
            ids_scores.append((row['id'], s))
    q90 = float(np.percentile(scores, 90))
    tier1 = {cid for cid, s in ids_scores if s >= q90}
    print(f'Tier 1 (P90 >= {q90:.2f}): {len(tier1)}')

    # Midpoint root
    print('Midpoint rooting...')
    tree.root_at_midpoint()
    tree.ladderize(reverse=True)

    print(f'Drawing ({args.theme} theme)...')
    draw_tree(tree, ref_category, berghia_class, tier1, args.output, args.dpi, args.theme)


if __name__ == '__main__':
    main()
