#!/usr/bin/env python3
"""
Rectangular ("vertical") renderings of the preliminary v2 GPCR phylogeny.

Two figures, both reusing the taxonomy colors + data-driven candidate tiers
from visualize_gpcr_tree.py (the circular renderer) so styling stays identical:

  full : whole tree as a rectangular phylogram (root left, 2,371 tips stacked
         vertically down the right), branch lengths kept with the same long-branch
         compression the circular figure uses, taxonomy color bar + Tier-1 stars.

  lse  : zoom into the large Berghia lineage-specific expansion (the 240-tip,
         ~96%-Berghia clade holding 42 of the 44 Tier-1 candidates), with the
         Tier-1 candidates annotated (rank #, id, score).

Usage:
  python3 scripts/visualize_gpcr_tree_rect.py \
    --tree preliminary/results/phylogenies/protein/v2/gpcrs.treefile \
    --ranking preliminary/results/ranking/ranked_candidates_sorted.csv \
    --outdir preliminary/results/phylogenies/protein/v2/figures
"""

import argparse
import csv
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as PhyloTree

# Reuse the established taxonomy map, group lookup, tier logic and clade-color
# helper from the circular renderer (same directory, importable).
import visualize_gpcr_tree as vgt
from visualize_gpcr_tree import (
    TAXON_COLORS, LEGEND_GROUPS, get_taxon_group, compute_candidate_tiers,
    _clade_color, validate_prefix_map,
)
from species_code_lookup import (
    AmbiguousSpeciesCodeError, build_code_species_map, species_for,
    berghia_display_name,
)


def is_berghia(name):
    return bool(name) and (name.startswith('Berste') or name.startswith('TRINITY_'))


def _ref_label(name, code_map):
    """Display label for a non-Berghia reference tip.

    Falls back to the raw leaf id when the species cannot be resolved, which
    includes the case where the code is ambiguous and ``species_for`` refuses
    to guess. An unresolved identifier is honest; a confident wrong binomial in
    a publication figure is not.
    """
    try:
        sp, locus = species_for(name, code_map)
    except AmbiguousSpeciesCodeError:
        return name
    if not sp:
        return name
    return f'{sp} ({locus})' if locus else sp


def load_ranking(ranking_csv):
    """Return (rank_of, score_of) keyed by candidate id (rank = 1-based row order)."""
    rank_of, score_of = {}, {}
    with open(ranking_csv) as f:
        for i, row in enumerate(csv.DictReader(f)):
            rank_of[row['id']] = i + 1
            score_of[row['id']] = float(row['rank_score'])
    return rank_of, score_of


def select_lse_clade(tree, tier1, min_berghia_frac=0.85, min_size=20):
    """Deterministically pick the large Berghia LSE.

    Among clades that are >= min_berghia_frac Berghia and >= min_size tips,
    take the one holding the most Tier-1 candidates, breaking ties toward the
    smallest such clade. On the v2 tree this is the 240-tip / ~96%-Berghia clade.
    """
    best_key, best_clade = None, None
    for cl in tree.get_nonterminals():
        leaves = cl.get_terminals()
        n = len(leaves)
        if n < min_size:
            continue
        b = sum(1 for l in leaves if is_berghia(l.name))
        if b / n < min_berghia_frac:
            continue
        nt = sum(1 for l in leaves if l.name in tier1)
        key = (nt, -n)  # most Tier-1, then smallest
        if best_key is None or key > best_key:
            best_key, best_clade = key, cl
    return best_clade


def _layout(tree, compress=True):
    """Rectangular coordinates: x = (compressed) root-depth, y = tip order.

    Returns (xn, y, n, clip, parent_map) where xn is normalized to [0, 1].
    """
    terms = tree.get_terminals()
    n = len(terms)
    y = {id(t): i for i, t in enumerate(terms)}
    for cl in tree.find_clades(order='postorder'):
        if not cl.is_terminal():
            ys = [y[id(c)] for c in cl.clades if id(c) in y]
            y[id(cl)] = sum(ys) / len(ys) if ys else 0.0

    depths = tree.depths()  # branch-length distance from root
    dvals = sorted(depths.values()) or [1.0]
    if compress:
        p95 = dvals[int(0.95 * len(dvals))]
        clip = max(p95 * 1.2, max(dvals) * 0.5)
    else:
        clip = max(dvals)
    clip = clip or 1.0
    xn = {id(cl): min(d, clip) / clip for cl, d in depths.items()}

    parent_map = {}
    for cl in tree.find_clades(order='preorder'):
        for ch in cl.clades:
            parent_map[id(ch)] = cl
    return xn, y, n, clip, parent_map


def _theme(theme):
    dark = theme == 'dark'
    TAXON_COLORS['Berghia'] = '#ffffff' if dark else '#2d2d2d'
    return {
        'dark': dark,
        'bg': '#000000' if dark else 'white',
        'text': '#e0e0e0' if dark else 'black',
        'branch_alpha': 0.85 if dark else 0.7,
        'root_marker': '#ffffff' if dark else 'black',
        'root_edge': '#000000' if dark else 'white',
        'star_edge': '#ffffff' if dark else 'black',
        'legend_bg': '#111111' if dark else 'white',
        'legend_edge': '#555555' if dark else '#cccccc',
        'tip_dot': '#888888' if dark else '#bbbbbb',
    }


def _scale_bar(ax, clip, theme, y_at, x0=0.02, fontsize=10):
    """Draw a substitutions/site scale bar sized to the compression."""
    nice = 0.5 if clip >= 1.0 else round(clip / 4, 2)
    w = nice / clip
    ax.plot([x0, x0 + w], [y_at, y_at], color=theme['text'], lw=2, solid_capstyle='butt')
    ax.text(x0 + w / 2, y_at, f'{nice:g} subs/site', ha='center', va='bottom',
            fontsize=fontsize, color=theme['text'])


def draw_full(tree, tier1, q90, output, dpi, theme_name):
    """Whole tree, rectangular phylogram, taxonomy color bar + Tier-1 stars."""
    th = _theme(theme_name)
    xn, y, n, clip, parent_map = _layout(tree, compress=True)

    fig_h = max(24.0, n * 0.018)
    fig, ax = plt.subplots(figsize=(18, fig_h))
    fig.patch.set_facecolor(th['bg'])
    ax.set_facecolor(th['bg'])

    color_cache = {}
    def cc(cl):
        k = id(cl)
        if k not in color_cache:
            color_cache[k] = _clade_color(cl)
        return color_cache[k]

    # branches (square phylogram)
    for cl in tree.find_clades(order='preorder'):
        col = cc(cl)
        if cl.clades:
            ys = [y[id(c)] for c in cl.clades]
            ax.plot([xn[id(cl)], xn[id(cl)]], [min(ys), max(ys)],
                    color=col, lw=0.25, alpha=th['branch_alpha'], solid_capstyle='round')
        p = parent_map.get(id(cl))
        if p is not None:
            ax.plot([xn[id(p)], xn[id(cl)]], [y[id(cl)], y[id(cl)]],
                    color=col, lw=0.25, alpha=th['branch_alpha'], solid_capstyle='round')

    # root marker
    ax.scatter([xn[id(tree.root)]], [y[id(tree.root)]], c=th['root_marker'], s=60,
               marker='D', edgecolors=th['root_edge'], linewidths=1.0, zorder=15)

    # taxonomy color bar (one row per tip) just past the tips, + Tier-1 stars
    bar_x0, bar_w = 1.03, 0.035
    for leaf in tree.get_terminals():
        name = leaf.name or ''
        group = get_taxon_group(name)
        yy = y[id(leaf)]
        ax.add_patch(mpatches.Rectangle((bar_x0, yy - 0.5), bar_w, 1.0,
                     facecolor=TAXON_COLORS.get(group, '#cccccc'), edgecolor='none',
                     alpha=0.95 if th['dark'] else 0.85))
        if name in tier1:
            ax.scatter([bar_x0 + bar_w + 0.02], [yy], c=TAXON_COLORS['Tier 1'], s=40,
                       marker='*', edgecolors=th['star_edge'], linewidths=0.5, zorder=10)

    ax.set_xlim(-0.02, 1.12)
    ax.set_ylim(n + 1, -2)  # tip 0 at top
    ax.axis('off')
    _scale_bar(ax, clip, th, y_at=n + 0.5)

    _legend(ax, tree, tier1, q90, th)
    ax.set_title('GPCR Phylogeny — Berghia stephanieae (rectangular)\n'
                 f'VT+I+R10 | {n} sequences | 1,000 UFBoot | Midpoint rooted, '
                 'ladderized | branch lengths compressed above p95',
                 fontsize=16, fontweight='bold', color=th['text'], pad=24)

    plt.savefig(output, dpi=dpi, bbox_inches='tight', facecolor=th['bg'])
    plt.close()
    print(f'Saved: {output}')


def draw_lse(clade, tier1, rank_of, score_of, code_map, output, dpi, theme_name,
             label_mode='all'):
    """Zoom of the Berghia LSE clade as a rectangular phylogram.

    label_mode controls labels + vertical compaction:
      'all'        - every tip labeled (Tier-1 bold, other Berghia gray, refs italic).
      'candidates' - label only Tier-1 candidates + non-Berghia refs; other Berghia
                     are unlabeled hairline rows, so the tree compacts vertically.
      'none'       - no text labels; keep the gold Tier-1 stars + taxon colors only,
                     tightly packed (smallest).
    """
    th = _theme(theme_name)
    berghia_gray = '#9a9a9a' if th['dark'] else '#555555'
    sub = PhyloTree(root=clade)
    sub.ladderize(reverse=True)
    xn, _y_uniform, n, clip, parent_map = _layout(sub, compress=False)

    terms = sub.get_terminals()

    def labeled(nm):
        if label_mode == 'none':
            return False
        if label_mode == 'candidates':
            return (nm in tier1) or (not is_berghia(nm))
        return True

    # Mode-specific row weights (drive compaction), figure size, marker/label sizes.
    if label_mode == 'none':
        weights = [1.0] * len(terms)
        row_in, fig_w, xlim_r = 0.04, 13.0, 1.10
        sz_star, sz_ref, sz_b, lw_branch = 48, 18, 4, 0.5
        font_t1 = font_b = font_ref = None
    elif label_mode == 'candidates':
        weights = [1.0 if labeled(l.name) else 0.10 for l in terms]
        row_in, fig_w, xlim_r = 0.26, 18.0, 2.25
        sz_star, sz_ref, sz_b, lw_branch = 100, 45, 11, 0.9
        font_t1, font_b, font_ref = 11.0, 11.0, 11.0
    else:  # 'all'
        weights = [1.0] * len(terms)
        row_in, fig_w, xlim_r = 0.21, 18.0, 2.25
        sz_star, sz_ref, sz_b, lw_branch = 110, 50, 13, 1.0
        font_t1, font_b, font_ref = 11.5, 11.0, 11.5

    # Variable-height tip slots (centered), internal nodes = mean of children.
    y = {}
    pos = 0.0
    for leaf, w in zip(terms, weights):
        y[id(leaf)] = pos + w / 2.0
        pos += w
    total = pos
    for cl in sub.find_clades(order='postorder'):
        if not cl.is_terminal():
            ys = [y[id(c)] for c in cl.clades if id(c) in y]
            y[id(cl)] = sum(ys) / len(ys) if ys else 0.0

    fig_h = max(8.0, total * row_in)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor(th['bg'])
    ax.set_facecolor(th['bg'])

    color_cache = {}
    def cc(cl):
        k = id(cl)
        if k not in color_cache:
            color_cache[k] = _clade_color(cl)
        return color_cache[k]

    for cl in sub.find_clades(order='preorder'):
        col = cc(cl)
        if cl.clades:
            ys = [y[id(c)] for c in cl.clades]
            ax.plot([xn[id(cl)], xn[id(cl)]], [min(ys), max(ys)],
                    color=col, lw=lw_branch, alpha=th['branch_alpha'], solid_capstyle='round')
        p = parent_map.get(id(cl))
        if p is not None:
            ax.plot([xn[id(p)], xn[id(cl)]], [y[id(cl)], y[id(cl)]],
                    color=col, lw=lw_branch, alpha=th['branch_alpha'], solid_capstyle='round')

    # Labeled tips get a leader to the aligned label column; others just a marker.
    label_x = 1.02
    for leaf in terms:
        name = leaf.name or ''
        yy, xx = y[id(leaf)], xn[id(leaf)]
        lab = labeled(name)
        if lab:
            ax.plot([xx + 0.004, label_x - 0.004], [yy, yy],
                    color=berghia_gray, lw=0.4, alpha=0.5, zorder=1)
        if name in tier1:
            ax.scatter([xx], [yy], c=TAXON_COLORS['Tier 1'], s=sz_star, marker='*',
                       edgecolors=th['star_edge'], linewidths=0.7, zorder=12)
            if lab:
                # Tier-1 ids come from the ranking CSV, which is not guaranteed
                # to hold only Berghia transcripts; label a non-Berghia tip as
                # what it actually is rather than as Berghia.
                label = (berghia_display_name(name) if is_berghia(name)
                         else _ref_label(name, code_map))
                ax.text(label_x, yy, f'#{rank_of[name]}  {label}',
                        va='center', ha='left', fontsize=font_t1, color=th['text'],
                        fontweight='bold', style='italic', zorder=12, clip_on=False)
        elif is_berghia(name):
            ax.scatter([xx], [yy], c=berghia_gray, s=sz_b, marker='o',
                       edgecolors='none', zorder=9)
            if lab:
                ax.text(label_x, yy, berghia_display_name(name), va='center', ha='left',
                        fontsize=font_b, color=berghia_gray, style='italic', zorder=9,
                        clip_on=False)
        else:
            group = get_taxon_group(name, code_map)
            col = TAXON_COLORS.get(group, '#cccccc')
            ax.scatter([xx], [yy], c=col, s=sz_ref, marker='o',
                       edgecolors=th['star_edge'], linewidths=0.4, zorder=11)
            if lab:
                label = _ref_label(name, code_map)
                ax.text(label_x, yy, label, va='center', ha='left', fontsize=font_ref,
                        color=col, style='italic', fontweight='bold', zorder=11, clip_on=False)

    ax.set_xlim(-0.02, xlim_r)
    ax.set_ylim(total + 2, -2)
    ax.axis('off')
    _scale_bar(ax, clip, th, y_at=total + 1.2, fontsize=14)

    # Legend: candidate markers, then the taxon-group colors used for the
    # non-Berghia reference tips (same scheme as the circular figure).
    from collections import defaultdict
    ref_group_counts = defaultdict(int)
    for leaf in sub.get_terminals():
        nm = leaf.name or ''
        if not is_berghia(nm):
            ref_group_counts[get_taxon_group(nm, code_map)] += 1

    handles = [
        plt.scatter([], [], c=TAXON_COLORS['Tier 1'], s=120, marker='*',
                    edgecolors=th['star_edge'], linewidths=0.6,
                    label='Tier-1 candidate'),
        mpatches.Patch(color=berghia_gray, label=r'$\mathit{B.\ stephanieae}$'),
    ]
    for group in LEGEND_GROUPS:
        if ref_group_counts.get(group):
            handles.append(mpatches.Patch(color=TAXON_COLORS[group], label=group))
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 1.0), fontsize=15,
              facecolor=th['legend_bg'], edgecolor=th['legend_edge'],
              labelcolor=th['text'], framealpha=0.95)
    ax.set_title(r'$\mathbfit{Berghia\ stephanieae}$ — large lineage-specific expansion (LSE)',
                 fontsize=21, fontweight='bold', color=th['text'], pad=20)

    plt.savefig(output, dpi=dpi, bbox_inches='tight', facecolor=th['bg'])
    plt.close()
    print(f'Saved: {output}  [{label_mode}: {len(terms)} tips, ~{fig_h:.0f} in tall]')


def _legend(ax, tree, tier1, q90, th):
    from collections import defaultdict
    counts = defaultdict(int)
    for leaf in tree.get_terminals():
        counts[get_taxon_group(leaf.name or '')] += 1
    patches = []
    for group in LEGEND_GROUPS:
        if counts.get(group, 0):
            patches.append(mpatches.Patch(color=TAXON_COLORS[group],
                                          label=f'{group} (n={counts[group]})'))
    patches.append(mpatches.Patch(color=TAXON_COLORS['Berghia'],
                                  label=f"B. stephanieae (n={counts.get('Berghia', 0)})"))
    patches.append(plt.scatter([], [], c=TAXON_COLORS['Tier 1'], s=100, marker='*',
                               edgecolors=th['star_edge'], linewidths=0.6,
                               label=f'Tier 1: score ≥ {q90:.1f} (P90, n={len(tier1)})'))
    patches.append(mlines.Line2D([], [], color=th['root_marker'], marker='D', linestyle='None',
                                 markersize=8, markerfacecolor=th['root_marker'],
                                 markeredgecolor=th['root_edge'], label='Root'))
    ax.legend(handles=patches, loc='upper left', bbox_to_anchor=(0.0, 1.0), fontsize=11,
              ncol=2, handlelength=1.5, facecolor=th['legend_bg'], edgecolor=th['legend_edge'],
              labelcolor=th['text'], framealpha=0.95)


def main():
    ap = argparse.ArgumentParser(description='Rectangular GPCR tree + Berghia LSE zoom')
    ap.add_argument('--tree', required=True)
    ap.add_argument('--ranking', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--figures', default='full,lse,lse_compact,lse_nolabels',
                    help='comma list: full, lse, lse_compact, lse_nolabels')
    ap.add_argument('--themes', default='light,dark', help='comma list: light,dark')
    ap.add_argument('--dpi-full', type=int, default=250)
    ap.add_argument('--dpi-lse', type=int, default=300)
    ap.add_argument('--skip-prefix-validation', action='store_true',
                    help='render even if PREFIX_TO_GROUP has drifted from the '
                         'reference proteomes (NOT for publication figures)')
    args = ap.parse_args()

    vgt._validate_prefix_map_or_die(args.skip_prefix_validation)

    os.makedirs(args.outdir, exist_ok=True)
    figures = [f.strip() for f in args.figures.split(',') if f.strip()]
    themes = [t.strip() for t in args.themes.split(',') if t.strip()]

    tier1, _tier2, q90, _gap = compute_candidate_tiers(args.ranking)
    rank_of, score_of = load_ranking(args.ranking)
    print(f'Tier 1 (P90, score >= {q90:.2f}): {len(tier1)} candidates')

    lse_specs = [('lse', 'all', 'berghia_lse_candidates'),
                 ('lse_compact', 'candidates', 'berghia_lse_compact'),
                 ('lse_nolabels', 'none', 'berghia_lse_nolabels')]
    want_lse = [s for s in lse_specs if s[0] in figures]

    code_map = build_code_species_map() if want_lse else {}
    if code_map:
        print(f'Resolved {len(code_map)} reference species codes')

    for theme in themes:
        # fresh root+ladderize per theme so coordinates are independent of prior mutation
        tree = Phylo.read(args.tree, 'newick')
        tree.root_at_midpoint()
        tree.ladderize(reverse=True)

        if 'full' in figures:
            for ext in ('png', 'pdf'):
                draw_full(tree, tier1, q90,
                          os.path.join(args.outdir, f'gpcrs_taxonomy_rect_{theme}.{ext}'),
                          args.dpi_full, theme)
        if want_lse:
            clade = select_lse_clade(tree, tier1)
            print(f'LSE clade: {len(clade.get_terminals())} tips')
            for _key, mode, stem in want_lse:
                for ext in ('png', 'pdf'):
                    draw_lse(clade, tier1, rank_of, score_of, code_map,
                             os.path.join(args.outdir, f'{stem}_{theme}.{ext}'),
                             args.dpi_lse, theme, label_mode=mode)


if __name__ == '__main__':
    main()
