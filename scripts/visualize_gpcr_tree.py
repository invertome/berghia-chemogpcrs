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
import glob
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from Bio import Phylo

PROJECT_ROOT = Path(__file__).resolve().parent.parent

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
    # Deliberately non-taxonomic neutrals. 'Ambiguous' is darker than 'Unknown'
    # so an unresolvable tip is visible rather than reading as a faint gap.
    'Ambiguous':        '#666666',
    'Unknown':          '#cccccc',
}

# Groups a legend enumerates, in order. 'Ambiguous'/'Unknown' are listed last and
# are only drawn when non-empty, so a figure never hides tips it could not place.
LEGEND_GROUPS = [
    'Gastropoda', 'Bivalvia', 'Cephalopoda', 'Annelida', 'Platyhelminthes',
    'Other Mollusca', 'Other Lophotrochozoa', 'Ambiguous', 'Unknown',
]

# Header-prefix codes claimed by proteomes in MORE THAN ONE taxonomic group.
# These can never be resolved from the code alone, so they are deliberately
# absent from PREFIX_TO_GROUP and get_taxon_group() returns 'Ambiguous' for
# them. Keep in sync with species_code_lookup.build_code_species_map(), which
# raises AmbiguousSpeciesCodeError for the same codes;
# validate_prefix_map() asserts both registers match the real proteomes.
#
# 'phau' (measured 2026-07): 1,078 headers in
# lse+one_to_one/gastropoda/109671_Physella_acuta.faa (Mollusca > Gastropoda)
# and 428 in .../other_lophotrochozoan_phyla/115415_Phoronis_australis.faa
# (Phoronida). The 's/^>phau_/>phaust_/' rename that would separate them is
# tracked separately and has NOT been applied; identifiers here are write-once,
# so the display layer stays honest instead of guessing.
AMBIGUOUS_PREFIXES = {
    'phau': ('Physella acuta', 'Phoronis australis'),
}

# Species prefix -> taxonomic group
PREFIX_TO_GROUP = {}
_groups = {
    'Gastropoda': [
        'alvmar', 'amba', 'anvo', 'aplcal', 'arvu', 'baar', 'baare', 'bigl', 'bipf',
        'bist', 'caun', 'chori', 'chsq', 'cobe', 'coco', 'cotr', 'dipe', 'drsu',
        'elch', 'elma', 'giae', 'gima', 'goco', 'gofi', 'goge', 'goku', 'gole',
        'hacra', 'hadi', 'hala', 'haru', 'logi', 'losc', 'lyst', 'metu', 'orid',
        # 'phau' (Physella acuta) is NOT listed: the same code is carried by
        # Phoronis australis. See AMBIGUOUS_PREFIXES.
        'pade', 'pape', 'pavu', 'phli', 'ploc', 'poca', 'raau', 'rave', 'trte',
    ],
    'Bivalvia': [
        'arma', 'arpu', 'atja', 'bapl', 'bofu', 'ceed', 'chfa', 'cobi', 'cofl', 'coku',
        'craan', 'craar', 'crahon', 'crgi', 'crpl', 'crvi', 'drpo', 'drro', 'frfr',
        # 'liant' = Lithophaga antillarum (772 seqs). Was 'lian2', which matched
        # no sequence anywhere, so those 772 rendered as 'Unknown'.
        'frwh', 'gate', 'hihi', 'hycu', 'lael', 'lifo', 'liant', 'lini', 'maho',
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
        'lute', 'mevu', 'owfu', 'paec', 'papa', 'pige', 'pldu', 'poma', 'ripis',
        'ripa', 'sinu', 'stbe', 'stli', 'tela', 'whpi',
    ],
    # 'scint' = Schistosoma intercalatum (71 seqs); was 'scin', which matched
    # nothing. 'ecol', 'meco', 'moex', 'pake', 'prxe', 'taas', 'tapi', 'taso'
    # were never mapped at all (361 seqs rendering as 'Unknown').
    'Platyhelminthes': [
        'atwi', 'cafo', 'clsi', 'dican', 'diden', 'dila', 'duja', 'eccap', 'eccan',
        'ecgr', 'ecmu', 'ecol', 'euni', 'fabu', 'fagi', 'fahe', 'gitig', 'gybu',
        'gysa', 'heame', 'hydi', 'hyta', 'macl', 'mahy', 'mali', 'meco', 'moex',
        'opfe', 'opvi', 'pahe', 'pake', 'pask', 'pawe', 'prcro', 'prfa', 'prxe',
        'rona', 'scbo', 'sccu', 'scgu', 'scha', 'scint', 'scja', 'scman', 'scmar',
        'scmat', 'scme', 'scro', 'scso', 'sctur', 'sper', 'sppr', 'taas', 'tacr',
        'tamu', 'tapi', 'tasa', 'taso', 'trre',
    ],
    'Other Mollusca': [
        'accr', 'acgra', 'epba', 'gato', 'gype', 'laant', 'mosw', 'move', 'wiar',
    ],
    # 'memmem' = Membranipora membranacea (105 seqs); was 'meme', which matched
    # nothing. 'phaust' (the Phoronis australis rename target) is NOT listed:
    # the rename has never been applied, so the real Phoronis code is still
    # 'phau' and is registered in AMBIGUOUS_PREFIXES instead. Listing 'phaust'
    # here is what previously made the map *look* fixed while every phoronid
    # leaf still fell through to the gastropod 'phau' entry.
    'Other Lophotrochozoa': [
        'bune', 'crmu', 'crpa', 'lian', 'liun', 'lilo', 'memmem', 'noge', 'pece',
        'phov', 'phpa', 'teme',
    ],
}
for group, prefixes in _groups.items():
    for p in prefixes:
        PREFIX_TO_GROUP[p] = group


def get_taxon_group(leaf_name, code_map=None):
    """Determine taxonomic group from leaf name.

    Returns ``'Ambiguous'`` when the leaf's code is claimed by proteomes in more
    than one taxonomic group, rather than silently resolving it to whichever
    group happened to be registered. ``code_map`` is an optional
    :class:`species_code_lookup.CodeSpeciesMap`; when supplied, collisions it
    discovered in the data are honoured on top of ``AMBIGUOUS_PREFIXES``, so a
    NEW collision is refused the first time it appears rather than at the next
    time somebody edits this file.

    Why 'Ambiguous' and not a raise: ``species_for`` raises because its caller
    can fall back to printing the raw leaf name, so refusing costs one label.
    This function is called once per tip and recursively for every internal
    clade in ``_clade_color``; raising would abort the whole figure over a
    handful of tips, destroying output that is correct for the other 99%.
    There is no honest fallback colour, so the collision gets its own visible
    category and its own legend count instead.
    """
    if not leaf_name:
        return 'Unknown'
    if leaf_name.startswith('Berste') or leaf_name.startswith('TRINITY_'):
        return 'Berghia'
    name = leaf_name[4:] if leaf_name.startswith('ref_') else leaf_name
    extra = getattr(code_map, 'ambiguous', None) or {}

    def _group(code):
        return 'Ambiguous' if (code in AMBIGUOUS_PREFIXES or code in extra) \
            else PREFIX_TO_GROUP.get(code)

    candidates = set(PREFIX_TO_GROUP) | set(AMBIGUOUS_PREFIXES) | set(extra)
    for p in sorted(candidates, key=len, reverse=True):
        if name.startswith(p + '_') or name == p:
            return _group(p)
    return _group(name.split('_')[0]) or 'Unknown'


# --- PREFIX_TO_GROUP drift validation -------------------------------------
#
# PREFIX_TO_GROUP is hand-maintained against auto-generated header prefixes, so
# it drifts silently: get_taxon_group()'s final lookup has a 'Unknown' default,
# which means a stale key and a missing key both render a figure without ever
# failing a run. The functions below derive the truth from the proteomes
# themselves and refuse a map that no longer matches.

# references/nath_et_al/<lse|one_to_one_ortholog>/<group_dir>/<taxid>_<Genus_species>.faa
# The directory is the authoritative taxonomic group (that is how the reference
# set is organised); the header prefix inside the file is the code.
DIR_TO_GROUP = {
    'annelida': 'Annelida',
    'bivalvia': 'Bivalvia',
    'cephalopoda': 'Cephalopoda',
    'gastropoda': 'Gastropoda',
    'other_lophotrochozoan_phyla': 'Other Lophotrochozoa',
    'other_molluscan_classes': 'Other Mollusca',
    'platyhelminthes': 'Platyhelminthes',
}


class PrefixMapDriftError(AssertionError):
    """PREFIX_TO_GROUP no longer matches the reference proteomes on disk."""


class ObservedCode:
    """What the proteomes actually say about one header-prefix code."""

    def __init__(self, group, count, ambiguous_over, files):
        self.group = group                      # None when ambiguous
        self.count = count                      # headers carrying this code
        self.ambiguous_over = ambiguous_over    # set of groups claiming it
        self.files = files                      # basenames of claiming proteomes

    def __repr__(self):
        return (f'ObservedCode(group={self.group!r}, count={self.count}, '
                f'ambiguous_over={sorted(self.ambiguous_over)})')


def scan_reference_codes(references_dir=None):
    """Scan reference proteomes -> ``{code: ObservedCode}``.

    Only the taxon-organised ``nath_et_al`` tree is in scope. Outgroup and
    Swiss-Prot references use pipe-delimited headers (``AcCRa|...``,
    ``O14804|aminergic|...``) that carry no species code and belong to no
    molluscan group, so 'Unknown' is the correct answer for them.
    """
    root = Path(references_dir) if references_dir else PROJECT_ROOT / 'references' / 'nath_et_al'
    per_group = defaultdict(Counter)   # code -> {group: n}
    files = defaultdict(set)           # code -> {basename}
    for faa in sorted(glob.glob(str(Path(root) / '**' / '*.faa'), recursive=True)):
        group = DIR_TO_GROUP.get(os.path.basename(os.path.dirname(faa)))
        if group is None:
            continue
        with open(faa) as fh:
            for line in fh:
                if line.startswith('>'):
                    code = line[1:].split()[0].split('_')[0]
                    per_group[code][group] += 1
                    files[code].add(os.path.basename(faa))
    observed = {}
    for code, counts in per_group.items():
        groups = set(counts)
        observed[code] = ObservedCode(
            group=next(iter(groups)) if len(groups) == 1 else None,
            count=sum(counts.values()),
            ambiguous_over=groups if len(groups) > 1 else set(),
            files=files[code])
    return observed


def validate_prefix_map(references_dir=None, prefix_to_group=None, ambiguous=None):
    """Raise :class:`PrefixMapDriftError` unless the map matches the proteomes.

    Checks, in the order they have historically broken:
      1. dead keys      — a mapped code matching ZERO sequences (half-applied rename)
      2. missing codes  — a real code with no entry (renders as 'Unknown')
      3. wrong group    — a mapped code whose proteome sits in another group
      4. collisions     — a code claimed by two groups must be declared ambiguous
                          and must NOT also sit in the group map

    Not run at import: it needs ``references/`` on disk, which unit tests and
    lightweight consumers do not have. It runs at the entry points that actually
    render figures, and in the test suite against the real reference directory.
    """
    p2g = PREFIX_TO_GROUP if prefix_to_group is None else prefix_to_group
    amb = set(AMBIGUOUS_PREFIXES if ambiguous is None else ambiguous)
    observed = scan_reference_codes(references_dir)
    real_amb = {c for c, o in observed.items() if o.ambiguous_over}
    problems = []

    dead = sorted(set(p2g) - set(observed))
    if dead:
        problems.append(
            f'{len(dead)} mapped code(s) match ZERO sequences (stale or '
            f'half-applied rename): {dead}')

    missing = sorted(set(observed) - set(p2g) - amb)
    if missing:
        n = sum(observed[c].count for c in missing)
        problems.append(
            f'{len(missing)} real header code(s) are unmapped and would render '
            f'as "Unknown" ({n} sequences): '
            + ', '.join(f'{c}({observed[c].count})' for c in missing))

    wrong = [(c, p2g[c], observed[c].group) for c in sorted(set(p2g) & set(observed))
             if observed[c].group is not None and p2g[c] != observed[c].group]
    if wrong:
        problems.append('mis-grouped code(s): ' + ', '.join(
            f'{c} mapped {m!r} but proteome is {o!r}' for c, m, o in wrong))

    undeclared = sorted(real_amb - amb)
    if undeclared:
        problems.append('code(s) claimed by >1 taxonomic group but not declared '
                        'ambiguous: ' + ', '.join(
                            f'{c} {sorted(observed[c].ambiguous_over)}' for c in undeclared))

    landmines = sorted(amb & set(p2g))
    if landmines:
        problems.append(f'ambiguous code(s) still present in the group map, where '
                        f'they resolve to one phylum: {landmines}')

    stale_amb = sorted(amb - real_amb)
    if stale_amb:
        problems.append(f'code(s) declared ambiguous but no longer collide in the '
                        f'data (rename applied? update the register): {stale_amb}')

    if problems:
        raise PrefixMapDriftError(
            'PREFIX_TO_GROUP is out of sync with the reference proteomes:\n  - '
            + '\n  - '.join(problems))
    return observed


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
    for group in LEGEND_GROUPS:
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


def _validate_prefix_map_or_die(skip=False):
    """Refuse to render a figure from a map that has drifted from the data.

    Runs at the entry point rather than at import: the check needs
    ``references/`` on disk (~0.1 s to scan), which unit tests and lightweight
    consumers of this module do not have. A missing reference directory is not
    an error here — only a directory that CONTRADICTS the map is.
    """
    if skip:
        print('WARNING: PREFIX_TO_GROUP validation skipped; taxon colours and '
              'legend counts are unverified.', file=sys.stderr)
        return
    if not (PROJECT_ROOT / 'references' / 'nath_et_al').is_dir():
        return
    validate_prefix_map()


def main():
    parser = argparse.ArgumentParser(description='Visualize GPCR tree with taxonomy coloring')
    parser.add_argument('--tree', required=True, help='Newick tree file')
    parser.add_argument('--ranking', required=True, help='Ranked candidates CSV')
    parser.add_argument('--output', required=True, help='Output PNG path')
    parser.add_argument('--theme', choices=['light', 'dark'], default='light')
    parser.add_argument('--dpi', type=int, default=300, help='Output DPI')
    parser.add_argument('--skip-prefix-validation', action='store_true',
                        help='render even if PREFIX_TO_GROUP has drifted from '
                             'the reference proteomes (NOT for publication figures)')
    args = parser.parse_args()

    _validate_prefix_map_or_die(args.skip_prefix_validation)

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
