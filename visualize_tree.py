#!/usr/bin/env python3
# visualize_tree.py
# Purpose: Visualize phylogenetic trees with multiple layouts and taxonomic coloring.
# Inputs: Tree file ($1), output prefix ($2)
# Outputs: Basic (${output_prefix}_basic.png), colored (${output_prefix}_colored.png), circular (${output_prefix}_circular.png) plots
# Logic: Plots basic tree, colors by taxid/type, and provides circular layout.
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces
from collections import defaultdict

tree_file = sys.argv[1]
output_prefix = sys.argv[2]

t = Tree(tree_file)

# Define color palette for different sequence types
COLORS = {
    'ref': '#1f77b4',      # Blue for reference sequences
    'berghia': '#2ca02c',  # Green for Berghia
    'taxid1': '#ff7f0e',   # Orange
    'taxid2': '#d62728',   # Red
    'default': '#7f7f7f'   # Gray for unknown
}


def get_leaf_color(leaf_name):
    """Determine color based on leaf name prefix."""
    if leaf_name.startswith('ref_'):
        return COLORS['ref']
    elif 'berghia' in leaf_name.lower():
        return COLORS['berghia']
    else:
        # Try to extract taxid prefix
        prefix = leaf_name.split('_')[0] if '_' in leaf_name else leaf_name
        if prefix in COLORS:
            return COLORS[prefix]
        return COLORS['default']


def get_leaf_category(leaf_name):
    """Categorize leaf for legend."""
    if leaf_name.startswith('ref_'):
        return 'Reference'
    elif 'berghia' in leaf_name.lower():
        return 'Berghia'
    else:
        return 'Other taxa'


# --- Basic tree plot ---
ts_basic = TreeStyle()
ts_basic.show_leaf_name = True
ts_basic.scale = 100
t.render(f"{output_prefix}_basic.png", w=800, units='px', tree_style=ts_basic)

# --- Colored tree plot (by taxon type) ---
ts_colored = TreeStyle()
ts_colored.show_leaf_name = True
ts_colored.scale = 100
ts_colored.legend_position = 1

# Count categories for legend
category_counts = defaultdict(int)

for leaf in t:
    # Create node style with color
    nstyle = NodeStyle()
    color = get_leaf_color(leaf.name)
    nstyle['fgcolor'] = color
    nstyle['size'] = 6
    leaf.set_style(nstyle)

    # Track categories
    category = get_leaf_category(leaf.name)
    category_counts[category] += 1

# Add legend
legend_colors = {
    'Reference': COLORS['ref'],
    'Berghia': COLORS['berghia'],
    'Other taxa': COLORS['default']
}

for category, color in legend_colors.items():
    if category_counts[category] > 0:
        ts_colored.legend.add_face(
            faces.CircleFace(5, color),
            column=0
        )
        ts_colored.legend.add_face(
            TextFace(f" {category} ({category_counts[category]})", fsize=10),
            column=1
        )

t.render(f"{output_prefix}_colored.png", w=1000, units='px', tree_style=ts_colored)

# --- Circular layout ---
ts_circular = TreeStyle()
ts_circular.mode = "c"
ts_circular.show_leaf_name = True
ts_circular.scale = 50

# Apply same coloring for circular view
for leaf in t:
    nstyle = NodeStyle()
    color = get_leaf_color(leaf.name)
    nstyle['fgcolor'] = color
    nstyle['size'] = 4
    leaf.set_style(nstyle)

t.render(f"{output_prefix}_circular.png", w=1000, units='px', tree_style=ts_circular)

print(f"Generated tree visualizations: {output_prefix}_basic.png, {output_prefix}_colored.png, {output_prefix}_circular.png")
