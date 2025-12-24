#!/usr/bin/env python3
# divergence_dating.py
# Purpose: Estimate divergence times for GPCR duplications using molecular clock methods.
# Inputs: Ultrametric tree, calibration points (optional)
# Outputs: Dated tree with divergence time estimates
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python divergence_dating.py <ultrametric_tree> <output_prefix> [calibration_file]
#
# Arguments:
#   ultrametric_tree   - Ultrametric tree from chronos (Newick format)
#   output_prefix      - Prefix for output files
#   calibration_file   - Optional: TSV with calibration points (node, min_age, max_age)

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Conditional ete3 import - may fail on Python 3.13+ due to removed cgi module
try:
    from ete3 import Tree, TreeStyle, NodeStyle, faces
    ETE3_AVAILABLE = True
except ImportError as e:
    print(f"Warning: ete3 not available ({e}). Tree visualization disabled.",
          file=sys.stderr)
    Tree = TreeStyle = NodeStyle = faces = None
    ETE3_AVAILABLE = False

# Geological time scale (Ma = million years ago)
GEOLOGICAL_PERIODS = {
    'Cambrian': (538.8, 485.4),
    'Ordovician': (485.4, 443.8),
    'Silurian': (443.8, 419.2),
    'Devonian': (419.2, 358.9),
    'Carboniferous': (358.9, 298.9),
    'Permian': (298.9, 251.9),
    'Triassic': (251.9, 201.4),
    'Jurassic': (201.4, 145.0),
    'Cretaceous': (145.0, 66.0),
    'Paleogene': (66.0, 23.0),
    'Neogene': (23.0, 2.6),
    'Quaternary': (2.6, 0)
}

# Known gastropod/mollusc calibration points (from fossil record)
DEFAULT_CALIBRATIONS = {
    'Gastropoda_root': {'min': 500, 'max': 540, 'source': 'Cambrian gastropod fossils'},
    'Heterobranchia_root': {'min': 200, 'max': 250, 'source': 'Jurassic heterobranch fossils'},
    'Nudibranchia_root': {'min': 100, 'max': 150, 'source': 'Cretaceous nudibranch estimates'},
    'Aeolidida_root': {'min': 50, 'max': 100, 'source': 'Cenozoic aeolid estimates'}
}


def load_calibrations(calibration_file: str) -> Dict:
    """
    Load calibration points from file.

    Expected format (TSV):
    node_name\tmin_age\tmax_age\tsource

    Args:
        calibration_file: Path to calibration file

    Returns:
        Dictionary of calibration points
    """
    calibrations = {}

    if not os.path.exists(calibration_file):
        return calibrations

    try:
        df = pd.read_csv(calibration_file, sep='\t')
        for _, row in df.iterrows():
            node_name = str(row.iloc[0])
            calibrations[node_name] = {
                'min': float(row.iloc[1]),
                'max': float(row.iloc[2]),
                'source': str(row.iloc[3]) if len(row) > 3 else 'user-defined'
            }
        print(f"Loaded {len(calibrations)} calibration points", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load calibrations: {e}", file=sys.stderr)

    return calibrations


def estimate_absolute_ages(tree: Tree,
                          calibrations: Dict,
                          root_age: float = 500.0) -> Tuple[Tree, Dict]:
    """
    Estimate absolute ages for tree nodes using calibration points.

    Uses linear interpolation based on relative branch lengths.

    Args:
        tree: Ultrametric tree (relative ages)
        calibrations: Calibration points with min/max ages
        root_age: Default root age if no calibration available (Ma)

    Returns:
        Tuple of (dated tree, node ages dictionary)
    """
    # Get tree height (root to tips in relative units)
    # Guard against empty tree
    leaves = list(tree.iter_leaves())
    if not leaves:
        print("Warning: Tree has no leaves, cannot estimate ages", file=sys.stderr)
        return tree, {}
    tree_height = max(tree.get_distance(leaf) for leaf in leaves)

    # Find applicable calibration points
    applied_calibrations = []

    for node_name, calib in calibrations.items():
        nodes = tree.search_nodes(name=node_name)
        if nodes:
            node = nodes[0]
            rel_depth = tree_height - node.get_distance(tree)
            applied_calibrations.append({
                'node': node,
                'rel_depth': rel_depth,
                'min_age': calib['min'],
                'max_age': calib['max'],
                'mean_age': (calib['min'] + calib['max']) / 2
            })

    # If no calibrations found, use default root age
    if not applied_calibrations:
        print(f"No calibration points matched. Using default root age: {root_age} Ma",
              file=sys.stderr)
        scale_factor = root_age / tree_height if tree_height > 0 else 1
    else:
        # Use regression to estimate scale factor
        rel_depths = [c['rel_depth'] for c in applied_calibrations]
        mean_ages = [c['mean_age'] for c in applied_calibrations]

        if len(applied_calibrations) >= 2:
            # Linear regression
            slope, intercept = np.polyfit(rel_depths, mean_ages, 1)
            scale_factor = slope
        else:
            # Single calibration point
            c = applied_calibrations[0]
            scale_factor = c['mean_age'] / c['rel_depth'] if c['rel_depth'] > 0 else root_age

    # Calculate absolute ages for all nodes
    node_ages = {}

    for node in tree.traverse():
        rel_depth = tree_height - node.get_distance(tree)
        abs_age = rel_depth * scale_factor

        node_ages[node.name if node.name else f"node_{id(node)}"] = {
            'relative_depth': rel_depth,
            'absolute_age_ma': abs_age,
            'geological_period': get_geological_period(abs_age),
            'is_leaf': node.is_leaf()
        }

        # Add age to node features
        node.add_feature('age_ma', abs_age)
        node.add_feature('geological_period', get_geological_period(abs_age))

    return tree, node_ages


def get_geological_period(age_ma: float) -> str:
    """Get geological period for an age in million years."""
    for period, (start, end) in GEOLOGICAL_PERIODS.items():
        if end <= age_ma <= start:
            return period
    if age_ma > 538.8:
        return 'Precambrian'
    return 'Recent'


def identify_duplication_ages(tree: Tree, node_ages: Dict) -> List[Dict]:
    """
    Identify and date duplication events in the tree.

    Duplications are identified as nodes where children share species.

    Args:
        tree: Dated tree
        node_ages: Dictionary of node ages

    Returns:
        List of duplication events with ages
    """
    duplications = []

    for node in tree.traverse():
        if node.is_leaf():
            continue

        # Check if this could be a duplication (children share taxa)
        child_taxa = []
        for child in node.children:
            taxa = set()
            for leaf in child:
                # Extract taxon from leaf name
                parts = leaf.name.split('_')
                if len(parts) >= 2:
                    taxa.add(parts[0])  # Assume first part is taxon
            child_taxa.append(taxa)

        # Check for overlap (shared taxa = likely duplication)
        if len(child_taxa) >= 2:
            overlap = child_taxa[0] & child_taxa[1]
            if overlap:
                node_name = node.name if node.name else f"node_{id(node)}"
                age_info = node_ages.get(node_name, {})

                duplications.append({
                    'node': node_name,
                    'age_ma': age_info.get('absolute_age_ma', 0),
                    'geological_period': age_info.get('geological_period', 'Unknown'),
                    'shared_taxa': list(overlap),
                    'n_descendants': len(list(node.iter_leaves())),
                    'support': getattr(node, 'support', None)
                })

    return duplications


def create_dated_tree_visualization(tree: Tree, output_file: str):
    """
    Create a visualization of the dated tree.

    Args:
        tree: Dated tree with age_ma feature
        output_file: Output file path (PDF/PNG)
    """
    try:
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = False
        ts.show_branch_support = True
        ts.scale = 50  # pixels per branch length unit

        # Color nodes by geological period
        period_colors = {
            'Cambrian': '#7fc97f',
            'Ordovician': '#beaed4',
            'Silurian': '#fdc086',
            'Devonian': '#ffff99',
            'Carboniferous': '#386cb0',
            'Permian': '#f0027f',
            'Triassic': '#bf5b17',
            'Jurassic': '#666666',
            'Cretaceous': '#a6cee3',
            'Paleogene': '#1f78b4',
            'Neogene': '#b2df8a',
            'Quaternary': '#33a02c',
            'Recent': '#fb9a99',
            'Precambrian': '#e31a1c'
        }

        for node in tree.traverse():
            nstyle = NodeStyle()

            if hasattr(node, 'geological_period'):
                period = node.geological_period
                color = period_colors.get(period, '#888888')
                nstyle['fgcolor'] = color

                if not node.is_leaf():
                    nstyle['size'] = 8
                    # Add age label
                    if hasattr(node, 'age_ma'):
                        age_face = faces.TextFace(f"{node.age_ma:.1f} Ma", fsize=6)
                        node.add_face(age_face, column=0, position="branch-top")

            node.set_style(nstyle)

        # Add title
        ts.title.add_face(faces.TextFace("Dated Phylogeny (Ma)", fsize=14), column=0)

        tree.render(output_file, tree_style=ts, w=1200)
        print(f"Tree visualization saved to: {output_file}", file=sys.stderr)

    except Exception as e:
        print(f"Warning: Could not create visualization: {e}", file=sys.stderr)


def main():
    """Main execution function."""
    if len(sys.argv) < 3:
        print("Usage: python divergence_dating.py <ultrametric_tree> <output_prefix> [calibration_file]")
        sys.exit(1)

    tree_file = sys.argv[1]
    output_prefix = sys.argv[2]
    calibration_file = sys.argv[3] if len(sys.argv) > 3 else None

    # Validate inputs
    if not os.path.exists(tree_file):
        print(f"Error: Tree file not found: {tree_file}", file=sys.stderr)
        sys.exit(1)

    # Load tree
    print(f"Loading tree: {tree_file}", file=sys.stderr)
    tree = Tree(tree_file, format=1)
    n_leaves = len(tree)
    print(f"Tree has {n_leaves} leaves", file=sys.stderr)

    # Load calibrations
    calibrations = DEFAULT_CALIBRATIONS.copy()
    if calibration_file:
        user_calibrations = load_calibrations(calibration_file)
        calibrations.update(user_calibrations)

    print(f"Using {len(calibrations)} calibration points", file=sys.stderr)

    # Estimate absolute ages
    print("\nEstimating divergence times...", file=sys.stderr)
    dated_tree, node_ages = estimate_absolute_ages(tree, calibrations)

    # Identify duplication events
    print("Identifying duplication events...", file=sys.stderr)
    duplications = identify_duplication_ages(dated_tree, node_ages)

    # Create output directory
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Write Outputs ---

    # 1. Dated tree in Newick format
    tree_output = f"{output_prefix}_dated.tre"
    dated_tree.write(outfile=tree_output, format=1)
    print(f"\nDated tree written to: {tree_output}", file=sys.stderr)

    # 2. Node ages table
    ages_df = pd.DataFrame([
        {'node': name, **info}
        for name, info in node_ages.items()
    ])
    ages_file = f"{output_prefix}_node_ages.tsv"
    ages_df.to_csv(ages_file, sep='\t', index=False)
    print(f"Node ages written to: {ages_file}", file=sys.stderr)

    # 3. Duplication events
    if duplications:
        dup_df = pd.DataFrame(duplications)
        dup_df = dup_df.sort_values('age_ma', ascending=False)
        dup_file = f"{output_prefix}_duplications.tsv"
        dup_df.to_csv(dup_file, sep='\t', index=False)
        print(f"Duplication events written to: {dup_file}", file=sys.stderr)

    # 4. Summary
    # Get age range (with guards for empty lists)
    leaf_ages = [info['absolute_age_ma'] for info in node_ages.values() if info.get('is_leaf', False)]
    internal_ages = [info['absolute_age_ma'] for info in node_ages.values() if not info.get('is_leaf', True)]

    summary = {
        'tree_file': tree_file,
        'n_leaves': n_leaves,
        'n_internal_nodes': len(internal_ages),
        'root_age_ma': max(internal_ages) if internal_ages else 0,
        'crown_age_ma': min(internal_ages) if internal_ages else 0,
        'n_duplications_detected': len(duplications),
        'duplication_age_range_ma': {
            'oldest': max(d['age_ma'] for d in duplications) if duplications else 0,
            'youngest': min(d['age_ma'] for d in duplications) if duplications else 0
        },
        'calibrations_used': list(calibrations.keys()),
        'geological_periods_spanned': list(set(info['geological_period'] for info in node_ages.values()))
    }

    summary_file = f"{output_prefix}_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Summary written to: {summary_file}", file=sys.stderr)

    # 5. Visualization (if possible)
    try:
        viz_file = f"{output_prefix}_dated_tree.pdf"
        create_dated_tree_visualization(dated_tree, viz_file)
    except Exception as e:
        print(f"Visualization skipped: {e}", file=sys.stderr)

    # Print summary
    print("\n=== Divergence Dating Summary ===", file=sys.stderr)
    print(f"Root age: {summary['root_age_ma']:.1f} Ma ({get_geological_period(summary['root_age_ma'])})",
          file=sys.stderr)
    print(f"Crown age: {summary['crown_age_ma']:.1f} Ma", file=sys.stderr)
    print(f"Duplications detected: {summary['n_duplications_detected']}", file=sys.stderr)

    if duplications:
        print("\nDuplication events by geological period:", file=sys.stderr)
        period_counts = {}
        for dup in duplications:
            period = dup['geological_period']
            period_counts[period] = period_counts.get(period, 0) + 1
        for period, count in sorted(period_counts.items(), key=lambda x: -x[1]):
            print(f"  {period}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
