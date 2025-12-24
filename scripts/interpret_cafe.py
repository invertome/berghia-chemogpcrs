#!/usr/bin/env python3
# interpret_cafe.py
# Purpose: Interpret CAFE5 gene family expansion results with taxonomic context.
# Classifies expansions as Aeolid-specific, Nudibranch-wide, Gastropod-wide, etc.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python interpret_cafe.py --cafe-output <dir> --species-tree <tree>
#                            --lse-levels <levels> --pvalue-threshold <float>
#                            --output <csv>

import os
import sys
import argparse
import re
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from ete3 import Tree


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Interpret CAFE5 expansions with taxonomic context'
    )
    parser.add_argument(
        '--cafe-output', required=True,
        help='Directory with CAFE5 output files'
    )
    parser.add_argument(
        '--species-tree', required=True,
        help='Species tree file (Newick format)'
    )
    parser.add_argument(
        '--lse-levels', required=True,
        help='Space-separated LSE levels (e.g., "Aeolids:tax1,tax2 Nudibranchs:tax3")'
    )
    parser.add_argument(
        '--pvalue-threshold', type=float, default=0.05,
        help='P-value threshold for significant expansions (default: 0.05)'
    )
    parser.add_argument(
        '--berghia-id', default='taxid_berghia',
        help='Berghia taxon ID in tree (default: taxid_berghia)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output CSV file'
    )
    return parser.parse_args()


def parse_lse_levels(lse_string: str) -> Dict[str, List[str]]:
    """
    Parse LSE level specification.

    Format: "LevelName:taxid1,taxid2 LevelName2:taxid3,taxid4"

    Returns dict: level_name -> list of taxids
    """
    levels = {}

    for level_spec in lse_string.split():
        if ':' in level_spec:
            name, taxids = level_spec.split(':', 1)
            levels[name] = [t.strip() for t in taxids.split(',')]

    return levels


def load_species_tree(tree_file: str) -> Optional[Tree]:
    """Load species tree from Newick file."""
    if not os.path.exists(tree_file):
        print(f"Warning: Species tree not found: {tree_file}", file=sys.stderr)
        return None

    try:
        tree = Tree(tree_file)
        return tree
    except Exception as e:
        print(f"Warning: Could not load species tree: {e}", file=sys.stderr)
        return None


def parse_cafe_output(cafe_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse CAFE5 output files.

    Looks for:
    - Base_change.tab / Gamma_change.tab (expansion/contraction counts)
    - Base_family_results.txt (p-values per family)
    - Base_branch_probabilities.tab (branch-specific p-values)

    Returns:
        (family_df, branch_df) - DataFrames with family and branch-level results
    """
    cafe_path = Path(cafe_dir)
    family_results = []
    branch_results = []

    if not cafe_path.exists():
        print(f"Warning: CAFE output directory not found: {cafe_dir}", file=sys.stderr)
        return pd.DataFrame(), pd.DataFrame()

    # Parse family results (p-values)
    for result_file in cafe_path.glob('*family_results*'):
        try:
            with open(result_file) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        # Format varies, try to extract family ID and p-value
                        family_id = parts[0]
                        try:
                            pvalue = float(parts[-1])  # P-value usually last column
                            family_results.append({
                                'orthogroup': family_id,
                                'family_pvalue': pvalue
                            })
                        except ValueError:
                            continue
        except Exception as e:
            print(f"Warning: Could not parse {result_file}: {e}", file=sys.stderr)

    # Parse change files (expansion/contraction per branch)
    for change_file in cafe_path.glob('*change.tab'):
        try:
            df = pd.read_csv(change_file, sep='\t')
            # Extract branch-specific changes
            for _, row in df.iterrows():
                family_id = row.iloc[0]  # First column is family ID
                for col in df.columns[1:]:
                    change = row[col]
                    if pd.notna(change) and change != 0:
                        branch_results.append({
                            'orthogroup': family_id,
                            'branch': col,
                            'change': int(change)
                        })
        except Exception as e:
            print(f"Warning: Could not parse {change_file}: {e}", file=sys.stderr)

    # Parse branch probabilities if available
    for prob_file in cafe_path.glob('*branch_probabilities*'):
        try:
            df = pd.read_csv(prob_file, sep='\t')
            for _, row in df.iterrows():
                family_id = row.iloc[0]
                for col in df.columns[1:]:
                    prob = row[col]
                    if pd.notna(prob):
                        # Find matching branch result and add p-value
                        for br in branch_results:
                            if br['orthogroup'] == family_id and br['branch'] == col:
                                br['branch_pvalue'] = float(prob)
                                break
        except Exception as e:
            print(f"Warning: Could not parse {prob_file}: {e}", file=sys.stderr)

    family_df = pd.DataFrame(family_results)
    branch_df = pd.DataFrame(branch_results)

    print(f"Parsed {len(family_df)} families, {len(branch_df)} branch events",
          file=sys.stderr)

    return family_df, branch_df


def determine_expansion_level(
    branch_name: str,
    species_tree: Optional[Tree],
    lse_levels: Dict[str, List[str]],
    berghia_id: str
) -> str:
    """
    Determine taxonomic level of expansion based on branch.

    Returns level name (e.g., "Aeolid-specific", "Nudibranch-wide").
    """
    if not species_tree:
        return "Unknown"

    # Try to find the branch in the tree
    # Branch names may be node names or leaf combinations

    # Check if branch name matches a known taxon
    for level_name, taxids in lse_levels.items():
        if branch_name in taxids:
            return f"{level_name}-specific"

    # Check if branch is an ancestor of specific groups
    try:
        # Find Berghia in tree
        berghia_node = None
        for node in species_tree.traverse():
            if berghia_id in node.name:
                berghia_node = node
                break

        if berghia_node:
            # Check each level
            for level_name, taxids in lse_levels.items():
                level_leaves = []
                for taxid in taxids:
                    for leaf in species_tree.get_leaves():
                        if taxid in leaf.name:
                            level_leaves.append(leaf)

                if level_leaves:
                    # Get ancestor of this level
                    try:
                        ancestor = species_tree.get_common_ancestor(level_leaves)
                        if branch_name in ancestor.name or ancestor.name in branch_name:
                            return f"{level_name}-wide"
                    except Exception:
                        continue

    except Exception as e:
        print(f"Warning: Could not determine level for branch {branch_name}: {e}",
              file=sys.stderr)

    # Check if it's a Berghia-specific expansion
    if berghia_id in branch_name:
        return "Berghia-specific"

    return "Other"


def calculate_expansion_metrics(
    family_df: pd.DataFrame,
    branch_df: pd.DataFrame,
    species_tree: Optional[Tree],
    lse_levels: Dict[str, List[str]],
    berghia_id: str,
    pvalue_threshold: float
) -> pd.DataFrame:
    """
    Calculate expansion metrics for each orthogroup.

    Returns DataFrame with expansion interpretation.
    """
    results = []

    # Group branch events by orthogroup
    for og in family_df['orthogroup'].unique():
        og_family = family_df[family_df['orthogroup'] == og]
        og_branches = branch_df[branch_df['orthogroup'] == og]

        family_pvalue = og_family['family_pvalue'].values[0] if not og_family.empty else 1.0

        if family_pvalue >= pvalue_threshold:
            continue  # Skip non-significant families

        # Find significant expansion branches
        expansions = og_branches[og_branches['change'] > 0]

        if expansions.empty:
            continue

        # Find the most significant expansion
        best_expansion = None
        best_pvalue = 1.0
        best_change = 0

        for _, exp in expansions.iterrows():
            branch_pval = exp.get('branch_pvalue', family_pvalue)
            if branch_pval < best_pvalue or (branch_pval == best_pvalue and exp['change'] > best_change):
                best_pvalue = branch_pval
                best_change = exp['change']
                best_expansion = exp

        if best_expansion is None:
            continue

        # Determine taxonomic level
        branch_name = best_expansion['branch']
        tax_level = determine_expansion_level(
            branch_name, species_tree, lse_levels, berghia_id
        )

        # Get Berghia copy count
        berghia_branches = og_branches[og_branches['branch'].str.contains(berghia_id, na=False)]
        berghia_copies = berghia_branches['change'].sum() if not berghia_branches.empty else 0

        # Calculate expansion fold (approximate)
        total_expansion = expansions['change'].sum()
        ancestral_estimate = max(1, total_expansion - best_change)  # Rough estimate
        expansion_fold = (ancestral_estimate + best_change) / max(1, ancestral_estimate)

        results.append({
            'orthogroup': og,
            'expansion_pvalue': best_pvalue,
            'expansion_branch': branch_name,
            'taxonomic_level': tax_level,
            'berghia_copies': int(berghia_copies),
            'ancestral_copies': int(ancestral_estimate),
            'expansion_fold': round(expansion_fold, 2),
            'total_expansion_events': len(expansions)
        })

    return pd.DataFrame(results)


def main():
    """Main execution function."""
    args = parse_arguments()

    # Parse LSE levels
    lse_levels = parse_lse_levels(args.lse_levels)
    print(f"LSE levels: {list(lse_levels.keys())}", file=sys.stderr)

    # Load species tree
    species_tree = load_species_tree(args.species_tree)

    # Parse CAFE output
    family_df, branch_df = parse_cafe_output(args.cafe_output)

    if family_df.empty:
        print("Warning: No CAFE results found", file=sys.stderr)
        pd.DataFrame(columns=[
            'orthogroup', 'expansion_pvalue', 'expansion_branch', 'taxonomic_level',
            'berghia_copies', 'ancestral_copies', 'expansion_fold', 'total_expansion_events'
        ]).to_csv(args.output, index=False)
        return

    # Calculate expansion metrics
    result_df = calculate_expansion_metrics(
        family_df, branch_df, species_tree, lse_levels,
        args.berghia_id, args.pvalue_threshold
    )

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write output
    result_df.to_csv(args.output, index=False)

    # Print summary
    print(f"\n=== CAFE Expansion Interpretation Summary ===", file=sys.stderr)
    print(f"Significant expansions: {len(result_df)}", file=sys.stderr)

    if not result_df.empty:
        print(f"By taxonomic level:", file=sys.stderr)
        for level, count in result_df['taxonomic_level'].value_counts().items():
            print(f"  {level}: {count}", file=sys.stderr)

    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
