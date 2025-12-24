#!/usr/bin/env python3
# coexpression_analysis.py
# Purpose: Calculate co-expression between GPCR candidates and G-proteins.
# Identifies GPCRs co-expressed with olfactory G-proteins (Golf, Gi, Go) in chemosensory tissues.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python coexpression_analysis.py --expression <csv> --gproteins <csv>
#                                   --candidates <fasta> --output <csv>

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
from scipy.stats import spearmanr
from Bio import SeqIO


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Calculate GPCR-Gprotein co-expression'
    )
    parser.add_argument(
        '--expression', required=True,
        help='Expression summary CSV from process_expression.py'
    )
    parser.add_argument(
        '--gproteins', required=True,
        help='Classified G-proteins CSV from classify_gproteins.py'
    )
    parser.add_argument(
        '--candidates', required=True,
        help='GPCR candidate FASTA file'
    )
    parser.add_argument(
        '--chemosensory-tissues', required=True,
        help='Comma-separated list of chemosensory tissues'
    )
    parser.add_argument(
        '--olfactory-classes', default='Golf,Gi,Go',
        help='Comma-separated olfactory G-protein classes (default: Golf,Gi,Go)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output CSV file'
    )
    return parser.parse_args()


def load_expression_data(expr_file: str) -> pd.DataFrame:
    """
    Load expression summary data.

    Returns DataFrame with gene expression per tissue.
    """
    if not os.path.exists(expr_file):
        return pd.DataFrame()

    try:
        df = pd.read_csv(expr_file)
        return df
    except Exception as e:
        print(f"Warning: Could not load expression data: {e}", file=sys.stderr)
        return pd.DataFrame()


def load_gproteins(gprotein_file: str) -> pd.DataFrame:
    """
    Load classified G-proteins.

    Returns DataFrame with gene_id, species, gprotein_class, subtype.
    """
    if not os.path.exists(gprotein_file):
        return pd.DataFrame()

    try:
        df = pd.read_csv(gprotein_file)
        return df
    except Exception as e:
        print(f"Warning: Could not load G-protein data: {e}", file=sys.stderr)
        return pd.DataFrame()


def get_candidate_ids(fasta_file: str) -> List[str]:
    """Extract candidate IDs from FASTA file."""
    ids = []
    if os.path.exists(fasta_file):
        for record in SeqIO.parse(fasta_file, 'fasta'):
            ids.append(record.id)
    return ids


def get_tissue_tpm_columns(expr_df: pd.DataFrame) -> List[str]:
    """Find TPM columns in expression DataFrame."""
    return [c for c in expr_df.columns if c.endswith('_tpm')]


def calculate_expression_correlation(
    gpcr_id: str,
    gprotein_id: str,
    expr_df: pd.DataFrame,
    tpm_cols: List[str]
) -> Tuple[float, int]:
    """
    Calculate Spearman correlation between GPCR and G-protein expression.

    Returns: (correlation, n_tissues)
    """
    gpcr_row = expr_df[expr_df['gene_id'] == gpcr_id]
    gprotein_row = expr_df[expr_df['gene_id'] == gprotein_id]

    if gpcr_row.empty or gprotein_row.empty:
        return 0.0, 0

    # Get expression values
    gpcr_vals = gpcr_row[tpm_cols].values.flatten()
    gprotein_vals = gprotein_row[tpm_cols].values.flatten()

    # Remove NaN pairs
    mask = ~(np.isnan(gpcr_vals) | np.isnan(gprotein_vals))
    gpcr_clean = gpcr_vals[mask]
    gprotein_clean = gprotein_vals[mask]

    if len(gpcr_clean) < 3:
        return 0.0, len(gpcr_clean)

    try:
        corr, pval = spearmanr(gpcr_clean, gprotein_clean)
        if np.isnan(corr):
            return 0.0, len(gpcr_clean)
        return corr, len(gpcr_clean)
    except Exception:
        return 0.0, 0


def check_coexpression_in_tissue(
    gpcr_id: str,
    gprotein_id: str,
    tissue: str,
    expr_df: pd.DataFrame,
    min_tpm: float = 1.0
) -> bool:
    """
    Check if GPCR and G-protein are co-expressed in a tissue.

    Returns True if both are expressed above threshold.
    """
    tpm_col = f"{tissue}_tpm"
    if tpm_col not in expr_df.columns:
        return False

    gpcr_row = expr_df[expr_df['gene_id'] == gpcr_id]
    gprotein_row = expr_df[expr_df['gene_id'] == gprotein_id]

    if gpcr_row.empty or gprotein_row.empty:
        return False

    gpcr_tpm = gpcr_row[tpm_col].values[0]
    gprotein_tpm = gprotein_row[tpm_col].values[0]

    return (
        pd.notna(gpcr_tpm) and gpcr_tpm >= min_tpm and
        pd.notna(gprotein_tpm) and gprotein_tpm >= min_tpm
    )


def calculate_coexpression_score(
    gpcr_id: str,
    gprotein_df: pd.DataFrame,
    expr_df: pd.DataFrame,
    chemo_tissues: List[str],
    olfactory_classes: List[str]
) -> Dict:
    """
    Calculate co-expression score for a GPCR candidate.

    Prioritizes co-expression with olfactory G-proteins in chemosensory tissues.

    Returns dict with best G-protein match and score.
    """
    tpm_cols = get_tissue_tpm_columns(expr_df)

    best_result = {
        'gpcr_id': gpcr_id,
        'best_coexpr_gprotein': '',
        'gprotein_class': '',
        'correlation': 0.0,
        'coexpr_tissue': '',
        'coexpr_score': 0.0
    }

    best_score = 0.0

    # Get species of GPCR (if in expression data)
    gpcr_row = expr_df[expr_df['gene_id'] == gpcr_id]
    gpcr_species = gpcr_row['species'].values[0] if not gpcr_row.empty and 'species' in gpcr_row.columns else None

    # Filter G-proteins to same species if possible
    if gpcr_species:
        species_gproteins = gprotein_df[gprotein_df['species'] == gpcr_species]
    else:
        species_gproteins = gprotein_df

    for _, gp_row in species_gproteins.iterrows():
        gprotein_id = gp_row['gene_id']
        gp_class = gp_row['gprotein_class']
        is_olfactory = gp_class in olfactory_classes

        # Calculate correlation
        corr, n_tissues = calculate_expression_correlation(
            gpcr_id, gprotein_id, expr_df, tpm_cols
        )

        if n_tissues < 3:
            continue

        # Base score from correlation
        score = max(0, corr)  # Only positive correlations

        # Check co-expression in chemosensory tissues
        coexpr_tissue = ''
        for tissue in chemo_tissues:
            if check_coexpression_in_tissue(gpcr_id, gprotein_id, tissue, expr_df):
                coexpr_tissue = tissue
                score *= 1.5  # 50% bonus for chemosensory co-expression
                break

        # Bonus for olfactory G-protein class
        if is_olfactory:
            score *= 1.5  # 50% bonus for olfactory G-proteins

        if score > best_score:
            best_score = score
            best_result = {
                'gpcr_id': gpcr_id,
                'best_coexpr_gprotein': gprotein_id,
                'gprotein_class': gp_class,
                'correlation': corr,
                'coexpr_tissue': coexpr_tissue,
                'coexpr_score': score
            }

    return best_result


def main():
    """Main execution function."""
    args = parse_arguments()

    # Parse tissue and class lists
    chemo_tissues = [t.strip() for t in args.chemosensory_tissues.split(',')]
    olfactory_classes = [c.strip() for c in args.olfactory_classes.split(',')]

    print(f"Chemosensory tissues: {chemo_tissues}", file=sys.stderr)
    print(f"Olfactory G-protein classes: {olfactory_classes}", file=sys.stderr)

    # Load data
    expr_df = load_expression_data(args.expression)
    gprotein_df = load_gproteins(args.gproteins)
    candidate_ids = get_candidate_ids(args.candidates)

    if expr_df.empty:
        print("Warning: No expression data available", file=sys.stderr)
        pd.DataFrame(columns=[
            'gpcr_id', 'best_coexpr_gprotein', 'gprotein_class',
            'correlation', 'coexpr_tissue', 'coexpr_score'
        ]).to_csv(args.output, index=False)
        return

    if gprotein_df.empty:
        print("Warning: No G-protein data available", file=sys.stderr)
        pd.DataFrame(columns=[
            'gpcr_id', 'best_coexpr_gprotein', 'gprotein_class',
            'correlation', 'coexpr_tissue', 'coexpr_score'
        ]).to_csv(args.output, index=False)
        return

    print(f"Loaded expression data for {len(expr_df)} genes", file=sys.stderr)
    print(f"Loaded {len(gprotein_df)} classified G-proteins", file=sys.stderr)
    print(f"Processing {len(candidate_ids)} GPCR candidates", file=sys.stderr)

    # Calculate co-expression for each candidate
    results = []
    for i, cand_id in enumerate(candidate_ids):
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{len(candidate_ids)} candidates...",
                  file=sys.stderr)

        result = calculate_coexpression_score(
            cand_id, gprotein_df, expr_df, chemo_tissues, olfactory_classes
        )
        results.append(result)

    result_df = pd.DataFrame(results)

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write output
    result_df.to_csv(args.output, index=False)

    # Print summary
    with_coexpr = (result_df['coexpr_score'] > 0).sum()
    with_olfactory = result_df[result_df['gprotein_class'].isin(olfactory_classes)]

    print(f"\n=== Co-expression Analysis Summary ===", file=sys.stderr)
    print(f"Candidates with G-protein co-expression: {with_coexpr}/{len(result_df)}",
          file=sys.stderr)
    print(f"Candidates co-expressed with olfactory G-proteins: {len(with_olfactory)}",
          file=sys.stderr)
    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
