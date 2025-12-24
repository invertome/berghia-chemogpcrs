#!/usr/bin/env python3
# process_expression.py
# Purpose: Process Salmon quant.sf files into per-gene, per-tissue TPM matrix with tau index.
# Calculates tissue specificity (tau) and chemosensory enrichment for ranking.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python process_expression.py --quant-dir <dir> --chemosensory-tissues <tissues>
#                                --other-tissues <tissues> --min-tpm <threshold>
#                                --output <output.csv>

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Process Salmon quant.sf files for expression analysis'
    )
    parser.add_argument(
        '--quant-dir', required=True,
        help='Directory with Salmon quant.sf outputs organized as: species/tissue/replicate/quant.sf'
    )
    parser.add_argument(
        '--chemosensory-tissues', required=True,
        help='Comma-separated list of chemosensory tissues (e.g., "rhinophore,oral_veil,tentacle")'
    )
    parser.add_argument(
        '--other-tissues', required=True,
        help='Comma-separated list of non-chemosensory tissues for tau calculation'
    )
    parser.add_argument(
        '--min-tpm', type=float, default=1.0,
        help='Minimum TPM for "expressed" (default: 1.0)'
    )
    parser.add_argument(
        '--tau-threshold', type=float, default=0.8,
        help='Tau index threshold for tissue-specific (default: 0.8)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output CSV file path'
    )
    return parser.parse_args()


def find_quant_files(quant_dir: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Find all quant.sf files in the directory structure.

    Expected structure: {quant_dir}/{species}/{tissue}/{replicate}/quant.sf

    Returns:
        Dict mapping species -> tissue -> list of quant.sf paths
    """
    quant_files = defaultdict(lambda: defaultdict(list))
    quant_path = Path(quant_dir)

    if not quant_path.exists():
        print(f"Warning: Quant directory not found: {quant_dir}", file=sys.stderr)
        return quant_files

    # Find all quant.sf files
    for qf in quant_path.rglob('quant.sf'):
        # Parse path components
        rel_path = qf.relative_to(quant_path)
        parts = rel_path.parts

        if len(parts) >= 3:
            # species/tissue/replicate/quant.sf
            species = parts[0]
            tissue = parts[1]
            quant_files[species][tissue].append(str(qf))
        elif len(parts) >= 2:
            # tissue/replicate/quant.sf (single species)
            species = "default"
            tissue = parts[0]
            quant_files[species][tissue].append(str(qf))

    n_files = sum(len(v) for d in quant_files.values() for v in d.values())
    print(f"Found {n_files} quant.sf files across {len(quant_files)} species",
          file=sys.stderr)

    return quant_files


def load_salmon_quant(quant_file: str) -> pd.DataFrame:
    """
    Load a single Salmon quant.sf file.

    Returns:
        DataFrame with columns: Name, TPM
    """
    try:
        df = pd.read_csv(quant_file, sep='\t')
        return df[['Name', 'TPM']]
    except Exception as e:
        print(f"Warning: Could not load {quant_file}: {e}", file=sys.stderr)
        return pd.DataFrame(columns=['Name', 'TPM'])


def aggregate_replicates(quant_files: List[str]) -> pd.DataFrame:
    """
    Aggregate TPM values across replicate quant.sf files.

    Takes the mean TPM across replicates.

    Returns:
        DataFrame with columns: gene_id, TPM
    """
    all_dfs = []
    for qf in quant_files:
        df = load_salmon_quant(qf)
        if not df.empty:
            all_dfs.append(df)

    if not all_dfs:
        return pd.DataFrame(columns=['gene_id', 'TPM'])

    # Concatenate and aggregate
    combined = pd.concat(all_dfs, ignore_index=True)
    aggregated = combined.groupby('Name')['TPM'].mean().reset_index()
    aggregated.columns = ['gene_id', 'TPM']

    return aggregated


def build_expression_matrix(quant_files: Dict[str, Dict[str, List[str]]]) -> pd.DataFrame:
    """
    Build a gene x tissue expression matrix.

    Returns:
        DataFrame with gene_id as index and tissue columns containing TPM values
    """
    all_data = []

    for species, tissues in quant_files.items():
        for tissue, files in tissues.items():
            tissue_df = aggregate_replicates(files)
            if not tissue_df.empty:
                tissue_df['species'] = species
                tissue_df['tissue'] = tissue
                all_data.append(tissue_df)

    if not all_data:
        return pd.DataFrame()

    # Combine all tissue data
    combined = pd.concat(all_data, ignore_index=True)

    # Pivot to create matrix
    matrix = combined.pivot_table(
        index=['gene_id', 'species'],
        columns='tissue',
        values='TPM',
        aggfunc='mean'
    ).reset_index()

    return matrix


def calculate_tau_index(tpm_row: pd.Series, tissue_cols: List[str]) -> float:
    """
    Calculate the tau index (tissue specificity) for a gene.

    Tau = sum(1 - x_i/max(x)) / (n - 1)

    Where x_i is expression in tissue i, n is number of tissues.
    Tau = 1 means perfectly tissue-specific.
    Tau = 0 means uniformly expressed.

    Args:
        tpm_row: Series of TPM values
        tissue_cols: List of tissue column names

    Returns:
        Tau index (0 to 1)
    """
    values = tpm_row[tissue_cols].dropna().values

    if len(values) < 2:
        return np.nan

    max_val = np.max(values)
    if max_val == 0:
        return 0.0  # Not expressed anywhere

    # Normalize and calculate tau
    normalized = values / max_val
    tau = np.sum(1 - normalized) / (len(values) - 1)

    return tau


def calculate_chemosensory_enrichment(tpm_row: pd.Series,
                                       chemo_tissues: List[str],
                                       other_tissues: List[str]) -> Tuple[float, float]:
    """
    Calculate chemosensory tissue enrichment.

    Returns fold-change of mean TPM in chemosensory vs non-chemosensory tissues.

    Args:
        tpm_row: Series of TPM values
        chemo_tissues: List of chemosensory tissue names
        other_tissues: List of non-chemosensory tissue names

    Returns:
        Tuple of (fold_enrichment, mean_chemosensory_tpm)
    """
    # Get available chemosensory tissues
    avail_chemo = [t for t in chemo_tissues if t in tpm_row.index and pd.notna(tpm_row[t])]
    avail_other = [t for t in other_tissues if t in tpm_row.index and pd.notna(tpm_row[t])]

    if not avail_chemo:
        return 0.0, 0.0

    mean_chemo = np.mean([tpm_row[t] for t in avail_chemo])

    if not avail_other:
        # No other tissues to compare, return chemosensory mean as enrichment
        return mean_chemo, mean_chemo

    mean_other = np.mean([tpm_row[t] for t in avail_other])

    # Calculate fold enrichment (add pseudocount to avoid division by zero)
    fold_enrichment = (mean_chemo + 0.1) / (mean_other + 0.1)

    return fold_enrichment, mean_chemo


def process_expression_data(quant_dir: str,
                           chemo_tissues: List[str],
                           other_tissues: List[str],
                           min_tpm: float,
                           tau_threshold: float) -> pd.DataFrame:
    """
    Main processing function to generate expression summary.

    Args:
        quant_dir: Directory with Salmon quant.sf files
        chemo_tissues: List of chemosensory tissue names
        other_tissues: List of non-chemosensory tissue names
        min_tpm: Minimum TPM threshold
        tau_threshold: Tau index threshold for tissue-specific

    Returns:
        DataFrame with expression summary
    """
    # Find and load quant files
    quant_files = find_quant_files(quant_dir)

    if not quant_files:
        print("No expression data found", file=sys.stderr)
        return pd.DataFrame()

    # Build expression matrix
    expr_matrix = build_expression_matrix(quant_files)

    if expr_matrix.empty:
        print("Could not build expression matrix", file=sys.stderr)
        return pd.DataFrame()

    print(f"Built expression matrix: {len(expr_matrix)} genes x {len(expr_matrix.columns)-2} tissues",
          file=sys.stderr)

    # Get tissue columns
    tissue_cols = [c for c in expr_matrix.columns if c not in ['gene_id', 'species']]
    all_tissues = list(set(chemo_tissues + other_tissues))
    available_tissues = [t for t in all_tissues if t in tissue_cols]

    print(f"Available tissues: {available_tissues}", file=sys.stderr)
    print(f"Chemosensory tissues found: {[t for t in chemo_tissues if t in tissue_cols]}",
          file=sys.stderr)

    # Calculate metrics for each gene
    results = []

    for idx, row in expr_matrix.iterrows():
        gene_id = row['gene_id']
        species = row['species']

        # Calculate tau index
        tau = calculate_tau_index(row, tissue_cols)

        # Calculate chemosensory enrichment
        enrichment, mean_chemo_tpm = calculate_chemosensory_enrichment(
            row, chemo_tissues, other_tissues
        )

        # Check if expressed in chemosensory tissues
        expressed_in_chemo = any(
            row.get(t, 0) >= min_tpm for t in chemo_tissues if t in row.index
        )

        # Check if tissue-specific to chemosensory tissues
        is_chemo_specific = False
        if tau >= tau_threshold and expressed_in_chemo:
            # Find which tissue has max expression
            max_tissue = None
            max_tpm = 0
            for t in tissue_cols:
                if pd.notna(row.get(t)) and row[t] > max_tpm:
                    max_tpm = row[t]
                    max_tissue = t
            if max_tissue in chemo_tissues:
                is_chemo_specific = True

        result = {
            'gene_id': gene_id,
            'species': species,
            'tau_index': tau,
            'chemosensory_enrichment': enrichment,
            'mean_chemosensory_tpm': mean_chemo_tpm,
            'expressed_in_chemosensory': expressed_in_chemo,
            'is_chemosensory_specific': is_chemo_specific
        }

        # Add individual tissue TPMs
        for t in tissue_cols:
            result[f'{t}_tpm'] = row.get(t, np.nan)

        results.append(result)

    return pd.DataFrame(results)


def main():
    """Main execution function."""
    args = parse_arguments()

    # Parse tissue lists
    chemo_tissues = [t.strip() for t in args.chemosensory_tissues.split(',')]
    other_tissues = [t.strip() for t in args.other_tissues.split(',')]

    print(f"Chemosensory tissues: {chemo_tissues}", file=sys.stderr)
    print(f"Other tissues: {other_tissues}", file=sys.stderr)

    # Process expression data
    result_df = process_expression_data(
        args.quant_dir,
        chemo_tissues,
        other_tissues,
        args.min_tpm,
        args.tau_threshold
    )

    if result_df.empty:
        print("Warning: No expression data processed", file=sys.stderr)
        # Create empty output with expected columns
        result_df = pd.DataFrame(columns=[
            'gene_id', 'species', 'tau_index', 'chemosensory_enrichment',
            'mean_chemosensory_tpm', 'expressed_in_chemosensory', 'is_chemosensory_specific'
        ])

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write output
    result_df.to_csv(args.output, index=False)

    # Print summary
    print(f"\n=== Expression Analysis Summary ===", file=sys.stderr)
    print(f"Total genes processed: {len(result_df)}", file=sys.stderr)
    print(f"Genes expressed in chemosensory tissues: {result_df['expressed_in_chemosensory'].sum()}",
          file=sys.stderr)
    print(f"Genes chemosensory-specific (tau >= {args.tau_threshold}): "
          f"{result_df['is_chemosensory_specific'].sum()}", file=sys.stderr)
    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
