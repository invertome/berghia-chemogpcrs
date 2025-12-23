#!/usr/bin/env python3
# expression_analysis.py
# Purpose: Analyze expression data from Salmon quant.sf files for GPCR candidates.
# Calculates tissue specificity (tau index) and identifies chemosensory-enriched candidates.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python expression_analysis.py <quant_dir> <candidates_file> <output_prefix> [tissue_annotations]
#
# Arguments:
#   quant_dir          - Directory containing Salmon output (subfolders with quant.sf)
#   candidates_file    - List of GPCR candidate IDs (one per line or CSV)
#   output_prefix      - Prefix for output files
#   tissue_annotations - Optional: TSV mapping sample names to tissue types

import pandas as pd
import numpy as np
import os
import sys
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# --- Configuration ---
MIN_TPM_THRESHOLD = float(os.getenv('MIN_TPM_THRESHOLD', 1.0))  # Minimum TPM for "expressed"
TAU_THRESHOLD = float(os.getenv('TAU_THRESHOLD', 0.8))  # Threshold for "tissue-specific"
CHEMOSENSORY_TISSUES = os.getenv('CHEMOSENSORY_TISSUES', 'rhinophore,oral_veil,tentacle,cephalic').split(',')


def parse_salmon_quants(quant_dir: Path) -> Tuple[pd.DataFrame, Dict]:
    """
    Parse all quant.sf files in a directory and return combined TPM matrix.

    Args:
        quant_dir: Directory containing sample subdirectories with quant.sf files

    Returns:
        Tuple of (expression_matrix, sample_stats):
            - expression_matrix: DataFrame with genes as rows and samples as columns (TPM values)
            - sample_stats: Dictionary with per-sample statistics
    """
    samples = []
    sample_stats = {}

    # Look for quant.sf files in subdirectories
    quant_files = list(quant_dir.glob("*/quant.sf"))

    if not quant_files:
        # Try direct quant.sf in directory
        quant_files = list(quant_dir.glob("quant.sf"))

    if not quant_files:
        print(f"Error: No quant.sf files found in {quant_dir}", file=sys.stderr)
        return pd.DataFrame(), {}

    print(f"Found {len(quant_files)} Salmon quant.sf files", file=sys.stderr)

    for sf in quant_files:
        # Get sample name from parent directory or filename
        if sf.parent.name != quant_dir.name:
            sample_name = sf.parent.name
        else:
            sample_name = sf.stem

        try:
            df = pd.read_csv(sf, sep="\t")

            # Store sample statistics
            sample_stats[sample_name] = {
                'total_transcripts': len(df),
                'total_reads': df['NumReads'].sum() if 'NumReads' in df.columns else 0,
                'expressed_genes': (df['TPM'] >= MIN_TPM_THRESHOLD).sum()
            }

            # Extract TPM column
            tpm_series = df.set_index("Name")["TPM"].rename(sample_name)
            samples.append(tpm_series)

        except Exception as e:
            print(f"Warning: Could not parse {sf}: {e}", file=sys.stderr)
            continue

    if not samples:
        return pd.DataFrame(), {}

    # Combine all samples into matrix
    expression_matrix = pd.concat(samples, axis=1)

    print(f"Expression matrix: {expression_matrix.shape[0]} genes x {expression_matrix.shape[1]} samples",
          file=sys.stderr)

    return expression_matrix, sample_stats


def calculate_tau_index(expression_row: pd.Series) -> float:
    """
    Calculate tissue specificity index (tau) for a gene.

    Tau ranges from 0 (broadly expressed) to 1 (tissue-specific).
    Formula: tau = sum(1 - x_i/x_max) / (n-1)

    Args:
        expression_row: TPM values across tissues

    Returns:
        Tau index (0-1)
    """
    # Filter out zeros and get non-zero values
    values = expression_row.values
    values = values[~np.isnan(values)]

    if len(values) <= 1:
        return np.nan

    max_val = np.max(values)
    if max_val == 0:
        return np.nan

    # Normalize values
    normalized = values / max_val

    # Calculate tau
    n = len(values)
    tau = np.sum(1 - normalized) / (n - 1)

    return tau


def calculate_specificity_scores(expression_matrix: pd.DataFrame,
                                  tissue_annotations: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Calculate expression metrics for all genes.

    Args:
        expression_matrix: TPM matrix (genes x samples)
        tissue_annotations: Optional mapping of sample names to tissue types

    Returns:
        DataFrame with expression metrics per gene
    """
    results = []

    for gene_id in expression_matrix.index:
        expr = expression_matrix.loc[gene_id]

        # Basic expression stats
        mean_tpm = expr.mean()
        max_tpm = expr.max()
        median_tpm = expr.median()
        std_tpm = expr.std()
        cv = std_tpm / mean_tpm if mean_tpm > 0 else np.nan

        # Count samples where gene is expressed
        n_expressed = (expr >= MIN_TPM_THRESHOLD).sum()
        pct_expressed = n_expressed / len(expr) * 100

        # Tau index (tissue specificity)
        tau = calculate_tau_index(expr)

        # Find top-expressing sample
        top_sample = expr.idxmax()
        top_tpm = expr.max()

        # Determine tissue type if annotations provided
        top_tissue = tissue_annotations.get(top_sample, 'unknown') if tissue_annotations else top_sample

        # Check if enriched in chemosensory tissues
        chemosensory_enriched = False
        chemosensory_fold = np.nan
        chemosensory_mean = 0

        if tissue_annotations:
            # Group samples by tissue type
            chemosensory_samples = [s for s, t in tissue_annotations.items()
                                    if any(ct in t.lower() for ct in CHEMOSENSORY_TISSUES)]
            other_samples = [s for s in expr.index if s not in chemosensory_samples]

            if chemosensory_samples and other_samples:
                chemosensory_mean = expr[chemosensory_samples].mean() if chemosensory_samples else 0
                other_mean = expr[other_samples].mean() if other_samples else 0

                if other_mean > 0:
                    chemosensory_fold = chemosensory_mean / other_mean
                    chemosensory_enriched = chemosensory_fold > 2  # 2-fold enrichment threshold

        results.append({
            'gene_id': gene_id,
            'mean_tpm': mean_tpm,
            'max_tpm': max_tpm,
            'median_tpm': median_tpm,
            'std_tpm': std_tpm,
            'cv': cv,
            'n_samples_expressed': n_expressed,
            'pct_samples_expressed': pct_expressed,
            'tau_index': tau,
            'is_tissue_specific': tau >= TAU_THRESHOLD if pd.notna(tau) else False,
            'top_sample': top_sample,
            'top_tissue': top_tissue,
            'top_tpm': top_tpm,
            'chemosensory_mean_tpm': chemosensory_mean,
            'chemosensory_fold_change': chemosensory_fold,
            'is_chemosensory_enriched': chemosensory_enriched
        })

    return pd.DataFrame(results)


def load_tissue_annotations(annotation_file: str) -> Dict[str, str]:
    """
    Load tissue annotations from a TSV file.

    Expected format:
    sample_name\ttissue_type

    Args:
        annotation_file: Path to annotation TSV

    Returns:
        Dictionary mapping sample names to tissue types
    """
    annotations = {}

    try:
        df = pd.read_csv(annotation_file, sep='\t', header=0)
        if len(df.columns) >= 2:
            for _, row in df.iterrows():
                annotations[str(row.iloc[0])] = str(row.iloc[1])
        print(f"Loaded tissue annotations for {len(annotations)} samples", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not load tissue annotations: {e}", file=sys.stderr)

    return annotations


def filter_candidates(expression_df: pd.DataFrame,
                      candidate_ids: List[str]) -> pd.DataFrame:
    """
    Filter expression data to include only GPCR candidates.

    Args:
        expression_df: Expression metrics DataFrame
        candidate_ids: List of candidate gene IDs

    Returns:
        Filtered DataFrame
    """
    # Create set for fast lookup
    candidate_set = set(candidate_ids)

    # Filter to candidates
    filtered = expression_df[expression_df['gene_id'].isin(candidate_set)]

    print(f"Filtered to {len(filtered)} / {len(candidate_ids)} candidates with expression data",
          file=sys.stderr)

    return filtered


def generate_expression_score(row: pd.Series) -> float:
    """
    Generate a composite expression score for ranking.

    Score incorporates:
    - Expression level (higher = better)
    - Tissue specificity (more specific = better for chemoreceptors)
    - Chemosensory enrichment (enriched = much better)

    Args:
        row: Series with expression metrics

    Returns:
        Composite score (0-1 range)
    """
    # Expression component (log-scaled)
    expr_component = np.log1p(row['mean_tpm']) / 10  # Normalize roughly to 0-1

    # Tissue specificity component
    tau_component = row['tau_index'] if pd.notna(row['tau_index']) else 0

    # Chemosensory enrichment bonus
    chemo_bonus = 0
    if row['is_chemosensory_enriched']:
        chemo_bonus = 0.3
    elif row['chemosensory_fold_change'] > 1.5:
        chemo_bonus = 0.1

    # Combine components
    score = 0.4 * min(expr_component, 1.0) + 0.4 * tau_component + 0.2 + chemo_bonus

    return min(score, 1.0)


def main():
    """Main execution function."""
    # Parse arguments
    if len(sys.argv) < 4:
        print("Usage: python expression_analysis.py <quant_dir> <candidates_file> <output_prefix> [tissue_annotations]")
        sys.exit(1)

    quant_dir = Path(sys.argv[1])
    candidates_file = sys.argv[2]
    output_prefix = sys.argv[3]
    tissue_annotation_file = sys.argv[4] if len(sys.argv) > 4 else None

    # Validate inputs
    if not quant_dir.exists():
        print(f"Error: Quant directory not found: {quant_dir}", file=sys.stderr)
        sys.exit(1)

    # Load tissue annotations if provided
    tissue_annotations = {}
    if tissue_annotation_file and os.path.exists(tissue_annotation_file):
        tissue_annotations = load_tissue_annotations(tissue_annotation_file)

    # Load candidate IDs
    candidate_ids = []
    try:
        with open(candidates_file) as f:
            for line in f:
                parts = line.strip().split(',')
                if parts:
                    candidate_ids.append(parts[0].strip())
        print(f"Loaded {len(candidate_ids)} candidate IDs", file=sys.stderr)
    except Exception as e:
        print(f"Error loading candidates file: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse Salmon quant files
    print("\nParsing Salmon quantification files...", file=sys.stderr)
    result = parse_salmon_quants(quant_dir)

    if isinstance(result, tuple):
        expression_matrix, sample_stats = result
    else:
        expression_matrix = result
        sample_stats = {}

    if expression_matrix.empty:
        print("Error: No expression data parsed", file=sys.stderr)
        sys.exit(1)

    # Calculate expression metrics
    print("\nCalculating expression metrics...", file=sys.stderr)
    expression_metrics = calculate_specificity_scores(expression_matrix, tissue_annotations)

    # Filter to candidates
    candidate_expression = filter_candidates(expression_metrics, candidate_ids)

    # Add expression scores for ranking
    candidate_expression['expression_score'] = candidate_expression.apply(
        generate_expression_score, axis=1
    )

    # Sort by expression score
    candidate_expression = candidate_expression.sort_values('expression_score', ascending=False)

    # --- Write Outputs ---
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Full expression matrix (for all genes)
    matrix_file = f"{output_prefix}_matrix.tsv"
    expression_matrix.to_csv(matrix_file, sep='\t')
    print(f"\nExpression matrix written to: {matrix_file}", file=sys.stderr)

    # 2. Expression metrics for candidates
    metrics_file = f"{output_prefix}_candidate_metrics.tsv"
    candidate_expression.to_csv(metrics_file, sep='\t', index=False)
    print(f"Candidate expression metrics written to: {metrics_file}", file=sys.stderr)

    # 3. Simple expression weights for ranking (id, weight format)
    weights_file = f"{output_prefix}_weights.csv"
    candidate_expression[['gene_id', 'expression_score']].rename(
        columns={'gene_id': 'id', 'expression_score': 'weight'}
    ).to_csv(weights_file, index=False, header=False)
    print(f"Expression weights written to: {weights_file}", file=sys.stderr)

    # 4. Summary statistics
    summary = {
        'total_genes': len(expression_matrix),
        'total_samples': expression_matrix.shape[1],
        'candidates_with_expression': len(candidate_expression),
        'tissue_specific_candidates': candidate_expression['is_tissue_specific'].sum(),
        'chemosensory_enriched_candidates': candidate_expression['is_chemosensory_enriched'].sum(),
        'min_tpm_threshold': MIN_TPM_THRESHOLD,
        'tau_threshold': TAU_THRESHOLD,
        'chemosensory_tissues': CHEMOSENSORY_TISSUES,
        'sample_stats': sample_stats
    }

    summary_file = f"{output_prefix}_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Summary written to: {summary_file}", file=sys.stderr)

    # Print summary
    print("\n=== Expression Analysis Summary ===", file=sys.stderr)
    print(f"  Total genes in matrix: {summary['total_genes']}", file=sys.stderr)
    print(f"  Total samples: {summary['total_samples']}", file=sys.stderr)
    print(f"  Candidates with expression: {summary['candidates_with_expression']}", file=sys.stderr)
    print(f"  Tissue-specific candidates (tau >= {TAU_THRESHOLD}): {summary['tissue_specific_candidates']}",
          file=sys.stderr)
    print(f"  Chemosensory-enriched candidates: {summary['chemosensory_enriched_candidates']}",
          file=sys.stderr)

    # Top candidates by expression score
    print("\n--- Top 10 Candidates by Expression Score ---", file=sys.stderr)
    top10 = candidate_expression.head(10)
    for _, row in top10.iterrows():
        chemo_str = " [CHEMOSENSORY]" if row['is_chemosensory_enriched'] else ""
        spec_str = " [SPECIFIC]" if row['is_tissue_specific'] else ""
        print(f"  {row['gene_id']}: score={row['expression_score']:.3f}, "
              f"tau={row['tau_index']:.2f}, max_TPM={row['max_tpm']:.1f}{chemo_str}{spec_str}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
