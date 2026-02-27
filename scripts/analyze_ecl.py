#!/usr/bin/env python3
# analyze_ecl.py
# Purpose: Analyze extracellular loop (ECL) divergence vs transmembrane (TM) conservation.
# High ECL divergence with conserved TM suggests ligand-binding diversification.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python analyze_ecl.py --alignments <file_or_glob> --deeptmhmm <dir>
#                         --min-ecl-length <int> --output <csv>

import os
import sys
import argparse
import glob
import re
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict, Counter
from Bio import AlignIO, SeqIO


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze ECL divergence vs TM conservation'
    )
    parser.add_argument(
        '--alignments', required=True,
        help='Alignment file path or glob pattern (e.g., "aligned_*.fasta")'
    )
    parser.add_argument(
        '--deeptmhmm', required=True,
        help='Directory with DeepTMHMM output files (predicted_topologies.3line and/or TMRs.gff3)'
    )
    parser.add_argument(
        '--min-ecl-length', type=int, default=5,
        help='Minimum ECL length to analyze (default: 5)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output CSV file'
    )
    return parser.parse_args()


def parse_deeptmhmm_output(deeptmhmm_dir: str) -> Dict[str, Dict]:
    """
    Parse DeepTMHMM predictions to get TM and ECL boundaries.

    Reads both 3-line format and GFF3 format. The 3-line format is preferred
    as it provides per-residue topology directly.

    Returns dict: gene_id -> {
        'tm_regions': [(start, end), ...],
        'ecl_regions': [(start, end), ...],  # Extracellular loops
        'icl_regions': [(start, end), ...],  # Intracellular loops
    }
    """
    predictions = {}
    deeptmhmm_path = Path(deeptmhmm_dir)

    if not deeptmhmm_path.exists():
        print(f"Warning: DeepTMHMM directory not found: {deeptmhmm_dir}", file=sys.stderr)
        return predictions

    # Look for 3line format files (both .3line extension and exact filename)
    threeline_files = list(deeptmhmm_path.glob('*.3line'))
    if not threeline_files:
        # Check for file named exactly "predicted_topologies.3line"
        candidate = deeptmhmm_path / 'predicted_topologies.3line'
        if candidate.exists():
            threeline_files = [candidate]

    for pred_file in threeline_files:
        try:
            with open(pred_file) as f:
                lines = f.readlines()
                i = 0
                while i < len(lines):
                    # DeepTMHMM 3-line format:
                    #   Line 1: >header [| TYPE]
                    #   Line 2: amino acid sequence
                    #   Line 3: topology labels (I/O/M/S/B/P)
                    if lines[i].startswith('>'):
                        gene_id = lines[i][1:].strip().split()[0]
                        if i + 2 < len(lines):
                            topology = lines[i + 2].strip()
                            regions = parse_topology_string(topology)
                            predictions[gene_id] = regions
                        i += 3
                    else:
                        i += 1
        except Exception as e:
            print(f"Warning: Could not parse {pred_file}: {e}", file=sys.stderr)

    # Also look for GFF3 format (as fallback if 3-line didn't provide predictions)
    for gff_file in deeptmhmm_path.glob('*.gff3'):
        try:
            gene_regions = defaultdict(lambda: {
                'tm_regions': [],
                'ecl_regions': [],
                'icl_regions': []
            })

            with open(gff_file) as f:
                for line in f:
                    if line.startswith('#') or line.startswith('//'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        gene_id = parts[0]
                        # DeepTMHMM GFF3: col1=feature, col2=start, col3=end
                        # Standard GFF3:  col2=feature, col3=start, col4=end
                        # Detect format: if col1 is a known feature name, use DeepTMHMM layout
                        if parts[1] in ('TMhelix', 'inside', 'outside', 'signal',
                                        'Beta sheet', 'periplasm'):
                            feature = parts[1]
                            start = int(parts[2])
                            end = int(parts[3])
                        elif len(parts) >= 5:
                            feature = parts[2]
                            start = int(parts[3])
                            end = int(parts[4])
                        else:
                            continue

                        if feature == 'TMhelix':
                            gene_regions[gene_id]['tm_regions'].append((start, end))
                        elif feature == 'outside':
                            gene_regions[gene_id]['ecl_regions'].append((start, end))
                        elif feature == 'inside':
                            gene_regions[gene_id]['icl_regions'].append((start, end))

            for gene_id, regions in gene_regions.items():
                if gene_id not in predictions:
                    predictions[gene_id] = regions
        except Exception as e:
            print(f"Warning: Could not parse {gff_file}: {e}", file=sys.stderr)

    print(f"Loaded TM predictions for {len(predictions)} genes", file=sys.stderr)
    return predictions


def parse_topology_string(topology: str) -> Dict:
    """
    Parse DeepTMHMM topology string.

    Format: string of i/o/M characters
    - i: inside (cytoplasmic)
    - o: outside (extracellular)
    - M: membrane (transmembrane)

    Returns dict with tm_regions, ecl_regions, icl_regions.
    Coordinates are 1-indexed, inclusive (matching DeepTMHMM convention).
    """
    regions = {
        'tm_regions': [],
        'ecl_regions': [],
        'icl_regions': []
    }

    current_type = None
    start = 0

    for i, char in enumerate(topology):
        if char.upper() == 'M':
            region_type = 'tm'
        elif char.upper() == 'O':
            region_type = 'ecl'
        elif char.upper() == 'I':
            region_type = 'icl'
        else:
            continue

        if region_type != current_type:
            if current_type is not None:
                regions[f'{current_type}_regions'].append((start + 1, i))
            current_type = region_type
            start = i

    # Add final region
    if current_type is not None:
        regions[f'{current_type}_regions'].append((start + 1, len(topology)))

    return regions


def calculate_column_conservation(alignment, positions: List[int]) -> float:
    """
    Calculate frequency-based conservation for specified alignment positions.

    Uses sum of squared amino acid frequencies per column (equivalent to
    expected pairwise identity). This is O(n) per column instead of O(n^2)
    for explicit pairwise comparison.

    Args:
        alignment: Bio.Align.MultipleSeqAlignment
        positions: List of column positions (0-indexed)

    Returns:
        Mean conservation score (0 to 1)
    """
    if not positions or len(alignment) < 2:
        return 0.0

    aln_length = alignment.get_alignment_length()
    valid_positions = [p for p in positions if 0 <= p < aln_length]
    if not valid_positions:
        return 0.0

    conservation_scores = []

    for pos in valid_positions:
        column = alignment[:, pos]
        # Count non-gap characters
        chars = [c.upper() for c in column if c not in ('-', 'X', '*')]
        n = len(chars)
        if n < 2:
            continue

        # Frequency-based conservation: sum(freq_i^2) = expected pairwise identity
        counts = Counter(chars)
        freq_sq_sum = sum((count / n) ** 2 for count in counts.values())
        conservation_scores.append(freq_sq_sum)

    return float(np.mean(conservation_scores)) if conservation_scores else 0.0


def find_sequence_in_alignment(gene_id: str, alignment) -> Optional[str]:
    """
    Find a sequence in the alignment by gene ID.

    Handles ID matching between DeepTMHMM predictions (original Berghia IDs
    like BersteEVm036785t1) and alignment IDs (which may be the same or
    ref_TAXID_N format for references).

    Returns the aligned sequence string, or None if not found.
    """
    for record in alignment:
        rec_id = record.id
        # Exact match
        if rec_id == gene_id:
            return str(record.seq)
        # Match ignoring trailing description
        if rec_id.split()[0] == gene_id:
            return str(record.seq)

    return None


def map_sequence_to_alignment(gene_id: str, alignment) -> Dict[int, int]:
    """
    Map sequence positions to alignment positions.

    Args:
        gene_id: sequence identifier to find in alignment
        alignment: Bio.Align.MultipleSeqAlignment

    Returns dict: seq_position (1-indexed) -> alignment_position (0-indexed)
    """
    aln_seq = find_sequence_in_alignment(gene_id, alignment)
    if not aln_seq:
        return {}

    mapping = {}
    seq_pos = 0

    for aln_pos, char in enumerate(aln_seq):
        if char != '-':
            seq_pos += 1
            mapping[seq_pos] = aln_pos

    return mapping


def get_orthogroup_name(aln_file: str, is_single_file: bool) -> str:
    """Extract orthogroup name from alignment file path."""
    stem = Path(aln_file).stem
    if is_single_file:
        # For the combined alignment, use a generic name
        return stem.replace('_aligned', '').replace('aligned_', '')
    return stem.replace('aligned_', '')


def analyze_ecl_divergence(
    alignment_files: List[str],
    tm_predictions: Dict[str, Dict],
    min_ecl_length: int,
    is_single_file: bool = False
) -> pd.DataFrame:
    """
    Analyze ECL vs TM conservation for each gene in alignments.

    Only genes with TM predictions are analyzed. In the combined alignment,
    this means only Berghia candidates (which have DeepTMHMM/KD predictions)
    get ECL scores — references serve as the conservation background.

    Returns DataFrame with divergence analysis per gene.
    """
    results = []

    for aln_file in alignment_files:
        try:
            alignment = AlignIO.read(aln_file, 'fasta')
            orthogroup = get_orthogroup_name(aln_file, is_single_file)
            n_seqs = len(alignment)

            print(f"  Analyzing {orthogroup} ({n_seqs} sequences, "
                  f"{alignment.get_alignment_length()} columns)...",
                  file=sys.stderr)

            # Count how many sequences have predictions
            genes_with_pred = [
                record.id for record in alignment
                if record.id in tm_predictions
            ]
            print(f"    {len(genes_with_pred)} sequences have TM predictions",
                  file=sys.stderr)

            for gene_id in genes_with_pred:
                pred = tm_predictions[gene_id]
                tm_regions = pred.get('tm_regions', [])
                ecl_regions = pred.get('ecl_regions', [])

                if not tm_regions or not ecl_regions:
                    continue

                # Map sequence to alignment
                pos_map = map_sequence_to_alignment(gene_id, alignment)

                if not pos_map:
                    continue

                # Get alignment positions for TM regions
                tm_positions = []
                for start, end in tm_regions:
                    for pos in range(start, end + 1):
                        if pos in pos_map:
                            tm_positions.append(pos_map[pos])

                # Get alignment positions for ECL regions (filtered by min length)
                ecl_positions = []
                ecl_lengths = []
                for start, end in ecl_regions:
                    ecl_len = end - start + 1
                    if ecl_len >= min_ecl_length:
                        ecl_lengths.append(ecl_len)
                        for pos in range(start, end + 1):
                            if pos in pos_map:
                                ecl_positions.append(pos_map[pos])

                if not tm_positions or not ecl_positions:
                    continue

                # Calculate conservation using frequency-based method
                tm_conservation = calculate_column_conservation(alignment, tm_positions)
                ecl_conservation = calculate_column_conservation(alignment, ecl_positions)

                # Calculate divergence (1 - conservation)
                tm_divergence = 1 - tm_conservation
                ecl_divergence = 1 - ecl_conservation

                # Calculate ratio
                if tm_divergence > 0.01:
                    ecl_tm_ratio = ecl_divergence / tm_divergence
                else:
                    # TM highly conserved, any ECL divergence is notable
                    ecl_tm_ratio = ecl_divergence * 10 if ecl_divergence > 0 else 0.0

                results.append({
                    'orthogroup': orthogroup,
                    'gene_id': gene_id,
                    'n_tm_regions': len(tm_regions),
                    'n_ecl_regions': len([l for l in ecl_lengths if l >= min_ecl_length]),
                    'mean_ecl_length': float(np.mean(ecl_lengths)) if ecl_lengths else 0,
                    'tm_conservation': tm_conservation,
                    'ecl_conservation': ecl_conservation,
                    'tm_divergence': tm_divergence,
                    'ecl_divergence': ecl_divergence,
                    'ecl_tm_ratio': ecl_tm_ratio,
                    'ecl_divergence_score': ecl_tm_ratio  # Alias for ranking
                })

        except Exception as e:
            print(f"Warning: Could not analyze {aln_file}: {e}", file=sys.stderr)

    return pd.DataFrame(results)


def main():
    """Main execution function."""
    args = parse_arguments()

    # Find alignment files — support both single file and glob pattern
    is_single_file = False
    if os.path.isfile(args.alignments):
        alignment_files = [args.alignments]
        is_single_file = True
    else:
        alignment_files = glob.glob(args.alignments)
        if not alignment_files:
            # Try as directory + pattern
            base_dir = os.path.dirname(args.alignments) or '.'
            pattern = os.path.basename(args.alignments)
            alignment_files = glob.glob(os.path.join(base_dir, pattern))

    if not alignment_files:
        print(f"Warning: No alignment files found matching {args.alignments}",
              file=sys.stderr)
        pd.DataFrame(columns=[
            'orthogroup', 'gene_id', 'n_tm_regions', 'n_ecl_regions',
            'mean_ecl_length', 'tm_conservation', 'ecl_conservation',
            'tm_divergence', 'ecl_divergence', 'ecl_tm_ratio', 'ecl_divergence_score'
        ]).to_csv(args.output, index=False)
        return

    print(f"Found {len(alignment_files)} alignment file(s)", file=sys.stderr)
    if is_single_file:
        print(f"  Single-file mode: {alignment_files[0]}", file=sys.stderr)

    # Parse DeepTMHMM predictions
    tm_predictions = parse_deeptmhmm_output(args.deeptmhmm)

    if not tm_predictions:
        print("Warning: No TM predictions found — cannot compute ECL divergence",
              file=sys.stderr)
        pd.DataFrame(columns=[
            'orthogroup', 'gene_id', 'n_tm_regions', 'n_ecl_regions',
            'mean_ecl_length', 'tm_conservation', 'ecl_conservation',
            'tm_divergence', 'ecl_divergence', 'ecl_tm_ratio', 'ecl_divergence_score'
        ]).to_csv(args.output, index=False)
        return

    # Analyze ECL divergence
    print("Analyzing ECL divergence...", file=sys.stderr)
    result_df = analyze_ecl_divergence(
        alignment_files, tm_predictions, args.min_ecl_length,
        is_single_file=is_single_file
    )

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write output
    result_df.to_csv(args.output, index=False)

    # Print summary
    print(f"\n=== ECL Divergence Analysis Summary ===", file=sys.stderr)
    print(f"Genes analyzed: {len(result_df)}", file=sys.stderr)

    if not result_df.empty:
        high_ratio = (result_df['ecl_tm_ratio'] >= 2.0).sum()
        print(f"Genes with high ECL/TM divergence ratio (>=2.0): {high_ratio}",
              file=sys.stderr)
        print(f"Mean ECL divergence: {result_df['ecl_divergence'].mean():.3f}",
              file=sys.stderr)
        print(f"Mean TM divergence: {result_df['tm_divergence'].mean():.3f}",
              file=sys.stderr)
        print(f"Mean ECL/TM ratio: {result_df['ecl_tm_ratio'].mean():.3f}",
              file=sys.stderr)

    print(f"Output written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
