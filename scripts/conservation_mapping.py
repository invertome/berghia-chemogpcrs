#!/usr/bin/env python3
# conservation_mapping.py
# Purpose: Calculate per-residue conservation scores and map to 3D structures.
# Identifies conserved functional motifs in GPCR candidates.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   python conservation_mapping.py <alignment> <output_prefix> [structure_dir] [reference_seq]
#
# Arguments:
#   alignment      - Multiple sequence alignment (FASTA format)
#   output_prefix  - Prefix for output files
#   structure_dir  - Optional: directory with PDB files for 3D mapping
#   reference_seq  - Optional: reference sequence ID for position numbering

import os
import sys
import json
import math
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

# Standard amino acid alphabet
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

# GPCR-specific conserved motifs
GPCR_MOTIFS = {
    'DRY': {
        'pattern': '[DE]R[YF]',
        'function': 'G-protein coupling',
        'region': 'TM3-ICL2 junction',
        'conservation_expected': 'high'
    },
    'NPxxY': {
        'pattern': 'NP..Y',
        'function': 'Receptor activation',
        'region': 'TM7',
        'conservation_expected': 'high'
    },
    'CWxP': {
        'pattern': 'CW.P',
        'function': 'Rotamer toggle switch',
        'region': 'TM6',
        'conservation_expected': 'high'
    },
    'ERC': {
        'pattern': 'ERC',
        'function': 'Mollusc-specific DRY variant',
        'region': 'TM3-ICL2',
        'conservation_expected': 'medium'
    },
    'TM3_Arg': {
        'pattern': 'R',
        'function': 'Conserved arginine (3.50)',
        'region': 'TM3',
        'conservation_expected': 'very high'
    }
}

# Amino acid groupings for conservation analysis
AA_GROUPS = {
    'hydrophobic': set('AILMFVWP'),
    'polar': set('STNQ'),
    'positive': set('KRH'),
    'negative': set('DE'),
    'aromatic': set('FYW'),
    'small': set('AGST'),
    'cysteine': set('C'),
    'proline': set('P'),
    'glycine': set('G')
}


def load_alignment(alignment_file: str) -> MultipleSeqAlignment:
    """Load multiple sequence alignment."""
    alignment = AlignIO.read(alignment_file, "fasta")
    print(f"Loaded alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} positions",
          file=sys.stderr)
    return alignment


def calculate_shannon_entropy(column: str) -> float:
    """
    Calculate Shannon entropy for an alignment column.

    Lower entropy = higher conservation.

    Args:
        column: String of amino acids at one alignment position

    Returns:
        Shannon entropy value (bits)
    """
    # Remove gaps
    column = column.replace('-', '').replace('X', '')

    if len(column) == 0:
        return float('nan')

    # Calculate frequencies
    freq = defaultdict(int)
    for aa in column:
        freq[aa] += 1

    total = len(column)
    entropy = 0.0

    for count in freq.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def calculate_conservation_score(column: str) -> float:
    """
    Calculate conservation score (0-1) for an alignment column.

    Uses inverse Shannon entropy normalized by maximum possible entropy.

    Args:
        column: String of amino acids at one alignment position

    Returns:
        Conservation score (0 = not conserved, 1 = fully conserved)
    """
    entropy = calculate_shannon_entropy(column)

    if math.isnan(entropy):
        return 0.0

    # Maximum entropy for 20 amino acids
    max_entropy = math.log2(20)

    # Convert to conservation score (inverse, normalized)
    conservation = 1.0 - (entropy / max_entropy)

    return max(0.0, min(1.0, conservation))


def calculate_gap_frequency(column: str) -> float:
    """Calculate the frequency of gaps in a column."""
    gap_count = column.count('-') + column.count('.')
    return gap_count / len(column) if len(column) > 0 else 1.0


def get_consensus_aa(column: str) -> Tuple[str, float]:
    """
    Get the consensus amino acid and its frequency.

    Args:
        column: String of amino acids at one alignment position

    Returns:
        Tuple of (consensus_aa, frequency)
    """
    # Remove gaps
    column_clean = column.replace('-', '').replace('X', '').replace('.', '')

    if len(column_clean) == 0:
        return '-', 0.0

    # Count amino acids
    counts = defaultdict(int)
    for aa in column_clean:
        counts[aa] += 1

    # Find most common
    consensus = max(counts, key=counts.get)
    frequency = counts[consensus] / len(column_clean)

    return consensus, frequency


def calculate_property_conservation(column: str) -> Dict[str, float]:
    """
    Calculate conservation of physicochemical properties.

    Args:
        column: String of amino acids at one alignment position

    Returns:
        Dictionary of property conservation scores
    """
    column_clean = column.replace('-', '').replace('X', '').replace('.', '')

    if len(column_clean) == 0:
        return {prop: 0.0 for prop in AA_GROUPS}

    property_scores = {}

    for prop, aa_set in AA_GROUPS.items():
        count = sum(1 for aa in column_clean if aa in aa_set)
        property_scores[prop] = count / len(column_clean)

    return property_scores


def analyze_alignment_conservation(alignment: MultipleSeqAlignment,
                                  reference_id: Optional[str] = None) -> pd.DataFrame:
    """
    Calculate conservation metrics for each alignment position.

    Args:
        alignment: Multiple sequence alignment
        reference_id: Optional reference sequence for position mapping

    Returns:
        DataFrame with conservation metrics per position
    """
    aln_length = alignment.get_alignment_length()

    # Find reference sequence if specified
    ref_seq = None
    ref_positions = {}

    if reference_id:
        for rec in alignment:
            if reference_id in rec.id:
                ref_seq = str(rec.seq)
                # Map alignment positions to reference positions
                ref_pos = 0
                for aln_pos, aa in enumerate(ref_seq):
                    if aa != '-':
                        ref_pos += 1
                        ref_positions[aln_pos] = ref_pos
                break

    results = []

    for pos in range(aln_length):
        # Get column
        column = ''.join(str(rec.seq)[pos] for rec in alignment)

        # Calculate metrics
        conservation = calculate_conservation_score(column)
        entropy = calculate_shannon_entropy(column)
        gap_freq = calculate_gap_frequency(column)
        consensus_aa, consensus_freq = get_consensus_aa(column)
        property_cons = calculate_property_conservation(column)

        # Reference position (if available)
        ref_position = ref_positions.get(pos, None)
        ref_aa = ref_seq[pos] if ref_seq else None

        results.append({
            'alignment_position': pos + 1,  # 1-indexed
            'reference_position': ref_position,
            'consensus_aa': consensus_aa,
            'consensus_frequency': consensus_freq,
            'reference_aa': ref_aa,
            'conservation_score': conservation,
            'shannon_entropy': entropy,
            'gap_frequency': gap_freq,
            'n_sequences': len(column) - column.count('-'),
            **{f'prop_{k}': v for k, v in property_cons.items()}
        })

    return pd.DataFrame(results)


def identify_conserved_motifs(conservation_df: pd.DataFrame,
                             alignment: MultipleSeqAlignment) -> List[Dict]:
    """
    Identify conserved GPCR motifs in the alignment.

    Args:
        conservation_df: Conservation metrics DataFrame
        alignment: Multiple sequence alignment

    Returns:
        List of identified motifs with positions and conservation
    """
    import re

    # Build consensus sequence
    consensus = ''.join(conservation_df['consensus_aa'].values)

    motifs_found = []

    for motif_name, motif_info in GPCR_MOTIFS.items():
        pattern = motif_info['pattern']

        for match in re.finditer(pattern, consensus):
            start = match.start()
            end = match.end()

            # Get conservation scores for motif positions
            motif_conservation = conservation_df.iloc[start:end]['conservation_score'].values

            motifs_found.append({
                'motif_name': motif_name,
                'pattern': pattern,
                'matched_sequence': match.group(),
                'alignment_start': start + 1,  # 1-indexed
                'alignment_end': end,
                'length': end - start,
                'mean_conservation': float(np.mean(motif_conservation)),
                'min_conservation': float(np.min(motif_conservation)),
                'function': motif_info['function'],
                'region': motif_info['region'],
                'expected_conservation': motif_info['conservation_expected']
            })

    return motifs_found


def identify_highly_conserved_regions(conservation_df: pd.DataFrame,
                                      window_size: int = 10,
                                      threshold: float = 0.7) -> List[Dict]:
    """
    Identify highly conserved regions using sliding window.

    Args:
        conservation_df: Conservation metrics DataFrame
        window_size: Size of sliding window
        threshold: Minimum mean conservation score

    Returns:
        List of highly conserved regions
    """
    scores = conservation_df['conservation_score'].values
    n_positions = len(scores)

    regions = []
    in_region = False
    region_start = 0

    for i in range(n_positions - window_size + 1):
        window_mean = np.mean(scores[i:i + window_size])

        if window_mean >= threshold and not in_region:
            in_region = True
            region_start = i
        elif window_mean < threshold and in_region:
            in_region = False
            regions.append({
                'start': region_start + 1,  # 1-indexed
                'end': i + window_size - 1,
                'length': i + window_size - 1 - region_start,
                'mean_conservation': float(np.mean(scores[region_start:i + window_size - 1]))
            })

    # Handle region at end
    if in_region:
        regions.append({
            'start': region_start + 1,
            'end': n_positions,
            'length': n_positions - region_start,
            'mean_conservation': float(np.mean(scores[region_start:]))
        })

    return regions


def map_conservation_to_structure(conservation_df: pd.DataFrame,
                                  pdb_file: str,
                                  output_file: str) -> bool:
    """
    Map conservation scores to PDB B-factor column for visualization.

    Args:
        conservation_df: Conservation metrics DataFrame
        pdb_file: Input PDB file
        output_file: Output PDB file with conservation in B-factor

    Returns:
        True if successful
    """
    try:
        from Bio.PDB import PDBParser, PDBIO

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)

        # Get conservation scores indexed by position
        scores = {int(row['reference_position']): row['conservation_score']
                 for _, row in conservation_df.iterrows()
                 if pd.notna(row['reference_position'])}

        # Map to B-factors
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_num = residue.get_id()[1]
                    cons_score = scores.get(res_num, 0.0)
                    # Scale to 0-100 for B-factor
                    b_factor = cons_score * 100

                    for atom in residue:
                        atom.set_bfactor(b_factor)

        # Write output
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)

        print(f"Conservation-mapped structure written to: {output_file}", file=sys.stderr)
        return True

    except Exception as e:
        print(f"Warning: Could not map conservation to structure: {e}", file=sys.stderr)
        return False


def generate_conservation_profile_plot(conservation_df: pd.DataFrame,
                                       output_file: str,
                                       motifs: List[Dict] = None):
    """
    Generate conservation profile plot.

    Args:
        conservation_df: Conservation metrics DataFrame
        output_file: Output plot file
        motifs: List of identified motifs to annotate
    """
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

        positions = conservation_df['alignment_position'].values
        conservation = conservation_df['conservation_score'].values
        gap_freq = conservation_df['gap_frequency'].values
        entropy = conservation_df['shannon_entropy'].values

        # Panel 1: Conservation score
        ax1 = axes[0]
        ax1.fill_between(positions, conservation, alpha=0.7, color='steelblue')
        ax1.axhline(y=0.7, color='red', linestyle='--', alpha=0.5, label='High conservation threshold')
        ax1.set_ylabel('Conservation Score')
        ax1.set_ylim(0, 1)
        ax1.legend(loc='upper right')
        ax1.set_title('Per-residue Conservation Profile')

        # Annotate motifs
        if motifs:
            for motif in motifs:
                ax1.axvspan(motif['alignment_start'], motif['alignment_end'],
                           alpha=0.3, color='green')
                ax1.annotate(motif['motif_name'],
                           xy=(motif['alignment_start'], 0.95),
                           fontsize=8, rotation=45)

        # Panel 2: Gap frequency
        ax2 = axes[1]
        ax2.fill_between(positions, gap_freq, alpha=0.7, color='orange')
        ax2.set_ylabel('Gap Frequency')
        ax2.set_ylim(0, 1)

        # Panel 3: Shannon entropy
        ax3 = axes[2]
        ax3.fill_between(positions, entropy, alpha=0.7, color='purple')
        ax3.set_ylabel('Shannon Entropy (bits)')
        ax3.set_xlabel('Alignment Position')

        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"Conservation profile plot saved to: {output_file}", file=sys.stderr)

    except Exception as e:
        print(f"Warning: Could not generate plot: {e}", file=sys.stderr)


def main():
    """Main execution function."""
    if len(sys.argv) < 3:
        print("Usage: python conservation_mapping.py <alignment> <output_prefix> [structure_dir] [reference_seq]")
        sys.exit(1)

    alignment_file = sys.argv[1]
    output_prefix = sys.argv[2]
    structure_dir = sys.argv[3] if len(sys.argv) > 3 else None
    reference_id = sys.argv[4] if len(sys.argv) > 4 else None

    # Validate inputs
    if not os.path.exists(alignment_file):
        print(f"Error: Alignment file not found: {alignment_file}", file=sys.stderr)
        sys.exit(1)

    # Load alignment
    print("Loading alignment...", file=sys.stderr)
    alignment = load_alignment(alignment_file)

    # If no reference specified, use first sequence
    if not reference_id:
        reference_id = alignment[0].id
        print(f"Using first sequence as reference: {reference_id}", file=sys.stderr)

    # Calculate conservation
    print("\nCalculating per-residue conservation...", file=sys.stderr)
    conservation_df = analyze_alignment_conservation(alignment, reference_id)

    # Identify conserved motifs
    print("Identifying conserved motifs...", file=sys.stderr)
    motifs = identify_conserved_motifs(conservation_df, alignment)

    # Identify highly conserved regions
    conserved_regions = identify_highly_conserved_regions(conservation_df)

    # Create output directory
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Write Outputs ---

    # 1. Full conservation table
    cons_file = f"{output_prefix}_conservation.tsv"
    conservation_df.to_csv(cons_file, sep='\t', index=False)
    print(f"\nConservation scores written to: {cons_file}", file=sys.stderr)

    # 2. Identified motifs
    if motifs:
        motifs_df = pd.DataFrame(motifs)
        motifs_file = f"{output_prefix}_motifs.tsv"
        motifs_df.to_csv(motifs_file, sep='\t', index=False)
        print(f"Identified motifs written to: {motifs_file}", file=sys.stderr)

    # 3. Conserved regions
    if conserved_regions:
        regions_df = pd.DataFrame(conserved_regions)
        regions_file = f"{output_prefix}_conserved_regions.tsv"
        regions_df.to_csv(regions_file, sep='\t', index=False)
        print(f"Conserved regions written to: {regions_file}", file=sys.stderr)

    # 4. Summary statistics
    summary = {
        'alignment_file': alignment_file,
        'n_sequences': len(alignment),
        'alignment_length': alignment.get_alignment_length(),
        'reference_sequence': reference_id,
        'mean_conservation': float(conservation_df['conservation_score'].mean()),
        'median_conservation': float(conservation_df['conservation_score'].median()),
        'highly_conserved_positions': int((conservation_df['conservation_score'] >= 0.9).sum()),
        'conserved_positions': int((conservation_df['conservation_score'] >= 0.7).sum()),
        'variable_positions': int((conservation_df['conservation_score'] < 0.3).sum()),
        'n_motifs_found': len(motifs),
        'n_conserved_regions': len(conserved_regions),
        'motifs_summary': [{'name': m['motif_name'], 'conservation': m['mean_conservation']}
                          for m in motifs]
    }

    summary_file = f"{output_prefix}_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    # 5. Conservation profile plot
    plot_file = f"{output_prefix}_profile.png"
    generate_conservation_profile_plot(conservation_df, plot_file, motifs)

    # 6. Map to structures if available
    if structure_dir and os.path.exists(structure_dir):
        print("\nMapping conservation to structures...", file=sys.stderr)
        structure_output_dir = Path(f"{output_prefix}_structures")
        structure_output_dir.mkdir(exist_ok=True)

        pdb_files = list(Path(structure_dir).glob("*.pdb"))
        for pdb_file in pdb_files[:10]:  # Limit to 10 structures
            output_pdb = structure_output_dir / f"{pdb_file.stem}_conservation.pdb"
            map_conservation_to_structure(conservation_df, str(pdb_file), str(output_pdb))

    # Print summary
    print("\n=== Conservation Analysis Summary ===", file=sys.stderr)
    print(f"Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} positions",
          file=sys.stderr)
    print(f"Mean conservation score: {summary['mean_conservation']:.3f}", file=sys.stderr)
    print(f"Highly conserved positions (≥0.9): {summary['highly_conserved_positions']}", file=sys.stderr)
    print(f"Conserved positions (≥0.7): {summary['conserved_positions']}", file=sys.stderr)
    print(f"Variable positions (<0.3): {summary['variable_positions']}", file=sys.stderr)

    if motifs:
        print(f"\nGPCR motifs identified: {len(motifs)}", file=sys.stderr)
        for motif in motifs:
            print(f"  {motif['motif_name']}: {motif['matched_sequence']} "
                  f"(pos {motif['alignment_start']}-{motif['alignment_end']}, "
                  f"conservation={motif['mean_conservation']:.2f})", file=sys.stderr)

    if conserved_regions:
        print(f"\nHighly conserved regions: {len(conserved_regions)}", file=sys.stderr)
        for region in conserved_regions[:5]:
            print(f"  Positions {region['start']}-{region['end']}: "
                  f"conservation={region['mean_conservation']:.2f}", file=sys.stderr)


if __name__ == "__main__":
    main()
