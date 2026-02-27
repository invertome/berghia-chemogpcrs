#!/usr/bin/env python3
"""
predict_tm_helices.py
Fast local transmembrane helix prediction using Kyte-Doolittle hydrophobicity.

Scans protein sequences with a sliding window to identify hydrophobic stretches
consistent with transmembrane alpha-helices. Designed as a fast pre-filter
(seconds on ~100K sequences) before slower methods like DeepTMHMM.

Output format matches DeepTMHMM prediction file:
    seqid\tpred_type\tconfidence\tSP\tn_tm_regions

Optional --topology-output writes DeepTMHMM-compatible .3line format:
    >header
    MEFSKNQSLM...                      (amino acid sequence)
    OOOOMMMMMMMMMMMMMMMMMMMMIIII...     (per-residue topology: I/O/M)

Usage:
    python predict_tm_helices.py input.fasta output_prediction [--min-tm 6]
    python predict_tm_helices.py input.fasta output_prediction --topology-output topo.3line
"""

import sys
import argparse


# Kyte-Doolittle hydrophobicity scale
KD_SCALE = {
    'A':  1.8, 'C':  2.5, 'D': -3.5, 'E': -3.5, 'F':  2.8,
    'G': -0.4, 'H': -3.2, 'I':  4.5, 'K': -3.9, 'L':  3.8,
    'M':  1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V':  4.2, 'W': -0.9, 'Y': -1.3,
}

# TM helix parameters
WINDOW_SIZE = 21       # Standard TM helix length
HYDRO_THRESHOLD = 1.6  # Mean hydrophobicity cutoff for TM segment
MIN_GAP = 8            # Minimum residues between TM helices (loop region)


def find_tm_segments(sequence):
    """
    Find putative transmembrane helix segments in a protein sequence.

    Uses sliding window hydrophobicity analysis with merging of overlapping
    predictions and enforcement of minimum loop lengths between helices.

    Returns list of (start, end) tuples (0-indexed, exclusive end).
    """
    seq = sequence.upper()
    n = len(seq)
    if n < WINDOW_SIZE:
        return []

    # Calculate hydrophobicity for each window position
    tm_positions = []
    for i in range(n - WINDOW_SIZE + 1):
        window = seq[i:i + WINDOW_SIZE]
        score = sum(KD_SCALE.get(aa, 0.0) for aa in window) / WINDOW_SIZE
        if score >= HYDRO_THRESHOLD:
            tm_positions.append(i)

    if not tm_positions:
        return []

    # Merge overlapping/adjacent TM windows into segments
    segments = []
    seg_start = tm_positions[0]
    seg_end = tm_positions[0] + WINDOW_SIZE

    for pos in tm_positions[1:]:
        if pos <= seg_end:
            seg_end = pos + WINDOW_SIZE
        else:
            segments.append((seg_start, seg_end))
            seg_start = pos
            seg_end = pos + WINDOW_SIZE
    segments.append((seg_start, seg_end))

    # Filter: enforce minimum gap between helices
    filtered = [segments[0]]
    for seg in segments[1:]:
        prev_end = filtered[-1][1]
        if seg[0] - prev_end >= MIN_GAP:
            filtered.append(seg)
        else:
            # Merge with previous if too close
            filtered[-1] = (filtered[-1][0], seg[1])

    return filtered


def count_tm_helices(sequence):
    """
    Count putative transmembrane helices in a protein sequence.

    Uses sliding window hydrophobicity analysis with merging of overlapping
    predictions and enforcement of minimum loop lengths between helices.
    """
    segments = find_tm_segments(sequence)

    # Count segments with reasonable TM helix length (15-45 residues)
    count = 0
    for start, end in segments:
        length = end - start
        if 15 <= length <= 45:
            count += 1
        elif length > 45:
            # Long hydrophobic stretch — likely multiple helices
            count += round(length / 25)

    return count


def build_topology_string(sequence, segments):
    """
    Build per-residue topology string from TM segment boundaries.

    Assumes standard Class A GPCR topology: N-terminus extracellular.
    - Before TM1: outside (o)
    - After TM1: inside (i)
    - After TM2: outside (o)
    - ... alternating

    Args:
        sequence: protein sequence string
        segments: list of (start, end) tuples from find_tm_segments (0-indexed, exclusive end)

    Returns:
        topology string of same length as sequence using i/M/o characters
    """
    n = len(sequence)
    topo = ['o'] * n  # Default: outside

    if not segments:
        return ''.join(topo)

    # Mark TM regions
    for start, end in segments:
        for i in range(start, min(end, n)):
            topo[i] = 'M'

    # Assign inside/outside to non-TM regions by alternating
    # GPCR N-terminus is extracellular (outside), so:
    # before TM1 = outside, after TM1 = inside, after TM2 = outside, etc.
    inside = False  # First loop (before TM1) is outside
    for i in range(n):
        if topo[i] == 'M':
            continue
        # Check if we just exited a TM region
        if i > 0 and topo[i - 1] == 'M':
            inside = not inside
        topo[i] = 'i' if inside else 'o'

    return ''.join(topo)


def parse_fasta(filepath):
    """Parse FASTA file yielding (header, sequence) tuples."""
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, ''.join(seq_parts)


def main():
    parser = argparse.ArgumentParser(
        description='Fast TM helix prediction via Kyte-Doolittle hydrophobicity'
    )
    parser.add_argument('fasta', help='Input protein FASTA file')
    parser.add_argument('output', help='Output prediction file (DeepTMHMM format)')
    parser.add_argument('--min-tm', type=int, default=0,
                        help='Only output sequences with >= N TM helices (0=all)')
    parser.add_argument('--topology-output',
                        help='Write DeepTMHMM-compatible .3line topology file')
    args = parser.parse_args()

    total = 0
    written = 0
    topo_entries = []

    with open(args.output, 'w') as out:
        for seq_id, sequence in parse_fasta(args.fasta):
            total += 1
            segments = find_tm_segments(sequence)

            # Count TM helices from segments
            tm_count = 0
            for start, end in segments:
                length = end - start
                if 15 <= length <= 45:
                    tm_count += 1
                elif length > 45:
                    tm_count += round(length / 25)

            if tm_count >= args.min_tm:
                pred_type = "TMH" if tm_count > 0 else "GLOB"
                # Confidence is approximate — hydrophobicity is less precise than DeepTMHMM
                conf = min(0.7 + 0.03 * tm_count, 0.95) if tm_count >= 6 else 0.5
                out.write(f"{seq_id}\tTMH\t{conf:.2f}\tSP\t{tm_count}\n")
                written += 1

            # Build topology string if requested
            if args.topology_output:
                topo_str = build_topology_string(sequence, segments)
                topo_entries.append((seq_id, sequence, topo_str))

            if total % 10000 == 0:
                print(f"  Processed {total} sequences...", file=sys.stderr)

    # Write topology file in DeepTMHMM 3-line format:
    #   >header
    #   amino_acid_sequence
    #   topology_labels (I/O/M characters)
    if args.topology_output and topo_entries:
        with open(args.topology_output, 'w') as f:
            for seq_id, sequence, topo_str in topo_entries:
                f.write(f">{seq_id}\n")
                f.write(f"{sequence}\n")
                f.write(f"{topo_str.upper()}\n")
        print(f"Wrote {len(topo_entries)} topology predictions to {args.topology_output}",
              file=sys.stderr)

    print(f"Processed {total} sequences, wrote {written} predictions", file=sys.stderr)


if __name__ == '__main__':
    main()
