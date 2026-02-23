#!/usr/bin/env python3
"""
predict_tm_helices.py
Fast local transmembrane helix prediction using Kyte-Doolittle hydrophobicity.

Scans protein sequences with a sliding window to identify hydrophobic stretches
consistent with transmembrane alpha-helices. Designed as a fast pre-filter
(seconds on ~100K sequences) before slower methods like DeepTMHMM.

Output format matches DeepTMHMM prediction file:
    seqid\tpred_type\tconfidence\tSP\tn_tm_regions

Usage:
    python predict_tm_helices.py input.fasta output_prediction [--min-tm 6]
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


def count_tm_helices(sequence):
    """
    Count putative transmembrane helices in a protein sequence.

    Uses sliding window hydrophobicity analysis with merging of overlapping
    predictions and enforcement of minimum loop lengths between helices.
    """
    seq = sequence.upper()
    n = len(seq)
    if n < WINDOW_SIZE:
        return 0

    # Calculate hydrophobicity for each window position
    tm_positions = []
    for i in range(n - WINDOW_SIZE + 1):
        window = seq[i:i + WINDOW_SIZE]
        score = sum(KD_SCALE.get(aa, 0.0) for aa in window) / WINDOW_SIZE
        if score >= HYDRO_THRESHOLD:
            tm_positions.append(i)

    if not tm_positions:
        return 0

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

    # Count segments with reasonable TM helix length (15-45 residues)
    count = 0
    for start, end in filtered:
        length = end - start
        if 15 <= length <= 45:
            count += 1
        elif length > 45:
            # Long hydrophobic stretch — likely multiple helices
            count += round(length / 25)

    return count


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
    args = parser.parse_args()

    total = 0
    written = 0

    with open(args.output, 'w') as out:
        for seq_id, sequence in parse_fasta(args.fasta):
            total += 1
            tm_count = count_tm_helices(sequence)
            if tm_count >= args.min_tm:
                pred_type = "TMH" if tm_count > 0 else "GLOB"
                # Confidence is approximate — hydrophobicity is less precise than DeepTMHMM
                conf = min(0.7 + 0.03 * tm_count, 0.95) if tm_count >= 6 else 0.5
                out.write(f"{seq_id}\tTMH\t{conf:.2f}\tSP\t{tm_count}\n")
                written += 1

            if total % 10000 == 0:
                print(f"  Processed {total} sequences...", file=sys.stderr)

    print(f"Processed {total} sequences, wrote {written} predictions", file=sys.stderr)


if __name__ == '__main__':
    main()
