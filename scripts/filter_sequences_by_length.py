#!/usr/bin/env python3
"""Filter FASTA sequences by length using statistical outlier detection.

Removes sequences outside acceptable length bounds:
- Lower bound: fixed minimum (default 250 aa, biological minimum for 6+ TM GPCRs)
- Upper bound: Tukey fence (Q3 + 1.5*IQR), floored at a configurable minimum
"""

import argparse
import sys
from pathlib import Path


def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def compute_upper_bound(lengths, method="tukey", floor=800):
    """Compute upper length bound from sequence length distribution."""
    s = sorted(lengths)
    n = len(s)
    q1 = s[n // 4]
    q3 = s[3 * n // 4]
    iqr = q3 - q1
    if method == "tukey":
        bound = q3 + 1.5 * iqr
    else:
        raise ValueError(f"Unknown method: {method}")
    return max(bound, floor)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output FASTA file (filtered)")
    parser.add_argument("--min-length", type=int, default=250,
                        help="Minimum sequence length (default: 250)")
    parser.add_argument("--max-length-method", default="tukey",
                        choices=["tukey"],
                        help="Method for upper bound (default: tukey)")
    parser.add_argument("--max-length-floor", type=int, default=800,
                        help="Minimum value for computed upper bound (default: 800)")
    parser.add_argument("--report", help="Path to write removal report TSV")
    args = parser.parse_args()

    # Read all sequences
    sequences = list(parse_fasta(args.input))
    if not sequences:
        print("Error: no sequences found in input", file=sys.stderr)
        sys.exit(1)

    lengths = [len(seq) for _, seq in sequences]
    upper = compute_upper_bound(lengths, args.max_length_method, args.max_length_floor)
    lower = args.min_length

    print(f"Length filter: min={lower}, max={upper:.0f} "
          f"(method={args.max_length_method}, floor={args.max_length_floor})",
          file=sys.stderr)

    kept = []
    removed = []
    for header, seq in sequences:
        L = len(seq)
        if L < lower:
            removed.append((header, L, f"too_short (<{lower})"))
        elif L > upper:
            removed.append((header, L, f"too_long (>{upper:.0f})"))
        else:
            kept.append((header, seq))

    # Write filtered output
    with open(args.output, "w") as f:
        for header, seq in kept:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

    # Write removal report
    if args.report:
        with open(args.report, "w") as f:
            f.write("sequence_id\tlength\treason\n")
            for sid, length, reason in removed:
                f.write(f"{sid}\t{length}\t{reason}\n")

    print(f"Kept {len(kept)}/{len(sequences)} sequences, "
          f"removed {len(removed)} ({len([r for r in removed if 'short' in r[2]])} short, "
          f"{len([r for r in removed if 'long' in r[2]])} long)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
