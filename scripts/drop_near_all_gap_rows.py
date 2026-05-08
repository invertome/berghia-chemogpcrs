#!/usr/bin/env python3
"""drop_near_all_gap_rows.py — defensive post-ClipKit alignment cleaner.

Bug context (commit 7db2a49, May 2026 classification tree build):
ClipKit's `kpic-smart-gap` mode trims columns based on alignment-wide
gappiness, which can leave a few sequences with ≥95 % gaps in the
output. IQ-TREE 3's auto-detection of sequence alphabet inspects the
first few rows; if those rows are near-all-gap it raises
"Unknown sequence type" and aborts. Even with `-st AA` forced, a row
of pure '-' carries no signal and just slows the run.

This script reads a FASTA alignment, drops every sequence whose gap
fraction is at or above the threshold (default 0.95), and writes the
remainder. It exits non-zero only on hard I/O errors; an alignment
where every row is dropped is reported but produces an empty (valid)
output file so downstream tooling can decide how to react.

Usage:
    python3 drop_near_all_gap_rows.py \\
        --input  trimmed.fa \\
        --output cleaned.fa \\
        [--gap-threshold 0.95]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def gap_fraction(seq: str) -> float:
    """Return the fraction of '-' characters in `seq`. Empty seq -> 1.0
    (treated as fully gappy so it gets dropped)."""
    if not seq:
        return 1.0
    return seq.count("-") / len(seq)


def filter_alignment(input_path: str, output_path: str,
                     gap_threshold: float = 0.95) -> tuple[int, int]:
    """Read FASTA at `input_path`, drop sequences with gap fraction
    ≥ `gap_threshold`, write the rest to `output_path`. Returns
    (n_kept, n_dropped)."""
    n_kept = 0
    n_dropped = 0
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    cur_id: str | None = None
    cur_lines: list[str] = []

    with open(input_path) as fin, open(output_path, "w") as fout:
        def _flush() -> None:
            nonlocal n_kept, n_dropped
            if cur_id is None:
                return
            seq = "".join(cur_lines)
            if gap_fraction(seq) < gap_threshold:
                fout.write(cur_id + "\n")
                fout.write(seq + "\n")
                n_kept += 1
            else:
                n_dropped += 1

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush()
                cur_id = line
                cur_lines = []
            else:
                cur_lines.append(line)
        _flush()

    return n_kept, n_dropped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--gap-threshold", type=float, default=0.95,
                    help="Drop sequences whose gap fraction is at least "
                         "this value (default 0.95). Must be in (0, 1].")
    args = ap.parse_args()

    if not (0 < args.gap_threshold <= 1.0):
        print(f"ERROR: --gap-threshold must be in (0, 1], got "
              f"{args.gap_threshold}", file=sys.stderr)
        return 2

    n_kept, n_dropped = filter_alignment(
        args.input, args.output, args.gap_threshold)
    print(f"[drop_gappy] kept {n_kept}, dropped {n_dropped} "
          f"(threshold gap_fraction >= {args.gap_threshold})",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
