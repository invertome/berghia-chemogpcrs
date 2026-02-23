#!/usr/bin/env python3
"""
split_mcl_to_fasta.py
Split MCL cluster output into per-orthogroup FASTA files.

MCL dump format: one cluster per line, tab-separated sequence IDs.
Reads the combined FASTA to extract sequences for each cluster.

Usage:
    python split_mcl_to_fasta.py clusters.mcl combined.fa output_dir \
        --min-seqs 3 --prefix OG_conserved
"""

import argparse
import os
import sys


def parse_fasta(filepath):
    """Parse FASTA file into dict of {header: sequence}."""
    seqs = {}
    header = None
    parts = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        seqs[header] = ''.join(parts)
    return seqs


def parse_mcl_clusters(filepath):
    """Parse MCL dump output (one cluster per line, tab-separated IDs)."""
    clusters = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            members = line.split('\t')
            if members:
                clusters.append(members)
    return clusters


def main():
    parser = argparse.ArgumentParser(
        description='Split MCL clusters into per-orthogroup FASTA files'
    )
    parser.add_argument('clusters', help='MCL cluster dump file')
    parser.add_argument('fasta', help='Combined FASTA file with all sequences')
    parser.add_argument('output_dir', help='Output directory for per-orthogroup FASTAs')
    parser.add_argument('--min-seqs', type=int, default=3,
                        help='Minimum sequences per orthogroup (default: 3)')
    parser.add_argument('--prefix', default='OG',
                        help='Prefix for orthogroup names (default: OG)')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Loading sequences from {args.fasta}...", file=sys.stderr)
    seqs = parse_fasta(args.fasta)
    print(f"Loaded {len(seqs)} sequences", file=sys.stderr)

    print(f"Parsing clusters from {args.clusters}...", file=sys.stderr)
    clusters = parse_mcl_clusters(args.clusters)
    print(f"Found {len(clusters)} clusters", file=sys.stderr)

    written = 0
    skipped = 0
    for i, members in enumerate(clusters):
        if len(members) < args.min_seqs:
            skipped += 1
            continue

        og_name = f"{args.prefix}_{i:05d}"
        outpath = os.path.join(args.output_dir, f"{og_name}.fa")

        with open(outpath, 'w') as out:
            found = 0
            for seq_id in members:
                if seq_id in seqs:
                    out.write(f">{seq_id}\n{seqs[seq_id]}\n")
                    found += 1
            if found < args.min_seqs:
                os.remove(outpath)
                skipped += 1
                continue

        written += 1

    print(f"Wrote {written} orthogroup files, skipped {skipped} "
          f"(< {args.min_seqs} sequences)", file=sys.stderr)


if __name__ == '__main__':
    main()
