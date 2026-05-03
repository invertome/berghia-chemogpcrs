#!/usr/bin/env python3
"""parse_jcvi_anchors.py — Convert JCVI MCscan anchor files to per-gene CSV.

Bead -e59. JCVI's `.anchors` file format:

    ###
    geneA1  geneB1  evalue
    geneA2  geneB2  evalue
    ...
    ###
    geneA9  geneB9  evalue
    ...

Each `###` block is one collinear segment between two scaffolds. We
aggregate per-gene-A counts: how many anchor blocks contain that gene,
how many anchor genes total are in those blocks (block-size proxy).

Output CSV columns:
    candidate_id            (gene name as it appears in column 1)
    n_anchor_blocks         number of distinct collinear blocks containing
                            this gene as gene-A
    total_anchor_genes      sum of block lengths (genes-A side) across all
                            blocks containing this gene
    target_species          --target-prefix value (provenance)

This becomes the canonical synteny input for `rank_candidates.py`,
replacing the legacy `synteny_ids.txt` line-count input that the audit
flagged (counted any minimap2 hit, not actual collinearity blocks).
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_anchors(anchor_path: str) -> list[list[tuple[str, str, float]]]:
    """Parse a JCVI .anchors file into a list of blocks.

    Each block is a list of (geneA, geneB, evalue) tuples.
    """
    blocks: list[list[tuple[str, str, float]]] = []
    cur: list[tuple[str, str, float]] = []
    with open(anchor_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if cur:
                    blocks.append(cur)
                    cur = []
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            geneA = parts[0].strip()
            geneB = parts[1].strip()
            try:
                evalue = float(parts[2]) if len(parts) >= 3 else float("nan")
            except ValueError:
                evalue = float("nan")
            cur.append((geneA, geneB, evalue))
    if cur:
        blocks.append(cur)
    return blocks


def per_gene_summary(blocks, target_species: str) -> pd.DataFrame:
    """Return DataFrame with one row per gene-A, summarising block memberships."""
    n_blocks: dict[str, int] = defaultdict(int)
    total_genes: dict[str, int] = defaultdict(int)
    for block in blocks:
        if not block:
            continue
        # All gene-A entries in this block see the same block size
        block_size = len(block)
        # Each gene appears at most once per block; deduplicate within-block
        seen_in_block: set[str] = set()
        for geneA, _geneB, _ev in block:
            if geneA in seen_in_block:
                continue
            seen_in_block.add(geneA)
            n_blocks[geneA] += 1
            total_genes[geneA] += block_size
    rows = [
        {
            "candidate_id": gid,
            "n_anchor_blocks": n_blocks[gid],
            "total_anchor_genes": total_genes[gid],
            "target_species": target_species,
        }
        for gid in sorted(n_blocks)
    ]
    return pd.DataFrame(rows, columns=[
        "candidate_id", "n_anchor_blocks", "total_anchor_genes", "target_species",
    ])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--anchors", required=True,
                    help="Input JCVI .anchors (or .lifted.anchors) file")
    ap.add_argument("--target-species", required=True,
                    help="Target species short name for the target_species column")
    ap.add_argument("--out", required=True, help="Output CSV path")
    args = ap.parse_args()

    blocks = parse_anchors(args.anchors)
    df = per_gene_summary(blocks, args.target_species)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    n_genes = len(df)
    n_blocks = len(blocks)
    print(f"Parsed {n_blocks} collinear blocks; {n_genes} gene-A entries; "
          f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
