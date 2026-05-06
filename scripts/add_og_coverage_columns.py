#!/usr/bin/env python3
"""add_og_coverage_columns.py — Augment ranked CSV with per-OG ref-CDS coverage.

Bead -ogc. Adds three columns to the ranked candidates CSV:

    og_n_ref_cds          number of OG members that have a CDS in the
                          merged reference CDS file (all_references_cds.fna)
    og_n_total            total number of OG members (across all species)
    og_dnds_reliability   'high' (>=10 ref CDS) | 'medium' (>=5) | 'low' (<5)

These are TRANSPARENCY columns — they don't change the ranking. They give
reviewers and wet-lab triagers context on whether each candidate's dN/dS
contribution to the score came from a well-sampled OG (high reliability)
or a degenerate one (low reliability).

Why this matters: bad-CDS rates are non-uniform across species (cate 81%,
aplcal 68%, arvu 58% miniprot bad-CDS rate per Unity B2 job 56821506).
The OGs where dN/dS would be most informative — chemoreceptor LSE
expansions in slug/snail genomes — are exactly the OGs where reference
CDS coverage is sparsest. Without the transparency annotation, reviewers
can't tell which candidates' rankings are robust vs noise-shaped.

Hooks into stage 07 alongside add_hcr_columns.py:

    rank_candidates.py -> ranked_candidates_sorted.csv
        |
    add_hcr_columns.py (cds_length_bp, paralog_min_identity, hcr_probe_friendly)
        |
    add_og_coverage_columns.py (og_n_ref_cds, og_n_total, og_dnds_reliability)
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


HIGH_THRESHOLD = 10
MEDIUM_THRESHOLD = 5


def reliability_flag(n_ref_cds: int) -> str:
    """Map reference-CDS count to a reliability bin.

    >=10 = high, 5-9 = medium, <5 = low.
    """
    if n_ref_cds >= HIGH_THRESHOLD:
        return "high"
    if n_ref_cds >= MEDIUM_THRESHOLD:
        return "medium"
    return "low"


def parse_cds_ids(cds_fasta_path: str) -> set[str]:
    """Return the set of sequence IDs (first whitespace token after '>')
    present in the merged reference CDS FASTA.

    Returns an empty set if the file doesn't exist (graceful — the caller
    decides whether that's an error or just missing context)."""
    ids: set[str] = set()
    if not cds_fasta_path or not os.path.exists(cds_fasta_path):
        return ids
    with open(cds_fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                tok = line[1:].split()[0] if len(line) > 1 else ""
                if tok:
                    ids.add(tok)
    return ids


def load_og_members(orthogroups_tsv_path: str) -> dict[str, list[str]]:
    """Parse OrthoFinder Orthogroups.tsv:

        Orthogroup<tab>species_1<tab>species_2<tab>...
        OG0000001<tab>gene_a, gene_b<tab>gene_c<tab>gene_d, gene_e
        ...

    Returns dict: og_id -> list of all member gene IDs (concatenated across
    species columns; empty cells skipped).
    """
    og_members: dict[str, list[str]] = {}
    if not orthogroups_tsv_path or not os.path.exists(orthogroups_tsv_path):
        return og_members
    with open(orthogroups_tsv_path) as f:
        next(f, None)  # discard header (species names)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            og_id = parts[0].strip()
            if not og_id:
                continue
            members: list[str] = []
            for cell in parts[1:]:
                cell = cell.strip()
                if not cell:
                    continue
                for g in cell.split(","):
                    g = g.strip()
                    if g:
                        members.append(g)
            og_members[og_id] = members
    return og_members


def add_coverage_columns(
    *,
    ranked_csv_path: str,
    cds_fasta_path: str,
    orthogroups_tsv_path: str,
    out_path: str,
    og_column: str = "orthogroup",
) -> int:
    """Read the ranked CSV, augment with og_n_ref_cds + og_n_total +
    og_dnds_reliability columns, write to out_path. Returns the number
    of rows annotated."""
    df = pd.read_csv(ranked_csv_path)
    cds_ids = parse_cds_ids(cds_fasta_path)
    og_members = load_og_members(orthogroups_tsv_path)

    if og_column not in df.columns:
        print(
            f"WARN: ranked CSV {ranked_csv_path!r} has no {og_column!r} column. "
            f"Coverage will be 0 / 'low' for all rows.",
            file=sys.stderr,
        )
        df["og_n_ref_cds"] = 0
        df["og_n_total"] = 0
        df["og_dnds_reliability"] = "low"
    else:
        # Pre-compute per-OG coverage so each row is a dict lookup, not a scan.
        og_cov: dict[str, int] = {}
        og_tot: dict[str, int] = {}
        for og_id, members in og_members.items():
            og_tot[og_id] = len(members)
            og_cov[og_id] = sum(1 for m in members if m in cds_ids)

        df["og_n_ref_cds"] = df[og_column].map(
            lambda og: og_cov.get(str(og), 0) if pd.notna(og) else 0
        ).astype(int)
        df["og_n_total"] = df[og_column].map(
            lambda og: og_tot.get(str(og), 0) if pd.notna(og) else 0
        ).astype(int)
        df["og_dnds_reliability"] = df["og_n_ref_cds"].apply(reliability_flag)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return len(df)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True, help="Input ranked CSV")
    ap.add_argument("--cds-fasta", required=True,
                    help="Path to merged reference CDS FASTA "
                         "(typically results/reference_sequences/cds/all_references_cds.fna)")
    ap.add_argument("--orthogroups-tsv", required=True,
                    help="Path to OrthoFinder Orthogroups.tsv")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--og-column", default="orthogroup",
                    help="Column name in ranked CSV holding the OG id "
                         "(default: 'orthogroup')")
    args = ap.parse_args()

    n_rows = add_coverage_columns(
        ranked_csv_path=args.ranked_csv,
        cds_fasta_path=args.cds_fasta,
        orthogroups_tsv_path=args.orthogroups_tsv,
        out_path=args.out,
        og_column=args.og_column,
    )
    print(f"Wrote {args.out} (annotated {n_rows} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
