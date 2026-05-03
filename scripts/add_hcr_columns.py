#!/usr/bin/env python3
"""add_hcr_columns.py — Augment ranked candidates CSV with HCR-friendliness columns.

Bead -xqz. Reads the ranked candidates CSV and adds:
  cds_length_bp           — from --cds-fasta (FASTA of CDS, e.g. GENOME_CDS)
  paralog_min_identity    — closest-paralog windowed % identity (from
                            --alignment FASTA, the trimmed protein MSA)
  hcr_probe_friendly      — boolean from _hcr_design_lib.is_hcr_probe_friendly

The actual HCR probe sequences are designed by Molecular Instruments'
service; this script just adds diagnostic columns so wet-lab staff can
prioritise candidates whose probe design is plausibly tractable.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _hcr_design_lib import (  # noqa: E402
    HCR_DEFAULT_MAX_PARALOG_IDENTITY,
    HCR_DEFAULT_MIN_CDS_LENGTH_BP,
    closest_paralog_identity,
    is_hcr_probe_friendly,
)


def fasta_records(path: str):
    """Minimal FASTA parser: yields (header_token, sequence)."""
    if not path:
        return
    h, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(parts)
                h = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if h is not None:
        yield h, "".join(parts)


def cds_lengths_from_fasta(cds_fasta: str) -> dict[str, int]:
    """Map gene/transcript ID -> CDS length in bp.

    For RefSeq CDS files the header is like:
      >lcl|NC_xxxxx.x_cds_XP_yyy.1_1234 [gene=...] [protein_id=XP_yyy.1] [...]
    We index by the first whitespace token AND by [gene=...] / [protein_id=...]
    extracted via simple substring scans, so the lookup matches whichever
    naming convention the rank_candidates 'id' column uses.
    """
    out: dict[str, int] = {}
    for header, seq in fasta_records(cds_fasta):
        ln = len(seq)
        out[header] = ln
        # Extract bracketed key=value pairs
        # (this avoids a regex dependency for a tiny grammar)
        i = 0
        while True:
            lb = header.find("[", i)
            if lb < 0:
                break
            rb = header.find("]", lb + 1)
            if rb < 0:
                break
            kv = header[lb + 1 : rb]
            i = rb + 1
            if "=" in kv:
                k, v = kv.split("=", 1)
                out.setdefault(v.strip(), ln)
                # Also index "gene=X" without the prefix for convenience
                out.setdefault(f"{k.strip()}={v.strip()}", ln)
    return out


def aligned_seqs_from_fasta(aln_fasta: str) -> dict[str, str]:
    return {h: s for h, s in fasta_records(aln_fasta)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--ranked-csv", required=True)
    ap.add_argument("--cds-fasta", default="",
                    help="Path to a CDS FASTA (e.g. $GENOME_CDS) for cds_length_bp")
    ap.add_argument("--alignment", default="",
                    help="Path to trimmed protein MSA for paralog distance")
    ap.add_argument("--window-aa", type=int, default=300)
    ap.add_argument("--step-aa", type=int, default=50)
    ap.add_argument("--max-paralog-identity", type=float,
                    default=HCR_DEFAULT_MAX_PARALOG_IDENTITY)
    ap.add_argument("--min-cds-length-bp", type=int,
                    default=HCR_DEFAULT_MIN_CDS_LENGTH_BP)
    ap.add_argument("--out", required=True,
                    help="Output CSV path (will overwrite if exists)")
    args = ap.parse_args()

    df = pd.read_csv(args.ranked_csv)
    if "id" not in df.columns:
        print("ERROR: ranked CSV missing 'id' column", file=sys.stderr)
        return 1

    cds_map = cds_lengths_from_fasta(args.cds_fasta) if args.cds_fasta else {}
    aligned = aligned_seqs_from_fasta(args.alignment) if args.alignment else {}

    cds_lengths = []
    paralog_idents = []
    paralog_ids = []
    friendly = []
    for _, row in df.iterrows():
        cid = str(row["id"])
        # CDS length lookup: try id, then any [gene=...] / [protein_id=...] alias
        cds_len = cds_map.get(cid)
        if cds_len is None:
            # Try suffix-based RefSeq protein id
            for k in (f"protein_id={cid}", f"gene={cid}"):
                cds_len = cds_map.get(k)
                if cds_len is not None:
                    break
        cds_lengths.append(cds_len if cds_len is not None else "")

        if aligned and cid in aligned:
            close_id, close_ident = closest_paralog_identity(
                cid, aligned, window_aa=args.window_aa, step_aa=args.step_aa,
            )
            paralog_idents.append(close_ident)
            paralog_ids.append(close_id or "")
        else:
            paralog_idents.append("")
            paralog_ids.append("")

        tandem_size = row.get("tandem_cluster_size")
        try:
            tandem_size = int(tandem_size) if pd.notna(tandem_size) else None
        except (TypeError, ValueError):
            tandem_size = None
        friendly.append(is_hcr_probe_friendly(
            cds_length_bp=cds_len,
            paralog_min_identity=(paralog_idents[-1]
                                  if paralog_idents[-1] != "" else None),
            tandem_cluster_size=tandem_size,
            max_paralog_identity=args.max_paralog_identity,
            min_cds_length_bp=args.min_cds_length_bp,
        ))

    df["cds_length_bp"] = cds_lengths
    df["paralog_min_identity"] = paralog_idents
    df["closest_paralog_id"] = paralog_ids
    df["hcr_probe_friendly"] = friendly

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    n_friendly = sum(1 for v in friendly if v)
    print(f"Wrote {len(df)} rows to {args.out}; {n_friendly} flagged "
          f"hcr_probe_friendly", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
