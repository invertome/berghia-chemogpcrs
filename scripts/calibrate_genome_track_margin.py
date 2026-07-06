#!/usr/bin/env python3
"""calibrate_genome_track_margin.py — empirical GENOME_TRACK_MIN_MARGIN calibration.

Ground-truthed TP/FP best-vs-second distinct-locus margin ROC. See
docs/plans/2026-07-05-genome-track-margin-calibration-design.md.
"""
from __future__ import annotations
import argparse, math, re, sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reconcile_candidates as rc

INF = float("inf")

def parse_cds_true_locus(header: str) -> tuple[str, int, int]:
    """(scaffold, min, max) from a RefSeq cds_from_genomic.fna header.

    scaffold from the ``lcl|<scaffold>_cds_<protid>_<n>`` id; span from every
    integer inside ``[location=...]`` (so join/complement/partial <,> all fold
    to the coordinate envelope — coarse locus overlap is all TP/FP needs)."""
    header = header.lstrip(">").strip()
    seqid = header.split()[0].split("|", 1)[-1]
    if "_cds_" not in seqid:
        raise ValueError(f"CDS id has no _cds_ scaffold prefix: {seqid}")
    scaffold = seqid.split("_cds_", 1)[0]
    m = re.search(r"\[location=([^\]]+)\]", header)
    if not m:
        raise ValueError(f"CDS header has no [location=]: {header[:80]}")
    coords = [int(x) for x in re.findall(r"\d+", m.group(1))]
    if not coords:
        raise ValueError(f"no coordinates in location: {m.group(1)}")
    return scaffold, min(coords), max(coords)
