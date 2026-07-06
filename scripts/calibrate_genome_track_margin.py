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


def distinct_locus_reps(placements: list[rc.Placement]) -> list[rc.Placement]:
    """One highest-identity rep per distinct locus, identity-descending —
    reconcile_candidates.py's own clustering, so 'distinct locus' matches the gate."""
    return sorted((min(c, key=rc._identity_rank) for c in rc._cluster_placements(placements)),
                  key=rc._identity_rank)

def _overlaps(p: rc.Placement, scaffold: str, start: int, end: int) -> bool:
    # tolerant inclusive overlap; paralog loci are kb apart so 1-base slop is irrelevant.
    return p.chrom == scaffold and p.start <= end and start <= p.end

def _margin(reps: list[rc.Placement]) -> float:
    return INF if len(reps) < 2 else reps[0].pct_identity - reps[1].pct_identity

@dataclass(frozen=True)
class MarginObs:
    tp: float | None
    fp: float | None
    fp_kind: str | None      # "natural" | "synthetic" | None

def query_margins(reps: list[rc.Placement], true_locus: tuple[str, int, int] | None) -> MarginObs:
    """TP margin (correct top locus) and the operative FP margin for one query.

    refseq-mode (true_locus given): top rep overlapping true -> TP; a paralog on
    top -> NATURAL FP; drop-true-locus among >=2 remaining -> SYNTHETIC FP.
    busco-mode (true_locus None): rep0 assumed true (single-copy)."""
    if not reps:
        return MarginObs(None, None, None)
    if true_locus is None:
        true_idx = 0
    else:
        sc, s, e = true_locus
        true_idx = next((i for i, r in enumerate(reps) if _overlaps(r, sc, s, e)), None)
    if true_idx is None:
        # true gene not mapped at all: top-vs-second is an absent-gene grab.
        fp = _margin(reps)
        return MarginObs(None, None if fp == INF else fp, None if fp == INF else "synthetic")
    if true_idx >= 1:
        return MarginObs(None, reps[0].pct_identity - reps[1].pct_identity, "natural")
    tp = _margin(reps)                    # true is rep0
    rest = reps[1:]                       # drop self
    fp = _margin(rest)
    return MarginObs(tp, None if fp == INF else fp, None if fp == INF else "synthetic")
