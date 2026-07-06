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


J_MIN_CLEAN = 0.5        # knee Youden-J below this ⇒ TP/FP overlap ⇒ advise splitting the tiers
GRID_STEP = 0.1

@dataclass(frozen=True)
class RocPoint:
    m: float
    retention: float     # P(tp margin >= m); inf counts as retained
    rejection: float     # P(fp margin < m)
    j: float             # retention + rejection - 1

@dataclass(frozen=True)
class MarginRecommendation:
    min_margin: float
    knee: RocPoint
    split_advisory: bool
    low_margin_threshold: float | None   # set only when split_advisory
    n_tp: int
    n_fp: int

def _grid(tp, fp):
    finite = [x for x in tp + fp if x != INF]
    hi = (max(finite) if finite else 1.0) + GRID_STEP
    n = int(hi / GRID_STEP) + 1
    return [round(i * GRID_STEP, 4) for i in range(n)]

def roc_points(tp, fp, grid=None) -> list[RocPoint]:
    grid = grid if grid is not None else _grid(tp, fp)
    n_tp, n_fp = len(tp), len(fp)
    out = []
    for m in grid:
        ret = (sum(1 for x in tp if x >= m) / n_tp) if n_tp else float("nan")
        rej = (sum(1 for x in fp if x < m) / n_fp) if n_fp else float("nan")
        out.append(RocPoint(m, ret, rej, ret + rej - 1))
    return out

def youden_knee(points) -> RocPoint:
    valid = [p for p in points if not math.isnan(p.j)]
    if not valid:
        raise ValueError("no valid ROC points (need both TP and FP samples)")
    return max(valid, key=lambda p: (p.j, -p.m))     # tie -> smaller margin (favor recall)

def recommend(tp, fp) -> MarginRecommendation:
    knee = youden_knee(roc_points(tp, fp))
    split = knee.j < J_MIN_CLEAN
    return MarginRecommendation(
        min_margin=round(knee.m, 2), knee=knee, split_advisory=split,
        low_margin_threshold=(round(_percentile([x for x in tp if x != INF] or [knee.m], 50.0), 2)
                              if split else None),
        n_tp=len(tp), n_fp=len(fp))

def _percentile(values, p):     # linear-interp percentile, stdlib (mirrors the id/cov calibrator)
    s = sorted(float(v) for v in values)
    if not s: raise ValueError("percentile of empty")
    if len(s) == 1: return s[0]
    k = (len(s) - 1) * (p / 100.0); lo, hi = math.floor(k), math.ceil(k)
    return s[int(k)] if lo == hi else s[lo] * (hi - k) + s[hi] * (k - lo)
