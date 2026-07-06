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


def load_cds_true_loci(cds_fasta: str) -> dict[str, tuple[str, int, int]]:
    loci = {}
    with open(cds_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                qid = line[1:].split()[0]
                loci[qid] = parse_cds_true_locus(line)
    return loci

def _placements_by_query(paf, gff):
    parsed = []
    if paf: parsed += rc.parse_minimap2_paf(paf)
    if gff: parsed += rc.parse_gmap_gff(gff)
    by_q: dict[str, list[rc.Placement]] = {}
    for p in parsed:
        by_q.setdefault(p.query, []).append(p)
    return by_q

def collect_margins(by_query, true_loci, restrict_ids=None):
    """-> (tp_list, fp_list). true_loci=None ⇒ busco-mode (rep0 truth)."""
    tp, fp = [], []
    ids = restrict_ids if restrict_ids is not None else list(by_query)
    for q in ids:
        pls = by_query.get(q)
        if not pls: continue
        reps = distinct_locus_reps(pls)
        tl = None if true_loci is None else true_loci.get(q)
        if true_loci is not None and tl is None: continue    # no truth for this id
        obs = query_margins(reps, tl)
        if obs.tp is not None: tp.append(obs.tp)
        if obs.fp is not None: fp.append(obs.fp)
    return tp, fp

def _synthetic():
    # 6 clean genes (TP margin ~8-11) + 4 absent-gene grabs (FP margin ~1-3);
    # FP ceiling 3.0 ⇒ Youden knee at 3.1 (clean separation, well inside 3<m<8).
    tp = [8.2, 9.1, 10.4, 8.8, 11.0, 9.6]
    fp = [1.1, 2.0, 3.0, 1.5]
    return tp, fp

def format_table(rec: MarginRecommendation) -> str:
    lines = ["# Margin-gate calibration",
             f"TP (correct) n={rec.n_tp}   FP (spurious) n={rec.n_fp}",
             f"knee: margin={rec.knee.m:g}  retention={rec.knee.retention:.3f}  "
             f"rejection={rec.knee.rejection:.3f}  J={rec.knee.j:.3f}",
             f"clean separation: {'NO (advise tier split)' if rec.split_advisory else 'YES'}",
             "", "## Recommended (write into config.sh after review)",
             f"GENOME_TRACK_MIN_MARGIN={rec.min_margin:g}"]
    if rec.split_advisory:
        lines.append(f"GENOME_TRACK_LOW_MARGIN={rec.low_margin_threshold:g}")
    return "\n".join(lines)

def write_recommendation(path, rec, dry_run=False):
    hdr = ["# DRY RUN synthetic — NOT for production"] if dry_run else []
    hdr += ["# scripts/calibrate_genome_track_margin.py",
            f"# knee J={rec.knee.j:.3f}; clean={not rec.split_advisory}",
            f"GENOME_TRACK_MIN_MARGIN={rec.min_margin:g}"]
    if rec.split_advisory:
        hdr.append(f"GENOME_TRACK_LOW_MARGIN={rec.low_margin_threshold:g}")
    Path(path).write_text("\n".join(hdr) + "\n")

def write_figure(path, tp, fp, rec):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"WARNING: matplotlib unavailable, skipping figure: {e}", file=sys.stderr); return
    pts = roc_points(tp, fp)
    ms = [p.m for p in pts]
    plt.figure(figsize=(6, 4))
    plt.plot(ms, [p.retention for p in pts], label="TP retention")
    plt.plot(ms, [p.rejection for p in pts], label="FP rejection")
    plt.axvline(rec.min_margin, ls="--", color="k", label=f"knee={rec.min_margin:g}")
    plt.xlabel("min_margin (%id points)"); plt.ylabel("fraction"); plt.legend(); plt.tight_layout()
    plt.savefig(path, dpi=150); plt.close()

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Calibrate GENOME_TRACK_MIN_MARGIN from a TP/FP margin ROC.")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--mode", choices=["refseq", "busco"], default="refseq")
    ap.add_argument("--minimap2-paf"); ap.add_argument("--gmap-gff")
    ap.add_argument("--refseq-cds", help="GPCR CDS FASTA (headers -> true loci) for refseq-mode")
    ap.add_argument("--busco-ids", help="restrict to these query ids for busco-mode")
    ap.add_argument("--out-recommendation"); ap.add_argument("--out-figure")
    a = ap.parse_args(argv)

    if a.dry_run:
        tp, fp = _synthetic()
    else:
        if not (a.minimap2_paf or a.gmap_gff):
            print("ERROR: need --minimap2-paf and/or --gmap-gff", file=sys.stderr); return 2
        by_q = _placements_by_query(a.minimap2_paf, a.gmap_gff)
        if a.mode == "refseq":
            if not a.refseq_cds:
                print("ERROR: refseq-mode needs --refseq-cds", file=sys.stderr); return 2
            tp, fp = collect_margins(by_q, load_cds_true_loci(a.refseq_cds))
        else:
            ids = [l.split()[0] for l in open(a.busco_ids) if l.strip() and not l.startswith("#")] if a.busco_ids else None
            tp, fp = collect_margins(by_q, None, restrict_ids=ids)

    if not tp or not fp:
        print(f"ERROR: insufficient data (TP={len(tp)}, FP={len(fp)})", file=sys.stderr); return 1
    rec = recommend(tp, fp)
    print(format_table(rec))
    if a.out_recommendation: write_recommendation(a.out_recommendation, rec, dry_run=a.dry_run)
    if a.out_figure: write_figure(a.out_figure, tp, fp, rec)
    return 0

if __name__ == "__main__":      # pragma: no cover
    raise SystemExit(main())
