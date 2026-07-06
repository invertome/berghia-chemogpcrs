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
    """One query's contribution to the ROC: its TP margin (``tp`` — the correct
    top-locus best-vs-second margin, ``inf`` for a lone locus, ``None`` when the
    top locus is not the true one) and its operative FP margin (``fp`` with
    ``fp_kind`` ``"natural"`` = a paralog genuinely outranked the true locus, or
    ``"synthetic"`` = the true-locus rep dropped to model an absent gene; both
    ``None`` when no spurious margin is available)."""
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
    """One threshold ``m`` on the margin ROC and its scores (see ``roc_points``):
    ``retention`` = TP kept, ``rejection`` = FP dropped, ``j`` = Youden-J
    (``retention + rejection - 1``)."""
    m: float
    retention: float     # P(tp margin >= m); inf counts as retained
    rejection: float     # P(fp margin < m)
    j: float             # retention + rejection - 1

@dataclass(frozen=True)
class MarginRecommendation:
    """The calibration verdict: the recommended ``min_margin``, the diagnostic
    ``knee`` RocPoint, and — when TP/FP overlap (``split_advisory``) — a
    ``low_margin_threshold`` for the advisory tier.

    ``min_margin`` is the knee ONLY on clean separation; on overlap it is the
    design-§5 keep-and-flag floor (p5 of the multi-locus TP) and the knee is
    retained as a DIAGNOSTIC of separability only (emitting it as a hard gate
    would drop most correct multi-locus placements — see ``recommend``).

    The knee is computed over the MULTI-LOCUS (finite-margin) TP population only
    (design §2: the margin gate acts only among >=2 competing loci; single-locus
    genes have margin=inf and are trivially retained, so including them would
    inflate retention and push the knee to an artifact). ``n_tp`` is the TOTAL TP
    count and ``n_tp_multi`` the multi-locus subset the knee is built on;
    ``full_pop_retention`` is the fraction of ALL TP (incl. the single-locus inf)
    retained at the recommended ``min_margin`` — reported as real-world CONTEXT.
    ``n_fp`` is the FP sample size."""
    min_margin: float
    knee: RocPoint
    split_advisory: bool
    low_margin_threshold: float | None   # set only when split_advisory
    n_tp: int
    n_fp: int
    n_tp_multi: int
    full_pop_retention: float

def _grid(tp, fp):
    """Threshold sweep grid: 0 .. (max finite margin + GRID_STEP) in GRID_STEP
    increments. ``round`` (not ``int``) on the point count so the grid always
    pads one full step BEYOND the largest data value — otherwise float
    truncation (e.g. 8.1/0.1 == 80.9999) clips the last point exactly at the max
    and the figure loses the tail of the curve (design §5 wants the full curve).
    ``inf`` margins are excluded from the span."""
    finite = [x for x in tp + fp if x != INF]
    hi = (max(finite) if finite else 1.0) + GRID_STEP
    n = round(hi / GRID_STEP) + 1
    return [round(i * GRID_STEP, 4) for i in range(n)]

def roc_points(tp, fp, grid=None) -> list[RocPoint]:
    """Evaluate TP-retention / FP-rejection / Youden-J at every grid threshold.

    For each ``m``: ``retention = P(tp margin >= m)`` (an ``inf`` margin — a
    single-locus gene with no competitor — is always retained) and
    ``rejection = P(fp margin < m)``. With no TP (resp. FP) samples the
    corresponding fraction is NaN (and ``j`` NaN), so ``youden_knee`` skips those
    points rather than dividing by zero."""
    grid = grid if grid is not None else _grid(tp, fp)
    n_tp, n_fp = len(tp), len(fp)
    out = []
    for m in grid:
        ret = (sum(1 for x in tp if x >= m) / n_tp) if n_tp else float("nan")
        rej = (sum(1 for x in fp if x < m) / n_fp) if n_fp else float("nan")
        out.append(RocPoint(m, ret, rej, ret + rej - 1))
    return out

def youden_knee(points) -> RocPoint:
    """The knee = the ROC point of maximum Youden-J (design §5). Ties break to
    the SMALLER margin (``-p.m``), i.e. favor recall / the more conservative
    gate. Raises if no point has a finite J (needs both TP and FP samples)."""
    valid = [p for p in points if not math.isnan(p.j)]
    if not valid:
        raise ValueError("no valid ROC points (need both TP and FP samples)")
    return max(valid, key=lambda p: (p.j, -p.m))     # tie -> smaller margin (favor recall)

def recommend(tp, fp) -> MarginRecommendation:
    """Turn TP/FP margin samples into a MarginRecommendation.

    The Youden knee is computed over the MULTI-LOCUS (finite-margin) TP
    population — ``tp_finite = [t for t in tp if t != inf]`` — vs ``fp`` (already
    all-finite by construction). Single-locus genes carry margin=inf and are
    trivially retained at any threshold, so including them inflates retention and
    pushes the knee to an artifact; design §2 defines the calibration on the
    >=2-paralog regime. Raises ``ValueError`` if there is no finite-margin TP or
    no FP (can't build a multi-locus ROC — surface it rather than fabricate one).

    The knee is a DIAGNOSTIC of separability, not always the recommendation:

    * CLEAN separation (knee J >= ``J_MIN_CLEAN``): the knee IS a usable hard
      gate, so ``min_margin = round(knee.m, 2)`` and ``low_margin_threshold`` is
      None (equal tiers).
    * OVERLAP (knee J < ``J_MIN_CLEAN``): the knee is NOT a usable hard gate — on
      overlapping TP/FP it lands where it would drop most correct multi-locus
      placements (and can invert the tiers). ``split_advisory`` is set and the
      recommendation switches to the design-§5 keep-and-flag posture: a LOW hard
      floor ``min_margin`` = p5 of the multi-locus TP (reject only near-ties) and
      a ``low_margin_threshold`` = p50 (flag placements below the typical
      multi-locus margin), guaranteed ``min_margin <= low_margin_threshold``. The
      tier split is gated on explicit user approval; the calibration only
      surfaces it.

    ``full_pop_retention`` — the fraction of ALL TP (incl. the single-locus inf)
    retained at the RECOMMENDED ``min_margin`` — is carried as real-world
    context."""
    tp_finite = [t for t in tp if t != INF]
    if not tp_finite:
        raise ValueError("no multi-locus (finite-margin) TP samples: the margin "
                         "gate is defined only on the >=2-paralog regime (design §2)")
    if not fp:
        raise ValueError("no FP margin samples: cannot build a multi-locus ROC")
    knee = youden_knee(roc_points(tp_finite, fp))       # DIAGNOSTIC of separability
    split = knee.j < J_MIN_CLEAN
    if split:
        # Overlap: decouple the recommendation from the (unusable-as-a-gate) knee
        # and emit the design-§5 keep-and-flag posture instead.
        rec_min_margin = round(_percentile(tp_finite, 5.0), 2)
        rec_low_margin = round(_percentile(tp_finite, 50.0), 2)
        if rec_min_margin > rec_low_margin:             # guarantee floor <= flag
            rec_min_margin = rec_low_margin
    else:
        rec_min_margin = round(knee.m, 2)
        rec_low_margin = None
    full_pop_retention = sum(1 for t in tp if t >= rec_min_margin) / len(tp)
    return MarginRecommendation(
        min_margin=rec_min_margin, knee=knee, split_advisory=split,
        low_margin_threshold=rec_low_margin,
        n_tp=len(tp), n_fp=len(fp), n_tp_multi=len(tp_finite),
        full_pop_retention=full_pop_retention)

def _percentile(values, p):     # linear-interp percentile, stdlib (mirrors the id/cov calibrator)
    s = sorted(float(v) for v in values)
    if not s: raise ValueError("percentile of empty")
    if len(s) == 1: return s[0]
    k = (len(s) - 1) * (p / 100.0); lo, hi = math.floor(k), math.ceil(k)
    return s[int(k)] if lo == hi else s[lo] * (hi - k) + s[hi] * (k - lo)


def load_cds_true_loci(cds_fasta: str) -> dict[str, tuple[str, int, int]]:
    """Map each RefSeq CDS record id (the full ``lcl|...`` header token) to its
    true genomic locus ``(scaffold, min, max)`` via ``parse_cds_true_locus`` —
    the refseq-mode ground truth. That id is the join key against the aligner
    output's query names, so it must match the PAF/GFF query ids exactly."""
    loci = {}
    with open(cds_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                qid = line[1:].split()[0]
                loci[qid] = parse_cds_true_locus(line)
    return loci

def _placements_by_query(paf, gff):
    """Parse the (optional) minimap2 PAF + GMAP GFF with reconcile_candidates'
    own parsers and bucket every placement by its query (transcript/CDS) id."""
    parsed = []
    if paf: parsed += rc.parse_minimap2_paf(paf)
    if gff: parsed += rc.parse_gmap_gff(gff)
    by_q: dict[str, list[rc.Placement]] = {}
    for p in parsed:
        by_q.setdefault(p.query, []).append(p)
    return by_q

def matched_fraction_warning(kind: str, matched: int, requested: int) -> str | None:
    """A stderr warning string when the requested ids only PARTIALLY join the
    alignment query names, else None. Mirrors the sibling id/cov calibrator
    (``calibrate_reconcile_thresholds.matched_fraction_warning``): a
    CDS-header ↔ PAF/GFF query-id format drift would otherwise silently compute a
    knee off a tiny (or empty) sample, so surface it and escalate the message
    when the matched fraction is very low."""
    if requested == 0 or matched >= requested:
        return None
    msg = f"WARNING: matched {matched} of {requested} requested {kind} ids"
    if matched == 0:
        return (msg + " — NONE matched; check the ids/CDS headers against the "
                "alignment query names (a format mismatch, e.g. version suffixes).")
    if matched / requested < 0.5:
        return (msg + " — LOW match fraction; likely an id-format mismatch "
                "between the CDS/ids and the PAF/GFF query names.")
    return msg

def collect_margins(by_query, true_loci, restrict_ids=None):
    """Collect the per-query TP and FP margins into ``(tp_list, fp_list)``.

    ``true_loci`` given ⇒ refseq-mode (header truth; requested = the CDS ids);
    ``true_loci=None`` ⇒ busco-mode (rep0 is truth; requested = ``restrict_ids``
    or every aligned query). Emits a stderr match-fraction WARNING (see
    ``matched_fraction_warning``) — surfacing the requested denominator — when
    the requested ids only partially join the alignment query names, so the knee
    is never computed off a silently-tiny sample after an id-format drift."""
    if restrict_ids is not None:
        requested = list(restrict_ids)
    elif true_loci is not None:
        requested = list(true_loci)
    else:
        requested = list(by_query)
    matched = sum(1 for q in requested if by_query.get(q))
    warning = matched_fraction_warning(
        "refseq" if true_loci is not None else "busco", matched, len(requested))
    if warning:
        print(warning, file=sys.stderr)

    tp, fp = [], []
    for q in requested:
        pls = by_query.get(q)
        if not pls:
            continue
        reps = distinct_locus_reps(pls)
        tl = None if true_loci is None else true_loci.get(q)
        if true_loci is not None and tl is None:
            continue    # no truth for this id
        obs = query_margins(reps, tl)
        if obs.tp is not None:
            tp.append(obs.tp)
        if obs.fp is not None:
            fp.append(obs.fp)
    return tp, fp

def _synthetic():
    # 6 clean genes (TP margin ~8-11) + 4 absent-gene grabs (FP margin ~1-3);
    # FP ceiling 3.0 ⇒ Youden knee at 3.1 (clean separation, well inside 3<m<8).
    tp = [8.2, 9.1, 10.4, 8.8, 11.0, 9.6]
    fp = [1.1, 2.0, 3.0, 1.5]
    return tp, fp

def format_table(rec: MarginRecommendation) -> str:
    """Render the recommendation as the aligned plain-text report printed to
    stdout: the TP breakdown (total / single-locus inf / multi-locus) + FP count,
    the MULTI-LOCUS knee (margin/retention/rejection/J), the full-population
    retention labeled as context, the clean-vs-overlap verdict, and the
    recommended ``GENOME_TRACK_MIN_MARGIN`` (plus ``GENOME_TRACK_LOW_MARGIN`` when
    a tier split is advised)."""
    n_single = rec.n_tp - rec.n_tp_multi
    lines = [
        "# Margin-gate calibration",
        f"TP (correct) n={rec.n_tp}  [single-locus/inf n={n_single}  "
        f"multi-locus n={rec.n_tp_multi}]   FP (spurious) n={rec.n_fp}",
        "",
        "knee is computed over the MULTI-LOCUS (finite-margin) population only "
        "(design §2: the margin gate acts only among >=2 competing loci; "
        "single-locus genes have margin=inf and are always retained).",
        f"multi-locus Youden knee (diagnostic): margin={rec.knee.m:g}  "
        f"retention={rec.knee.retention:.3f}  rejection={rec.knee.rejection:.3f}  "
        f"J={rec.knee.j:.3f}",
        f"full-population retention at the recommended min_margin={rec.min_margin:g} "
        f"(CONTEXT, incl. single-locus inf): {rec.full_pop_retention:.3f}",
        f"clean separation: {'NO (advise tier split)' if rec.split_advisory else 'YES'}"]
    if rec.split_advisory:
        lines.append(
            "  -> TP/FP margins OVERLAP: the Youden knee is NOT a usable hard "
            "gate (it would drop most correct multi-locus placements). "
            "Recommending the design-§5 keep-and-flag posture below.")
    lines += ["", "## Recommended (write into config.sh after review)",
              f"GENOME_TRACK_MIN_MARGIN={rec.min_margin:g}"]
    if rec.split_advisory:
        lines.append(f"GENOME_TRACK_LOW_MARGIN={rec.low_margin_threshold:g}")
    return "\n".join(lines)

def write_recommendation(path, rec, dry_run=False):
    """Write the recommendation as a sourceable shell file
    (``GENOME_TRACK_MIN_MARGIN=`` [+ ``GENOME_TRACK_LOW_MARGIN=`` when split],
    with commented provenance). ``dry_run`` stamps a ``NOT for production`` first
    line so a synthetic demo artifact can never be sourced into a real run
    (mirrors the sibling calibrator)."""
    hdr = ["# DRY RUN synthetic — NOT for production"] if dry_run else []
    hdr += ["# scripts/calibrate_genome_track_margin.py",
            f"# multi-locus Youden knee (diagnostic) J={rec.knee.j:.3f} at "
            f"margin={rec.knee.m:g}; clean={not rec.split_advisory}",
            f"# n_tp_multi={rec.n_tp_multi} of n_tp={rec.n_tp} "
            f"(single-locus/inf={rec.n_tp - rec.n_tp_multi}) vs n_fp={rec.n_fp} (design §2)"]
    if rec.split_advisory:
        hdr.append("# TP/FP margins OVERLAP: the recommendation is the design-§5 "
                   "keep-and-flag posture (low hard floor + median flag), NOT the knee.")
    hdr += [f"# full-population retention at min_margin (context)={rec.full_pop_retention:.3f}",
            f"GENOME_TRACK_MIN_MARGIN={rec.min_margin:g}"]
    if rec.split_advisory:
        hdr.append(f"GENOME_TRACK_LOW_MARGIN={rec.low_margin_threshold:g}")
    Path(path).write_text("\n".join(hdr) + "\n")

def write_figure(path, tp, fp, rec):
    """Write the margin-ROC PNG: TP-retention and FP-rejection vs ``m`` with the
    knee marked (design §5). matplotlib is imported lazily with the ``Agg``
    backend; if it is unavailable the figure is skipped with a stderr warning
    rather than crashing the calibration run."""
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"WARNING: matplotlib unavailable, skipping figure: {e}", file=sys.stderr); return
    tp_finite = [t for t in tp if t != INF]      # design §2: the gate acts only here
    pts = roc_points(tp_finite, fp)
    ms = [p.m for p in pts]
    plt.figure(figsize=(6, 4))
    plt.plot(ms, [p.retention for p in pts], label="multi-locus TP retention")
    plt.plot(ms, [p.rejection for p in pts], label="FP rejection")
    # The diagnostic Youden knee stays visible (light gray); the black lines are
    # the RECOMMENDED gate (min_margin) + the low_margin flag tier when split.
    plt.axvline(rec.knee.m, color="0.7", lw=0.8, label=f"Youden knee={rec.knee.m:g}")
    plt.axvline(rec.min_margin, ls="--", color="k", label=f"min_margin={rec.min_margin:g}")
    if rec.split_advisory:
        plt.axvline(rec.low_margin_threshold, ls=":", color="k",
                    label=f"low_margin flag={rec.low_margin_threshold:g}")
    plt.xlabel("min_margin (%id points)"); plt.ylabel("fraction")
    plt.title(f"multi-locus margin ROC (n_tp_multi={rec.n_tp_multi}, n_fp={rec.n_fp})")
    plt.legend(); plt.tight_layout()
    plt.savefig(path, dpi=150); plt.close()

def _build_arg_parser() -> argparse.ArgumentParser:
    """Build the argparse parser (factored out for testability, mirroring the
    sibling calibrator)."""
    p = argparse.ArgumentParser(
        prog="calibrate_genome_track_margin.py",
        description="Calibrate GENOME_TRACK_MIN_MARGIN from a ground-truthed "
                    "TP/FP best-vs-second distinct-locus margin ROC.")
    p.add_argument("--dry-run", action="store_true",
                   help="run the whole computation on a synthetic worked example "
                        "(no real data / Unity) and print the table")
    p.add_argument("--mode", choices=["refseq", "busco"], default="refseq",
                   help="refseq: true loci from --refseq-cds headers; busco: "
                        "rep0 is truth (single-copy), restricted to --busco-ids "
                        "(default: refseq)")
    p.add_argument("--minimap2-paf", help="transcript/CDS -> genome minimap2 PAF")
    p.add_argument("--gmap-gff", help="transcript/CDS -> genome GMAP gff3_gene")
    p.add_argument("--refseq-cds",
                   help="GPCR CDS FASTA whose headers give the true loci "
                        "(required for refseq-mode)")
    p.add_argument("--busco-ids",
                   help="file of BUSCO single-copy query ids, one per line, to "
                        "restrict busco-mode to")
    p.add_argument("--out-recommendation",
                   help="path to write the sourceable GENOME_TRACK_MIN_MARGIN "
                        "(+ GENOME_TRACK_LOW_MARGIN when split) recommendation")
    p.add_argument("--out-figure",
                   help="path to write the margin-ROC PNG (matplotlib)")
    return p


def main(argv=None) -> int:
    """CLI: build the TP/FP margin samples (synthetic ``--dry-run``, or from the
    aligner outputs in refseq/busco mode), compute the Youden-knee
    recommendation, print the table, and optionally write the recommendation file
    + ROC figure. Returns a process exit code (0 ok; 1 insufficient data; 2 bad
    arguments)."""
    a = _build_arg_parser().parse_args(argv)

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

    # A margin ROC needs a MULTI-LOCUS (finite-margin) TP sample AND an FP
    # sample (design §2); single-locus inf margins alone cannot place a knee.
    # Pre-validate here so the CLI exits cleanly rather than letting
    # recommend()'s ValueError surface as a traceback.
    tp_finite = [t for t in tp if t != INF]
    if not tp_finite or not fp:
        print(f"ERROR: insufficient multi-locus data for a margin ROC "
              f"(multi-locus TP={len(tp_finite)}, FP={len(fp)}; total TP={len(tp)}). "
              f"The margin gate is defined on the >=2-paralog regime (design §2).",
              file=sys.stderr)
        return 1
    rec = recommend(tp, fp)
    print(format_table(rec))
    if a.out_recommendation: write_recommendation(a.out_recommendation, rec, dry_run=a.dry_run)
    if a.out_figure: write_figure(a.out_figure, tp, fp, rec)
    return 0

if __name__ == "__main__":      # pragma: no cover
    raise SystemExit(main())
