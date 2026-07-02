#!/usr/bin/env python3
"""calibrate_reconcile_thresholds.py — genome-track %match gate calibration.

Finalizes the ``reconcile_candidates.py`` accept bar (``--min-id`` /
``--min-cov``, wired through config.sh's ``GENOME_TRACK_MIN_ID`` /
``GENOME_TRACK_MIN_COV``) from real Berghia data instead of the design's
working hypothesis (95 / 90).

Why a stringent bar (design §4/§9): a Berghia transcript is aligned to the
Berghia genome — its OWN species — so a true same-gene match is near
identical. Anything materially below that is a paralog cross-map, a novel
gene, or a real absence. The calibration makes the bar defensible:

  * true-same-gene distribution — the transcript->genome %id/%cov of KNOWN
    single-copy genes (the BUSCO single-copy set). Each such transcript has
    exactly one genomic copy, so its best placement IS its own gene; the
    distribution's floor (5th percentile) is how low a genuine same-gene
    match can go.
  * paralog-family %id ceiling — for the paralog families (the chemoreceptor
    candidate transcripts), each transcript's NON-self placements are
    cross-maps onto sibling paralog loci; their %id ceiling (a high
    percentile) is how high a paralog cross-map can reach.

A clean bar sits BELOW the true-same-gene floor and ABOVE the paralog ceiling.
The script prints the separation table and a recommended ``(min_id, min_cov)``;
the operator writes the finalized values into config.sh. (This task leaves the
config at the working hypothesis; the real run finalizes it.)

REAL-mode inputs (pre-computed, so no aligner is re-run here and no alignment
math is reimplemented — the parsers are reused from reconcile_candidates.py):

    --minimap2-paf / --gmap-gff   transcript->genome alignments of BOTH the
                                  BUSCO single-copy transcripts and the
                                  paralog-family (chemoreceptor candidate)
                                  transcripts, produced exactly as stage 02c
                                  does, e.g.:
                                    minimap2 -x splice <genome> <transcripts> > aln.paf
                                    gmap -d <db> -f gff3_gene <transcripts> > aln.gff3
    --busco-single-copy-ids FILE  one transcript id per line — the true-
                                  same-gene set (from the stage-03a BUSCO run's
                                  single_copy_busco_sequences).
    --paralog-ids FILE            one transcript id per line — the paralog
                                  families (the stage-02 chemoreceptor
                                  candidates). Optional; without it the paralog
                                  ceiling is not estimated.

``--dry-run`` runs the whole computation on a synthetic worked example (no
Unity, no real data) and prints the separation table — a self-contained
demonstration of the method.

Design spec: docs/plans/2026-06-30-berghia-genome-track-reconciliation-design.md (§4, §9)
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from dataclasses import dataclass

import reconcile_candidates as rc

# The design §4 working hypothesis (mirrors config.sh's GENOME_TRACK_MIN_ID /
# _MIN_COV defaults). Used as the fallback recommendation when the data show
# NO clean gap between paralogs and true same-gene matches — the calibrator
# never invents a bar it cannot justify.
WORKING_MIN_ID = 95.0
WORKING_MIN_COV = 90.0
# Percentile used to summarize each distribution: the 5th percentile is the
# true-same-gene FLOOR (the bar must be below it), and the 95th percentile is
# the paralog CEILING (the bar must be above it).
FLOOR_PCT = 5.0
DEFAULT_PARALOG_PCT = 95.0
# Minimum gap (percentage points) between the paralog ceiling and the true-
# same-gene floor for the separation to count as CLEAN. Below this the two
# distributions effectively touch: the calibrator refuses to recommend a
# tightened bar (which would also risk rounding the recommendation onto a
# boundary), and falls back to the working hypothesis.
MIN_SEPARATION_GAP = 0.1


@dataclass(frozen=True)
class Distribution:
    """Three summary statistics of a %id or %cov distribution."""
    minimum: float
    median: float
    p05: float


@dataclass(frozen=True)
class SeparationResult:
    """The calibration verdict: the two distributions, the paralog ceiling, and
    the recommended accept bar (see module docstring for the decision logic)."""
    id_dist: Distribution
    cov_dist: Distribution
    paralog_ceiling: float | None
    paralog_pct: float
    n_same_gene: int
    n_paralog: int
    clean_separation: bool
    rec_min_id: float
    rec_min_cov: float
    separation_margin: float | None   # id_dist.p05 - paralog_ceiling (None if no paralog data)


# ---- statistics -------------------------------------------------------------

def _percentile(values, p: float) -> float:
    """Linear-interpolated percentile (the numpy default method), stdlib only.

    ``p`` is in [0, 100]. Raises on empty input (a percentile of nothing is
    undefined — fail loud rather than return a fabricated number)."""
    if not values:
        raise ValueError("percentile of an empty sequence is undefined")
    s = sorted(float(v) for v in values)
    if len(s) == 1:
        return s[0]
    k = (len(s) - 1) * (p / 100.0)
    lo = math.floor(k)
    hi = math.ceil(k)
    if lo == hi:
        return s[int(k)]
    return s[lo] * (hi - k) + s[hi] * (k - lo)


def _distribution(values) -> Distribution:
    """Summarize a distribution as (min, median, 5th-percentile)."""
    return Distribution(minimum=min(values),
                        median=_percentile(values, 50.0),
                        p05=_percentile(values, FLOOR_PCT))


# ---- distinct-locus reduction (self vs paralog cross-maps) -----------------

def _distinct_locus_reps(placements: list[rc.Placement]) -> list[rc.Placement]:
    """Reduce one transcript's placements to one representative per DISTINCT
    genomic locus, highest identity first.

    Delegates to reconcile_candidates.py's OWN locus clustering
    (``rc._cluster_placements`` — connected components of same-chrom/strand
    overlap, which is TRANSITIVE), so "distinct locus" means exactly what the
    reconciliation means: an A–B–C overlap chain is ONE locus even though A and
    C are disjoint. A greedy pairwise scan would instead keep C as a spurious
    second locus and inflate the paralog ceiling with a within-gene alignment —
    the very quantity this script calibrates the bar against.

    Each cluster is represented by its highest-identity placement
    (``rc._identity_rank`` orders identity-descending), and the reps are sorted
    the same way: ``reps[0]`` is the self (own-gene) locus, ``reps[1:]`` are the
    paralog cross-maps."""
    return sorted((min(cluster, key=rc._identity_rank)
                   for cluster in rc._cluster_placements(placements)),
                  key=rc._identity_rank)


def same_gene_reps(placements_by_query: dict[str, list[rc.Placement]],
                   ids) -> list[rc.Placement]:
    """The self (best-locus) placement of each BUSCO single-copy transcript —
    the true-same-gene sample. Ids absent from the alignments are skipped."""
    reps: list[rc.Placement] = []
    for q in ids:
        pls = placements_by_query.get(q)
        if not pls:
            continue
        reps.append(_distinct_locus_reps(pls)[0])
    return reps


def paralog_identities(placements_by_query: dict[str, list[rc.Placement]],
                       ids) -> list[float]:
    """The identities of every NON-self (paralog cross-map) placement across the
    paralog-family transcripts — the paralog sample. A transcript with a single
    locus (no sibling cross-map) contributes nothing."""
    out: list[float] = []
    for q in ids:
        pls = placements_by_query.get(q)
        if not pls:
            continue
        out.extend(r.pct_identity for r in _distinct_locus_reps(pls)[1:])
    return out


# ---- separation + recommendation -------------------------------------------

def compute_separation(same_reps: list[rc.Placement],
                       paralog_ids: list[float],
                       paralog_pct: float = DEFAULT_PARALOG_PCT,
                       working_min_id: float = WORKING_MIN_ID,
                       working_min_cov: float = WORKING_MIN_COV) -> SeparationResult:
    """Turn the two samples into distributions + a recommended accept bar.

    A clean bar exists iff the paralog ceiling (``paralog_pct``-th percentile of
    the cross-map identities) is strictly below the true-same-gene identity
    floor (5th percentile). When it is, ``min_id`` is placed at their midpoint
    and ``min_cov`` just below the true-same-gene coverage floor. When it is NOT
    (paralogs reach into the same-gene range, or no paralog data was supplied),
    the recommendation falls back to the working hypothesis and
    ``clean_separation`` is False so the operator sees the calibration could not
    tighten the bar."""
    if not same_reps:
        raise ValueError("no BUSCO single-copy self-placements: cannot calibrate")

    id_dist = _distribution([r.pct_identity for r in same_reps])
    cov_dist = _distribution([r.pct_coverage for r in same_reps])

    if paralog_ids:
        ceiling: float | None = _percentile(paralog_ids, paralog_pct)
        margin: float | None = id_dist.p05 - ceiling
        # A clean bar needs a REAL gap, not a hairline touch (< MIN_SEPARATION_GAP
        # is effectively no separation, and would risk the recommendation below
        # rounding onto the ceiling).
        clean = margin >= MIN_SEPARATION_GAP
    else:
        ceiling = None
        margin = None
        clean = False

    if clean:
        rec_min_id = round((ceiling + id_dist.p05) / 2.0, 1)
        # Rounding to 1 dp can land ON a boundary; force the bar STRICTLY inside
        # (ceiling < min_id < floor) via the exact midpoint, which always is.
        if not (ceiling < rec_min_id < id_dist.p05):
            rec_min_id = (ceiling + id_dist.p05) / 2.0
        rec_min_cov = float(math.floor(cov_dist.p05))
    else:
        rec_min_id = working_min_id
        rec_min_cov = working_min_cov

    return SeparationResult(
        id_dist=id_dist, cov_dist=cov_dist, paralog_ceiling=ceiling,
        paralog_pct=paralog_pct, n_same_gene=len(same_reps),
        n_paralog=len(paralog_ids), clean_separation=clean,
        rec_min_id=rec_min_id, rec_min_cov=rec_min_cov,
        separation_margin=margin)


def format_table(result: SeparationResult, title: str = "Separation table") -> str:
    """Render the separation table + recommendation as aligned plain text."""
    r = result
    ceil_s = "n/a (no paralog data)" if r.paralog_ceiling is None else f"{r.paralog_ceiling:g}"
    margin_s = "n/a" if r.separation_margin is None else f"{r.separation_margin:g}"
    lines = [
        f"# {title}",
        "",
        f"true-same-gene set (BUSCO single-copy): n = {r.n_same_gene}",
        f"paralog cross-map set:                  n = {r.n_paralog}",
        "",
        "metric                                    value",
        f"true-same-gene identity  min             {r.id_dist.minimum:g}",
        f"true-same-gene identity  median          {r.id_dist.median:g}",
        f"true-same-gene identity  p05 (FLOOR)      {r.id_dist.p05:g}",
        f"true-same-gene coverage  min             {r.cov_dist.minimum:g}",
        f"true-same-gene coverage  median          {r.cov_dist.median:g}",
        f"true-same-gene coverage  p05 (FLOOR)      {r.cov_dist.p05:g}",
        f"paralog identity p{r.paralog_pct:g} (CEILING)          {ceil_s}",
        f"separation margin (floor - ceiling)      {margin_s}",
        "",
        f"clean separation: {'YES' if r.clean_separation else 'NO'}",
    ]
    if not r.clean_separation:
        lines.append("  -> NO clean gap; falling back to the working hypothesis "
                     "(paralogs reach the same-gene range, or no paralog data).")
    lines += [
        "",
        "## Recommended accept bar (write into config.sh)",
        f"GENOME_TRACK_MIN_ID={r.rec_min_id:g}",
        f"GENOME_TRACK_MIN_COV={r.rec_min_cov:g}",
        "",
        "min_id sits between the paralog ceiling and the true-same-gene floor; "
        "min_cov sits at/below the true-same-gene coverage floor.",
    ]
    return "\n".join(lines)


# ---- synthetic worked example (--dry-run) ----------------------------------

def synthetic_inputs():
    """A synthetic worked example: BUSCO single-copy transcripts (one clean
    self-hit each, ~98-100% id) and paralog-family transcripts (a near-perfect
    self-hit plus ~88-93% cross-maps to two sibling loci). Returns
    ``(placements_by_query, busco_ids, paralog_ids)`` — the same shape the REAL
    parsers produce, so the dry-run exercises the true computation path."""
    placements: dict[str, list[rc.Placement]] = {}

    busco_specs = [(99.8, 99.0), (99.5, 98.0), (99.1, 97.0), (98.7, 96.0),
                   (98.2, 95.0), (100.0, 99.0), (99.0, 98.0), (97.9, 94.0),
                   (98.5, 96.0), (99.3, 97.0)]
    busco_ids = []
    for i, (pid, cov) in enumerate(busco_specs):
        q = f"BUSCO_{i}"
        base = 1000 + i * 20000
        placements[q] = [rc.Placement(q, "chr1", base, base + 1500, "+",
                                      pid, cov, "minimap2", 60.0)]
        busco_ids.append(q)

    # self id, (cross-map id, cross-map id) at two sibling loci.
    para_specs = [(98.0, (90.0, 88.0)), (97.6, (91.0, 89.0)),
                  (98.5, (92.0, 90.0)), (99.0, (89.0, 87.0)),
                  (97.8, (93.0, 91.0)), (98.2, (90.0, 88.0))]
    paralog_ids = []
    for j, (self_id, (c1, c2)) in enumerate(para_specs):
        q = f"PARA_{j}"
        base = 1_000_000 + j * 20000
        placements[q] = [
            rc.Placement(q, "chr2", base, base + 1200, "+", self_id, 96.0, "minimap2", 60.0),
            rc.Placement(q, "chr2", base + 4000, base + 5200, "+", c1, 95.0, "minimap2", 60.0),
            rc.Placement(q, "chr2", base + 8000, base + 9200, "+", c2, 95.0, "minimap2", 60.0),
        ]
        paralog_ids.append(q)

    return placements, busco_ids, paralog_ids


# ---- CLI -------------------------------------------------------------------

def _read_ids(path: str) -> list[str]:
    """One id per line; blank and ``#``-commented lines skipped; the first
    whitespace token of each line is taken (tolerating trailing columns)."""
    ids: list[str] = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.append(s.split()[0])
    return ids


def matched_fraction_warning(kind: str, matched: int, requested: int) -> str | None:
    """A stderr warning when the requested ids only PARTIALLY match the
    alignment query names, else None. A format mismatch between the ids file and
    the PAF/GFF query names would otherwise silently compute a confident bar off
    a tiny (or empty) distribution — so surface it, and escalate the message
    when the matched fraction is very low."""
    if requested == 0 or matched >= requested:
        return None
    msg = f"WARNING: matched {matched} of {requested} requested {kind} ids"
    if matched == 0:
        return (msg + " — NONE matched; check the ids file against the alignment "
                "query names (a format mismatch, e.g. version suffixes).")
    if matched / requested < 0.5:
        return (msg + " — LOW match fraction; likely an id-format mismatch "
                "between the ids file and the PAF/GFF query names.")
    return msg


def _placements_by_query(minimap2_paf, gmap_gff) -> dict[str, list[rc.Placement]]:
    """Parse the (optional) alignment files with reconcile_candidates.py's own
    parsers and bucket every placement by its transcript id."""
    parsed: list[rc.Placement] = []
    if minimap2_paf:
        parsed += rc.parse_minimap2_paf(minimap2_paf)
    if gmap_gff:
        parsed += rc.parse_gmap_gff(gmap_gff)
    by_query: dict[str, list[rc.Placement]] = defaultdict(list)
    for p in parsed:
        by_query[p.query].append(p)
    return by_query


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="calibrate_reconcile_thresholds.py",
        description="Calibrate the genome-track %match accept bar (min_id / "
                    "min_cov) from the Berghia BUSCO single-copy set + paralog "
                    "cross-maps. Prints a separation table + recommendation.")
    p.add_argument("--dry-run", action="store_true",
                   help="run the whole computation on a synthetic worked "
                        "example (no real data / Unity) and print the table")
    p.add_argument("--minimap2-paf", help="transcript->genome minimap2 PAF")
    p.add_argument("--gmap-gff", help="transcript->genome GMAP gff3_gene")
    p.add_argument("--busco-single-copy-ids",
                   help="file of BUSCO single-copy transcript ids (true "
                        "same-gene set), one per line")
    p.add_argument("--paralog-ids",
                   help="file of paralog-family (chemoreceptor candidate) "
                        "transcript ids, one per line (for the paralog ceiling)")
    p.add_argument("--paralog-pct", type=float, default=DEFAULT_PARALOG_PCT,
                   help="percentile for the paralog identity ceiling "
                        f"(default {DEFAULT_PARALOG_PCT:g})")
    p.add_argument("--out",
                   help="optional path to write the recommended cutoffs "
                        "(sourceable GENOME_TRACK_MIN_ID / _MIN_COV lines)")
    return p


def _write_recommendation(path: str, result: SeparationResult,
                          dry_run: bool = False) -> None:
    """Write the recommended cutoffs as sourceable shell, citing the numbers.

    In ``--dry-run`` the values come from the SYNTHETIC example, so the file is
    stamped ``NOT for production`` on its first line to keep a demo artifact from
    ever being sourced into a real run."""
    ceil_s = "n/a" if result.paralog_ceiling is None else f"{result.paralog_ceiling:g}"
    header = ["# DRY RUN synthetic — NOT for production"] if dry_run else []
    header += [
        "# Recommended genome-track %match gate, from "
        "scripts/calibrate_reconcile_thresholds.py",
        f"# true-same-gene identity floor (p05) = {result.id_dist.p05:g}; "
        f"paralog identity ceiling (p{result.paralog_pct:g}) = {ceil_s}; "
        f"clean_separation = {result.clean_separation}",
        f"GENOME_TRACK_MIN_ID={result.rec_min_id:g}",
        f"GENOME_TRACK_MIN_COV={result.rec_min_cov:g}",
    ]
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(header) + "\n")


def main(argv=None) -> int:
    args = _build_arg_parser().parse_args(argv)

    if args.dry_run:
        placements, busco_ids, paralog_ids = synthetic_inputs()
        result = compute_separation(same_gene_reps(placements, busco_ids),
                                    paralog_identities(placements, paralog_ids),
                                    paralog_pct=args.paralog_pct)
        print(format_table(result, title="DRY RUN — synthetic worked example"))
        if args.out:
            _write_recommendation(args.out, result, dry_run=True)
        return 0

    if not args.busco_single_copy_ids:
        print("ERROR: --busco-single-copy-ids is required in real mode "
              "(or use --dry-run).", file=sys.stderr)
        return 2
    if not (args.minimap2_paf or args.gmap_gff):
        print("ERROR: at least one of --minimap2-paf / --gmap-gff is required "
              "in real mode (or use --dry-run).", file=sys.stderr)
        return 2

    by_query = _placements_by_query(args.minimap2_paf, args.gmap_gff)
    busco_ids = _read_ids(args.busco_single_copy_ids)
    para_ids = _read_ids(args.paralog_ids) if args.paralog_ids else []

    reps = same_gene_reps(by_query, busco_ids)
    # len(reps) == the number of requested BUSCO ids that matched an alignment
    # (one self-rep per matched id) — surface a partial/empty match.
    w = matched_fraction_warning("same-gene", len(reps), len(busco_ids))
    if w:
        print(w, file=sys.stderr)
    if not reps:
        print("ERROR: none of the BUSCO single-copy ids were found in the "
              "alignment files; cannot calibrate.", file=sys.stderr)
        return 1
    para = paralog_identities(by_query, para_ids)
    matched_para = sum(1 for q in para_ids if by_query.get(q))
    w = matched_fraction_warning("paralog", matched_para, len(para_ids))
    if w:
        print(w, file=sys.stderr)

    result = compute_separation(reps, para, paralog_pct=args.paralog_pct)
    print(format_table(result))
    if args.out:
        _write_recommendation(args.out, result)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
