#!/usr/bin/env python3
"""reconcile_candidates.py — genome-track candidate reconciliation.

Merges the genome-derived (Berghia RefSeq protein) chemoreceptor
candidate track with the transcriptome-derived track into one
non-redundant, genome-anchored, provenance-labeled gene set.

Task 1 (this pass) implements only the first stage of the pipeline:

    parse placements -> %match gate

i.e. the ``Placement`` model, the four aligner-output parsers
(minimap2 PAF, GMAP GFF3, miniprot GFF3, BLASTp outfmt6+qcovs), and
the absolute + best-vs-second-margin gate functions. Locus grouping,
isoform collapse, representative selection, the placement cascade,
and the output writers are added in later tasks against this same
module.

All functions here are pure (stdlib only) beyond reading the given
input file paths — no network access, no subprocess calls.

Design spec: docs/plans/2026-06-30-berghia-genome-track-reconciliation-design.md
"""
from __future__ import annotations

import os
from dataclasses import dataclass

# miniprot's native GFF3 mRNA row has no coverage-equivalent attribute
# (the `Target=<id> <start> <end>` range alone doesn't reveal the full
# query protein length, so %coverage isn't derivable from that row
# alone). When no coverage attribute is present at all, pct_coverage
# falls back to this sentinel. It is deliberately below any realistic
# min_cov gate value, so a miniprot placement without corroborating
# coverage information conservatively fails pass_gate() and the
# placement cascade (Task 4) falls through to the next method rather
# than trusting an unverified placement.
MINIPROT_COVERAGE_SENTINEL = 0.0


@dataclass(frozen=True)
class Placement:
    """A single candidate-to-genome placement from one aligner.

    Frozen (immutable + hashable): a Placement is a value record, never
    mutated after construction, and Task 2's locus dedup needs it
    hashable — matching the repo's value-record convention.
    """
    query: str
    chrom: str
    start: int
    end: int
    strand: str
    pct_identity: float
    pct_coverage: float
    method: str
    score: float


# ---- shared helpers --------------------------------------------------------

def _parse_gff_attributes(field: str) -> dict[str, str]:
    """Parse a GFF3 column-9 attribute string ("key=value;key=value")
    into a dict. Tokens without '=' are ignored rather than raising."""
    attrs: dict[str, str] = {}
    for tok in field.split(";"):
        tok = tok.strip()
        if not tok or "=" not in tok:
            continue
        key, _, value = tok.partition("=")
        attrs[key.strip()] = value.strip()
    return attrs


def _attr(attrs: dict[str, str], *names: str) -> str | None:
    """Look up the first present key among `names` (case-sensitive
    first, then a case-insensitive fallback pass), to tolerate minor
    capitalization drift across tool versions without guessing at
    unrelated synonyms."""
    for name in names:
        if name in attrs:
            return attrs[name]
    lowered = {k.lower(): v for k, v in attrs.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


# ---- parse_minimap2_paf -----------------------------------------------------

def parse_minimap2_paf(path: str) -> list[Placement]:
    """Parse minimap2 PAF output into Placements.

    Column layout (0-based indices): qname[0] qlen[1] qstart[2] qend[3]
    strand[4] tname[5] tlen[6] tstart[7] tend[8] nmatch[9] alen[10]
    mapq[11] ...(optional SAM-like tags).

    identity = 100 * nmatch / alen; coverage = 100 * (qend-qstart) / qlen.
    """
    if not os.path.exists(path):
        return []

    out: list[Placement] = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            try:
                qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
                tstart, tend = int(f[7]), int(f[8])
                nmatch, alen = int(f[9]), int(f[10])
                mapq = float(f[11])
            except ValueError:
                continue
            pct_identity = 100.0 * nmatch / alen if alen else 0.0
            pct_coverage = 100.0 * (qe - qs) / qlen if qlen else 0.0
            out.append(Placement(query=f[0], chrom=f[5], start=tstart, end=tend,
                                  strand=f[4], pct_identity=pct_identity,
                                  pct_coverage=pct_coverage, method="minimap2",
                                  score=mapq))
    return out


# ---- parse_gmap_gff ---------------------------------------------------------

def parse_gmap_gff(path: str) -> list[Placement]:
    """Parse GMAP GFF3 (gff3_gene format) output into Placements.

    GMAP writes `coverage=` and `identity=` (both already percentages)
    as attributes on the row that summarizes an alignment path — in
    standard `gff3_gene` output that is the mRNA row; the parent gene
    row does not carry them. Rather than special-case a feature-type
    preference, this parser accepts either a `gene` or `mRNA` row and
    only emits a Placement for rows that actually carry both
    attributes — so gene summary rows (and exon/CDS detail rows) are
    skipped naturally, without assuming a fixed row order.
    """
    if not os.path.exists(path):
        return []

    out: list[Placement] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            feature = f[2]
            if feature not in ("mRNA", "gene"):
                continue
            attrs = _parse_gff_attributes(f[8])
            identity_s = _attr(attrs, "identity")
            coverage_s = _attr(attrs, "coverage")
            if identity_s is None or coverage_s is None:
                continue
            try:
                start, end = int(f[3]), int(f[4])
                pct_identity = float(identity_s)
                pct_coverage = float(coverage_s)
            except ValueError:
                continue
            # Fail closed on the join key: use the bare transcript id from
            # Name= or Target=; never fall back to ID= (which is
            # `<id>.pathN.mrnaN` — a silently-wrong key that would break
            # Task 4 concordance). No usable id -> skip the row (the
            # cascade then falls through to miniprot).
            query = _attr(attrs, "Name") or _attr(attrs, "Target")
            if query is None:
                continue
            query = query.split()[0]  # Target="<id> start end [strand]" -> id
            score_s = _attr(attrs, "matches")
            try:
                score = float(score_s) if score_s is not None else 0.0
            except ValueError:
                score = 0.0
            out.append(Placement(query=query, chrom=f[0], start=start, end=end,
                                  strand=f[6], pct_identity=pct_identity,
                                  pct_coverage=pct_coverage, method="gmap",
                                  score=score))
    return out


# ---- parse_miniprot_gff ------------------------------------------------------

def parse_miniprot_gff(path: str) -> list[Placement]:
    """Parse miniprot GFF3 output into Placements.

    miniprot emits `Identity=` (a 0-1 fraction, converted to percent
    here) and `Target=<protein_id> <start> <end> [strand]` on its
    mRNA row; column 6 carries a numeric alignment score. It has no
    native coverage-equivalent attribute, so pct_coverage uses an
    explicit `Coverage=`/`coverage=` attribute when present (not part
    of stock miniprot output, but tolerated for forward-compat), else
    the documented MINIPROT_COVERAGE_SENTINEL (see module docstring)
    so unverified placements conservatively fail the %match gate.
    """
    if not os.path.exists(path):
        return []

    out: list[Placement] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            if f[2] != "mRNA":
                continue
            attrs = _parse_gff_attributes(f[8])
            identity_s = _attr(attrs, "Identity")
            target = _attr(attrs, "Target")
            if identity_s is None or not target:
                continue
            try:
                start, end = int(f[3]), int(f[4])
                pct_identity = 100.0 * float(identity_s)
            except ValueError:
                continue
            query = target.split()[0]
            coverage_s = _attr(attrs, "Coverage", "coverage")
            try:
                pct_coverage = (float(coverage_s) if coverage_s is not None
                                 else MINIPROT_COVERAGE_SENTINEL)
            except ValueError:
                pct_coverage = MINIPROT_COVERAGE_SENTINEL
            try:
                score = float(f[5])
            except ValueError:
                score = 0.0
            out.append(Placement(query=query, chrom=f[0], start=start, end=end,
                                  strand=f[6], pct_identity=pct_identity,
                                  pct_coverage=pct_coverage, method="miniprot",
                                  score=score))
    return out


# ---- parse_blastp_tab --------------------------------------------------------

def parse_blastp_tab(path: str, coords: dict[str, tuple[str, int, int, str]]
                      ) -> list[Placement]:
    """Parse BLASTp outfmt6 (+ qcovs) output into Placements.

    Columns: qseqid sseqid pident length mismatch gapopen qstart qend
    sstart send evalue bitscore qcovs (13 tab-separated fields).

    `coords` maps a subject RefSeq protein id to its gene locus
    `(chrom, start, end, strand)`. A subject id absent from `coords`
    is skipped (conservative — no genomic coordinates to place at).
    """
    if not os.path.exists(path):
        return []

    out: list[Placement] = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            qseqid, sseqid = f[0], f[1]
            locus = coords.get(sseqid)
            if locus is None:
                continue
            try:
                pident = float(f[2])
                bitscore = float(f[11])
                qcovs = float(f[12])
            except ValueError:
                continue
            chrom, start, end, strand = locus
            out.append(Placement(query=qseqid, chrom=chrom, start=start, end=end,
                                  strand=strand, pct_identity=pident,
                                  pct_coverage=qcovs, method="blastp",
                                  score=bitscore))
    return out


# ---- %match gate -------------------------------------------------------------

def pass_gate(p: Placement | None, min_id: float, min_cov: float) -> bool:
    """True iff `p` clears both the identity and coverage bars."""
    return p is not None and p.pct_identity >= min_id and p.pct_coverage >= min_cov


def best_and_margin(placements: list[Placement | None]
                     ) -> tuple[Placement | None, float]:
    """Return the best placement (by pct_identity) and its margin over
    the runner-up. margin = inf when there is exactly one placement;
    (None, 0.0) when there are none."""
    ps = sorted([p for p in placements if p], key=lambda p: p.pct_identity, reverse=True)
    if not ps:
        return None, 0.0
    if len(ps) == 1:
        return ps[0], float("inf")
    return ps[0], ps[0].pct_identity - ps[1].pct_identity
