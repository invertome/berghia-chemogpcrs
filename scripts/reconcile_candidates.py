#!/usr/bin/env python3
"""reconcile_candidates.py — genome-track candidate reconciliation.

Merges the genome-derived (Berghia RefSeq protein) chemoreceptor
candidate track with the transcriptome-derived track into one
non-redundant, genome-anchored, provenance-labeled gene set.

Task 1 (this pass) implements only the first stage of the pipeline:

    parse placements -> %match gate

i.e. the ``Placement`` model, the four aligner-output parsers
(minimap2 PAF, GMAP GFF3, miniprot GFF3, BLASTp outfmt6+qcovs), and
the absolute + best-vs-second-margin gate functions. Task 2 adds
locus grouping + isoform collapse (``PlacedModel``, ``Locus``,
``group_into_loci``). Representative selection, the placement
cascade, and the output writers are added in later tasks against
this same module.

All functions here are stdlib-only with no network access. The one
optional subprocess is ``group_into_loci``'s ``bedtools merge`` fast
path, used only when the binary is available and always backed by an
equivalent pure-Python interval sweep (so the tests need no binary).

Design spec: docs/plans/2026-06-30-berghia-genome-track-reconciliation-design.md
"""
from __future__ import annotations

import os
import shutil
import subprocess
from collections import defaultdict
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


# ---- locus grouping + isoform collapse (Task 2) ------------------------------

BEDTOOLS_DEFAULT = "bedtools"


@dataclass(frozen=True)
class PlacedModel:
    """A candidate gene model positioned on the genome.

    Distinct from ``Placement`` (which normalizes one aligner hit): a
    ``PlacedModel`` is a candidate gene model with settled genomic
    coordinates, tagged with the track it came from — ``source`` is a
    free-form string whose downstream values are ``"genome"`` and
    ``"txome"``. These are the units grouped into loci. Frozen +
    hashable, matching the module's value-record convention; ``source``
    must be carried through grouping unchanged (Task 3 reads it to
    label provenance).

    ``complete`` (True = complete ORF) and ``length`` (protein length in
    aa) are the two axes Task 3's representative selection ranks on. Both
    are appended with defaults so every Task-2 construction (which omits
    them) stays valid and unchanged.
    """
    query: str
    chrom: str
    start: int
    end: int
    strand: str
    source: str
    complete: bool = False
    length: int = 0


@dataclass(frozen=True)
class Locus:
    """One locus (= gene): a same-chromosome, same-strand overlap cluster
    of ``PlacedModel``s. ``start``/``end`` span ``min(member starts)`` /
    ``max(member ends)``; ``members`` is the tuple of models in the locus,
    sorted by ``query`` for deterministic output. A tuple (not a list)
    keeps ``Locus`` genuinely immutable and hashable — matching the
    ``Placement``/``PlacedModel`` value-record contract, so Tasks 3-4 can
    use loci in sets/dicts.
    """
    chrom: str
    start: int
    end: int
    strand: str
    members: tuple[PlacedModel, ...]


def group_into_loci(records: list[PlacedModel]) -> list[Locus]:
    """Group placed candidate models into loci (= genes) by same-strand
    genomic overlap, collapsing isoforms.

    Two ``PlacedModel``s share a locus iff they are on the same
    chromosome and strand and share >0 bases under half-open
    ``[start, end)`` semantics (``a.start < b.end and b.start < a.end``).
    Grouping is transitive (connected components), so a chain of pairwise
    overlaps forms one locus, while book-ended intervals (touching at a
    boundary, 0 shared bases) stay separate.

    Uses ``bedtools merge`` when the binary named by ``$BEDTOOLS``
    (default ``"bedtools"``) is on PATH, else an equivalent pure-Python
    interval sweep; both paths return identical loci. Output is
    deterministic: loci sorted by ``(chrom, start, strand)`` and each
    locus's ``members`` sorted by ``query``.
    """
    if not records:
        return []
    binary = os.environ.get("BEDTOOLS", BEDTOOLS_DEFAULT)
    if shutil.which(binary) is not None:
        clusters = _clusters_via_bedtools(records, binary)
    else:
        clusters = _clusters_via_sweep(records)
    loci = [_locus_from_members(members) for members in clusters]
    loci.sort(key=lambda locus: (locus.chrom, locus.start, locus.strand))
    return loci


def _locus_from_members(members: list[PlacedModel]) -> Locus:
    """Build a Locus from a non-empty cluster of same-chrom/same-strand
    models: span ``min(start)``/``max(end)`` with members sorted by query.
    chrom/strand are taken from the first member (uniform by construction).
    """
    return Locus(chrom=members[0].chrom,
                 start=min(m.start for m in members),
                 end=max(m.end for m in members),
                 strand=members[0].strand,
                 members=tuple(sorted(members, key=lambda m: m.query)))


def _clusters_via_sweep(records: list[PlacedModel]) -> list[list[PlacedModel]]:
    """Pure-Python clustering: within each ``(chrom, strand)`` group, sort
    by start and sweep, folding in any model that overlaps the running
    merged span (``model.start < cluster_end``). Half-open semantics — a
    model whose start equals the cluster end (book-ended) opens a new
    cluster. Equivalent to connected components of the overlap graph
    because the swept span stays contiguous (each folded-in model starts
    before the current right edge)."""
    by_key: dict[tuple[str, str], list[PlacedModel]] = defaultdict(list)
    for r in records:
        by_key[(r.chrom, r.strand)].append(r)

    clusters: list[list[PlacedModel]] = []
    for key in by_key:
        ordered = sorted(by_key[key], key=lambda r: (r.start, r.end))
        cluster: list[PlacedModel] = []
        cluster_end = 0
        for r in ordered:
            if cluster and r.start < cluster_end:
                cluster.append(r)
                cluster_end = max(cluster_end, r.end)
            else:
                if cluster:
                    clusters.append(cluster)
                cluster = [r]
                cluster_end = r.end
        if cluster:
            clusters.append(cluster)
    return clusters


def _clusters_via_bedtools(records: list[PlacedModel], binary: str
                            ) -> list[list[PlacedModel]]:
    """bedtools clustering: emit each model as a 6-column BED row keyed by
    its index (column 4), sorted by ``(chrom, start, end)``, and run
    ``bedtools merge -s -d -1 -c 4 -o collapse``. ``-s`` merges within a
    strand; ``-d -1`` requires >=1 bp of overlap so book-ended features
    are NOT merged (bedtools' default ``-d 0`` would merge them, breaking
    the gap guard); ``-c 4 -o collapse`` lists each merged interval's
    member indices, which map back to the original records."""
    indexed = sorted(enumerate(records),
                     key=lambda t: (t[1].chrom, t[1].start, t[1].end))
    bed = "".join(f"{r.chrom}\t{r.start}\t{r.end}\t{i}\t.\t{r.strand}\n"
                  for i, r in indexed)
    proc = subprocess.run(
        [binary, "merge", "-s", "-d", "-1", "-c", "4", "-o", "collapse", "-i", "-"],
        input=bed, capture_output=True, text=True, check=True)
    return _clusters_from_merge_output(proc.stdout, records)


def _clusters_from_merge_output(output: str, records: list[PlacedModel]
                                 ) -> list[list[PlacedModel]]:
    """Map ``bedtools merge -c 4 -o collapse`` output back to clusters.
    The collapsed index list is the last tab field of each line (``-s``
    may prepend a strand column, so read from the right); each non-blank
    line becomes one cluster of the referenced records."""
    clusters: list[list[PlacedModel]] = []
    for line in output.splitlines():
        if not line.strip():
            continue
        idx_field = line.split("\t")[-1]
        clusters.append([records[int(x)] for x in idx_field.split(",")])
    return clusters


# ---- representative + provenance + QC flags (Task 3) -------------------------

# Placeholder QC thresholds. These are NOT final scientific values — they
# are CALIBRATED IN TASK 8 against the positive-control set and BUSCO
# single-copy genes. They are surfaced as qc_flags() parameters so Task 8
# can retune them without touching any call site; the defaults below only
# make the function usable in isolation.
LOW_MARGIN_THRESHOLD_DEFAULT = 2.0        # best-vs-second %identity margin (percentage points)
SOURCE_DISAGREEMENT_FRAC_DEFAULT = 0.20   # >20% length divergence between source reps (design §5)


def _representative(members) -> PlacedModel:
    """Pick the representative model from a non-empty group of members.

    Preference order (design §5): a complete ORF before a partial one
    (so a complete SHORT model beats a partial LONG one); then greatest
    protein ``length``; then the genome-anchored (RefSeq) model over a
    transcriptome one; then the smallest ``query`` id for full
    determinism. Expressed as the ``min`` under a composite key whose
    every component sorts the preferred model first.
    """
    return min(members, key=lambda m: (not m.complete, -m.length,
                                       m.source != "genome", m.query))


def pick_representative(locus: Locus) -> PlacedModel:
    """Representative ``PlacedModel`` for a locus (= gene). See
    ``_representative`` for the preference order."""
    return _representative(locus.members)


def provenance(locus: Locus) -> str:
    """Label a locus ``both`` / ``genome_only`` / ``transcriptome_only``
    from the source tracks present among its members: ``both`` iff a
    genome (RefSeq) model AND a transcriptome model share the locus; an
    all-txome group (including the unplaced bucket) is
    ``transcriptome_only``."""
    sources = {m.source for m in locus.members}
    if "genome" in sources and "txome" in sources:
        return "both"
    if "genome" in sources:
        return "genome_only"
    return "transcriptome_only"


def _gene_interval(gene) -> tuple[str, int, int, str]:
    """Normalize a RefSeq gene interval to ``(chrom, start, end, strand)``,
    accepting either a 4-tuple or any object exposing
    ``.chrom/.start/.end/.strand`` (e.g. a ``Locus`` or ``PlacedModel``)."""
    if isinstance(gene, tuple):
        return gene
    return (gene.chrom, gene.start, gene.end, gene.strand)


def _member_overlaps_gene(m: PlacedModel, gene: tuple[str, int, int, str]) -> bool:
    """Same-chromosome, same-strand, half-open ``[start, end)`` overlap
    (the Task-2 overlap contract) between a member and a gene interval."""
    gchrom, gstart, gend, gstrand = gene
    return (m.chrom == gchrom and m.strand == gstrand
            and m.start < gend and gstart < m.end)


def _is_chimeric(locus: Locus, refseq_genes) -> bool:
    """True iff any single member overlaps >=2 DISTINCT RefSeq genes — the
    fused-model failure mode. Gene intervals are de-duplicated (a set of
    normalized tuples) so an identical interval listed twice can't fake a
    chimera."""
    genes = {_gene_interval(g) for g in refseq_genes}
    for m in locus.members:
        hits = sum(1 for gene in genes if _member_overlaps_gene(m, gene))
        if hits >= 2:
            return True
    return False


def _source_disagreement(locus: Locus, frac: float) -> bool:
    """True iff the genome-side and txome-side representatives differ in
    protein ``length`` by more than ``frac`` (relative to the longer of
    the two). The caller only invokes this on ``both`` provenance, so each
    side is non-empty; a defensive guard still returns False if either
    side is missing or neither carries a length."""
    genome = [m for m in locus.members if m.source == "genome"]
    txome = [m for m in locus.members if m.source == "txome"]
    if not genome or not txome:
        return False
    g_rep = _representative(genome)
    t_rep = _representative(txome)
    longer = max(g_rep.length, t_rep.length)
    if longer <= 0:
        return False    # no length information -> cannot judge disagreement
    return abs(g_rep.length - t_rep.length) / longer > frac


def qc_flags(locus: Locus, *, refseq_genes=None, margin: float | None = None,
             low_margin_threshold: float = LOW_MARGIN_THRESHOLD_DEFAULT,
             multi_mapping: bool = False,
             source_disagreement_frac: float = SOURCE_DISAGREEMENT_FRAC_DEFAULT
             ) -> tuple[str, ...]:
    """Deterministic (sorted) tuple of QC flags for a locus.

    Each flag isolates one failure mode; the caller supplies exactly the
    context each needs (design §8):

    * ``chimeric`` — a member overlaps >=2 distinct ``refseq_genes``
      (a list of ``Locus`` / ``(chrom, start, end, strand)`` intervals);
      skipped entirely when ``refseq_genes`` is None.
    * ``partial_only`` — no member has a complete ORF.
    * ``source_disagreement`` — a ``both`` locus whose genome-side and
      txome-side representatives differ in length by more than
      ``source_disagreement_frac``.
    * ``low_margin`` — the best-vs-second placement %identity ``margin``
      is below ``low_margin_threshold`` (ignored when ``margin`` is None,
      i.e. a single unambiguous placement).
    * ``multi_mapping`` — caller-supplied: the transcript had >=2
      gate-passing loci (ambiguous placement; only Task 4 can know this,
      so it is passed in rather than re-derived here).

    ``low_margin_threshold`` (and the margin scale generally) is a
    Task-8-calibrated placeholder — see the module-level constants.
    """
    flags: set[str] = set()
    if refseq_genes is not None and _is_chimeric(locus, refseq_genes):
        flags.add("chimeric")
    if not any(m.complete for m in locus.members):
        flags.add("partial_only")
    if provenance(locus) == "both" and _source_disagreement(locus, source_disagreement_frac):
        flags.add("source_disagreement")
    if margin is not None and margin < low_margin_threshold:
        flags.add("low_margin")
    if multi_mapping:
        flags.add("multi_mapping")
    return tuple(sorted(flags))
