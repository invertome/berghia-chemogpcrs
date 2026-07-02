"""Tests for scripts/reconcile_candidates.py — Task 3.

Genome-track + candidate reconciliation feature (bead: genome-track
stage 02c). Task 3 covers only representative selection, provenance
labeling, and QC flags over an already-grouped ``Locus``. The placement
cascade, ``reconcile()``, and the output writers belong to later tasks
and are out of scope here.

Rules under test (design doc §5, §8):
  * representative = complete-ORF first, then longest, tie -> the
    genome-anchored (RefSeq) model, final tie -> smallest query id
    (deterministic);
  * provenance = both / genome_only / transcriptome_only from the source
    labels present among a locus's members;
  * QC flags (sorted, deterministic) drawn from chimeric / partial_only /
    source_disagreement / low_margin / multi_mapping.

Loci are built directly here (not via ``group_into_loci``) so each test
controls its members precisely. The numeric QC thresholds are
placeholders calibrated in Task 8, so the threshold-sensitive tests pass
explicit values rather than leaning on the module defaults.
"""
from __future__ import annotations

import reconcile_candidates as rc


# ---- builders -------------------------------------------------------------

def _pm(query, *, source="txome", complete=False, length=0,
        chrom="chr1", start=100, end=500, strand="+"):
    return rc.PlacedModel(query, chrom, start, end, strand,
                          source=source, complete=complete, length=length)


def _locus(members, *, chrom="chr1", strand="+"):
    starts = [m.start for m in members]
    ends = [m.end for m in members]
    return rc.Locus(chrom, min(starts), max(ends), strand, tuple(members))


# ---- PlacedModel: complete + length fields --------------------------------

def test_placed_model_complete_and_length_default_false_and_zero():
    # Task-2 call sites omit the new fields; they must default so every
    # existing construction (and Task 2's tests) stays valid.
    pm = rc.PlacedModel("t1", "chr1", 100, 500, "+", source="txome")
    assert pm.complete is False
    assert pm.length == 0


def test_placed_model_accepts_complete_and_length():
    pm = rc.PlacedModel("g1", "chr1", 100, 500, "+", source="genome",
                        complete=True, length=312)
    assert pm.complete is True
    assert pm.length == 312


# ---- pick_representative --------------------------------------------------

def test_pick_representative_complete_short_beats_partial_long():
    # a COMPLETE but short ORF beats a PARTIAL but longer one:
    # completeness dominates length in the preference order.
    short_complete = _pm("c", complete=True, length=100)
    long_partial = _pm("p", complete=False, length=900)
    locus = _locus([short_complete, long_partial])
    assert rc.pick_representative(locus) is short_complete


def test_pick_representative_longest_among_completes():
    # ONLY length may select the winner: the longer model has the
    # lexically LARGER query, so a query-only tiebreak would (wrongly)
    # pick the short one. Kills the "remove -m.length" mutant, under which
    # both models are complete+txome and the min falls to query "a".
    small = _pm("a", complete=True, length=200)
    big = _pm("z", complete=True, length=400)
    locus = _locus([small, big])
    assert rc.pick_representative(locus) is big


def test_pick_representative_genome_wins_tie_over_txome():
    # equal completeness + length -> the genome-anchored (RefSeq) model
    # wins. The genome id is lexically LARGER, so a query-only tiebreak
    # would (wrongly) pick the txome one: this proves source outranks
    # query.
    genome = _pm("z_genome", source="genome", complete=True, length=300)
    txome = _pm("a_txome", source="txome", complete=True, length=300)
    locus = _locus([genome, txome])
    assert rc.pick_representative(locus) is genome


def test_pick_representative_final_tiebreak_smallest_query():
    # identical completeness + length + source -> deterministic by the
    # smallest query id, independent of member input order.
    a = _pm("aaa", source="genome", complete=True, length=300)
    b = _pm("bbb", source="genome", complete=True, length=300)
    locus = _locus([b, a])   # deliberately out of order
    assert rc.pick_representative(locus) is a


# ---- provenance -----------------------------------------------------------

def test_provenance_both_when_genome_and_txome_present():
    locus = _locus([_pm("g", source="genome"), _pm("t", source="txome")])
    assert rc.provenance(locus) == "both"


def test_provenance_genome_only():
    locus = _locus([_pm("g1", source="genome"), _pm("g2", source="genome")])
    assert rc.provenance(locus) == "genome_only"


def test_provenance_transcriptome_only():
    # incl. the unplaced case: an all-txome group.
    locus = _locus([_pm("t1", source="txome"), _pm("t2", source="txome")])
    assert rc.provenance(locus) == "transcriptome_only"


# ---- qc_flags: chimeric ---------------------------------------------------

def test_qc_chimeric_member_spanning_two_refseq_genes():
    # one fused txome model [100,900) straddles two adjacent same-strand
    # RefSeq genes -> the chimeric failure mode.
    fused = _pm("fused", start=100, end=900, complete=True, length=400)
    locus = _locus([fused])
    refseq = [("chr1", 100, 400, "+"), ("chr1", 500, 800, "+")]
    assert "chimeric" in rc.qc_flags(locus, refseq_genes=refseq)


def test_qc_not_chimeric_when_member_hits_one_gene():
    clean = _pm("clean", start=100, end=400, complete=True, length=300)
    locus = _locus([clean])
    refseq = [("chr1", 100, 400, "+"), ("chr1", 500, 800, "+")]
    assert "chimeric" not in rc.qc_flags(locus, refseq_genes=refseq)


def test_qc_chimeric_ignores_opposite_strand_gene():
    # member overlaps one same-strand gene + one opposite-strand gene ->
    # only the same-strand one counts, so NOT chimeric (strand-aware, per
    # the Task-2 overlap contract).
    m = _pm("m", start=100, end=900, strand="+", complete=True, length=400)
    locus = _locus([m], strand="+")
    refseq = [("chr1", 100, 400, "+"), ("chr1", 500, 800, "-")]
    assert "chimeric" not in rc.qc_flags(locus, refseq_genes=refseq)


def test_qc_chimeric_accepts_locus_objects_as_genes():
    # refseq_genes may be Locus/PlacedModel-like objects, not just
    # 4-tuples: a fused member straddling two same-strand gene loci is
    # still chimeric.
    fused = _pm("fused", start=100, end=900, complete=True, length=400)
    locus = _locus([fused])
    genes = [rc.Locus("chr1", 100, 400, "+", ()),
             rc.Locus("chr1", 500, 800, "+", ())]
    assert "chimeric" in rc.qc_flags(locus, refseq_genes=genes)


def test_qc_chimeric_not_flagged_when_refseq_genes_none():
    fused = _pm("fused", start=100, end=900, complete=True, length=400)
    locus = _locus([fused])
    assert "chimeric" not in rc.qc_flags(locus, refseq_genes=None)


def test_qc_chimeric_via_two_distinct_genes_with_a_duplicate_listed():
    # refseq_genes lists gene A twice plus a genuinely distinct gene B; the
    # fused member straddles A and B, so it IS chimeric via the 2 DISTINCT
    # genes (the duplicate is irrelevant).
    fused = _pm("fused", start=100, end=900, complete=True, length=400)
    locus = _locus([fused])
    a = ("chr1", 100, 400, "+")
    b = ("chr1", 500, 800, "+")
    assert "chimeric" in rc.qc_flags(locus, refseq_genes=[a, a, b])


def test_qc_not_chimeric_when_same_gene_interval_listed_twice():
    # the ONLY overlapped gene is listed twice; de-dup collapses it to one,
    # so the member overlaps just 1 DISTINCT gene -> NOT chimeric. Kills the
    # "drop the set() de-dup" mutant, under which the duplicate is counted
    # twice and fakes a 2-gene chimera.
    m = _pm("m", start=100, end=400, complete=True, length=300)
    locus = _locus([m])
    a = ("chr1", 100, 400, "+")
    assert "chimeric" not in rc.qc_flags(locus, refseq_genes=[a, a])


def test_qc_chimeric_respects_chromosome_of_gene():
    # member on chr1 overlaps a chr1 gene (1 distinct) and, by coordinate
    # only, a chr2 gene. The chromosome check must reject the chr2 gene, so
    # the member overlaps 1 distinct gene -> NOT chimeric. Kills the
    # "drop m.chrom == gchrom" mutant, which would count the chr2 gene and
    # flag a spurious chimera.
    m = _pm("m", chrom="chr1", start=100, end=900, complete=True, length=400)
    locus = _locus([m], chrom="chr1")
    same_chrom = ("chr1", 100, 400, "+")
    other_chrom = ("chr2", 200, 800, "+")   # overlaps m's coords, wrong chrom
    assert "chimeric" not in rc.qc_flags(locus, refseq_genes=[same_chrom, other_chrom])


def test_qc_not_chimeric_bookended_gene_half_open_boundary():
    # member [100,500) overlaps gene A [100,300) (1 distinct) and is
    # book-ended with gene B [500,900): B.start == member.end, so 0 shared
    # bases under half-open [start, end). B must NOT count -> 1 distinct
    # gene -> NOT chimeric. Kills a `gstart < m.end` -> `<=` boundary
    # mutant, which would treat the book-end as an overlap (the Task-2
    # gap-guard semantics).
    m = _pm("m", start=100, end=500, complete=True, length=300)
    locus = _locus([m])
    a = ("chr1", 100, 300, "+")
    bookend = ("chr1", 500, 900, "+")   # starts exactly at m.end
    assert "chimeric" not in rc.qc_flags(locus, refseq_genes=[a, bookend])


# ---- qc_flags: partial_only -----------------------------------------------

def test_qc_partial_only_when_no_complete_member():
    locus = _locus([_pm("a", complete=False), _pm("b", complete=False)])
    assert "partial_only" in rc.qc_flags(locus)


def test_qc_not_partial_only_when_any_member_complete():
    locus = _locus([_pm("a", complete=False), _pm("b", complete=True)])
    assert "partial_only" not in rc.qc_flags(locus)


# ---- qc_flags: source_disagreement ----------------------------------------

def test_qc_source_disagreement_above_20pct():
    # both-provenance locus; genome rep 100 aa vs txome rep 79 aa ->
    # 21/100 = 0.21 > 0.20 -> flagged.
    g = _pm("g", source="genome", complete=True, length=100)
    t = _pm("t", source="txome", complete=True, length=79)
    locus = _locus([g, t])
    assert "source_disagreement" in rc.qc_flags(locus)


def test_qc_no_source_disagreement_within_20pct():
    # 20/100 = 0.20 is NOT strictly > 0.20 -> not flagged (boundary).
    g = _pm("g", source="genome", complete=True, length=100)
    t = _pm("t", source="txome", complete=True, length=80)
    locus = _locus([g, t])
    assert "source_disagreement" not in rc.qc_flags(locus)


def test_qc_source_disagreement_requires_both_provenance():
    # a genome_only locus with wildly different member lengths is NOT a
    # source disagreement (there is no txome side to disagree with).
    g1 = _pm("g1", source="genome", complete=True, length=100)
    g2 = _pm("g2", source="genome", complete=True, length=10)
    locus = _locus([g1, g2])
    assert "source_disagreement" not in rc.qc_flags(locus)


def test_qc_source_disagreement_uses_per_source_representatives():
    # >=2 members per source, and the true representative is NOT the first
    # member on either side (each side's partial is longer but loses to the
    # complete ORF). The disagreement must come from the REAL per-source
    # reps -- genome rep 850 vs txome rep 100 -> 750/850 = 0.88 > 0.20 ->
    # flagged. Kills BOTH the `_representative(genome)->genome[0]` mutant
    # (genome[0]=100 vs txome rep 100 -> 0.0, no fire) AND the
    # `_representative(txome)->txome[0]` mutant (genome rep 850 vs
    # txome[0]=900 -> 0.056, no fire).
    genome = [_pm("g1", source="genome", complete=False, length=100),
              _pm("g2", source="genome", complete=True, length=850)]
    txome = [_pm("t1", source="txome", complete=False, length=900),
             _pm("t2", source="txome", complete=True, length=100)]
    locus = _locus([*genome, *txome])
    assert "source_disagreement" in rc.qc_flags(locus)


def test_qc_source_disagreement_flags_when_txome_rep_is_longer():
    # txome rep (850) is LONGER than genome rep (100): dropping abs() makes
    # (100 - 850)/850 negative, which is not > 0.20, so the flag would not
    # fire. Pins abs() -- the existing >20% test has the genome side longer
    # (a positive difference), so it survives an abs() deletion.
    g = _pm("g", source="genome", complete=True, length=100)
    t = _pm("t", source="txome", complete=True, length=850)
    locus = _locus([g, t])
    assert "source_disagreement" in rc.qc_flags(locus)


def test_qc_source_disagreement_zero_length_reps_no_error():
    # both-provenance locus whose two reps both have length 0 (unknown /
    # missing length): the `longer <= 0` guard must short-circuit, giving
    # no source_disagreement and, crucially, no ZeroDivisionError. Removing
    # the guard makes qc_flags raise on the 0/0 division.
    g = _pm("g", source="genome", complete=True, length=0)
    t = _pm("t", source="txome", complete=True, length=0)
    locus = _locus([g, t])
    flags = rc.qc_flags(locus)
    assert "source_disagreement" not in flags


# ---- qc_flags: low_margin -------------------------------------------------

def test_qc_low_margin_below_threshold():
    locus = _locus([_pm("a", complete=True, length=300)])
    flags = rc.qc_flags(locus, margin=1.0, low_margin_threshold=2.0)
    assert "low_margin" in flags


def test_qc_no_low_margin_at_or_above_threshold():
    locus = _locus([_pm("a", complete=True, length=300)])
    assert "low_margin" not in rc.qc_flags(locus, margin=2.0, low_margin_threshold=2.0)
    assert "low_margin" not in rc.qc_flags(locus, margin=5.0, low_margin_threshold=2.0)


def test_qc_low_margin_none_margin_not_flagged():
    # margin is None for a single (unambiguous) placement: the absence of
    # a runner-up is NOT a low-margin failure.
    locus = _locus([_pm("a", complete=True, length=300)])
    assert "low_margin" not in rc.qc_flags(locus, margin=None, low_margin_threshold=2.0)


# ---- qc_flags: multi_mapping ----------------------------------------------

def test_qc_multi_mapping_passthrough_true():
    locus = _locus([_pm("a", complete=True, length=300)])
    assert "multi_mapping" in rc.qc_flags(locus, multi_mapping=True)


def test_qc_multi_mapping_false_not_flagged():
    locus = _locus([_pm("a", complete=True, length=300)])
    assert "multi_mapping" not in rc.qc_flags(locus, multi_mapping=False)


# ---- qc_flags: determinism ------------------------------------------------

def test_qc_flags_sorted_and_deterministic():
    # a partial fused txome model, ambiguous placement, thin margin ->
    # several flags at once; output is a sorted tuple.
    fused = _pm("fused", start=100, end=900, complete=False, length=400)
    locus = _locus([fused])
    refseq = [("chr1", 100, 400, "+"), ("chr1", 500, 800, "+")]
    flags = rc.qc_flags(locus, refseq_genes=refseq, margin=0.5,
                        low_margin_threshold=2.0, multi_mapping=True)
    assert flags == ("chimeric", "low_margin", "multi_mapping", "partial_only")
    assert list(flags) == sorted(flags)


def test_qc_flags_empty_when_clean():
    # complete both-locus, matched lengths, wide margin, unambiguous,
    # single-gene members -> no flags.
    g = _pm("g", source="genome", start=100, end=400, complete=True, length=300)
    t = _pm("t", source="txome", start=100, end=400, complete=True, length=300)
    locus = _locus([g, t])
    refseq = [("chr1", 100, 400, "+")]
    flags = rc.qc_flags(locus, refseq_genes=refseq, margin=50.0,
                        low_margin_threshold=2.0, multi_mapping=False)
    assert flags == ()


def test_qc_flags_returns_tuple():
    locus = _locus([_pm("a", complete=True, length=300)])
    assert isinstance(rc.qc_flags(locus), tuple)
