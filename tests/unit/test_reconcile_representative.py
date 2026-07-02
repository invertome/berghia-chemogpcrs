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
    small = _pm("s", complete=True, length=200)
    big = _pm("b", complete=True, length=400)
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
