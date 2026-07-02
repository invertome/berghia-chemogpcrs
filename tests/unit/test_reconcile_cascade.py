"""Tests for scripts/reconcile_candidates.py — Task 4 (placement cascade
+ reconcile entry point).

Genome-track + candidate reconciliation feature (bead: genome-track stage
02c). Task 4 wires Tasks 1-3 together: the dual-aligner concordance check,
the placement cascade (concordant minimap2+gmap -> miniprot -> RBH ->
unplaced), and the top-level ``reconcile()`` that emits one
``ReconciledGene`` per gene.

Design (design doc §3-§5, §8):
  * cascade STOPS at the first confident hit and NEVER drops a transcript
    (unplaced is a labeled outcome, not a deletion);
  * concordance REQUIRES two independent spliced aligners agreeing — a
    lone minimap2 hit (no gmap agreement) is NOT confident and falls
    through to miniprot;
  * the %match gate is absolute pct_identity AND pct_coverage; the
    best-vs-second margin is surfaced (a QC signal), not a hard filter, so
    near-identical LSE paralogs are placed-and-flagged rather than dropped.

The numeric thresholds are Task-8 placeholders, so threshold-sensitive
tests pass explicit ``Thresholds`` rather than leaning on defaults.
"""
from __future__ import annotations

from dataclasses import fields

import reconcile_candidates as rc


# ---- builders -------------------------------------------------------------

def _pl(query="q", chrom="chr1", start=1000, end=2000, strand="+",
        pid=98.0, cov=95.0, method="minimap2", score=60.0):
    return rc.Placement(query, chrom, start, end, strand, pid, cov, method, score)


THR = rc.Thresholds()   # 95 / 90 / 2 / 2 defaults


# ---- PlacedModel: new n_tm field ------------------------------------------

def test_placed_model_n_tm_defaults_zero():
    # Task 2/3 call sites omit n_tm; it must default so their constructions
    # (and tests) stay valid.
    pm = rc.PlacedModel("t1", "chr1", 100, 500, "+", source="txome")
    assert pm.n_tm == 0


def test_placed_model_accepts_n_tm():
    pm = rc.PlacedModel("g1", "chr1", 100, 500, "+", source="genome",
                        complete=True, length=312, n_tm=7)
    assert pm.n_tm == 7


# ---- Candidate dataclass --------------------------------------------------

def test_candidate_defaults():
    c = rc.Candidate("BersteEVm1")
    assert c.query == "BersteEVm1"
    assert c.complete is False
    assert c.length == 0
    assert c.n_tm == 0


def test_candidate_carries_all_axes():
    c = rc.Candidate("BersteEVm1", complete=True, length=340, n_tm=7)
    assert c.complete is True and c.length == 340 and c.n_tm == 7


# ---- Thresholds dataclass -------------------------------------------------

def test_thresholds_defaults():
    t = rc.Thresholds()
    assert t.min_id == 95.0
    assert t.min_cov == 90.0
    assert t.min_margin == 2.0
    assert t.low_margin_threshold == 2.0


# ---- concordant -----------------------------------------------------------

def test_concordant_overlapping_and_gated_returns_higher_identity():
    # both pass the gate and overlap -> the agreed placement is the
    # higher-identity of the two (here gmap at 98 beats minimap2 at 96).
    mm = _pl(chrom="chr1", start=1000, end=2000, pid=96.0, cov=95.0, method="minimap2")
    gm = _pl(chrom="chr1", start=1010, end=1990, pid=98.0, cov=94.0, method="gmap")
    agreed = rc.concordant(mm, gm, THR.min_id, THR.min_cov)
    assert agreed is not None
    assert agreed.method == "gmap" and agreed.pct_identity == 98.0


def test_concordant_returns_minimap2_when_it_is_higher():
    # symmetric: minimap2 higher -> the agreed placement is minimap2's.
    mm = _pl(pid=99.0, cov=95.0, method="minimap2")
    gm = _pl(start=1010, end=1990, pid=97.0, cov=94.0, method="gmap")
    agreed = rc.concordant(mm, gm, THR.min_id, THR.min_cov)
    assert agreed.method == "minimap2" and agreed.pct_identity == 99.0


def test_concordant_none_when_loci_disjoint():
    # both pass the gate but the two aligners place at DIFFERENT loci
    # (no overlap) -> not concordant.
    mm = _pl(chrom="chr1", start=1000, end=2000, pid=98.0, cov=95.0, method="minimap2")
    gm = _pl(chrom="chr1", start=5000, end=6000, pid=98.0, cov=95.0, method="gmap")
    assert rc.concordant(mm, gm, THR.min_id, THR.min_cov) is None


def test_concordant_none_when_opposite_strand():
    # same coordinates, opposite strand -> not the same locus (strand-aware
    # overlap, per the Task-2 contract) -> not concordant.
    mm = _pl(chrom="chr1", start=1000, end=2000, strand="+", method="minimap2")
    gm = _pl(chrom="chr1", start=1000, end=2000, strand="-", method="gmap")
    assert rc.concordant(mm, gm, THR.min_id, THR.min_cov) is None


def test_concordant_none_when_one_fails_gate():
    # they overlap, but gmap is below the identity bar -> not concordant
    # (both must independently clear the gate).
    mm = _pl(pid=98.0, cov=95.0, method="minimap2")
    gm = _pl(start=1010, end=1990, pid=80.0, cov=95.0, method="gmap")
    assert rc.concordant(mm, gm, THR.min_id, THR.min_cov) is None


def test_concordant_none_when_bookended_zero_overlap():
    # [1000,2000) and [2000,3000) touch at 2000 with 0 shared bases under
    # half-open semantics -> not the same locus -> not concordant.
    mm = _pl(chrom="chr1", start=1000, end=2000, method="minimap2")
    gm = _pl(chrom="chr1", start=2000, end=3000, method="gmap")
    assert rc.concordant(mm, gm, THR.min_id, THR.min_cov) is None


# ---- place_transcript: cascade ordering -----------------------------------

def test_minimap2_alone_without_gmap_falls_through_to_miniprot():
    # THE key concordance rule: a lone minimap2 hit (no gmap to agree with)
    # is NOT a confident placement. Even a stellar minimap2 alignment must
    # fall through to the protein (miniprot) step.
    mm = [_pl(pid=99.0, cov=99.0, method="minimap2")]
    miniprot = [_pl(chrom="chr1", start=1000, end=2000, pid=96.0, cov=92.0, method="miniprot")]
    res = rc.place_transcript("q", mm, [], miniprot, [], THR)
    assert res.method == "miniprot"
    assert res.confidence == "medium"
    assert res.placement.method == "miniprot"


def test_concordant_pair_yields_high_confidence():
    mm = [_pl(pid=98.0, cov=95.0, method="minimap2")]
    gmap = [_pl(start=1010, end=1990, pid=97.0, cov=94.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.method == "minimap2+gmap_concordant"
    assert res.confidence == "high"
    assert res.placement is not None


def test_miniprot_used_when_concordance_fails_but_gate_passes():
    # minimap2 and gmap disagree on locus (no concordance); miniprot clears
    # the gate -> medium-confidence protein placement.
    mm = [_pl(chrom="chr1", start=1000, end=2000, method="minimap2")]
    gmap = [_pl(chrom="chr1", start=5000, end=6000, method="gmap")]   # different locus
    miniprot = [_pl(chrom="chr1", start=1000, end=2000, pid=96.0, cov=91.0, method="miniprot")]
    res = rc.place_transcript("q", mm, gmap, miniprot, [], THR)
    assert res.method == "miniprot" and res.confidence == "medium"


def test_rbh_used_when_miniprot_fails_gate():
    # miniprot is present but sub-threshold; RBH clears the gate -> the
    # low-confidence homology fallback fires.
    miniprot = [_pl(pid=90.0, cov=92.0, method="miniprot")]   # below id bar
    rbh = [_pl(chrom="chr1", start=1000, end=2000, pid=96.0, cov=93.0, method="blastp")]
    res = rc.place_transcript("q", [], [], miniprot, rbh, THR)
    assert res.method == "rbh" and res.confidence == "low"
    assert res.placement.method == "blastp"


def test_all_methods_fail_yields_unplaced():
    res = rc.place_transcript("q", [], [], [], [], THR)
    assert res.placement is None
    assert res.method == "unplaced" and res.confidence == "none"
    assert res.margin is None and res.multi_mapping is False


def test_placements_present_but_all_below_gate_is_unplaced():
    # every hit is sub-threshold -> unplaced (never silently accepted).
    mm = [_pl(pid=80.0, cov=99.0, method="minimap2")]
    gmap = [_pl(start=1010, end=1990, pid=80.0, cov=99.0, method="gmap")]
    miniprot = [_pl(pid=94.0, cov=99.0, method="miniprot")]
    res = rc.place_transcript("q", mm, gmap, miniprot, [], THR)
    assert res.method == "unplaced"


def test_concordant_step_ignores_below_gate_pair():
    # minimap2+gmap overlap but BOTH are below the coverage bar; concordance
    # is gated, so this does NOT count as a confident placement.
    mm = [_pl(pid=99.0, cov=50.0, method="minimap2")]
    gmap = [_pl(start=1010, end=1990, pid=99.0, cov=50.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.method == "unplaced"


# ---- place_transcript: margin + multi_mapping -----------------------------

def test_single_confident_placement_has_infinite_margin_not_multi():
    mm = [_pl(pid=98.0, cov=95.0, method="minimap2")]
    gmap = [_pl(start=1010, end=1990, pid=97.0, cov=94.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.margin == float("inf")
    assert res.multi_mapping is False


def test_duplicate_concordant_agreements_stay_single_locus():
    # one minimap2 hit concordant with TWO overlapping gmap hits at the SAME
    # locus -> both pairs agree on the SAME (higher-identity minimap2)
    # placement. Distinct-locus clustering keeps this ONE locus: margin inf,
    # not multi_mapping.
    mm = [_pl(chrom="chr1", start=1000, end=2000, pid=99.0, cov=95.0, method="minimap2")]
    gmap = [_pl(chrom="chr1", start=1010, end=1990, pid=97.0, cov=94.0, method="gmap"),
            _pl(chrom="chr1", start=1020, end=1980, pid=96.0, cov=94.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.margin == float("inf")
    assert res.multi_mapping is False


# ---- place_transcript: §4 distinct-locus margin GATE ----------------------

def test_two_distinct_loci_margin_above_min_is_confident_multi_mapping():
    # two distinct concordant loci (chr1 id98, chr5 id95); the top locus
    # beats the runner-up by 3.0 >= min_margin 2.0 -> CONFIDENT at the top
    # (chr1) locus, and still multi_mapping (>=2 distinct loci). Both aligners
    # clear the id gate at BOTH loci (chr5 gmap id95, not below the bar).
    mm = [_pl(chrom="chr1", start=1000, end=2000, pid=98.0, cov=95.0, method="minimap2"),
          _pl(chrom="chr5", start=1000, end=2000, pid=95.0, cov=95.0, method="minimap2")]
    gmap = [_pl(chrom="chr1", start=1010, end=1990, pid=97.0, cov=94.0, method="gmap"),
            _pl(chrom="chr5", start=1005, end=1995, pid=95.0, cov=94.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.method == "minimap2+gmap_concordant" and res.confidence == "high"
    assert res.placement.chrom == "chr1"          # the top-identity locus
    assert res.multi_mapping is True
    assert abs(res.margin - 3.0) < 1e-9           # 98 - 95 across distinct loci


def test_two_distinct_loci_margin_below_min_is_ambiguous_falls_through():
    # two near-tie concordant loci (chr1 id98, chr5 id97); margin 1.0 <
    # min_margin 2.0 -> the concordance step is AMBIGUOUS (the §4 hard gate
    # rejects it) and the cascade falls through to miniprot. This test FAILS
    # if the margin gate is removed (concordance would win at chr1 instead).
    mm = [_pl(chrom="chr1", start=1000, end=2000, pid=98.0, cov=95.0, method="minimap2"),
          _pl(chrom="chr5", start=1000, end=2000, pid=97.0, cov=95.0, method="minimap2")]
    gmap = [_pl(chrom="chr1", start=1010, end=1990, pid=97.0, cov=94.0, method="gmap"),
            _pl(chrom="chr5", start=1005, end=1995, pid=96.0, cov=94.0, method="gmap")]
    miniprot = [_pl(chrom="chr1", start=1000, end=2000, pid=96.0, cov=91.0, method="miniprot")]
    res = rc.place_transcript("q", mm, gmap, miniprot, [], THR)
    assert res.method == "miniprot" and res.confidence == "medium"


def test_two_distinct_loci_margin_at_min_is_confident_boundary():
    # margin exactly == min_margin (2.0) clears the gate (>= is inclusive).
    mm = [_pl(chrom="chr1", start=1000, end=2000, pid=98.0, cov=95.0, method="minimap2"),
          _pl(chrom="chr5", start=1000, end=2000, pid=96.0, cov=95.0, method="minimap2")]
    gmap = [_pl(chrom="chr1", start=1010, end=1990, pid=97.0, cov=94.0, method="gmap"),
            _pl(chrom="chr5", start=1005, end=1995, pid=95.0, cov=94.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.method == "minimap2+gmap_concordant"
    assert abs(res.margin - 2.0) < 1e-9


def test_fully_ambiguous_cascade_is_unplaced_never_dropped():
    # concordance ambiguous (2 loci, margin 1.0 < 2.0) and NO protein /
    # homology fallback -> unplaced (retained, not dropped).
    mm = [_pl(chrom="chr1", start=1000, end=2000, pid=98.0, cov=95.0, method="minimap2"),
          _pl(chrom="chr5", start=1000, end=2000, pid=97.0, cov=95.0, method="minimap2")]
    gmap = [_pl(chrom="chr1", start=1010, end=1990, pid=97.0, cov=94.0, method="gmap"),
            _pl(chrom="chr5", start=1005, end=1995, pid=96.0, cov=94.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.placement is None and res.method == "unplaced"


def test_single_locus_multiple_within_locus_alignments_is_confident_inf_margin():
    # THE spec-review case: a lone minimap2 id96 hit agrees with TWO
    # higher-identity gmap hits (id97, id98) at ONE locus -> TWO agreed
    # placements at the SAME locus. Distinct-locus clustering collapses them
    # to one locus (margin inf), so the transcript is placed CONFIDENTLY at
    # the best within-locus alignment (id98), NOT spuriously unplaced by a
    # within-locus "margin" of 1.0.
    mm = [_pl(chrom="chr1", start=1000, end=2000, pid=96.0, cov=95.0, method="minimap2")]
    gmap = [_pl(chrom="chr1", start=1010, end=1990, pid=97.0, cov=94.0, method="gmap"),
            _pl(chrom="chr1", start=1020, end=1980, pid=98.0, cov=94.0, method="gmap")]
    res = rc.place_transcript("q", mm, gmap, [], [], THR)
    assert res.method == "minimap2+gmap_concordant" and res.confidence == "high"
    assert res.margin == float("inf")
    assert res.multi_mapping is False
    assert res.placement.pct_identity == 98.0     # best within the single locus


def test_miniprot_ambiguous_two_loci_falls_through_to_rbh():
    # the margin gate applies at EVERY cascade step, not just concordance:
    # miniprot lands two near-tie distinct loci (margin 1.0 < 2.0) ->
    # ambiguous -> RBH resolves it unambiguously.
    miniprot = [_pl(chrom="chr1", start=1000, end=2000, pid=98.0, cov=95.0, method="miniprot"),
                _pl(chrom="chr5", start=1000, end=2000, pid=97.0, cov=95.0, method="miniprot")]
    rbh = [_pl(chrom="chr8", start=1000, end=2000, pid=96.0, cov=93.0, method="blastp")]
    res = rc.place_transcript("q", [], [], miniprot, rbh, THR)
    assert res.method == "rbh" and res.confidence == "low"


# ---- PlacementResult shape ------------------------------------------------

def test_placement_result_is_five_field_namedtuple():
    res = rc.place_transcript("q", [], [], [], [], THR)
    assert res._fields == ("placement", "method", "confidence", "margin", "multi_mapping")


# ---- ReconciledGene schema (pins §8 for the Task-5 writers) ---------------

def test_reconciled_gene_has_exact_section8_schema():
    names = tuple(f.name for f in fields(rc.ReconciledGene))
    assert names == (
        "gene_id", "representative_id", "representative_source", "provenance",
        "genome_locus", "placement_method", "placement_confidence",
        "pct_identity", "pct_coverage", "best_vs_second_margin",
        "transcriptome_isoform_ids", "genome_isoform_ids", "n_tm",
        "completeness", "qc_flags",
    )


# ---- reconcile(): the full four-outcome fixture ---------------------------

def _run_fixture():
    """One RefSeq+txome overlap (both), one lone RefSeq gene (genome_only),
    one txome placed where no RefSeq gene sits (transcriptome_only), and one
    txome candidate whose only hit is sub-threshold (unplaced). Inputs are
    deliberately unsorted to exercise deterministic ordering."""
    g_both = rc.PlacedModel("XP_both", "chr1", 1000, 2000, "+",
                            source="genome", complete=True, length=350, n_tm=7)
    g_only = rc.PlacedModel("XP_only", "chr2", 5000, 6000, "+",
                            source="genome", complete=True, length=300, n_tm=6)
    refseq_candidates = [g_only, g_both]                       # unsorted
    refseq_loci = [("chr2", 5000, 6000, "+"), ("chr1", 1000, 2000, "+")]

    c_both = rc.Candidate("Bst_both", complete=True, length=340, n_tm=7)
    c_tx = rc.Candidate("Bst_txonly", complete=True, length=320, n_tm=8)
    c_unpl = rc.Candidate("Bst_unplaced", complete=True, length=200, n_tm=5)
    txome = [c_unpl, c_both, c_tx]                             # unsorted

    placements = {
        "Bst_both": {
            "minimap2": [_pl("Bst_both", "chr1", 1000, 2000, "+", 98.0, 95.0, "minimap2")],
            "gmap":     [_pl("Bst_both", "chr1", 1010, 1990, "+", 97.0, 94.0, "gmap")],
        },
        "Bst_txonly": {
            "minimap2": [_pl("Bst_txonly", "chr3", 8000, 9000, "+", 99.0, 96.0, "minimap2")],
            "gmap":     [_pl("Bst_txonly", "chr3", 8005, 8995, "+", 98.0, 95.0, "gmap")],
        },
        "Bst_unplaced": {
            # a real but sub-threshold minimap2 hit, no gmap: proves an
            # unplaced candidate that DID have placements is still emitted.
            "minimap2": [_pl("Bst_unplaced", "chr9", 1, 1000, "+", 80.0, 99.0, "minimap2")],
        },
    }
    return rc.reconcile(txome, refseq_candidates, placements, refseq_loci, rc.Thresholds())


def _by_id():
    return {g.gene_id: g for g in _run_fixture()}


def test_reconcile_emits_exactly_four_genes():
    assert len(_run_fixture()) == 4


def test_reconcile_both_gene():
    g = _by_id()["chr1:1000-2000:+"]
    assert g.provenance == "both"
    assert g.genome_locus == "chr1:1000-2000:+"
    # genome model (350 aa) outlengths the txome model (340 aa) -> RefSeq rep.
    assert g.representative_id == "XP_both"
    assert g.representative_source == "genome"
    # placement_* describe the txome evidence's route in (concordant).
    assert g.placement_method == "minimap2+gmap_concordant"
    assert g.placement_confidence == "high"
    assert g.pct_identity == 98.0 and g.pct_coverage == 95.0
    # single distinct locus -> distinct-locus margin is inf in the output.
    assert g.best_vs_second_margin == float("inf")
    assert g.completeness == "complete"
    assert g.n_tm == 7
    assert g.transcriptome_isoform_ids == ("Bst_both",)
    assert g.genome_isoform_ids == ("XP_both",)
    assert g.qc_flags == ()


def test_reconcile_genome_only_gene_is_refseq_native():
    g = _by_id()["chr2:5000-6000:+"]
    assert g.provenance == "genome_only"
    assert g.representative_id == "XP_only" and g.representative_source == "genome"
    assert g.placement_method == "refseq_native"
    assert g.placement_confidence == "high"
    # no txome placement -> the placement metrics are absent, not fabricated.
    assert g.pct_identity is None and g.pct_coverage is None
    assert g.best_vs_second_margin is None
    assert g.transcriptome_isoform_ids == ()
    assert g.genome_isoform_ids == ("XP_only",)
    assert g.n_tm == 6
    assert g.qc_flags == ()


def test_reconcile_transcriptome_only_placed_gene():
    g = _by_id()["chr3:8000-9000:+"]
    assert g.provenance == "transcriptome_only"
    assert g.representative_id == "Bst_txonly"
    assert g.representative_source == "txome"
    assert g.placement_method == "minimap2+gmap_concordant"
    assert g.placement_confidence == "high"
    assert g.pct_identity == 99.0 and g.pct_coverage == 96.0
    # n_tm travels from the Candidate through the built PlacedModel to the row.
    assert g.n_tm == 8
    assert g.genome_isoform_ids == ()


def test_reconcile_unplaced_never_dropped():
    g = _by_id()["unplaced:Bst_unplaced"]
    assert g.provenance == "transcriptome_only"
    assert g.genome_locus == "unplaced"
    assert g.placement_method == "unplaced"
    assert g.placement_confidence == "none"
    assert g.representative_id == "Bst_unplaced"
    assert g.representative_source == "txome"
    assert g.completeness == "complete"
    assert g.pct_identity is None and g.pct_coverage is None
    assert "unplaced" in g.qc_flags


def test_reconcile_output_is_sorted_by_gene_id_and_deterministic():
    ids_a = [g.gene_id for g in _run_fixture()]
    ids_b = [g.gene_id for g in _run_fixture()]
    assert ids_a == ids_b                          # deterministic
    assert ids_a == sorted(ids_a)                  # stable, gene_id order
    assert ids_a == ["chr1:1000-2000:+", "chr2:5000-6000:+",
                     "chr3:8000-9000:+", "unplaced:Bst_unplaced"]


def test_reconcile_candidate_with_no_placement_entry_is_unplaced():
    # a candidate entirely absent from the placements map -> unplaced,
    # still emitted (never dropped).
    c = rc.Candidate("Bst_orphan", complete=False, length=150, n_tm=4)
    genes = rc.reconcile([c], [], {}, [], rc.Thresholds())
    assert len(genes) == 1
    g = genes[0]
    assert g.gene_id == "unplaced:Bst_orphan"
    assert g.placement_method == "unplaced"
    assert g.completeness == "partial"
    assert "unplaced" in g.qc_flags


def test_reconcile_multi_locus_confident_carries_finite_distinct_locus_margin():
    # a transcript confidently placed at its top of TWO distinct loci
    # (chr1 id98 beats chr7 id95 by 3.0 >= min_margin): the output row is at
    # the top locus, carries the FINITE distinct-locus margin, and is flagged
    # multi_mapping. Pins that the §4 distinct-locus margin flows to output.
    c = rc.Candidate("Bst_multi", complete=True, length=300, n_tm=6)
    placements = {"Bst_multi": {
        "minimap2": [_pl("Bst_multi", "chr1", 1000, 2000, "+", 98.0, 95.0, "minimap2"),
                     _pl("Bst_multi", "chr7", 1000, 2000, "+", 95.0, 95.0, "minimap2")],
        "gmap":     [_pl("Bst_multi", "chr1", 1010, 1990, "+", 97.0, 94.0, "gmap"),
                     _pl("Bst_multi", "chr7", 1005, 1995, "+", 95.0, 94.0, "gmap")],
    }}
    genes = rc.reconcile([c], [], placements, [], rc.Thresholds())
    assert len(genes) == 1
    g = genes[0]
    assert g.genome_locus == "chr1:1000-2000:+"          # the top locus
    assert abs(g.best_vs_second_margin - 3.0) < 1e-9     # 98 - 95 distinct-locus
    assert "multi_mapping" in g.qc_flags
