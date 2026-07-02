"""Positive-control scientific validation for scripts/reconcile_candidates.py.

Genome-track + candidate reconciliation feature (bead: genome-track stage
02c), Task 8. These are the CORE scientific-correctness guards for the
reconciliation design: they pin the two claims the whole feature rests on,
using synthetic fixtures fed to the public ``reconcile()`` entry point.

Design claims under test (design doc §3-§5, §9):

  1. TANDEM PARALOG ARRAY — the signature chemoreceptor LSE signal. Several
     paralogous genes ~90% identical to EACH OTHER sit at DISTINCT genomic
     loci. Each transcript aligns near-perfectly to its OWN locus but only
     ~90% to its siblings. The STRINGENT %match gate (design §4) drops the
     ~90% cross-maps, so every transcript is placed at its own locus and the
     distinct loci stay SEPARATE genes. Collapsing them would erase exactly
     the expansion this pipeline exists to find.

  2. ISOFORM / 1:1 COLLAPSE — a housekeeping GPCR present in BOTH tracks
     (a RefSeq model + one-or-more transcriptome isoforms) that map to ONE
     locus collapses to a SINGLE ``both`` gene, with every source isoform id
     retained (locus == gene; isoforms collapse, provenance is preserved).

  3. GATE BOUNDARY — a placement materially below the identity bar
     (94% under min_id=95) is NOT the same gene: it is rejected, the RefSeq
     gene stands alone, and the transcript falls to ``unplaced`` rather than
     silently merging. The contrast case (96% >= 95) DOES merge — proving the
     bar, not an accident of the fixture, drives the separation.

Unlike a red-green unit test of new code, these VALIDATE the already-
implemented ``reconcile`` against the design's scientific contract: they
pass on the correct implementation and would fail if a future change
collapsed paralogs, split isoforms, or loosened the gate. Threshold-
sensitive cases pass explicit ``Thresholds`` (the numeric defaults are
themselves Task-8-calibrated placeholders).
"""
from __future__ import annotations

import reconcile_candidates as rc


# ---- builders -------------------------------------------------------------

def _pl(query, chrom, start, end, strand="+", pid=98.0, cov=95.0,
        method="minimap2", score=60.0):
    """A single aligner placement (matches the cascade tests' helper)."""
    return rc.Placement(query, chrom, start, end, strand, pid, cov, method, score)


def _concordant_bundle(query, chrom, start, end, strand="+", pid=98.0, cov=95.0):
    """A minimap2+gmap pair that AGREES at one locus (both gate-passing), the
    high-confidence cascade route. gmap is offset slightly and a touch lower
    identity, mirroring two independent spliced aligners on the same gene."""
    return {
        "minimap2": [_pl(query, chrom, start, end, strand, pid, cov, "minimap2")],
        "gmap": [_pl(query, chrom, start + 10, end - 10, strand,
                     pid - 1.0, cov - 1.0, "gmap")],
    }


# ============================================================================
# 1. TANDEM PARALOG ARRAY — paralogs at distinct loci stay separate genes.
# ============================================================================

def test_tandem_paralog_array_yields_three_separate_genes():
    """Three ~90%-identical paralogs in a tandem array (three distinct loci on
    one scaffold). Each transcript maps near-perfectly (98%) to its OWN locus
    AND cross-maps at 90% (concordantly!) to its two siblings. The stringent
    gate (min_id=95) drops the 90% cross-maps, so each transcript is placed at
    its own locus and the three loci remain THREE genes — never collapsed."""
    thr = rc.Thresholds(min_id=95.0, min_cov=90.0, min_margin=2.0)
    loci = {
        "A": ("chr1", 1000, 2000),
        "B": ("chr1", 3000, 4000),   # non-overlapping with A (2000 < 3000)
        "C": ("chr1", 5000, 6000),   # non-overlapping with B (4000 < 5000)
    }
    txome = [rc.Candidate("cand_A", complete=True, length=330, n_tm=7),
             rc.Candidate("cand_B", complete=True, length=328, n_tm=7),
             rc.Candidate("cand_C", complete=True, length=331, n_tm=7)]

    placements = {}
    for tag, cand in (("A", "cand_A"), ("B", "cand_B"), ("C", "cand_C")):
        own_chrom, own_s, own_e = loci[tag]
        mm = [_pl(cand, own_chrom, own_s, own_e, "+", 98.0, 95.0, "minimap2")]
        gm = [_pl(cand, own_chrom, own_s + 10, own_e - 10, "+", 97.0, 94.0, "gmap")]
        # Sub-gate (90%) concordant cross-maps to the OTHER two paralog loci —
        # these must be filtered out by the gate, not collapse the array.
        for other in loci:
            if other == tag:
                continue
            oc, os_, oe = loci[other]
            mm.append(_pl(cand, oc, os_, oe, "+", 90.0, 95.0, "minimap2"))
            gm.append(_pl(cand, oc, os_ + 10, oe - 10, "+", 90.0, 94.0, "gmap"))
        placements[cand] = {"minimap2": mm, "gmap": gm}

    genes = rc.reconcile(txome, [], placements, [], thr)

    # The critical assertion: 3 paralogs -> 3 genes (NOT collapsed to 1).
    assert len(genes) == 3
    assert {g.gene_id for g in genes} == {
        "chr1:1000-2000:+", "chr1:3000-4000:+", "chr1:5000-6000:+"}
    # Each is its own transcriptome_only gene at a distinct locus with its own
    # representative — no two paralogs folded into one gene.
    assert {g.representative_id for g in genes} == {"cand_A", "cand_B", "cand_C"}
    assert all(g.provenance == "transcriptome_only" for g in genes)
    assert all(g.placement_method == "minimap2+gmap_concordant" for g in genes)
    # Each gene carries exactly ONE transcriptome isoform (no sibling merged in).
    assert all(len(g.transcriptome_isoform_ids) == 1 for g in genes)


def test_tandem_paralogs_are_not_collapsed_by_locus_grouping():
    """Lower-level guard on the same claim: three same-strand, non-overlapping
    placed paralog models group into THREE loci, never one. This isolates the
    ``group_into_loci`` half of the separation (independent of the cascade)."""
    models = [
        rc.PlacedModel("cand_A", "chr1", 1000, 2000, "+", source="txome"),
        rc.PlacedModel("cand_B", "chr1", 3000, 4000, "+", source="txome"),
        rc.PlacedModel("cand_C", "chr1", 5000, 6000, "+", source="txome"),
    ]
    assert len(rc.group_into_loci(models)) == 3


# ============================================================================
# 2. ISOFORM / 1:1 COLLAPSE — one gene present in both tracks -> one `both`.
# ============================================================================

def test_housekeeping_gpcr_txome_and_refseq_collapse_to_one_both_gene():
    """A housekeeping GPCR present in BOTH tracks: one RefSeq model and one
    transcriptome candidate that place at the SAME locus (isoforms of one
    gene) collapse to a SINGLE ``both`` gene, retaining both source ids."""
    thr = rc.Thresholds(min_id=95.0, min_cov=90.0, min_margin=2.0)
    refseq = [rc.PlacedModel("XP_hk", "chr4", 10000, 12000, "+",
                             source="genome", complete=True, length=360, n_tm=7)]
    refseq_loci = [("chr4", 10000, 12000, "+")]
    txome = [rc.Candidate("Bst_hk", complete=True, length=352, n_tm=7)]
    placements = {"Bst_hk": _concordant_bundle("Bst_hk", "chr4", 10000, 12000)}

    genes = rc.reconcile(txome, refseq, placements, refseq_loci, thr)

    assert len(genes) == 1
    g = genes[0]
    assert g.provenance == "both"
    assert g.genome_locus == "chr4:10000-12000:+"
    # Both source ids are retained (never dropped in the collapse).
    assert g.genome_isoform_ids == ("XP_hk",)
    assert g.transcriptome_isoform_ids == ("Bst_hk",)
    # RefSeq (360 aa) outlengths the txome model (352 aa) -> genome rep.
    assert g.representative_source == "genome"


def test_multiple_txome_isoforms_collapse_at_one_locus():
    """Two transcriptome isoforms of one gene, plus its RefSeq model, all at a
    single locus -> one ``both`` gene with BOTH txome isoform ids merged (locus
    == gene: isoforms collapse rather than proliferate genes)."""
    thr = rc.Thresholds(min_id=95.0, min_cov=90.0, min_margin=2.0)
    refseq = [rc.PlacedModel("XP_iso", "chr2", 500, 2500, "+",
                             source="genome", complete=True, length=400, n_tm=7)]
    refseq_loci = [("chr2", 500, 2500, "+")]
    txome = [rc.Candidate("Bst_iso1", complete=True, length=390, n_tm=7),
             rc.Candidate("Bst_iso2", complete=False, length=210, n_tm=5)]
    placements = {
        "Bst_iso1": _concordant_bundle("Bst_iso1", "chr2", 500, 2500),
        # Overlapping the same locus (500 < 2000 and 600 < 2500) -> same gene.
        "Bst_iso2": _concordant_bundle("Bst_iso2", "chr2", 600, 2000),
    }

    genes = rc.reconcile(txome, refseq, placements, refseq_loci, thr)

    assert len(genes) == 1
    g = genes[0]
    assert g.provenance == "both"
    assert g.transcriptome_isoform_ids == ("Bst_iso1", "Bst_iso2")
    assert g.genome_isoform_ids == ("XP_iso",)


# ============================================================================
# 3. GATE BOUNDARY — the stringent %id bar separates a near-miss as a
#    DIFFERENT gene, and the contrast case merges (the bar drives it).
# ============================================================================

def test_subthreshold_identity_txome_not_merged_into_refseq_gene():
    """A transcript placing at a RefSeq gene's locus at 94% identity, under
    min_id=95, is NOT the same gene: it fails the gate, the RefSeq gene stands
    alone (``genome_only``) and the transcript is retained as ``unplaced`` —
    two genes, NEVER a spurious ``both`` merge."""
    thr = rc.Thresholds(min_id=95.0, min_cov=90.0, min_margin=2.0)
    refseq = [rc.PlacedModel("XP_g", "chr1", 1000, 2000, "+",
                             source="genome", complete=True, length=350, n_tm=7)]
    refseq_loci = [("chr1", 1000, 2000, "+")]
    txome = [rc.Candidate("Bst_q", complete=True, length=340, n_tm=7)]
    # 94% id concordant at the RefSeq locus — one point below the bar.
    placements = {"Bst_q": _concordant_bundle("Bst_q", "chr1", 1000, 2000, pid=94.0)}

    genes = rc.reconcile(txome, refseq, placements, refseq_loci, thr)

    assert len(genes) == 2
    by_id = {g.gene_id: g for g in genes}
    assert by_id["chr1:1000-2000:+"].provenance == "genome_only"
    assert by_id["chr1:1000-2000:+"].transcriptome_isoform_ids == ()
    assert by_id["unplaced:Bst_q"].provenance == "transcriptome_only"
    assert by_id["unplaced:Bst_q"].placement_method == "unplaced"
    assert "unplaced" in by_id["unplaced:Bst_q"].qc_flags
    # No gene is `both`: the near-miss did not merge.
    assert all(g.provenance != "both" for g in genes)


def test_above_gate_identity_txome_merges_into_both():
    """Contrast to the boundary case: the SAME locus and candidate, but at 96%
    identity (>= min_id=95), DOES merge into one ``both`` gene. Together with
    the 94% case this proves the identity BAR drives the separation, not the
    fixture geometry."""
    thr = rc.Thresholds(min_id=95.0, min_cov=90.0, min_margin=2.0)
    refseq = [rc.PlacedModel("XP_g", "chr1", 1000, 2000, "+",
                             source="genome", complete=True, length=350, n_tm=7)]
    refseq_loci = [("chr1", 1000, 2000, "+")]
    txome = [rc.Candidate("Bst_q", complete=True, length=340, n_tm=7)]
    placements = {"Bst_q": _concordant_bundle("Bst_q", "chr1", 1000, 2000, pid=96.0)}

    genes = rc.reconcile(txome, refseq, placements, refseq_loci, thr)

    assert len(genes) == 1
    g = genes[0]
    assert g.provenance == "both"
    assert g.genome_isoform_ids == ("XP_g",)
    assert g.transcriptome_isoform_ids == ("Bst_q",)
