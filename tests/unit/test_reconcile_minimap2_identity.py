"""Tests for the intron-aware minimap2 identity fix (bead kfqz, fix F1).

Under ``minimap2 -x splice`` the PAF alignment-block length ``alen`` (col 11)
spans the genomic footprint INCLUDING introns, so for a multi-exon transcript
``alen >> nmatch`` and the legacy ``100 * nmatch / alen`` identity collapses far
below the 95 gate — while GMAP reports a true per-base identity and passes. That
mismatch (not a real locus disagreement) is why minimap2+GMAP concordance
produced almost nothing and ~300 transcripts were wrongly left unplaced.

The fix computes the TRUE per-base identity from minimap2's own intron-aware
signals: the ``de:f``/``dv:f`` divergence tag (identity = 100*(1-div)) or, when
absent, the ``cs:Z`` string (matches / (matches+subs+indels), the ``~`` intron
op excluded). A record with neither source FAILS LOUD rather than silently
falling back to ``nmatch/alen`` or 0 — that silent fallback IS the bug.

Empirically verified against minimap2 2.30: ``-x splice --cs`` emits ``de:f``
(gap-compressed per-base divergence, NOT ``dv:f`` — ``dv:f`` appears only
without base alignment and is coarse) plus the intron-aware ``cs`` string.
"""
from __future__ import annotations

import pytest

import reconcile_candidates as rc


# ---- the bug: multi-exon alen spans introns, nmatch/alen collapses ----------

def test_paf_multiexon_identity_from_divergence_not_alen(tmp_path):
    """A 3-exon transcript: nmatch=894 over alen=2400 (alen spans two ~750 bp
    introns). The legacy nmatch/alen formula gives 37.25% and fails the 95
    gate; the de:f divergence tag gives the true 99.33% and passes."""
    paf = tmp_path / "m.paf"
    # q qlen qs qe strand t tlen ts te nmatch alen mapq  tags...
    paf.write_text(
        "tx1\t900\t0\t900\t+\tchr1\t100000\t1000\t3400\t894\t2400\t60\t"
        "tp:A:P\tde:f:0.0067\tcg:Z:300M750N300M750N300M\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert 100.0 * 894 / 2400 < 95.0             # the bug: nmatch/alen = 37.25%
    assert abs(p.pct_identity - 99.33) < 0.01    # 100 * (1 - 0.0067)
    assert p.pct_coverage == 100.0               # 100 * (900-0)/900 (unaffected)
    assert rc.pass_gate(p, 95.0, 90.0) is True   # now clears the gate


def test_paf_multiexon_identity_from_cs_excludes_intron(tmp_path):
    """cs fallback (no de:f/dv:f tag): the ``~`` intron op is skipped, so a
    two-exon perfect alignment reads 100% even though alen (col 11) spans the
    intron and nmatch/alen would give 30%."""
    paf = tmp_path / "m.paf"
    paf.write_text(
        "tx2\t300\t0\t300\t+\tchr1\t100000\t1000\t2000\t300\t1000\t60\t"
        "tp:A:P\tcs:Z::150~gt700ag:150\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert 100.0 * 300 / 1000 < 95.0             # the bug: nmatch/alen = 30%
    assert p.pct_identity == 100.0               # 300 matches / 300 aligned
    assert rc.pass_gate(p, 95.0, 90.0) is True


def test_paf_cs_arithmetic_matches_over_aligned_columns(tmp_path):
    """Pins the exact cs formula: matches / (matches + subs + ins + del),
    introns excluded. cs :100*ac+gg:50-c:48 -> 198 matches, 1 sub, 2 ins,
    1 del -> 198/202."""
    paf = tmp_path / "m.paf"
    paf.write_text(
        "tx3\t200\t0\t200\t+\tchr1\t5000\t100\t400\t198\t300\t60\t"
        "cs:Z::100*ac+gg:50-c:48\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert abs(p.pct_identity - 100.0 * 198 / 202) < 1e-9


# ---- no regression on a gap-free single-exon record -------------------------

def test_paf_single_exon_divergence_agrees_with_nmatch_over_alen(tmp_path):
    """Single-exon alignment has no intron in alen, so the de:f identity and
    the legacy nmatch/alen value agree — the fix does not move single-exon
    identities."""
    paf = tmp_path / "m.paf"
    paf.write_text(
        "tx1\t300\t0\t300\t+\tchr1\t10000\t500\t800\t297\t300\t60\t"
        "tp:A:P\tde:f:0.01\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert abs(p.pct_identity - 99.0) < 1e-6            # 100 * (1 - 0.01)
    assert abs(p.pct_identity - 100.0 * 297 / 300) < 1e-6
    assert p.pct_coverage == 100.0


def test_paf_dv_tag_accepted_when_de_absent(tmp_path):
    """Without base alignment minimap2 emits dv:f (not de:f); it is still an
    intron-aware divergence and is accepted as the identity source."""
    paf = tmp_path / "m.paf"
    paf.write_text(
        "tx1\t300\t0\t300\t+\tchr1\t10000\t500\t800\t297\t300\t60\t"
        "tp:A:P\tdv:f:0.02\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert abs(p.pct_identity - 98.0) < 1e-6            # 100 * (1 - 0.02)


# ---- FAIL LOUD: no intron-aware identity source -----------------------------

def test_paf_no_identity_source_raises_naming_record(tmp_path):
    """A record with neither a de:f/dv:f divergence tag nor a cs:Z tag has no
    intron-aware identity source. The parser MUST NOT silently fall back to
    nmatch/alen (the bug) or to 0 — it raises, naming the record."""
    paf = tmp_path / "m.paf"
    paf.write_text(
        "tx1\t900\t0\t900\t+\tchr1\t100000\t1000\t3400\t894\t2400\t60\ttp:A:P\n")
    with pytest.raises(ValueError, match="tx1"):
        rc.parse_minimap2_paf(str(paf))


def test_paf_no_identity_source_does_not_return_nmatch_over_alen(tmp_path):
    """Guards the specific silent-failure mode: the raise fires INSTEAD of
    returning the buggy 37.25% (nmatch/alen) placement."""
    paf = tmp_path / "m.paf"
    paf.write_text(
        "tx1\t900\t0\t900\t+\tchr1\t100000\t1000\t3400\t894\t2400\t60\ttp:A:P\n")
    with pytest.raises(ValueError):
        rc.parse_minimap2_paf(str(paf))


# ---- end-to-end: multi-exon transcript is now placed concordantly -----------

def test_concordant_multiexon_placed_high_confidence(tmp_path):
    """The whole point of bead kfqz: a multi-exon transcript that minimap2 and
    GMAP both place at the same locus is now concordant, because the minimap2
    identity is intron-aware (de:f -> 99.33%) and clears the 95 gate. Under the
    old nmatch/alen math the minimap2 side scored 37% and the pair was never
    concordant, so the transcript was wrongly left unplaced."""
    paf = tmp_path / "m.paf"
    paf.write_text(
        "tx1\t900\t0\t900\t+\tchr1\t100000\t1000\t3400\t894\t2400\t60\t"
        "tp:A:P\tde:f:0.0067\tcg:Z:300M750N300M750N300M\n")
    mm = rc.parse_minimap2_paf(str(paf))
    gmap = [rc.Placement("tx1", "chr1", 1000, 3400, "+", 99.0, 100.0, "gmap", 298)]

    # New (fixed) identity -> the pair is concordant, one distinct locus.
    agreed = rc._concordant_placements(mm, gmap, 95.0, 90.0)
    assert len(agreed) == 1

    res = rc.place_transcript("tx1", mm, gmap, [], [], rc.Thresholds())
    assert res.placement is not None
    assert res.method == "minimap2+gmap_concordant"
    assert res.confidence == "high"

    # Contrast: the OLD nmatch/alen identity (37.25%) fails the gate, so the
    # same pair would NOT have been concordant (the transcript stayed unplaced).
    old = [rc.Placement("tx1", "chr1", 1000, 3400, "+", 100.0 * 894 / 2400,
                        100.0, "minimap2", 60)]
    assert rc._concordant_placements(old, gmap, 95.0, 90.0) == []
