"""Tests for the paralog-discrimination tightening in
``scripts/recover_cds_from_assemblies`` (bead -8st, 2026-05-08).

The previous CDS-recovery logic accepted any miniprot hit whose
back-translated identity to the query exceeded 50 %. For chemoreceptor
LSE expansions where paralogs share 80–95 % identity, that threshold is
useless: a wrong-paralog miniprot match clears it easily, and the wrong
paralog's CDS gets attached to the query protein's identity. dN/dS for
that "gene" then mixes one paralog's protein with another paralog's CDS
and produces a meaningless omega.

Three layers of tightening are exercised here:

    L1 (alignment-based identity, raised threshold):
        ``_alignment_identity`` returns identity over a global pairwise
        alignment, not zip-compared positions. ``CDS_PROTEIN_IDENTITY_MIN``
        defaults to 0.95 — paralog matches at 0.85 are now rejected as
        ``low_identity``.

    L2 (ambiguity margin):
        When a query has multiple miniprot hits, the best must beat the
        second-best by at least ``CDS_AMBIGUITY_MARGIN`` (default 0.05)
        to be assigned. Otherwise the assignment is ``ambiguous_top_hits``.

    L3 (genomic-locus dedup):
        Two query proteins whose accepted hits land at the same genomic
        locus (overlap ≥ ``CDS_OVERLAP_FRACTION``) cannot both be right.
        Only the higher-identity one keeps the CDS; the other's hit is
        recorded as ``paralog_overlap_rejected``.
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

import pytest

# Defer import until we can patch env vars per test (the thresholds are
# captured at module import time).
_RECOVER_MOD_NAME = "recover_cds_from_assemblies"


def _fresh_module(env: dict[str, str] | None = None):
    """Reload the recover_cds module under a controlled environment.

    Module-level constants (CDS_PROTEIN_IDENTITY_MIN, etc.) are read from
    os.environ at import time, so we mutate the environment before
    `importlib.reload` to test threshold customisation cleanly.
    """
    for k in ("CDS_PROTEIN_IDENTITY_MIN", "CDS_AMBIGUITY_MARGIN",
              "CDS_OVERLAP_FRACTION"):
        os.environ.pop(k, None)
    if env:
        os.environ.update(env)
    if _RECOVER_MOD_NAME in sys.modules:
        return importlib.reload(sys.modules[_RECOVER_MOD_NAME])
    return importlib.import_module(_RECOVER_MOD_NAME)


# ---- L1: _alignment_identity --------------------------------------------

def test_alignment_identity_perfect_match() -> None:
    mod = _fresh_module()
    assert mod._alignment_identity("MKILFV", "MKILFV") == pytest.approx(1.0)


def test_alignment_identity_paralog_pair() -> None:
    """Two GPCR-like sequences differing at ~15 % of positions should
    score around 0.85 — exactly the regime where the previous 0.5
    threshold leaks paralog confusions but a 0.95 threshold catches them."""
    mod = _fresh_module()
    a = "MKLLVALSILVTAQAQDLNTPLSALLAFTNRSPCTEYWDQYTVRSGCLAQGRLATFTLPL"
    # 9 substitutions / 60 = 15 % divergence
    b = "MKVLVALSILVTAQAQDLDTPLSALLPFTNRSPCREYWDLYTVRSACLAQGNLATFTLPV"
    ident = mod._alignment_identity(a, b)
    assert 0.80 <= ident <= 0.92, f"expected ~0.85, got {ident:.3f}"


def test_alignment_identity_handles_empty_inputs() -> None:
    mod = _fresh_module()
    assert mod._alignment_identity("", "MKIL") == 0.0
    assert mod._alignment_identity("MKIL", "") == 0.0
    assert mod._alignment_identity("", "") == 0.0


def test_alignment_identity_rejects_paralog_at_default_threshold() -> None:
    """The 0.95 default threshold must reject a 0.85-identity paralog
    pair — that's the load-bearing chemoreceptor-discrimination behavior."""
    mod = _fresh_module()
    a = "MKLLVALSILVTAQAQDLNTPLSALLAFTNRSPCTEYWDQYTVRSGCLAQGRLATFTLPL"
    b = "MKVLVALSILVTAQAQDLDTPLSALLPFTNRSPCREYWDLYTVRSACLAQGNLATFTLPV"
    assert mod._alignment_identity(a, b) < mod.CDS_PROTEIN_IDENTITY_MIN


# ---- L2: ambiguity margin (via direct constant + integration) -----------

def test_ambiguity_margin_default_is_five_percent() -> None:
    mod = _fresh_module()
    assert mod.CDS_AMBIGUITY_MARGIN == pytest.approx(0.05)


def test_ambiguity_margin_overridable_via_env() -> None:
    mod = _fresh_module({"CDS_AMBIGUITY_MARGIN": "0.10"})
    assert mod.CDS_AMBIGUITY_MARGIN == pytest.approx(0.10)


# ---- L3: _intervals_overlap ---------------------------------------------

def test_intervals_overlap_disjoint_returns_false() -> None:
    mod = _fresh_module()
    a = [(100, 200), (300, 400)]
    b = [(1000, 1100)]
    assert mod._intervals_overlap(a, b) is False


def test_intervals_overlap_adjacent_no_overlap() -> None:
    """Touching but not overlapping intervals do not register."""
    mod = _fresh_module()
    a = [(100, 200)]
    b = [(201, 300)]
    assert mod._intervals_overlap(a, b) is False


def test_intervals_overlap_substantial_overlap() -> None:
    """When the smaller span is mostly contained in the larger, the
    overlap fraction exceeds the default 0.5 and the function returns
    True — this is the paralog-locus collision signature."""
    mod = _fresh_module()
    a = [(100, 1100)]    # 1001 bp
    b = [(500, 1100)]    # 601 bp, fully contained in a
    assert mod._intervals_overlap(a, b) is True


def test_intervals_overlap_partial_below_threshold() -> None:
    """A 30 % overlap of the smaller does not count under the default
    0.5 threshold (different paralogs can have neighboring loci with
    a small genuine overlap at exon boundaries)."""
    mod = _fresh_module()
    a = [(100, 1000)]    # 901 bp
    b = [(900, 1100)]    # 201 bp; overlap = 101 bp = 50.2 % of b
    # 50.2 % is just over the threshold, so True — this is the sensitive
    # boundary. Test the just-below case explicitly.
    a2 = [(100, 1000)]
    b2 = [(951, 1100)]   # 150 bp; overlap = 50 bp = 33.3 %
    assert mod._intervals_overlap(a2, b2) is False


def test_intervals_overlap_handles_empty_input() -> None:
    mod = _fresh_module()
    assert mod._intervals_overlap([], [(100, 200)]) is False
    assert mod._intervals_overlap([(100, 200)], []) is False


def test_intervals_overlap_threshold_overridable() -> None:
    mod = _fresh_module({"CDS_OVERLAP_FRACTION": "0.10"})
    a = [(100, 1000)]
    b = [(951, 1100)]    # 33.3 % of smaller — passes 0.10 but not 0.50
    assert mod._intervals_overlap(a, b) is True


# ---- Threshold defaults are appropriately strict ------------------------

def test_protein_identity_default_is_strict() -> None:
    """Document the default — 0.95 is the load-bearing setting for the
    chemoreceptor LSE-expansion use case. If this ever loosens we want
    a test failure to force a deliberate decision."""
    mod = _fresh_module()
    assert mod.CDS_PROTEIN_IDENTITY_MIN == pytest.approx(0.95), (
        "CDS_PROTEIN_IDENTITY_MIN default must be 0.95; loosening it "
        "would re-enable the paralog-confusion bug from bead -8st."
    )


def test_protein_identity_overridable_via_env() -> None:
    mod = _fresh_module({"CDS_PROTEIN_IDENTITY_MIN": "0.80"})
    assert mod.CDS_PROTEIN_IDENTITY_MIN == pytest.approx(0.80)


# Cleanup: restore default-importable module state for other tests.
@pytest.fixture(autouse=True, scope="module")
def _restore_module_after_tests():
    yield
    _fresh_module()
