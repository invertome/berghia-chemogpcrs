"""Tests for scripts/or_microswitch.py -- the cheap OR-microswitch structural
flag from Ballesteros-Weinstein (BW) numbered residues (Task 6 of the ML/PLM
chemoreceptor ranking plan, docs/plans/2026-07-01-ml-plm-chemoreceptor-
ranking.md).

The upstream BW-numbering assignment (structurally aligning an AlphaFold
model to a numbered class-A GPCR template) is OUT OF SCOPE here -- a Unity/
upstream step. This module is a pure, cheap classifier over an
already-produced ``{bw_position: residue}`` map for one candidate.
"""
from or_microswitch import or_microswitch_flag


# --- flag = 1: altered toggle + conserved motif anchor ----------------------

def test_altered_toggle_tyrosine_with_conserved_anchor_is_flagged():
    assert or_microswitch_flag({"6.48": "Y", "6.50": "P"}) == 1


def test_altered_toggle_phenylalanine_variant_is_also_flagged():
    assert or_microswitch_flag({"6.48": "F", "6.50": "P"}) == 1


# --- flag = 0: BW residues present, but the pattern isn't met ---------------

def test_canonical_tryptophan_toggle_is_not_flagged():
    assert or_microswitch_flag({"6.48": "W", "6.50": "P"}) == 0


def test_altered_toggle_without_conserved_anchor_is_not_flagged():
    # 6.48 substituted, but the CWxP anchor proline is ALSO gone -- not the
    # documented altered-motif signature, just an unrelated substitution.
    assert or_microswitch_flag({"6.48": "Y", "6.50": "A"}) == 0


def test_neither_position_altered_is_not_flagged():
    assert or_microswitch_flag({"6.48": "W", "6.50": "A"}) == 0


# --- flag = None: insufficient BW-numbered residues -------------------------

def test_missing_toggle_position_is_none():
    assert or_microswitch_flag({"6.50": "P"}) is None


def test_missing_anchor_position_is_none():
    assert or_microswitch_flag({"6.48": "Y"}) is None


def test_empty_map_is_none():
    assert or_microswitch_flag({}) is None


def test_none_input_is_none():
    assert or_microswitch_flag(None) is None


# --- robustness --------------------------------------------------------------

def test_lowercase_residue_codes_are_handled():
    assert or_microswitch_flag({"6.48": "y", "6.50": "p"}) == 1


def test_extra_unrelated_positions_are_ignored():
    assert or_microswitch_flag({"6.48": "Y", "6.50": "P", "3.50": "R"}) == 1
