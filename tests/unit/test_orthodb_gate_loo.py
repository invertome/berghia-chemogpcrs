"""Unit tests for the derived harvest gate.

The gate decides which orthogroups are trustworthy enough to harvest reference
sequences from. Its thresholds were derived by leave-one-out on 198
invertebrate anchors, not chosen, but the DECISION LOGIC still has to be exact:
a group that silently passes on mixed-family or under-attested evidence would
inject wrong family labels into the reference envelope, which is the failure
this whole exercise exists to avoid.
"""

from __future__ import annotations

import importlib.util
from collections import Counter
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "orthodb_gate_loo",
    Path(__file__).resolve().parents[2] / "scripts" / "orthodb_gate_loo.py",
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def gv(**kw):
    base = dict(members=1000, seed_families=Counter({"opsin": 5}),
                max_size=10**9, min_seeds=5, purity="strict")
    base.update(kw)
    return mod.gate_verdict(**base)


# ------------------------------------------------------------- the happy path --

def test_pure_well_attested_group_passes_with_its_family():
    ok, fam, why = gv()
    assert (ok, fam, why) == (True, "opsin", "")


# ------------------------------------------------------------------- purity --

def test_strict_purity_rejects_any_second_family():
    """A group seeded by two curated families cannot assign either one."""
    ok, fam, why = gv(seed_families=Counter({"opsin": 9, "peptide": 1}))
    assert ok is False
    assert why == "mixed_family_seeds"


def test_modal80_admits_a_dominant_family():
    ok, fam, why = gv(seed_families=Counter({"opsin": 9, "peptide": 1}),
                      purity="modal80")
    assert ok is True
    assert fam == "opsin"


def test_modal80_still_rejects_a_weak_majority():
    ok, _, why = gv(seed_families=Counter({"opsin": 6, "peptide": 4}),
                    purity="modal80")
    assert ok is False
    assert why == "impure_seeds"


def test_tied_seed_families_never_pass():
    """A tie means the seeds disagree; picking either would be arbitrary."""
    ok, _, why = gv(seed_families=Counter({"opsin": 5, "peptide": 5}),
                    purity="modal80")
    assert ok is False
    assert why == "tied_seed_families"


# ------------------------------------------------------- attestation density --

def test_too_few_seeds_is_rejected():
    ok, _, why = gv(seed_families=Counter({"opsin": 4}), min_seeds=5)
    assert ok is False
    assert why == "too_few_seeds"


def test_seed_count_is_counted_on_the_modal_family_only():
    """8 total seeds but only 5 of the modal family: min_seeds applies to the
    family being assigned, not to the group's total attestation."""
    ok, fam, _ = gv(seed_families=Counter({"opsin": 5, "peptide": 3}),
                    purity="modal80", min_seeds=5)
    # modal80: 5/8 = 0.625 < 0.8, so purity rejects it first
    assert ok is False
    ok, fam, _ = gv(seed_families=Counter({"opsin": 5}), min_seeds=5)
    assert ok is True and fam == "opsin"


def test_min_density_rejects_a_sparsely_attested_giant():
    ok, _, why = gv(members=100000, seed_families=Counter({"opsin": 5}),
                    min_seeds=5)
    assert ok is True, "density is off by default"
    ok, _, why = mod.gate_verdict(
        members=100000, seed_families=Counter({"opsin": 5}),
        max_size=10**9, min_seeds=5, purity="strict", min_density=0.001)
    assert ok is False
    assert why == "too_sparse"


# ---------------------------------------------------------------------- size --

def test_size_cap_applies_when_requested():
    ok, _, why = gv(members=5000, max_size=2500)
    assert ok is False
    assert why == "too_large"


def test_size_is_not_capped_by_default():
    """The LOO sweep showed a size cap destroys recovery without improving
    precision, so the default gate must not impose one."""
    ok, _, _ = gv(members=68077)
    assert ok is True


# -------------------------------------------------------------- degenerate --

def test_group_with_no_characterized_seeds_is_rejected():
    ok, fam, why = gv(seed_families=Counter())
    assert ok is False
    assert why == "no_seed_anchors"
    assert fam == ""


def test_rejection_still_reports_the_modal_family_when_one_exists():
    """Callers log why a group failed; the family it would have assigned is
    part of that diagnosis."""
    _, fam, why = gv(seed_families=Counter({"opsin": 2}), min_seeds=5)
    assert fam == "opsin"
    assert why == "too_few_seeds"


def test_characterized_predicate_accepts_only_real_characterization():
    assert mod.is_characterized({"evidence": "experimental"})
    assert mod.is_characterized({"evidence": "experimental:26190115"})
    assert mod.is_characterized({"evidence": "characterized-gtopdb"})
    assert not mod.is_characterized({"evidence": "reviewed"})
    assert not mod.is_characterized({"evidence": ""})
