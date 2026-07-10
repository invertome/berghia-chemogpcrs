"""dN/dS reliability in the PRODUCTION scorer (_rank_candidates_lib.calculate_fair_rank_score).

Bead -8st / berghia-chemogpcrs-0rg. The per-OG dN/dS reliability weight
(min(1, n_ref_cds / DNDS_RELIABILITY_FULL)) must scale the positive AND purifying
axes in BOTH the score numerator and the completeness denominator (total_weight),
so an under-supported omega estimate cleanly falls out and the candidate is judged
on its other axes via fair-scoring — the option-(ii) "clean axis removal" design
documented in rank_candidates.py:133-136. It must NOT be treated as missing data
(which would drag completeness down and double-penalize a reference-poor OG).
"""
from __future__ import annotations

import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts")
sys.path.insert(0, os.path.normpath(SCRIPTS))

import _rank_candidates_lib as lib  # noqa: E402

WEIGHTS = {"phylo": 2.0, "positive": 1.0, "purifying": 1.0, "synteny": 3.0}


def test_reliability_1_is_noop():
    scores = {"phylo": 0.9, "positive": 0.8, "purifying": 0.0, "synteny": None}
    base = lib.calculate_fair_rank_score(scores, WEIGHTS)
    rel = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=1.0)
    assert rel == pytest.approx(base)


def test_reliability_0_equals_axis_removed_entirely():
    # Option (ii): dnds_reliability=0 must equal removing the dN/dS axes from
    # BOTH scores and weights (clean removal) — NOT setting them to None.
    scores = {"phylo": 0.9, "positive": 0.8, "purifying": 0.7, "synteny": None}
    rw0 = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=0.0)
    removed = lib.calculate_fair_rank_score(
        {"phylo": 0.9, "synteny": None}, {"phylo": 2.0, "synteny": 3.0})
    assert rw0 == pytest.approx(removed)


def test_reliability_0_not_same_as_missing_data():
    # Distinguishes option (ii) from option (i): treating unreliable dN/dS as
    # missing data (None, total_weight unchanged) would give a LOWER score than
    # clean removal, because completeness would be dragged down. Phylo-dominant
    # weights keep both completeness values above the 0.4 floor so the
    # difference is visible.
    W = {"phylo": 5.0, "positive": 1.0, "purifying": 1.0}
    rw0 = lib.calculate_fair_rank_score(
        {"phylo": 0.9, "positive": 0.8, "purifying": 0.7}, W, dnds_reliability=0.0)
    as_missing = lib.calculate_fair_rank_score(
        {"phylo": 0.9, "positive": None, "purifying": None}, W)
    assert rw0 > as_missing


def test_reliability_half_is_between():
    scores = {"phylo": 0.2, "positive": 0.9, "purifying": 0.0, "synteny": None}
    s0 = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=0.0)
    s5 = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=0.5)
    s1 = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=1.0)
    # positive is the strong signal; more reliability -> more of it counts.
    assert s0 < s5 < s1


def test_reliability_scales_both_positive_and_purifying():
    # A purifying-only row must also respond to reliability (both dN/dS axes
    # share the same aBSREL omega foundation).
    scores = {"phylo": 0.0, "positive": 0.0, "purifying": 0.9, "synteny": None}
    s_full = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=1.0)
    s_half = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=0.5)
    assert s_half < s_full


def test_diagnostics_completeness_excludes_unreliable_dnds():
    # With dN/dS removed (rw=0), evidence_completeness is computed over the
    # remaining axes' weights only.
    scores = {"phylo": 0.9, "positive": 0.8, "purifying": 0.7, "synteny": None}
    out = lib.calculate_fair_rank_score(scores, WEIGHTS, dnds_reliability=0.0,
                                        return_diagnostics=True)
    # Only phylo (weight 2) has data; effective total = phylo(2)+synteny(3)=5.
    assert out["available_weight"] == pytest.approx(2.0)
    assert out["evidence_completeness"] == pytest.approx(max(0.4, 2.0 / 5.0))
