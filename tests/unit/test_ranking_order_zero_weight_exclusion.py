"""A signal weighted 0 must not vote under RANK_METHOD=rankagg.

config.sh sets ``PURIFYING_WEIGHT=0`` deliberately (bead -ea9: chemoreceptor
discovery rewards diversifying selection on extracellular loops, NOT
whole-gene purifying selection -- a conserved housekeeping GPCR should not
rank highly for the chemoreceptor question). Under the weighted scorer that
zero removes the axis entirely.

RRA is weight-free BY DESIGN and should stay so. But `purifying` is an
unconditional entry in SIGNAL_SPEC, so under the production default
RANK_METHOD=rankagg it votes at FULL strength -- the zero weight's EXCLUSION
semantics are silently discarded. Measured at realistic scale (790
candidates, 12 signals): a conserved GPCR ranked 1st on purifying and around
the 65th percentile everywhere else lands at position 52 instead of 71, a
19-place lift from an axis meant to contribute nothing.

The mechanism here is an EXPLICIT exclusion set, not a reinterpretation of
weights as RRA inputs: zero-weight signals are named, the exclusion is
visible in the emitted provenance, and RANKAGG_EXCLUDED_SIGNALS can override
it in either direction.

Same shape, second half: `phylo` and `lse_divergence` are ungated (flag None) in
SIGNAL_SPEC but gate on `has_phylo_data` in the weighted scorer
(rank_candidates.py:1513-1514), so bead o98's phylo-absence fix is live on
only one of the two ranking paths.
"""
from __future__ import annotations

import pandas as pd
import pytest

import rank_aggregation as ra


# --------------------------------------------------------------------------
# excluded_signals_from_weights: the explicit, visible derivation
# --------------------------------------------------------------------------

def test_zero_weight_signals_become_the_exclusion_set() -> None:
    weights = {"phylo": 2, "purifying": 0, "positive": 1, "synteny": 2}
    assert ra.excluded_signals_from_weights(weights) == {"purifying"}


def test_no_zero_weights_excludes_nothing() -> None:
    assert ra.excluded_signals_from_weights(
        {"phylo": 2, "positive": 1}) == set()


def test_several_zero_weights_all_excluded() -> None:
    assert ra.excluded_signals_from_weights(
        {"phylo": 0, "purifying": 0, "positive": 1}) == {"phylo", "purifying"}


def test_env_override_replaces_the_derived_set(monkeypatch) -> None:
    """The policy stays reversible without editing weights."""
    monkeypatch.setenv("RANKAGG_EXCLUDED_SIGNALS", "synteny")
    assert ra.excluded_signals_from_weights(
        {"purifying": 0, "synteny": 2}) == {"synteny"}


def test_env_override_can_disable_exclusion_entirely(monkeypatch) -> None:
    """Disabling the exclusion now needs the explicit `none` token.

    This used to be spelled as the empty string. That made the override a
    trapdoor (bead wtwi): `export RANKAGG_EXCLUDED_SIGNALS=` produces the same
    "" that a deliberate "exclude nothing" produced, so an accidental empty
    assignment silently restored every zero-weighted signal to full strength.
    The capability is unchanged; only the spelling is, so that the deliberate
    value can no longer be produced by accident. See
    tests/unit/test_small_rankagg_excluded_signals_trapdoor.py.
    """
    monkeypatch.setenv("RANKAGG_EXCLUDED_SIGNALS", "none")
    assert ra.excluded_signals_from_weights({"purifying": 0}) == set()


def test_empty_env_override_falls_back_to_the_derivation(monkeypatch) -> None:
    """The trapdoor itself: an empty value states no policy, so it decides nothing."""
    monkeypatch.setenv("RANKAGG_EXCLUDED_SIGNALS", "")
    assert ra.excluded_signals_from_weights({"purifying": 0}) == {"purifying"}


# --------------------------------------------------------------------------
# the exclusion must actually remove the vote
# --------------------------------------------------------------------------

N = 100


def _df():
    """`conserved` tops purifying and sits ~65th percentile on every other axis.

    Sized at N=100 rather than a handful because RRA's ``min(rho * m, 1.0)``
    correction saturates on small inputs -- every candidate would hit the cap
    and tie, and the ordering would fall through to id-string order, masking
    the effect under test. (That saturation is itself an open P0; see the
    session report.)
    """
    rows = []
    for i in range(N):
        rows.append({
            "id": f"c{i:03d}",
            "phylo_score_norm": 1.0 - i / N,
            "purifying_score_norm": i / N,
            "positive_score_norm": 1.0 - i / N,
            "lse_divergence_score_norm": 1.0 - i / N,
            "has_phylo_data": True,
            "has_dnds_data": True,
        })
    rows.append({
        "id": "conserved",
        "phylo_score_norm": 0.35,
        "purifying_score_norm": 5.0,   # best on the zero-weight axis
        "positive_score_norm": 0.35,
        "lse_divergence_score_norm": 0.35,
        "has_phylo_data": True,
        "has_dnds_data": True,
    })
    return pd.DataFrame(rows)


def test_excluded_signal_is_absent_from_the_ranklists() -> None:
    lists = ra.build_ranklists_from_df(_df(), excluded={"purifying"})
    assert "purifying" not in lists
    assert "positive" in lists, "only the named signal may be dropped"


def test_exclusion_removes_the_zero_weight_lift() -> None:
    """The conserved GPCR must not be lifted by the axis meant to be silent."""
    df = _df()
    with_vote = ra.aggregate(ra.build_ranklists_from_df(df), method="rra")
    without = ra.aggregate(
        ra.build_ranklists_from_df(df, excluded={"purifying"}), method="rra")
    assert with_vote.index("conserved") < without.index("conserved"), (
        "excluding the zero-weight purifying axis should DEMOTE the "
        "purifying-topping candidate"
    )


def test_default_still_includes_every_signal() -> None:
    """No `excluded` argument = today's behaviour, exactly."""
    df = _df()
    assert set(ra.build_ranklists_from_df(df)) == set(
        ra.build_ranklists_from_df(df, excluded=set()))


def test_rerank_output_honours_the_exclusion(monkeypatch) -> None:
    df = _df().assign(final_rank=range(1, N + 2))
    kept = ra.rerank_output(df, "rankagg", excluded={"purifying"})
    voted = ra.rerank_output(df, "rankagg")
    assert list(kept["id"]).index("conserved") > list(voted["id"]).index("conserved")


# --------------------------------------------------------------------------
# phylo / lse_divergence must gate on has_phylo_data on BOTH paths (bead o98)
# --------------------------------------------------------------------------

def test_phylo_and_lse_divergence_gate_on_has_phylo_data() -> None:
    spec = {entry[0]: entry[1] for entry in ra.SIGNAL_SPEC}
    assert spec["phylo"] == "has_phylo_data", (
        "phylo votes unconditionally under rankagg while the weighted scorer "
        "gates it on has_phylo_data (bead o98)"
    )
    assert spec["lse_divergence"] == "has_phylo_data"


def test_selection_axes_gate_on_has_dnds_data() -> None:
    spec = {entry[0]: entry[1] for entry in ra.SIGNAL_SPEC}
    assert spec["purifying"] == "has_dnds_data"
    assert spec["positive"] == "has_dnds_data"


def test_candidate_outside_the_tree_casts_no_phylo_vote() -> None:
    df = pd.DataFrame([
        {"id": "in_tree", "phylo_score_norm": 0.9,
         "lse_divergence_score_norm": 0.9, "has_phylo_data": True},
        {"id": "off_tree", "phylo_score_norm": 0.0,
         "lse_divergence_score_norm": 0.0, "has_phylo_data": False},
    ])
    lists = ra.build_ranklists_from_df(df)
    assert "off_tree" not in lists.get("phylo", {}), (
        "an off-tree candidate must have no phylo signal, not a present 0.0"
    )
    assert "off_tree" not in lists.get("lse_divergence", {})
    assert "in_tree" in lists["phylo"]


def test_missing_flag_column_keeps_the_base_signals_voting() -> None:
    """Legacy CSVs with no has_phylo_data column must not lose every signal.

    Gated signals whose flag column is absent are skipped entirely, so making
    phylo/lse_divergence/purifying/positive gated would silently empty the whole
    ranklist for any pre-flag CSV. They must fall back to voting.
    """
    df = pd.DataFrame([
        {"id": "a", "phylo_score_norm": 0.9, "purifying_score_norm": 0.5},
        {"id": "b", "phylo_score_norm": 0.1, "purifying_score_norm": 0.2},
    ])
    lists = ra.build_ranklists_from_df(df)
    assert "phylo" in lists and len(lists["phylo"]) == 2
    assert "purifying" in lists and len(lists["purifying"]) == 2
