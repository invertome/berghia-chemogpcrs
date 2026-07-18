"""Tests for the per-OG dN/dS-axis reliability handling (bead -8st, 2026-05-08;
application point moved 2026-07-18).

The per-OG ``dnds_reliability_weight`` (in [0, 1]) down-weights the dN/dS signal
for under-supported OGs (chemoreceptor LSE expansions where most paralogs lack a
reference CDS), where an aBSREL/BUSTED estimate is underpowered. That reliability
used to be applied as a MULTIPLIER *inside* ``calculate_rank_score``. It has now
moved UPSTREAM: ``dnds_reliability.reliability_shrink`` shrinks the
``positive_score_norm``/``purifying_score_norm`` columns toward the cohort NEUTRAL
(median) once, on the dataframe BOTH ranking paths read (the weighted composite
and the production rank-aggregation ranklist). So an underpowered candidate
contributes a neutral value (never a penalizing low score), and
``calculate_rank_score`` itself is now reliability-AGNOSTIC (re-applying the
weight here would double-count it). These tests assert that moved behavior.
"""
from __future__ import annotations

import importlib
import os
import sys
from typing import Any

import pandas as pd
import pytest

from dnds_reliability import reliability_shrink   # conftest adds scripts/ to path


def _import_rank_candidates_lib():
    """Load rank_candidates.py as a module without running its top-level
    side effects (which require a populated results dir + CLI args).

    rank_candidates.py is a script, so we extract `calculate_rank_score`
    via importlib.util after reading the source — this is the same trick
    test_ranking_lib.py uses.
    """
    return _exec_calc_rank_score_only()


def _exec_calc_rank_score_only() -> Any:
    """Compile rank_candidates.py up to and including ``calculate_rank_score``,
    then return the function. Avoids running the script's CLI arg parsing
    and downstream loaders.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.normpath(os.path.join(here, '..', '..'))
    src_path = os.path.join(repo, 'scripts', 'rank_candidates.py')
    with open(src_path) as f:
        src = f.read()
    # Find the `def calculate_rank_score` start and the next top-level def
    # after it; execute just that fragment in a controlled namespace.
    start = src.find('def calculate_rank_score(')
    assert start != -1, "calculate_rank_score not found"
    # End at the next top-level def
    end_marker = '\n\ndef '
    end = src.find(end_marker, start + 1)
    fragment = src[start:end] if end != -1 else src[start:]
    ns: dict[str, Any] = {'pd': pd, 'np': __import__('numpy')}
    exec(fragment, ns)
    return ns['calculate_rank_score']


# ---- Plumbing fixtures ---------------------------------------------------

@pytest.fixture(scope='module')
def calc():
    return _exec_calc_rank_score_only()


def _row(positive: float = 0.0, purifying: float = 0.0,
         phylo: float = 0.0, lse: float = 0.0,
         dnds_rw: float = 1.0) -> dict[str, Any]:
    """Build a minimal scoring row. Synteny / expression / etc. are all
    marked unavailable so the multiplier's effect is isolated."""
    return {
        'phylo_score_norm': phylo,
        'purifying_score_norm': purifying,
        'positive_score_norm': positive,
        'lse_depth_score_norm': lse,
        'dnds_reliability_weight': dnds_rw,
        'has_synteny_data': False,
        'has_expression_data': False,
        'has_chemosensory_expr_data': False,
        'has_gprotein_data': False,
        'has_ecl_data': False,
        'has_expansion_data': False,
        'has_og_confidence_data': False,
    }


WEIGHTS = {
    'phylo': 2.0, 'purifying': 1.0, 'positive': 1.0, 'lse_depth': 1.0,
    'synteny': 3.0, 'expr': 1.0, 'chemosensory_expr': 3.0,
    'gprotein_coexpr': 2.0, 'ecl_divergence': 2.0, 'expansion': 1.0,
    'og_confidence': 1.0,
}


# ---- Reliability now shrinks the norm columns UPSTREAM (not in the scorer) --

def test_full_reliability_matches_baseline(calc) -> None:
    """At dnds_rw = 1.0 the score must equal the plain baseline (identity)."""
    df = pd.DataFrame([_row(positive=0.8, dnds_rw=1.0)])
    out = calc(df, WEIGHTS).iloc[0]
    # Active axes: phylo (0)*2 + purifying (0)*1 + positive (0.8)*1 + lse (0)*1
    # total_weight = 2 + 1 + 1 + 1 = 5; max_possible_weight = sum(WEIGHTS).
    expected = (0.0 + 0.0 + 0.8 + 0.0) / 5.0 * sum(WEIGHTS.values())
    assert out == pytest.approx(expected)


def test_calculate_rank_score_is_reliability_agnostic(calc) -> None:
    """Reliability has moved UPSTREAM into the norm columns, so
    ``calculate_rank_score`` must no longer down-weight dN/dS by
    ``dnds_reliability_weight``. Varying it on an otherwise identical row must
    NOT change the score — re-applying it here would double-count the weight
    (once in the already-shrunk norm value, once in the scorer)."""
    s_full = calc(pd.DataFrame([_row(positive=0.8, dnds_rw=1.0)]), WEIGHTS).iloc[0]
    s_half = calc(pd.DataFrame([_row(positive=0.8, dnds_rw=0.5)]), WEIGHTS).iloc[0]
    s_zero = calc(pd.DataFrame([_row(positive=0.8, dnds_rw=0.0)]), WEIGHTS).iloc[0]
    assert s_full == pytest.approx(s_half) == pytest.approx(s_zero)


def test_underpowered_positive_shrinks_to_neutral_not_bottom(calc) -> None:
    """The reliability EFFECT is preserved, only its application point moved:
    an underpowered (reliability approx 0) candidate whose RAW
    positive_score_norm is LOW is shrunk to the cohort NEUTRAL (median) before
    the scorer reads it, so it contributes a neutral value rather than a
    penalizing low one. After the shrink it must score identically to a
    candidate that genuinely sits at the median positive signal — not at the
    bottom — while a fully-powered candidate is unchanged."""
    raw = {"poor": 0.1, "mid": 0.5, "rich": 0.9}   # median 0.5
    rel = {"poor": 0.0, "mid": 1.0, "rich": 1.0}
    shrunk = reliability_shrink(raw, rel)
    assert shrunk["poor"] == pytest.approx(0.5)          # lifted to neutral
    assert shrunk["rich"] == pytest.approx(0.9)          # fully-powered unchanged
    poor = calc(pd.DataFrame([_row(positive=shrunk["poor"])]), WEIGHTS).iloc[0]
    genuine_mid = calc(pd.DataFrame([_row(positive=0.5)]), WEIGHTS).iloc[0]
    bottom = calc(pd.DataFrame([_row(positive=0.1)]), WEIGHTS).iloc[0]
    assert poor == pytest.approx(genuine_mid)
    assert poor > bottom


def test_shrink_applies_to_purifying_axis_too(calc) -> None:
    """Both norm columns are shrunk in the pre-write block — positive AND
    purifying share the same aBSREL omega foundation. An underpowered purifying
    signal is likewise lifted to neutral, so a purifying-only candidate scores
    the same as one genuinely at the median rather than being penalized."""
    raw = {"a": 0.9, "b": 0.1}   # median 0.5
    shrunk = reliability_shrink(raw, {"a": 1.0, "b": 0.0})
    assert shrunk["b"] == pytest.approx(0.5)   # underpowered purifying -> neutral
    neutralized = calc(pd.DataFrame([_row(purifying=shrunk["b"])]), WEIGHTS).iloc[0]
    genuine_mid = calc(pd.DataFrame([_row(purifying=0.5)]), WEIGHTS).iloc[0]
    assert neutralized == pytest.approx(genuine_mid)


def test_missing_reliability_defaults_to_full_identity() -> None:
    """A candidate absent from the reliability map defaults to full reliability
    (w = 1.0), so its norm value round-trips unchanged through the shrink — the
    no-op that protects pre-fix / legacy data from being silently neutralized."""
    out = reliability_shrink({"a": 0.9, "b": 0.1}, {"a": 1.0})   # b absent
    assert out["b"] == pytest.approx(0.1)
    assert out["a"] == pytest.approx(0.9)


def test_tandem_cluster_axis_included_in_sensitivity(calc) -> None:
    """Bead 0rg (MEDIUM-4): the sensitivity scorer (calculate_rank_score) must
    include the tandem-cluster axis — the field's signature chemoreceptor signal
    — not silently ignore it. Previously the sensitivity perturbed the tandem
    weight while the score never used it, so weight_importance for tandem read
    ~0 (misleading for the highest-weight axis)."""
    row = _row(positive=0.0, dnds_rw=1.0)
    row['has_tandem_cluster_data'] = True
    row['tandem_cluster_score_norm'] = 1.0
    df = pd.DataFrame([row])
    w = {**WEIGHTS, 'tandem_cluster': 2.5}
    out = calc(df, w).iloc[0]
    assert out > 0, "tandem-cluster contribution must reach the sensitivity score"
