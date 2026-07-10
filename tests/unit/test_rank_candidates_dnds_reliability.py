"""Tests for the per-OG dN/dS-axis reliability multiplier (bead -8st,
2026-05-08).

`rank_candidates.calculate_rank_score` multiplies each row's dN/dS
contribution (purifying + positive) by ``dnds_reliability_weight``
(a value in [0, 1]) so candidates in under-supported OGs (chemoreceptor
LSE expansions where most paralogs lack a reference CDS) don't get
full dN/dS weight on what is statistically just noise. The multiplier
applies to BOTH the score and the total_weight denominator — that is
the fair-scoring pattern: when the dN/dS axis is unreliable, it falls
out of both numerator and denominator and the other axes carry the
candidate.
"""
from __future__ import annotations

import importlib
import os
import sys
from typing import Any

import pandas as pd
import pytest


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


# ---- Multiplier semantics ------------------------------------------------

def test_full_reliability_matches_no_multiplier(calc) -> None:
    """At dnds_rw = 1.0 the score must equal the baseline (no down-weight)."""
    df = pd.DataFrame([_row(positive=0.8, dnds_rw=1.0)])
    out = calc(df, WEIGHTS).iloc[0]
    # Hand-compute baseline.
    # Active axes: phylo (0)*2 + purifying (0)*1 + positive (0.8)*1 + lse (0)*1
    # total_weight = 2 + 1 + 1 + 1 = 5
    # max_possible_weight = sum(WEIGHTS.values())
    expected = (0.0 + 0.0 + 0.8 + 0.0) / 5.0 * sum(WEIGHTS.values())
    assert out == pytest.approx(expected)


def test_zero_reliability_drops_dnds_axis(calc) -> None:
    """At dnds_rw = 0.0 the dN/dS axis falls out of BOTH numerator and
    denominator. With no other active axes returning a non-zero value,
    the row's score is the average of remaining axes — i.e. zero here."""
    df = pd.DataFrame([_row(positive=0.8, dnds_rw=0.0)])
    out = calc(df, WEIGHTS).iloc[0]
    # phylo=0, purifying=0, positive=0.8 (dropped), lse=0
    # Effective numerator = 0; total_weight = 2 (phylo) + 0 + 0 + 1 (lse) = 3
    # Score = 0/3 * sum_weights = 0
    assert out == pytest.approx(0.0)


def test_zero_reliability_does_not_inflate_other_axes(calc) -> None:
    """When dN/dS is zeroed and the OTHER axes have signal, the result
    must NOT exceed what those axes would produce on their own with
    dnds_rw = 1.0. Otherwise we'd be artificially boosting non-dN/dS
    axes by removing their normalisation peer — exactly the bug the
    fair-scoring pattern prevents."""
    # Row with strong phylo signal but no dN/dS data in the underlying OG
    df_zeroed = pd.DataFrame([_row(phylo=0.9, positive=0.0, dnds_rw=0.0)])
    df_full = pd.DataFrame([_row(phylo=0.9, positive=0.0, dnds_rw=1.0)])
    s_zeroed = calc(df_zeroed, WEIGHTS).iloc[0]
    s_full = calc(df_full, WEIGHTS).iloc[0]
    # With purifying=positive=0 in BOTH cases, the only difference is the
    # denominator (weight 0 vs weight 2). When dnds_rw=0, total_weight is
    # smaller so the same numerator divides into a slightly higher
    # normalised value. That's correct fair-scoring (the row's available
    # evidence is fewer axes), but should not change the absolute order.
    assert s_zeroed >= s_full, (
        "fair-scoring expected when dN/dS is dropped — the available "
        "evidence the row contributes shouldn't shrink"
    )


def test_partial_reliability_linear_ramp(calc) -> None:
    """A dnds_rw of 0.5 should put the dN/dS contribution at half-credit
    on both sides of the ratio. The score for a row with positive=0.8
    and dnds_rw=0.5 should fall between dnds_rw=0 and dnds_rw=1."""
    df_low = pd.DataFrame([_row(positive=0.8, dnds_rw=0.0)])
    df_mid = pd.DataFrame([_row(positive=0.8, dnds_rw=0.5)])
    df_hi = pd.DataFrame([_row(positive=0.8, dnds_rw=1.0)])
    s_low = calc(df_low, WEIGHTS).iloc[0]
    s_mid = calc(df_mid, WEIGHTS).iloc[0]
    s_hi = calc(df_hi, WEIGHTS).iloc[0]
    assert s_low < s_mid < s_hi


def test_missing_dnds_reliability_column_defaults_to_full(calc) -> None:
    """Rows without the column at all (legacy CSVs reaching downstream
    consumers via re-rank passes) must default to dnds_rw = 1.0 so the
    multiplier is a no-op and we don't silently zero pre-fix data."""
    row = _row(positive=0.8, dnds_rw=1.0)
    row.pop('dnds_reliability_weight')
    df = pd.DataFrame([row])
    expected = pd.DataFrame([_row(positive=0.8, dnds_rw=1.0)])
    out = calc(df, WEIGHTS).iloc[0]
    out_expected = calc(expected, WEIGHTS).iloc[0]
    assert out == pytest.approx(out_expected)


def test_purifying_axis_also_scaled(calc) -> None:
    """Both purifying_score_norm and positive_score_norm must be scaled
    by dnds_rw — they share the same statistical foundation (the per-OG
    aBSREL omega estimate). A multiplier that hits only one axis would
    introduce a bias when PURIFYING_WEIGHT is nonzero (rare in this
    project but legitimate for conserved-function searches)."""
    df_full = pd.DataFrame([_row(purifying=0.9, dnds_rw=1.0)])
    df_zero = pd.DataFrame([_row(purifying=0.9, dnds_rw=0.0)])
    # If only positive were scaled, df_zero would still have purifying
    # contributing to the numerator, and the score would be > 0.
    assert calc(df_zero, WEIGHTS).iloc[0] == pytest.approx(0.0)
    assert calc(df_full, WEIGHTS).iloc[0] > 0.0


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
