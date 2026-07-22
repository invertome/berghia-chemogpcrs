"""Bead od2f: the downstream RRA consumers must not let a zeroed/quarantined
signal vote.

rank_candidates.py (the canonical ranker) excludes zero-weight and
orthology-quarantined signals from the aggregation:

    excluded = excluded_signals_from_weights(SCORING_WEIGHTS) | QUARANTINED_AXES

but it runs the whole pipeline on import, so it cannot be imported. The six
downstream consumers -- select_probe_batch, rank_confidence, permutation_null,
shortlist_impact, compare_ranking_methods, audit_rra_correlation_sensitivity --
each read a ranked CSV produced under the SAME environment. Before this fix they
called ``build_ranklists_from_df`` WITHOUT ``excluded=``, so a signal the
production shortlist deliberately silenced (e.g. PURIFYING_WEIGHT=0, or the
orthology-quarantined axes) still voted at full strength in their analyses.

``rank_aggregation.production_excluded_signals`` rebuilds the identical set from
the environment; the consumers now thread it into every ranklist build. This
module checks that reconstruction and, for the two consumers named in the bead
(rank_confidence, permutation_null), that the excluded signal really is absent
from the ranklists they produce.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

# conftest.py puts scripts/ on sys.path
import rank_aggregation as ra
import permutation_null
import rank_confidence


# ---------------------------------------------------------------------------
# production_excluded_signals: faithful reconstruction of the canonical set
# ---------------------------------------------------------------------------

def test_key_parity_with_signal_spec() -> None:
    """The weight-env spec must cover exactly the 13 non-channel SIGNAL_SPEC
    signals, so a newly-added scored signal cannot silently escape the
    exclusion derivation."""
    spec_keys = {e[0] for e in ra.SIGNAL_SPEC if len(e) == 2}
    env_keys = {key for key, _, _ in ra._PRODUCTION_WEIGHT_ENV}
    assert spec_keys == env_keys, sorted(spec_keys ^ env_keys)


def test_zero_weight_axis_is_excluded_without_quarantine(monkeypatch) -> None:
    """With orthology trusted, only the zero-weight axes are excluded."""
    monkeypatch.setenv("ORTHOLOGY_SOURCE_TRUSTED", "1")
    monkeypatch.setenv("PURIFYING_WEIGHT", "0")
    monkeypatch.delenv("RANKAGG_EXCLUDED_SIGNALS", raising=False)
    assert ra.production_excluded_signals() == {"purifying"}


def test_nothing_excluded_when_all_weighted_and_trusted(monkeypatch) -> None:
    monkeypatch.setenv("ORTHOLOGY_SOURCE_TRUSTED", "1")
    monkeypatch.setenv("PURIFYING_WEIGHT", "2")
    monkeypatch.delenv("RANKAGG_EXCLUDED_SIGNALS", raising=False)
    assert ra.production_excluded_signals() == set()


def test_quarantine_adds_the_orthogroup_axes(monkeypatch) -> None:
    """Untrusted orthology quarantines the orthogroup-derived axes on top of
    any zero-weight exclusion."""
    monkeypatch.setenv("ORTHOLOGY_SOURCE_TRUSTED", "0")
    monkeypatch.setenv("PURIFYING_WEIGHT", "0")
    monkeypatch.delenv("RANKAGG_EXCLUDED_SIGNALS", raising=False)
    assert ra.production_excluded_signals() == {
        "purifying", "positive", "expansion", "og_confidence"}


def test_env_override_is_honoured(monkeypatch) -> None:
    """RANKAGG_EXCLUDED_SIGNALS still overrides the weight derivation, and the
    quarantine is unioned on top (matching rank_candidates.py)."""
    monkeypatch.setenv("ORTHOLOGY_SOURCE_TRUSTED", "1")
    monkeypatch.setenv("RANKAGG_EXCLUDED_SIGNALS", "synteny")
    assert ra.production_excluded_signals() == {"synteny"}


# ---------------------------------------------------------------------------
# the two consumers named in the bead must drop the excluded signal
# ---------------------------------------------------------------------------

def _df_with_purifying_vote() -> pd.DataFrame:
    """A df where `purifying` and `positive` WOULD vote if not excluded:
    varying score columns plus has_dnds_data=True. `phylo` is the kept control."""
    ids = [f"c{i}" for i in range(6)]
    n = len(ids)
    lin = [1.0 - i / n for i in range(n)]
    inv = [i / n for i in range(n)]
    return pd.DataFrame({
        "id": ids,
        "phylo_score": lin,
        "has_phylo_data": [True] * n,
        "purifying_score": inv,
        "positive_score": lin,
        "has_dnds_data": [True] * n,
    })


def _spy(real):
    captured = {}

    def spy(df, *args, **kwargs):
        result = real(df, *args, **kwargs)
        captured["excluded"] = kwargs.get("excluded")
        captured["ranklists"] = result
        return result

    return spy, captured


@pytest.fixture(autouse=True)
def _production_default_env(monkeypatch):
    """Pin the production default (orthology untrusted, purifying weight 0) so
    the two integration tests below run against a known exclusion set."""
    monkeypatch.setenv("ORTHOLOGY_SOURCE_TRUSTED", "0")
    monkeypatch.setenv("PURIFYING_WEIGHT", "0")
    monkeypatch.delenv("RANKAGG_EXCLUDED_SIGNALS", raising=False)


def test_fixture_would_vote_purifying_without_exclusion() -> None:
    """Precondition (must-accept side): the fixture DOES produce a purifying
    vote when nothing is excluded, so its absence below is caused by the fix,
    not by the fixture lacking the column."""
    lists = ra.build_ranklists_from_df(_df_with_purifying_vote())
    assert "purifying" in lists and "positive" in lists
    assert "phylo" in lists


def test_permutation_null_excludes_the_zeroed_signal(tmp_path, monkeypatch) -> None:
    real = ra.build_ranklists_from_df
    spy, captured = _spy(real)
    # permutation_null.main does `from rank_aggregation import ...` at call time,
    # so patching the rank_aggregation attribute is what the local import binds.
    monkeypatch.setattr(ra, "build_ranklists_from_df", spy)

    csv_path = tmp_path / "ranked.csv"
    _df_with_purifying_vote().to_csv(csv_path, index=False)
    out_path = tmp_path / "null.json"
    permutation_null.main([
        "--ranked-csv", str(csv_path), "--out", str(out_path),
        "--k", "2", "--n-perm", "20"])

    assert "purifying" in captured["excluded"], (
        "permutation_null must pass the production exclusion set")
    assert "purifying" not in captured["ranklists"], (
        "the zeroed/quarantined purifying signal still voted in the null")
    assert "positive" not in captured["ranklists"]
    assert "phylo" in captured["ranklists"], "the exclusion must be surgical"
    assert out_path.exists() and json.loads(out_path.read_text())["k"] == 2


def test_rank_confidence_excludes_the_zeroed_signal(tmp_path, monkeypatch) -> None:
    real = ra.build_ranklists_from_df
    spy, captured = _spy(real)
    # rank_confidence imports the name at module scope.
    monkeypatch.setattr(rank_confidence, "build_ranklists_from_df", spy)

    csv_path = tmp_path / "ranked.csv"
    _df_with_purifying_vote().to_csv(csv_path, index=False)
    rank_confidence.annotate_ranked_csv(str(csv_path), k=2, n_boot=20)

    assert "purifying" in captured["excluded"], (
        "rank_confidence must pass the production exclusion set")
    assert "purifying" not in captured["ranklists"], (
        "the zeroed/quarantined purifying signal still voted in the rank CIs")
    assert "positive" not in captured["ranklists"]
    assert "phylo" in captured["ranklists"], "the exclusion must be surgical"
