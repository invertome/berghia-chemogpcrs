"""``RANKAGG_EXCLUDED_SIGNALS=""`` must not disarm the exclusion derivation.

The defect (bead wtwi)
----------------------
``excluded_signals_from_weights`` branched on ``os.getenv(...) is not None``::

    override = os.getenv(env_var)
    if override is not None:
        return {n.strip() for n in override.split(",") if n.strip()}
    return {name for name, w in weights.items() if float(w) == 0.0}

``export RANKAGG_EXCLUDED_SIGNALS=`` sets the variable to ``""``, which is not
None, so the override branch is taken and returns the EMPTY SET. The
weight-derived exclusion is skipped entirely and every zero-weighted signal
votes at full strength -- precisely the failure the function exists to prevent.

The production case is ``PURIFYING_WEIGHT=0``. Under the default
``RANK_METHOD=rankagg`` a stray empty assignment silently restores whole-gene
purifying selection as a full-strength voter, which is the wrong signal for
chemoreceptor discovery and reorders the shortlist.

What makes it a trapdoor rather than an ordinary bug: the INTENDED value
("exclude nothing") and the ACCIDENTAL value (an empty export, an unset shell
variable interpolated into a wrapper, a CI default that resolved to nothing)
are the SAME STRING. No amount of care at the call site can tell them apart,
so the ambiguity has to be removed from the contract itself.

The contract now
----------------
  * unset                -> derive exclusions from the weights (safe default)
  * ``""`` / whitespace  -> NO intent expressed; fall back to the derivation,
                            and say so on stderr
  * ``"none"``           -> explicitly exclude nothing
  * ``"a,b"``            -> exclude exactly those

This is the cross-cutting invariant applied to configuration: an empty string
is not a measurement of intent, so it must not be turned into a consequential
decision that happens to look deliberate.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import rank_aggregation as ra  # noqa: E402

ENV = "RANKAGG_EXCLUDED_SIGNALS"

# The production weighting this protects: PURIFYING_WEIGHT=0 means "retire this
# axis", not "weight it lightly".
PROD_WEIGHTS = {"phylo": 2.0, "purifying": 0.0, "positive": 1.5, "synteny": 1.0}


# --------------------------------------------------------------------------
# the trapdoor
# --------------------------------------------------------------------------
def test_empty_env_does_not_silently_re_enable_a_zero_weighted_signal(monkeypatch):
    """The defect, stated directly: an empty export must not un-exclude."""
    monkeypatch.setenv(ENV, "")
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == {"purifying"}, (
        "an accidental `export RANKAGG_EXCLUDED_SIGNALS=` disarmed the "
        "weight-derived exclusion; purifying would vote at full strength"
    )


def test_whitespace_only_env_is_treated_the_same_as_empty(monkeypatch):
    """`export VAR=' '` is the same accident with a space in it."""
    monkeypatch.setenv(ENV, "   ")
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == {"purifying"}


def test_comma_only_env_is_treated_the_same_as_empty(monkeypatch):
    """`export VAR=,` -- a list that lost its items -- is also no intent."""
    monkeypatch.setenv(ENV, ",")
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == {"purifying"}


def test_empty_env_is_announced_not_silent(monkeypatch, capsys):
    """Falling back is a decision, so it must be visible."""
    monkeypatch.setenv(ENV, "")
    ra.excluded_signals_from_weights(PROD_WEIGHTS)
    err = capsys.readouterr().err
    assert ENV in err and "purifying" in err, (
        f"the fallback was silent; operators cannot see it happened:\n{err}"
    )


# --------------------------------------------------------------------------
# "exclude nothing" is still reachable -- but only on purpose
# --------------------------------------------------------------------------
def test_explicit_none_token_excludes_nothing(monkeypatch):
    monkeypatch.setenv(ENV, "none")
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == set()


def test_none_token_is_case_and_space_insensitive(monkeypatch):
    monkeypatch.setenv(ENV, "  NONE  ")
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == set()


def test_none_mixed_with_real_signals_is_a_hard_error(monkeypatch):
    """"exclude nothing, and also exclude synteny" is contradictory."""
    monkeypatch.setenv(ENV, "none,synteny")
    with pytest.raises(ValueError, match="none"):
        ra.excluded_signals_from_weights(PROD_WEIGHTS)


# --------------------------------------------------------------------------
# everything else is unchanged
# --------------------------------------------------------------------------
def test_unset_still_derives_from_weights(monkeypatch):
    monkeypatch.delenv(ENV, raising=False)
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == {"purifying"}


def test_unset_with_no_zero_weights_excludes_nothing(monkeypatch):
    monkeypatch.delenv(ENV, raising=False)
    assert ra.excluded_signals_from_weights({"phylo": 2.0, "synteny": 1.0}) == set()


def test_explicit_list_still_overrides_verbatim(monkeypatch):
    """An explicit list replaces the derivation, it does not union with it."""
    monkeypatch.setenv(ENV, "synteny")
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == {"synteny"}


def test_explicit_list_tolerates_spacing(monkeypatch):
    monkeypatch.setenv(ENV, " synteny , positive ")
    assert ra.excluded_signals_from_weights(PROD_WEIGHTS) == {"synteny", "positive"}


# --------------------------------------------------------------------------
# the exclusion actually reaches the aggregator
# --------------------------------------------------------------------------
def test_excluded_signal_casts_no_vote_in_build_ranklists(monkeypatch):
    """End-to-end: the derived exclusion must remove the ranklist, not just the name.

    Grounded in the column shape rank_candidates.py actually writes -- the raw
    ``<signal>_score`` columns plus their ``has_*_data`` flags -- rather than a
    schema invented for the test.
    """
    pd = pytest.importorskip("pandas")
    monkeypatch.setenv(ENV, "")
    df = pd.DataFrame({
        "id": ["cand_1", "cand_2", "cand_3"],
        "phylo_score": [0.9, 0.5, 0.1],
        "has_phylo_data": [True, True, True],
        "purifying_score": [0.1, 0.5, 0.9],
        "has_purifying_data": [True, True, True],
    })
    excluded = ra.excluded_signals_from_weights(PROD_WEIGHTS)
    ranklists = ra.build_ranklists_from_df(df, excluded=excluded)

    assert "purifying" not in ranklists, (
        "purifying still votes despite PURIFYING_WEIGHT=0 -- the empty-string "
        "trapdoor reaches all the way into the aggregation"
    )
    assert "phylo" in ranklists, "the exclusion must be surgical, not total"
