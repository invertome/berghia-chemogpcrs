"""A constant (zero-variance) signal must WARN but keep voting (bead zv1j).

A signal column whose non-null values are all equal cannot rank anything -- it
gives every candidate that votes on it the same normalized rank -- yet it still
casts a vote in the RRA fusion and so subtly distorts the aggregate. The
approved behaviour (WARN-AND-KEEP) is to emit a clear warning naming the
offending signal while leaving the ranklist exactly as before: visibility with
zero change to the computed ranking.

Fail-loud (reject/drop the signal) was deliberately rejected here because it
breaks the many existing constant-column fixtures; that migration is deferred.

The module's warning convention is ``print("WARNING: ...", file=sys.stderr)``
(see ``excluded_signals_from_weights``), so the warning is captured on stderr.
"""
from __future__ import annotations

import pandas as pd

import rank_aggregation as ra


def test_constant_signal_warns_and_still_votes(capsys):
    # 'phylo' is genuinely constant across all three candidates -> it cannot
    # discriminate, but WARN-AND-KEEP means it stays in the ranklists.
    df = pd.DataFrame(
        {
            "id": ["x", "y", "z"],
            "phylo_score": [0.5, 0.5, 0.5],   # constant -> must warn
            "positive_score": [0.9, 0.5, 0.1],  # varies -> must not warn
        }
    )
    rl = ra.build_ranklists_from_df(df)

    err = capsys.readouterr().err
    assert "phylo" in err, f"the constant signal was not named in a warning:\n{err}"
    assert "WARNING" in err

    # behaviour unchanged: the constant signal is still present and still holds
    # every candidate's value.
    assert "phylo" in rl
    assert set(rl["phylo"]) == {"x", "y", "z"}
    assert set(rl["phylo"].values()) == {0.5}


def test_varying_signal_does_not_warn(capsys):
    df = pd.DataFrame(
        {
            "id": ["x", "y", "z"],
            "phylo_score": [0.9, 0.5, 0.1],
            "positive_score": [0.2, 0.4, 0.6],
        }
    )
    rl = ra.build_ranklists_from_df(df)

    err = capsys.readouterr().err
    assert "constant" not in err.lower()
    assert "phylo" not in err and "positive" not in err
    # both signals vote, none flagged
    assert set(rl["phylo"]) == {"x", "y", "z"}


def test_single_candidate_frame_does_not_warn(capsys):
    # One candidate -> every signal is trivially "constant" but there is nothing
    # to rank, so warning would be noise. The guard is len(votes) > 1.
    df = pd.DataFrame(
        {
            "id": ["only"],
            "phylo_score": [0.5],
            "positive_score": [0.5],
        }
    )
    rl = ra.build_ranklists_from_df(df)

    err = capsys.readouterr().err
    assert err == "", f"a single-candidate frame must not warn:\n{err}"
    assert set(rl["phylo"]) == {"only"}


def test_entirely_nan_signal_does_not_warn(capsys):
    # A signal with no data at all is dropped (empty ranklist) and is a
    # different, already-handled case -- it must NOT trip the constant warning.
    df = pd.DataFrame(
        {
            "id": ["x", "y"],
            "phylo_score": [0.9, 0.1],
            "synteny_score": [float("nan"), float("nan")],
            "has_synteny_data": [True, True],
        }
    )
    rl = ra.build_ranklists_from_df(df)

    err = capsys.readouterr().err
    assert "synteny" not in err, f"an all-NaN signal must not warn:\n{err}"
    assert "synteny" not in rl
