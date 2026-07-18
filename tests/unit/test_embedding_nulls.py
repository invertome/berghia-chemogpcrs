"""Unit tests for embedding_nulls (A4 residue-shuffle embedding null, v4bs.5)."""
from __future__ import annotations

import collections

import numpy as np

from embedding_nulls import empirical_pvalue, main, novelty_collapse, shuffle_residues


def _rng(seed=0):
    return np.random.default_rng(seed)


def test_shuffle_residues_frac_zero_is_identity():
    seq = "MKTAYIAKQR"
    assert shuffle_residues(seq, 0.0, _rng()) == seq


def test_shuffle_residues_frac_one_preserves_multiset_and_shuffles():
    # Long, non-uniform sequence so a full permutation almost surely reorders it.
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
    out = shuffle_residues(seq, 1.0, _rng(7))
    assert len(out) == len(seq)
    assert collections.Counter(out) == collections.Counter(seq)
    # frac=1.0 permutes every position; on a long non-uniform seq it reorders.
    assert out != seq


def test_shuffle_residues_is_seed_deterministic_and_uses_rng():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
    # Same seed + same args -> identical output.
    assert shuffle_residues(seq, 0.5, _rng(3)) == shuffle_residues(seq, 0.5, _rng(3))
    # Randomness flows through the passed rng: two seeds can differ.
    assert shuffle_residues(seq, 0.5, _rng(1)) != shuffle_residues(seq, 0.5, _rng(2))


def test_shuffle_residues_intermediate_frac_preserves_multiset():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGT"
    out = shuffle_residues(seq, 0.4, _rng(11))
    assert len(out) == len(seq)
    assert collections.Counter(out) == collections.Counter(seq)


def test_empirical_pvalue_observed_above_all_nulls():
    # Right-tail: observed strictly greater than every null -> 1/(m+1).
    nulls = [0.1, 0.2, 0.3, 0.4]
    assert empirical_pvalue(0.9, nulls) == 1.0 / (len(nulls) + 1)


def test_empirical_pvalue_observed_below_all_nulls():
    # observed <= every null -> all m nulls counted -> (m+1)/(m+1) == 1.0.
    nulls = [0.5, 0.6, 0.7]
    assert empirical_pvalue(0.1, nulls) == 1.0


def test_empirical_pvalue_hand_computed_case():
    # nulls >= observed(0.5): {0.5, 0.6, 0.9} -> b=3, m=5 -> (3+1)/(5+1).
    nulls = [0.1, 0.4, 0.5, 0.6, 0.9]
    assert empirical_pvalue(0.5, nulls) == (3 + 1) / (5 + 1)


def test_novelty_collapse_observed_far_above_nulls():
    nulls = [0.10, 0.12, 0.09, 0.11, 0.10]
    d = novelty_collapse(1.5, nulls)
    assert d["excess"] > 5.0           # many std above the null mean
    assert d["collapsed"] == False     # signal stands above shuffle noise
    assert d["null_mean"] == np.mean(nulls)


def test_novelty_collapse_observed_within_null_cloud():
    nulls = [0.40, 0.55, 0.48, 0.60, 0.45]
    d = novelty_collapse(0.30, nulls)  # below the null mean
    assert d["collapsed"] == True      # shuffled null reproduces/exceeds observed
    assert d["excess"] < 0.0


def test_novelty_collapse_zero_null_std_does_not_raise():
    nulls = [0.5, 0.5, 0.5]
    # observed below the (degenerate) null mean -> no positive excess, collapsed.
    d_below = novelty_collapse(0.2, nulls)
    assert d_below["null_std"] == 0.0
    assert d_below["excess"] == 0.0
    assert d_below["collapsed"] == True
    # observed above a zero-variance null -> infinite excess, not collapsed.
    d_above = novelty_collapse(0.8, nulls)
    assert d_above["excess"] == float("inf")
    assert d_above["collapsed"] == False


def _write_tsv(path, rows):
    path.write_text(
        "candidate_id\tnovelty\n"
        + "".join(f"{cid}\t{nov}\n" for cid, nov in rows)
    )


def test_main_writes_pvalue_and_collapse_tsv(tmp_path):
    import pandas as pd

    obs = tmp_path / "observed.tsv"
    _write_tsv(obs, [("c1", 0.9), ("c2", 0.1)])
    r1 = tmp_path / "null_r1.tsv"
    _write_tsv(r1, [("c1", 0.1), ("c2", 0.4)])
    r2 = tmp_path / "null_r2.tsv"
    _write_tsv(r2, [("c1", 0.2), ("c2", 0.5)])
    r3 = tmp_path / "null_r3.tsv"
    _write_tsv(r3, [("c1", 0.3), ("c2", 0.6)])
    out = tmp_path / "out.tsv"

    main([
        "--observed-tsv", str(obs),
        "--null-tsv", str(r1), str(r2), str(r3),
        "--out", str(out),
    ])

    df = pd.read_csv(out, sep="\t").set_index("id")
    assert list(df.columns) == ["novelty", "null_mean", "p_value", "collapsed"]
    # c1: nulls [0.1,0.2,0.3], observed 0.9 -> b=0 -> p=1/4; mean 0.2; not collapsed.
    assert df.loc["c1", "p_value"] == 1.0 / 4.0
    assert df.loc["c1", "null_mean"] == 0.2
    assert not bool(df.loc["c1", "collapsed"])
    # c2: nulls [0.4,0.5,0.6], observed 0.1 -> b=3 -> p=1.0; mean 0.5; collapsed.
    assert df.loc["c2", "p_value"] == 1.0
    assert df.loc["c2", "null_mean"] == 0.5
    assert bool(df.loc["c2", "collapsed"])
