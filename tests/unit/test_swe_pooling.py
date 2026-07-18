"""Unit tests for swe_pooling (A5, epic v4bs.6) — Sliced-Wasserstein pooling.

Each test encodes one acceptance criterion from the A5 spec: fixed output
shape, permutation invariance (SET function), seed determinism,
length-robustness ordering, unit-norm slice directions, and the D=1 /
single-residue edge cases. Imports are by bare module name (conftest puts
scripts/ on sys.path).
"""
from __future__ import annotations

import numpy as np

from swe_pooling import random_unit_directions, swe_pool


def test_output_shape_is_slices_times_quantiles():
    rng = np.random.default_rng(0)
    vecs = rng.standard_normal((40, 12))
    out = swe_pool(vecs, n_slices=16, n_quantiles=5, seed=0)
    assert out.shape == (16 * 5,)


def test_permutation_invariance_is_a_set_function():
    rng = np.random.default_rng(1)
    vecs = rng.standard_normal((30, 7))
    perm = rng.permutation(vecs.shape[0])
    shuffled = vecs[perm]
    a = swe_pool(vecs, n_slices=8, n_quantiles=4, seed=5)
    b = swe_pool(shuffled, n_slices=8, n_quantiles=4, seed=5)
    assert np.array_equal(a, b)


def test_determinism_same_seed_identical_diff_seed_differs():
    rng = np.random.default_rng(2)
    vecs = rng.standard_normal((50, 9))
    a = swe_pool(vecs, n_slices=16, n_quantiles=6, seed=11)
    a_again = swe_pool(vecs, n_slices=16, n_quantiles=6, seed=11)
    b = swe_pool(vecs, n_slices=16, n_quantiles=6, seed=12)
    assert np.array_equal(a, a_again)
    assert not np.allclose(a, b)


def test_length_robustness_same_dist_closer_than_shifted_dist():
    # Same underlying Gaussian, different L (50 vs 500) -> close pooled vectors;
    # a mean-shifted distribution -> clearly larger pooled difference.
    rng = np.random.default_rng(7)
    dim = 6
    big = rng.standard_normal((500, dim))
    small = big[:50]
    shifted = big + 3.0  # different distribution (shifted mean)

    p_small = swe_pool(small, n_slices=48, n_quantiles=8, seed=0)
    p_big = swe_pool(big, n_slices=48, n_quantiles=8, seed=0)
    p_shift = swe_pool(shifted, n_slices=48, n_quantiles=8, seed=0)

    same_dist_diff = float(np.mean(np.abs(p_small - p_big)))
    diff_dist_diff = float(np.mean(np.abs(p_big - p_shift)))

    # length change moves the pooled vector far less than a real shift does
    assert same_dist_diff < diff_dist_diff
    assert same_dist_diff < 0.5 * diff_dist_diff


def test_slice_directions_are_unit_vectors():
    rng = np.random.default_rng(0)
    dirs = random_unit_directions(32, 7, rng)
    assert dirs.shape == (32, 7)
    norms = np.linalg.norm(dirs, axis=1)
    assert np.allclose(norms, 1.0)


def test_dim_one_input_does_not_raise():
    rng = np.random.default_rng(3)
    vecs = rng.standard_normal((25, 1))
    out = swe_pool(vecs, n_slices=8, n_quantiles=4, seed=0)
    assert out.shape == (32,)
    assert np.all(np.isfinite(out))


def test_single_residue_input_quantiles_are_that_value():
    # One residue: every quantile of a singleton equals that value, so for
    # D=1 each slice's block is all +/- the projected scalar; no raise.
    vecs = np.array([[0.5]])
    out = swe_pool(vecs, n_slices=4, n_quantiles=3, seed=0)
    assert out.shape == (12,)
    # each length-n_quantiles block is constant (all quantiles equal)
    blocks = out.reshape(4, 3)
    assert np.allclose(blocks, blocks[:, :1])
    assert np.allclose(np.abs(out), 0.5)
