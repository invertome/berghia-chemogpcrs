"""Tests for the Mahalanobis scorer added to scripts/embedding_evidence.py.

Phase-0 optimization (bead cw3): the exclusion/novelty signal moves from raw
cosine (anisotropy-fragile) to tied-covariance Mahalanobis distance
(preprocessing-invariant), with multi-prototype family centroids. Emits
`S_known` (nearest-prototype family assignment) and `S_novel` (plain
Mahalanobis min-distance to the known families — the positive novelty axis).
"""
from __future__ import annotations

import numpy as np

from embedding_evidence import (
    family_prototypes,
    mahalanobis_assignment,
    mahalanobis_channel,
    mahalanobis_sq,
    novelty_score,
    shrinkage_precision,
)


def test_mahalanobis_sq_is_the_precision_quadratic_form():
    # (x-mu)^T P (x-mu) with a diagonal precision = weighted squared distance.
    precision = np.array([[2.0, 0.0], [0.0, 8.0]])
    x = np.array([1.0, 1.0])
    mu = np.array([0.0, 0.0])
    # 2*1^2 + 8*1^2 = 10
    assert np.isclose(mahalanobis_sq(x, mu, precision), 10.0)


def test_mahalanobis_sq_zero_at_the_centroid():
    precision = np.array([[3.0, 1.0], [1.0, 5.0]])
    mu = np.array([2.0, -1.0])
    assert np.isclose(mahalanobis_sq(mu, mu, precision), 0.0)


def test_shrinkage_precision_normalizes_anisotropic_spread():
    # A cloud wide along dim 0 (sd~10), narrow along dim 1 (sd~1).
    rng = np.random.RandomState(0)
    X = rng.randn(4000, 2) * np.array([10.0, 1.0])
    precision = shrinkage_precision(X)
    origin = np.zeros(2)
    # A 10-unit step along the wide axis and a 1-unit step along the narrow axis
    # are both ~1 sd. Mahalanobis rescales each axis by its spread, so both land
    # at ~1 (within a small factor) — whereas raw Euclidean would rate them 100x
    # apart (100 vs 1). This is the whole point: distance in the cloud's own units.
    d_wide = mahalanobis_sq(np.array([10.0, 0.0]), origin, precision)
    d_narrow = mahalanobis_sq(np.array([0.0, 1.0]), origin, precision)
    assert 0.5 < d_wide / d_narrow < 2.0
    maha_ratio = max(d_wide, d_narrow) / min(d_wide, d_narrow)
    assert maha_ratio < 0.1 * (100.0 / 1.0)  # far more balanced than Euclidean


def test_shrinkage_precision_is_usable_when_n_less_than_d():
    # 5 samples, 10 dims: the raw sample covariance is singular; shrinkage keeps
    # the precision finite and positive-definite (our real case: 206 refs, 960-d).
    rng = np.random.RandomState(1)
    X = rng.randn(5, 10)
    precision = shrinkage_precision(X)
    assert precision.shape == (10, 10)
    assert np.all(np.isfinite(precision))
    assert mahalanobis_sq(np.ones(10), np.zeros(10), precision) > 0


def test_family_prototypes_captures_a_bimodal_family():
    # A family with two well-separated modes (+x and -x) -> prototypes reach both,
    # unlike a single mean which would land near the origin between the modes.
    rng = np.random.RandomState(0)
    mode_a = rng.randn(20, 4) + np.array([10.0, 0, 0, 0])
    mode_b = rng.randn(20, 4) + np.array([-10.0, 0, 0, 0])
    emb = {f"a{i}": v for i, v in enumerate(mode_a)}
    emb.update({f"b{i}": v for i, v in enumerate(mode_b)})
    labels = {k: "fam" for k in emb}
    protos = family_prototypes(emb, labels, k=3)
    P = protos["fam"]
    assert P.shape[1] == 4
    assert P[:, 0].max() > 5.0 and P[:, 0].min() < -5.0  # both modes represented


def test_family_prototypes_singleton_family_falls_back_to_mean():
    emb = {"x": np.array([1.0, 2.0, 3.0])}
    labels = {"x": "solo"}
    protos = family_prototypes(emb, labels, k=3)
    assert protos["solo"].shape[0] == 1          # too few to cluster -> one proto
    assert np.allclose(protos["solo"][0], [1.0, 2.0, 3.0])


def test_mahalanobis_assignment_picks_the_nearest_family():
    prototypes = {"A": np.array([[0.0, 0.0]]), "B": np.array([[10.0, 10.0]])}
    precision = np.eye(2)
    family, dist = mahalanobis_assignment(np.array([0.5, 0.5]), prototypes, precision)
    assert family == "A"
    assert dist >= 0


def test_mahalanobis_assignment_uses_nearest_prototype_within_a_family():
    # Family A is bimodal (2 prototypes); a query near A's SECOND prototype must
    # still assign to A, not to the single-prototype family B.
    prototypes = {
        "A": np.array([[0.0, 0.0], [10.0, 0.0]]),
        "B": np.array([[0.0, 10.0]]),
    }
    precision = np.eye(2)
    family, _ = mahalanobis_assignment(np.array([9.5, 0.0]), prototypes, precision)
    assert family == "A"


def test_novelty_score_plain_ranks_a_divergent_point_higher():
    # S_novel (plain): far from the only known family => more novel.
    prototypes = {"A": np.array([[0.0, 0.0]])}
    precision = np.eye(2)
    near = novelty_score(np.array([0.3, 0.0]), prototypes, precision)
    far = novelty_score(np.array([5.0, 0.0]), prototypes, precision)
    assert far > near


def test_novelty_score_relative_subtracts_the_background():
    # Relative Mahalanobis (config option): min-distance minus distance to the
    # background centroid (Ren et al. 2021) — the conservative variant.
    prototypes = {"A": np.array([[0.0, 0.0]])}
    precision = np.eye(2)
    bg = np.array([2.0, 0.0])
    x = np.array([5.0, 0.0])
    plain = novelty_score(x, prototypes, precision)
    rel = novelty_score(x, prototypes, precision, background_mean=bg)
    assert np.isclose(rel, plain - mahalanobis_sq(x, bg, precision))
    assert rel < plain


def _two_family_reference(rng, dim=8):
    ref, labels = {}, {}
    for i in range(12):
        ref[f"A{i}"] = rng.randn(dim) + np.eye(dim)[0] * 5.0; labels[f"A{i}"] = "A"
        ref[f"B{i}"] = rng.randn(dim) - np.eye(dim)[0] * 5.0; labels[f"B{i}"] = "B"
    return ref, labels


def test_mahalanobis_channel_emits_novelty_and_family_per_candidate():
    rng = np.random.RandomState(0)
    ref, labels = _two_family_reference(rng)
    candidates = {
        "like_A": np.eye(8)[0] * 5.0,        # sits in family A
        "novel": np.eye(8)[1] * 50.0,        # far from both families
    }
    channel = mahalanobis_channel(candidates, ref, labels, k=3)
    assert set(channel) == {"like_A", "novel"}
    # the divergent candidate scores higher novelty (the positive axis)
    assert channel["novel"]["emb_novelty"] > channel["like_A"]["emb_novelty"]
    # the known-like candidate is annotated with its family
    assert channel["like_A"]["emb_nonchemo_family"] == "A"


def test_mahalanobis_channel_key_contract_never_a_chemoreceptor_score():
    rng = np.random.RandomState(1)
    ref, labels = _two_family_reference(rng)
    channel = mahalanobis_channel({"c": np.eye(8)[0]}, ref, labels, k=3)
    entry = channel["c"]
    assert set(entry) == {
        "emb_nonchemo_family", "emb_novelty", "has_emb_data", "emb_leakage_flag"
    }
    assert entry["has_emb_data"] is True
    assert entry["emb_leakage_flag"] is True
    assert not any("chemoreceptor" in key.lower() for key in entry)


def test_mahalanobis_channel_relative_flag_changes_the_novelty():
    rng = np.random.RandomState(2)
    ref, labels = _two_family_reference(rng)
    cand = {"c": np.eye(8)[1] * 20.0}
    plain = mahalanobis_channel(cand, ref, labels, relative=False)["c"]["emb_novelty"]
    rel = mahalanobis_channel(cand, ref, labels, relative=True)["c"]["emb_novelty"]
    assert plain != rel
