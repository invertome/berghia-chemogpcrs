"""Tests for scripts/embedding_evidence.py.

Rank Task 5 — ESM-C embedding-evidence channel. Council HARD RULE: PLM
embeddings are used ONLY for exclusion (distance to known NON-chemoreceptor
family centroids) and recall (broad class-A-GPCR-ness), NEVER as a positive
"similarity to a known chemoreceptor" score — a real lineage-specific
expansion is *supposed* to look unlike characterized genes, so positive
resemblance is anti-correlated with the target. Every score also carries an
explicit homology-leakage-uncertainty flag (Berghia is ~absent from PLM
pretraining corpora).

Coverage:
    - load_embeddings: reads a {id: vector} dict from .npz; missing file -> {}
    - centroid: L2-normalized mean over present ids; absent ids ignored;
      empty selection -> zero vector
    - nearest_family: (family, cosine_sim) of the max-cosine centroid,
      deterministic tie-break, safe on zero vectors / empty centroid dicts
    - embedding_channel: emits exclusion (emb_nonchemo_*) + recall
      (emb_classA_sim) signals + has_emb_data/emb_leakage_flag, and is
      guarded to NEVER emit a chemoreceptor-similarity key
"""
from __future__ import annotations

import math

import numpy as np

from embedding_evidence import (
    centroid,
    embedding_channel,
    load_embeddings,
    nearest_family,
)


# ---- load_embeddings -----------------------------------------------------

def test_load_embeddings_missing_file_returns_empty_dict(tmp_path):
    assert load_embeddings(str(tmp_path / "missing.npz")) == {}


def test_load_embeddings_reads_per_id_arrays(tmp_path):
    npz_path = tmp_path / "emb.npz"
    np.savez(npz_path, gene_1=np.array([1.0, 2.0]), gene_2=np.array([3.0, 4.0]))
    result = load_embeddings(str(npz_path))
    assert set(result) == {"gene_1", "gene_2"}
    assert np.array_equal(result["gene_1"], np.array([1.0, 2.0]))
    assert np.array_equal(result["gene_2"], np.array([3.0, 4.0]))


# ---- centroid --------------------------------------------------------------

def test_centroid_l2_normalizes_the_mean():
    embeddings = {"a": np.array([3.0, 4.0]), "b": np.array([3.0, 4.0])}
    result = centroid(embeddings, ["a", "b"])
    assert np.allclose(result, [0.6, 0.8])
    assert math.isclose(float(np.linalg.norm(result)), 1.0, rel_tol=1e-9)


def test_centroid_ignores_ids_absent_from_embeddings():
    embeddings = {"a": np.array([1.0, 0.0]), "b": np.array([0.0, 1.0])}
    # "c" is requested but has no embedding -> must not affect the mean
    result = centroid(embeddings, ["a", "b", "c"])
    assert np.allclose(result, [1 / math.sqrt(2), 1 / math.sqrt(2)])


def test_centroid_empty_selection_returns_zero_vector_with_inferred_dim():
    embeddings = {"a": np.array([1.0, 2.0, 3.0])}
    result = centroid(embeddings, [])
    assert np.array_equal(result, np.zeros(3))


def test_centroid_fully_empty_returns_zero_length_vector():
    result = centroid({}, [])
    assert result.size == 0


# ---- nearest_family ---------------------------------------------------------

def test_nearest_family_picks_higher_cosine():
    vec = np.array([1.0, 0.0])
    family_centroids = {
        "bioamine": np.array([0.0, 1.0]),  # orthogonal -> cos = 0
        "peptide": np.array([1.0, 0.0]),   # identical direction -> cos = 1
    }
    family, sim = nearest_family(vec, family_centroids)
    assert family == "peptide"
    assert math.isclose(sim, 1.0, rel_tol=1e-9)


def test_nearest_family_ties_break_alphabetically_for_determinism():
    vec = np.array([1.0, 1.0])
    family_centroids = {
        "zeta": np.array([2.0, 2.0]),   # same direction as vec -> cos = 1
        "alpha": np.array([1.0, 1.0]),  # same direction as vec -> cos = 1
    }
    family, sim = nearest_family(vec, family_centroids)
    assert family == "alpha"
    assert math.isclose(sim, 1.0, rel_tol=1e-9)


def test_nearest_family_empty_centroids_returns_none_and_zero():
    family, sim = nearest_family(np.array([1.0, 0.0]), {})
    assert family is None
    assert sim == 0.0


def test_nearest_family_zero_vector_returns_zero_similarity_not_nan():
    vec = np.zeros(2)
    family_centroids = {"bioamine": np.array([1.0, 0.0])}
    family, sim = nearest_family(vec, family_centroids)
    assert family == "bioamine"
    assert sim == 0.0
    assert not math.isnan(sim)


# ---- embedding_channel -------------------------------------------------------

_EXPECTED_KEYS = {
    "emb_nonchemo_sim",
    "emb_nonchemo_family",
    "emb_classA_sim",
    "has_emb_data",
    "emb_leakage_flag",
}


def test_embedding_channel_has_required_keys():
    embeddings = {"cand_1": np.array([1.0, 0.0])}
    nonchemo_centroids = {
        "bioamine": np.array([1.0, 0.0]),
        "peptide": np.array([0.0, 1.0]),
    }
    classA_centroid = np.array([1.0, 0.0])
    result = embedding_channel(embeddings, nonchemo_centroids, classA_centroid)
    assert set(result["cand_1"]) == _EXPECTED_KEYS


def test_embedding_channel_picks_nearest_nonchemo_family():
    embeddings = {"cand_1": np.array([0.0, 1.0])}
    nonchemo_centroids = {
        "bioamine": np.array([1.0, 0.0]),
        "peptide": np.array([0.0, 1.0]),
    }
    classA_centroid = np.array([1.0, 1.0])
    result = embedding_channel(embeddings, nonchemo_centroids, classA_centroid)
    assert result["cand_1"]["emb_nonchemo_family"] == "peptide"
    assert math.isclose(result["cand_1"]["emb_nonchemo_sim"], 1.0, rel_tol=1e-9)


def test_embedding_channel_classA_sim_matches_manual_cosine():
    embeddings = {"cand_1": np.array([1.0, 1.0])}
    nonchemo_centroids = {"bioamine": np.array([1.0, 0.0])}
    classA_centroid = np.array([1.0, 1.0])
    result = embedding_channel(embeddings, nonchemo_centroids, classA_centroid)
    assert math.isclose(result["cand_1"]["emb_classA_sim"], 1.0, rel_tol=1e-9)


def test_embedding_channel_has_emb_data_and_leakage_flag_are_true():
    embeddings = {"cand_1": np.array([1.0, 0.0])}
    result = embedding_channel(
        embeddings, {"bioamine": np.array([1.0, 0.0])}, np.array([1.0, 0.0])
    )
    assert result["cand_1"]["has_emb_data"] is True
    assert result["cand_1"]["emb_leakage_flag"] is True


def test_embedding_channel_similarity_values_are_python_floats():
    embeddings = {"cand_1": np.array([1.0, 0.0])}
    result = embedding_channel(
        embeddings, {"bioamine": np.array([1.0, 0.0])}, np.array([1.0, 0.0])
    )
    assert isinstance(result["cand_1"]["emb_nonchemo_sim"], float)
    assert isinstance(result["cand_1"]["emb_classA_sim"], float)


def test_embedding_channel_absent_candidates_are_simply_absent():
    """A candidate with no embedding gets no entry at all — downstream
    joins are documented to treat a missing id as has_emb_data=False;
    this pure function never fabricates a score for data it doesn't have."""
    embeddings = {"cand_1": np.array([1.0, 0.0])}
    result = embedding_channel(
        embeddings, {"bioamine": np.array([1.0, 0.0])}, np.array([1.0, 0.0])
    )
    assert "cand_2" not in result
    assert set(result) == {"cand_1"}


def test_embedding_channel_never_produces_a_chemoreceptor_similarity_key():
    """HARD RULE guard (council review): positive similarity-to-known-
    chemoreceptor is banned outright. No emitted key may reference a
    chemoreceptor similarity, for any candidate."""
    embeddings = {
        "cand_1": np.array([1.0, 0.0]),
        "cand_2": np.array([0.0, 1.0]),
    }
    nonchemo_centroids = {
        "bioamine": np.array([1.0, 0.0]),
        "peptide": np.array([0.0, 1.0]),
    }
    classA_centroid = np.array([0.5, 0.5])
    result = embedding_channel(embeddings, nonchemo_centroids, classA_centroid)
    assert result  # sanity: the guard below must not vacuously pass
    for candidate_id, fields in result.items():
        for key in fields:
            assert "chemoreceptor_sim" not in key.lower(), (
                f"banned positive-similarity key {key!r} found for {candidate_id}"
            )
        # belt-and-braces: the key set is EXACTLY the documented set, so no
        # stray field (chemoreceptor-similarity or otherwise) can sneak in
        assert set(fields) == _EXPECTED_KEYS
