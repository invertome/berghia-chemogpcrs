#!/usr/bin/env python3
"""embedding_evidence.py — ESM-C embedding-evidence channel for candidate ranking.

HARD RULE (council review, bead berghia-chemogpcrs-875 / Rank Task 5): protein
language model (PLM) embeddings are used ONLY for two council-sanctioned
purposes:

    1. EXCLUSION — cosine similarity to centroids of known NON-chemoreceptor
       GPCR families (bioamine, peptide, opsin, ...). A confident match
       corroborates the existing classifier's non-chemoreceptor call.
    2. RECALL — cosine similarity to a broad "class-A GPCR-ness" centroid,
       used only to confirm a candidate looks like a class-A GPCR at all.

This module MUST NEVER emit a positive "similarity to a known chemoreceptor"
score. A real lineage-specific-expansion (LSE) chemoreceptor in Berghia is
*expected* to look unlike any characterized gene — that is the definition of
a novel expansion — so a positive resemblance-to-known-chemoreceptors signal
is anti-correlated with the very candidates the user wants to surface. There
is deliberately no "chemoreceptor centroid" argument anywhere in this module,
and `embedding_channel` never emits a key referencing one.

Every score this module emits also carries an explicit `emb_leakage_flag`:
Berghia (and nudibranchs generally) are ~absent from PLM pretraining
corpora, so any embedding-based judgment about the Berghia proteome carries
unusually high homology-leakage uncertainty relative to, say, a mammalian
protein. The flag is a blanket caveat attached to every scored candidate,
not a per-sequence computed value.

This module is a PURE scorer: it consumes precomputed embeddings (a `.npz`
produced on Unity by `scripts/unity/run_esmc_embeddings.sh`) plus reference
family centroids. It never imports torch or esm and never runs model
inference — that only happens in the Unity wrapper script, so these
functions stay importable and unit-testable without a GPU or those deps.

Typical usage (library):
    embeddings = load_embeddings("candidates_esmc300m.npz")
    reference = load_embeddings("reference_esmc300m.npz")
    nonchemo_centroids = {
        "bioamine": centroid(reference, bioamine_reference_ids),
        "peptide": centroid(reference, peptide_reference_ids),
        ...
    }
    classA_centroid = centroid(reference, classA_reference_ids)
    channel = embedding_channel(embeddings, nonchemo_centroids, classA_centroid)
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

import numpy as np


def load_embeddings(npz_path: str) -> Dict[str, np.ndarray]:
    """Load a `{id: vector}` dict from an .npz produced by the ESM-C wrapper.

    Returns {} if the file doesn't exist (e.g. the embedding step hasn't run
    yet, or ran for a different candidate set) — callers treat that as "no
    embedding evidence available" and set has_emb_data=False, never a crash.
    """
    if not os.path.exists(npz_path):
        return {}
    with np.load(npz_path) as data:
        return {key: data[key] for key in data.files}


def centroid(embeddings: Dict[str, np.ndarray], ids: List[str]) -> np.ndarray:
    """L2-normalized mean vector over the `ids` present in `embeddings`.

    ids absent from `embeddings` are silently ignored (a reference id list
    may legitimately span a bigger dataset than what actually got embedded).
    If none of `ids` resolve to an embedding, returns a zero vector — its
    dimensionality is inferred from any embedding already in `embeddings`,
    or length 0 if `embeddings` itself is empty.
    """
    vectors = [np.asarray(embeddings[i], dtype=float) for i in ids if i in embeddings]
    if not vectors:
        dim = len(next(iter(embeddings.values()))) if embeddings else 0
        return np.zeros(dim, dtype=float)
    mean_vec = np.mean(vectors, axis=0)
    norm = np.linalg.norm(mean_vec)
    return mean_vec if norm == 0 else mean_vec / norm


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine similarity; 0.0 (never NaN) when either vector has zero norm."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return float(np.dot(a, b) / (norm_a * norm_b))


def nearest_family(
    vec: np.ndarray, family_centroids: Dict[str, np.ndarray]
) -> Tuple[Optional[str], float]:
    """(family, cosine_sim) of the max-cosine reference family centroid.

    `family_centroids` maps a family label (e.g. "bioamine", "peptide") to
    its reference centroid vector — NEVER a chemoreceptor centroid; see the
    module-level HARD RULE. Ties are broken alphabetically by family name so
    the result is deterministic. An empty `family_centroids` has nothing to
    compare against, so it returns (None, 0.0) rather than raising.
    """
    if not family_centroids:
        return None, 0.0
    families = sorted(family_centroids)
    sims = {f: _cosine_similarity(vec, family_centroids[f]) for f in families}
    best_family = max(families, key=lambda f: sims[f])
    return best_family, sims[best_family]


def embedding_channel(
    embeddings: Dict[str, np.ndarray],
    nonchemo_centroids: Dict[str, np.ndarray],
    classA_centroid: np.ndarray,
) -> Dict[str, Dict[str, object]]:
    """Per-candidate exclusion + recall signals. NEVER a chemoreceptor-similarity score.

    Only ids present in `embeddings` get an entry — a candidate with no
    embedding is simply absent from the returned dict. Downstream joins
    (rank_candidates.py) are documented to treat a missing id as
    has_emb_data=False; this function never fabricates a score for data it
    doesn't have.

    Returned per-id dict has EXACTLY these keys:
        emb_nonchemo_sim    - cosine sim to the nearest non-chemoreceptor
                              family centroid (EXCLUSION signal)
        emb_nonchemo_family - which non-chemoreceptor family that was
        emb_classA_sim      - cosine sim to the broad class-A-GPCR centroid
                              (RECALL signal)
        has_emb_data        - always True for entries in this dict
        emb_leakage_flag    - always True; Berghia is ~absent from PLM
                              pretraining corpora, so every embedding
                              judgment here is leakage-uncertain
    """
    channel: Dict[str, Dict[str, object]] = {}
    for candidate_id, vec in embeddings.items():
        family, nonchemo_sim = nearest_family(vec, nonchemo_centroids)
        classA_sim = _cosine_similarity(vec, classA_centroid)
        channel[candidate_id] = {
            "emb_nonchemo_sim": float(nonchemo_sim),
            "emb_nonchemo_family": family,
            "emb_classA_sim": float(classA_sim),
            "has_emb_data": True,
            "emb_leakage_flag": True,
        }
    return channel
