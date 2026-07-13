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


def mahalanobis_sq(x: np.ndarray, mu: np.ndarray, precision: np.ndarray) -> float:
    """Squared Mahalanobis distance ``(x-mu)ᵀ P (x-mu)`` for a precision matrix P.

    Measures distance in the data's own spread (P = Σ⁻¹ whitens the cloud), so it
    is affine-invariant — mean-centering/rescaling the embeddings cannot change
    it, which is why it is immune to the ESM-C anisotropy that made raw cosine
    preprocessing-fragile.
    """
    d = np.asarray(x, dtype=float) - np.asarray(mu, dtype=float)
    return float(d @ np.asarray(precision, dtype=float) @ d)


def shrinkage_precision(matrix: np.ndarray) -> np.ndarray:
    """Precision matrix (Σ⁻¹) of ``matrix``'s rows via Ledoit-Wolf shrinkage.

    Shrinkage is essential here: the reference set has far fewer sequences
    (~206) than embedding dimensions (960/1152), so the raw sample covariance is
    singular and its naive inverse unusable. Ledoit-Wolf (2004) shrinks toward a
    scaled identity, yielding a well-conditioned, invertible Σ for a stable tied
    covariance. sklearn is imported lazily so the pure-cosine functions above
    stay importable without it.
    """
    from sklearn.covariance import LedoitWolf

    cov = LedoitWolf().fit(np.asarray(matrix, dtype=float)).covariance_
    return np.linalg.pinv(cov)


def family_prototypes(
    embeddings: Dict[str, np.ndarray],
    labels: Dict[str, str],
    k: int = 3,
) -> Dict[str, np.ndarray]:
    """{family: (n_proto, dim) prototype matrix} via k-means within each family.

    Multi-prototype centroids fix the multi-modality problem: families like
    aminergic/peptide are several distinct receptor subtypes, so a single mean
    lands in empty space between the modes. k-means sub-prototypes track each
    mode. A family with fewer than `k` (or <3) members can't be clustered, so it
    falls back to a single mean prototype. `labels` maps id -> family; ids absent
    from `embeddings` are ignored. sklearn is imported lazily.
    """
    from collections import defaultdict

    from sklearn.cluster import KMeans

    by_family: Dict[str, list] = defaultdict(list)
    for seq_id, family in labels.items():
        if seq_id in embeddings:
            by_family[family].append(np.asarray(embeddings[seq_id], dtype=float))

    prototypes: Dict[str, np.ndarray] = {}
    for family, vectors in by_family.items():
        matrix = np.array(vectors)
        if len(matrix) < max(k, 3):
            prototypes[family] = matrix.mean(axis=0, keepdims=True)
        else:
            km = KMeans(n_clusters=k, n_init=5, random_state=0).fit(matrix)
            prototypes[family] = km.cluster_centers_
    return prototypes


def mahalanobis_assignment(
    vec: np.ndarray,
    prototypes: Dict[str, np.ndarray],
    precision: np.ndarray,
) -> Tuple[Optional[str], float]:
    """(family, min squared-Mahalanobis distance) of the nearest prototype across
    ALL families' prototypes.

    `prototypes` maps family -> (n_proto, dim) matrix (from `family_prototypes`).
    The returned family is `S_known` (which known non-chemoreceptor family the
    candidate most resembles); the returned distance is the raw novelty quantity
    (small = sits inside a known family, large = divergent from all of them).
    Ties break alphabetically by family. Empty `prototypes` -> (None, inf).
    """
    best_family: Optional[str] = None
    best_dist = float("inf")
    for family in sorted(prototypes):
        protos = np.atleast_2d(prototypes[family])
        dist = min(mahalanobis_sq(vec, p, precision) for p in protos)
        if dist < best_dist:
            best_dist, best_family = dist, family
    return best_family, best_dist


def novelty_score(
    vec: np.ndarray,
    prototypes: Dict[str, np.ndarray],
    precision: np.ndarray,
    background_mean: Optional[np.ndarray] = None,
) -> float:
    """`S_novel`: how divergent a candidate is from ALL known families (high =
    novel/divergent — the positive ranking axis the user wants surfaced).

    Default = PLAIN Mahalanobis: the min squared-Mahalanobis distance to any known
    family prototype. If `background_mean` is given, RELATIVE Mahalanobis (Ren et
    al. 2021): that distance minus the distance to a single background centroid —
    the conservative config option that discounts generic outlierness (kept
    available; plain is the default per the near-OOD result + design intent).
    """
    _, min_dist = mahalanobis_assignment(vec, prototypes, precision)
    if background_mean is None:
        return min_dist
    return min_dist - mahalanobis_sq(vec, background_mean, precision)


def mahalanobis_channel(
    candidate_embeddings: Dict[str, np.ndarray],
    reference_embeddings: Dict[str, np.ndarray],
    reference_labels: Dict[str, str],
    k: int = 3,
    relative: bool = False,
) -> Dict[str, Dict[str, object]]:
    """Per-candidate Mahalanobis exclusion/novelty channel (Phase-0 scorer).

    Builds, from the labelled reference set: multi-prototype family centroids, a
    tied (within-family-pooled) Ledoit-Wolf precision, and a global background
    mean. Then per candidate emits:
        emb_nonchemo_family - nearest known non-chemoreceptor family (annotation)
        emb_novelty         - S_novel: divergence from ALL known families
                              (plain Mahalanobis; relative if `relative=True`).
                              HIGH = novel/divergent — the POSITIVE ranking axis.
        has_emb_data        - always True for entries present here
        emb_leakage_flag    - always True (Berghia ~absent from PLM corpora)

    Like the cosine channel, this NEVER emits a positive similarity-to-a-known-
    chemoreceptor score: novelty is distance FROM known non-chemoreceptor
    families, and a divergent candidate is exactly what we want to rank up.
    Candidates absent from `candidate_embeddings` simply don't appear.
    """
    ref_ids = [i for i in reference_labels if i in reference_embeddings]
    prototypes = family_prototypes(reference_embeddings, reference_labels, k)

    # Tied covariance: pool each family's mean-centered vectors, then shrink.
    by_family: Dict[str, list] = {}
    for seq_id in ref_ids:
        by_family.setdefault(reference_labels[seq_id], []).append(
            np.asarray(reference_embeddings[seq_id], dtype=float)
        )
    centered = np.vstack(
        [np.array(v) - np.array(v).mean(axis=0) for v in by_family.values()]
    )
    precision = shrinkage_precision(centered)
    background_mean = np.array([reference_embeddings[i] for i in ref_ids],
                               dtype=float).mean(axis=0)

    channel: Dict[str, Dict[str, object]] = {}
    for candidate_id, vec in candidate_embeddings.items():
        vec = np.asarray(vec, dtype=float)
        family, _ = mahalanobis_assignment(vec, prototypes, precision)
        nov = novelty_score(
            vec, prototypes, precision,
            background_mean=background_mean if relative else None,
        )
        channel[candidate_id] = {
            "emb_nonchemo_family": family,
            "emb_novelty": float(nov),
            "has_emb_data": True,
            "emb_leakage_flag": True,
        }
    return channel


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
