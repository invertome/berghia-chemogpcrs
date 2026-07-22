#!/usr/bin/env python3
"""taar_embedding_placement.py — should the 36 vertebrate TAAR2-9 be ADDED to
the embedding reference anchors (as their own prototype)?

Answer it from embedding evidence, not homology intuition. The decision rule
(spec docs/plans/2026-07-22-taar-embedding-placement.md, bead h0y0) is a
conjunction:

    ADD (as own prototype)  iff
      (Q1) the 36 TAAR are statistically INDISTINGUISHABLE from the aminergic
           anchors in the production tied-covariance Mahalanobis metric, AND
      (Q2) adding them does NOT suppress candidate novelty beyond a random-36
           null.
    else HOLD.

This module is a PURE numpy/scipy/sklearn analysis layer. It REUSES the
production scorer scripts/embedding_evidence.py for every distance/novelty
computation (load_embeddings, family_prototypes, shrinkage_precision,
mahalanobis_sq, mahalanobis_assignment, novelty_score, mahalanobis_channel) —
it never reimplements the scorer. Everything here is testable LOCALLY with
synthetic npz + synthetic embeddings; no GPU, no Unity, no torch.

Fail-loud discipline: a missing id-prefix, an empty set, a degenerate anchor
family, or a zero-overlap TSV join RAISES rather than emitting a
plausible-but-wrong number. Coverage (n per set) is reported on every run.

Statistics implemented (name -> what it computes):
  bh_fdr                     Benjamini-Hochberg adjusted p-values.
  conformal_pvalues          split-conformal set-membership p vs aminergic
                             (answer of record for Q1).
  rbf_mmd2 / energy_distance two-sample statistics in the whitened metric.
  two_sample_tests           MMD + energy with a permutation null.
  family_probe_posteriors    calibrated multiclass family probe -> P(aminergic).
  knn_class_purity           family composition of each TAAR's k nearest anchors.
  lda_illustration_axis      shrinkage-LDA aminergic<->probe axis (VIZ ONLY).
  novelty_leave_in_out       Q2 exact leave-in/leave-out novelty (prod scorer).
  nearest_anchor_localization does any candidate's nearest prototype become TAAR?
  rank_agreement             Wilcoxon signed-rank + Kendall tau-b + top-k Jaccard.
  random36_null              Monte-Carlo "add 36 random anchors" null for Q2.
  robustness_novelty         novelty under Mahalanobis / kNN-dist / LOF (+ RelMaha).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# scripts/ on sys.path so the production scorer imports whether we are run as a
# CLI (scripts/analysis/…) or imported as analysis.taar_embedding_placement.
_SCRIPTS = Path(__file__).resolve().parent.parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from embedding_evidence import (  # noqa: E402  (reused production scorer)
    family_prototypes,
    load_embeddings,
    mahalanobis_assignment,
    mahalanobis_channel,
    mahalanobis_sq,
    novelty_score,
    shrinkage_precision,
)

ANCHOR_PREFIX = "ANCHOR_"
PROBE_PREFIX = "PROBE_"
TAAR_PREFIX = "TAAR_"
CANDIDATE_PREFIX = "BersteEVm"
AMINERGIC = "aminergic"


@dataclass
class Config:
    k: int = 3                      # k-means prototypes per family
    n_perm: int = 2000              # two-sample permutation replicates
    n_random: int = 500             # random-36 Monte-Carlo draws
    knn: int = 15                   # neighbours for purity / kNN novelty / LOF
    alpha: float = 0.05             # FDR / significance level
    seed: int = 0
    split_frac: float = 0.5         # conformal train/calibration split
    min_aminergic: int = 4          # below this, conformal is undefined -> fail
    topk: Tuple[int, ...] = (20, 50)
    calib: str = "sigmoid"          # temperature-like; "isotonic" also allowed


# --------------------------------------------------------------------------- #
# sets, coverage, family join                                                  #
# --------------------------------------------------------------------------- #
def split_by_prefix(embeddings: Dict[str, np.ndarray]) -> Dict[str, Dict[str, np.ndarray]]:
    """Partition {id: vec} into anchor/candidate/probe/taar by id prefix.

    Fail loud on an id that matches no known prefix — an unclassified id is a
    silent-corruption hazard (it would vanish from every downstream set).
    """
    out: Dict[str, Dict[str, np.ndarray]] = {
        "anchor": {}, "candidate": {}, "probe": {}, "taar": {}
    }
    unknown: List[str] = []
    for sid, vec in embeddings.items():
        if sid.startswith(ANCHOR_PREFIX):
            out["anchor"][sid] = vec
        elif sid.startswith(TAAR_PREFIX):
            out["taar"][sid] = vec
        elif sid.startswith(PROBE_PREFIX):
            out["probe"][sid] = vec
        elif sid.startswith(CANDIDATE_PREFIX):
            out["candidate"][sid] = vec
        else:
            unknown.append(sid)
    if unknown:
        raise ValueError(
            f"{len(unknown)} id(s) match no known prefix "
            f"({ANCHOR_PREFIX}/{CANDIDATE_PREFIX}/{PROBE_PREFIX}/{TAAR_PREFIX}); "
            f"first few: {unknown[:5]}"
        )
    return out


def map_anchor_families(anchor_ids: Sequence[str], anchor_labels_tsv: str) -> Dict[str, str]:
    """{ANCHOR_<class>_<tier>_<acc>: family} joined to the TSV by accession.

    Accession = the id token(s) after the 3rd underscore (UniProt accessions
    never contain '_'). Fail loud if the join overlap is empty or if any anchor
    id cannot be mapped — a silent partial join would drop real anchors from the
    reference geometry.
    """
    acc2fam: Dict[str, str] = {}
    with open(anchor_labels_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_acc = header.index("accession")
            i_fam = header.index("family")
        except ValueError as exc:
            raise ValueError(
                f"{anchor_labels_tsv} must have 'accession' and 'family' columns; "
                f"got {header}"
            ) from exc
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_acc, i_fam):
                continue
            acc2fam[parts[i_acc]] = parts[i_fam]

    fam: Dict[str, str] = {}
    unmapped: List[str] = []
    for aid in anchor_ids:
        toks = aid.split("_", 3)
        acc = toks[3] if len(toks) == 4 else aid
        if acc in acc2fam:
            fam[aid] = acc2fam[acc]
        else:
            unmapped.append(aid)
    if not fam:
        raise ValueError(
            f"zero-overlap join: none of {len(list(anchor_ids))} anchor ids "
            f"matched an accession in {anchor_labels_tsv}"
        )
    if unmapped:
        raise ValueError(
            f"{len(unmapped)} anchor id(s) unmapped to a family (e.g. "
            f"{unmapped[:5]}); refusing a silent partial join"
        )
    return fam


def validate_sets(
    sets: Dict[str, Dict[str, np.ndarray]],
    anchor_families: Dict[str, str],
    cfg: Config,
) -> None:
    """Raise on any degeneracy that would make a metric meaningless."""
    for name, d in sets.items():
        if not d:
            raise ValueError(f"empty set: {name} has 0 members")
    missing = [a for a in sets["anchor"] if a not in anchor_families]
    if missing:
        raise ValueError(
            f"{len(missing)} anchor(s) have no family label (e.g. {missing[:5]})"
        )
    fam_counts: Dict[str, int] = {}
    for a in sets["anchor"]:
        fam_counts[anchor_families[a]] = fam_counts.get(anchor_families[a], 0) + 1
    if AMINERGIC not in fam_counts:
        raise ValueError(f"no '{AMINERGIC}' anchors — Q1 pole is undefined")
    if fam_counts[AMINERGIC] < cfg.min_aminergic:
        raise ValueError(
            f"too few '{AMINERGIC}' anchors ({fam_counts[AMINERGIC]} < "
            f"{cfg.min_aminergic}); conformal calibration undefined "
            f"(single-member family)"
        )


def _coverage(sets, anchor_families) -> Dict[str, object]:
    fam_counts: Dict[str, int] = {}
    for a in sets["anchor"]:
        fam_counts[anchor_families[a]] = fam_counts.get(anchor_families[a], 0) + 1
    return {
        "anchor": len(sets["anchor"]),
        "candidate": len(sets["candidate"]),
        "probe": len(sets["probe"]),
        "taar": len(sets["taar"]),
        "anchor_families": fam_counts,
    }


# --------------------------------------------------------------------------- #
# geometry: tied precision + whitening (production metric)                     #
# --------------------------------------------------------------------------- #
def tied_precision(embeddings: Dict[str, np.ndarray], labels: Dict[str, str]) -> np.ndarray:
    """Tied (within-family-pooled) Ledoit-Wolf precision — identical construction
    to embedding_evidence.mahalanobis_channel, so the whitened metric here IS the
    production metric."""
    by_family: Dict[str, list] = {}
    for sid, family in labels.items():
        if sid in embeddings:
            by_family.setdefault(family, []).append(
                np.asarray(embeddings[sid], dtype=float)
            )
    centered = np.vstack(
        [np.array(v) - np.array(v).mean(axis=0) for v in by_family.values()]
    )
    return shrinkage_precision(centered)


def whitening_matrix(precision: np.ndarray) -> np.ndarray:
    """Matrix L with rows whitened by ``X @ L`` so Euclidean distance in the
    whitened space equals tied-cov Mahalanobis (precision = L Lᵀ)."""
    P = np.asarray(precision, dtype=float)
    try:
        return np.linalg.cholesky(P)
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh((P + P.T) / 2.0)
        w = np.clip(w, 0.0, None)
        return V @ np.diag(np.sqrt(w)) @ V.T


def _stack(embeddings: Dict[str, np.ndarray], ids: Sequence[str]) -> np.ndarray:
    return np.array([np.asarray(embeddings[i], dtype=float) for i in ids])


# --------------------------------------------------------------------------- #
# generic statistics                                                          #
# --------------------------------------------------------------------------- #
def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (returned in the input order)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(1, n + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def rbf_mmd2(Kxx: np.ndarray, Kyy: np.ndarray, Kxy: np.ndarray) -> float:
    """Biased MMD² V-statistic from precomputed kernel blocks."""
    return float(Kxx.mean() + Kyy.mean() - 2.0 * Kxy.mean())


def energy_distance(Dxx: np.ndarray, Dyy: np.ndarray, Dxy: np.ndarray) -> float:
    """Székely-Rizzo energy distance from precomputed distance blocks."""
    return float(2.0 * Dxy.mean() - Dxx.mean() - Dyy.mean())


def two_sample_tests(
    A: np.ndarray, B: np.ndarray, n_perm: int = 2000, seed: int = 0
) -> Dict[str, float]:
    """MMD (RBF, median-heuristic bandwidth) and energy distance between two
    (already whitened) samples, each with a label-permutation null.

    Returns mmd_stat, mmd_p, energy_stat, energy_p. p = (1 + #{perm >= obs}) /
    (n_perm + 1). Distances are Euclidean in the whitened space = Mahalanobis.
    """
    from scipy.spatial.distance import cdist

    A = np.atleast_2d(np.asarray(A, dtype=float))
    B = np.atleast_2d(np.asarray(B, dtype=float))
    m, n = len(A), len(B)
    Z = np.vstack([A, B])
    D = cdist(Z, Z)                       # pooled pairwise Euclidean distances
    tri = D[np.triu_indices_from(D, k=1)]
    sigma = np.median(tri[tri > 0]) if np.any(tri > 0) else 1.0
    gamma = 1.0 / (2.0 * sigma * sigma)
    K = np.exp(-(D ** 2) * gamma)

    def _mmd(idx_a, idx_b):
        return rbf_mmd2(K[np.ix_(idx_a, idx_a)], K[np.ix_(idx_b, idx_b)],
                        K[np.ix_(idx_a, idx_b)])

    def _energy(idx_a, idx_b):
        return energy_distance(D[np.ix_(idx_a, idx_a)], D[np.ix_(idx_b, idx_b)],
                               D[np.ix_(idx_a, idx_b)])

    a0, b0 = np.arange(m), np.arange(m, m + n)
    mmd_obs, energy_obs = _mmd(a0, b0), _energy(a0, b0)

    rng = np.random.RandomState(seed)
    mmd_ge = energy_ge = 0
    all_idx = np.arange(m + n)
    for _ in range(n_perm):
        perm = rng.permutation(all_idx)
        ia, ib = perm[:m], perm[m:]
        if _mmd(ia, ib) >= mmd_obs:
            mmd_ge += 1
        if _energy(ia, ib) >= energy_obs:
            energy_ge += 1
    return {
        "mmd_stat": mmd_obs,
        "mmd_p": (1 + mmd_ge) / (n_perm + 1),
        "energy_stat": energy_obs,
        "energy_p": (1 + energy_ge) / (n_perm + 1),
    }


# --------------------------------------------------------------------------- #
# Q1 — placement                                                              #
# --------------------------------------------------------------------------- #
def conformal_pvalues(
    taar_white: np.ndarray, aminergic_white: np.ndarray, cfg: Config
) -> List[float]:
    """Split-conformal set-membership p-value of each TAAR vs the aminergic
    anchors. Conformity (nonconformity) score = min squared distance to
    aminergic prototypes fit on a disjoint train split; per-TAAR
    p = (1 + #{cal >= s_taar}) / (n_cal + 1). Small p -> non-aminergic.
    """
    rng = np.random.RandomState(cfg.seed)
    idx = rng.permutation(len(aminergic_white))
    n_tr = max(cfg.k, int(round(cfg.split_frac * len(idx))))
    n_tr = min(n_tr, len(idx) - 1)        # keep >=1 calibration point
    train = aminergic_white[idx[:n_tr]]
    calib = aminergic_white[idx[n_tr:]]

    proto = family_prototypes(
        {str(i): v for i, v in enumerate(train)},
        {str(i): AMINERGIC for i in range(len(train))},
        k=cfg.k,
    )[AMINERGIC]

    def score(x):  # min squared Euclidean = tied-cov Mahalanobis (pre-whitened)
        return float(min(np.sum((x - p) ** 2) for p in np.atleast_2d(proto)))

    cal_scores = np.array([score(x) for x in calib])
    n_cal = len(cal_scores)
    pvals = []
    for x in np.atleast_2d(taar_white):
        s = score(x)
        pvals.append((1 + int(np.sum(cal_scores >= s))) / (n_cal + 1))
    return pvals


def family_probe_posteriors(
    anchor_white: np.ndarray,
    anchor_fams: Sequence[str],
    taar_white: np.ndarray,
    cfg: Config,
) -> Tuple[List[Dict[str, float]], List[str], str]:
    """Calibrated multiclass family probe -> per-TAAR posterior over families
    (including P(aminergic)). Returns (posteriors, class_order, calib_method).
    Fail-soft on calibration (SUPPORT metric): falls back to an uncalibrated
    multinomial fit, labelled as such."""
    from collections import Counter

    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.linear_model import LogisticRegression

    classes = sorted(set(anchor_fams))
    y = np.asarray(anchor_fams)
    base = LogisticRegression(max_iter=2000, C=1.0)
    min_class = min(Counter(anchor_fams).values())
    method = "uncalibrated"
    clf = None
    if min_class >= 2:
        cv = int(min(3, min_class))
        try:
            clf = CalibratedClassifierCV(base, method=cfg.calib, cv=cv)
            clf.fit(anchor_white, y)
            method = cfg.calib
        except Exception:
            # Deliberate fail-soft on a SUPPORT metric only: if calibration
            # cannot fit (tiny/degenerate class in a CV fold) we downgrade to an
            # uncalibrated multinomial fit. This is TRACEABLE — `method` stays
            # "uncalibrated" and is written to every probe row — never a silent
            # substitution, and this probe never feeds the verdict.
            clf = None
    if clf is None:
        clf = base.fit(anchor_white, y)
    order = list(clf.classes_)
    proba = clf.predict_proba(np.atleast_2d(taar_white))
    posteriors = [{order[j]: float(row[j]) for j in range(len(order))} for row in proba]
    predicted = [order[int(np.argmax(row))] for row in proba]
    return posteriors, predicted, method, order


def knn_class_purity(
    anchor_white: np.ndarray,
    anchor_fams: Sequence[str],
    taar_white: np.ndarray,
    k: int,
) -> List[Dict[str, object]]:
    """Family composition of each TAAR's k nearest anchors (whitened metric)."""
    from collections import Counter

    from sklearn.neighbors import NearestNeighbors

    k = min(k, len(anchor_white))
    nn = NearestNeighbors(n_neighbors=k).fit(anchor_white)
    _, ind = nn.kneighbors(np.atleast_2d(taar_white))
    fams = np.asarray(anchor_fams)
    rows = []
    for neigh in ind:
        comp = Counter(fams[neigh])
        majority = comp.most_common(1)[0][0]
        rows.append({
            "k": k,
            "aminergic_fraction": comp.get(AMINERGIC, 0) / k,
            "majority_family": majority,
            "composition": dict(comp),
        })
    return rows


def lda_illustration_axis(
    aminergic_white: np.ndarray,
    probe_white: np.ndarray,
    project: Dict[str, np.ndarray],
) -> Dict[str, List[float]]:
    """Shrinkage-LDA aminergic<->probe 1-D axis; project the named point sets.
    ILLUSTRATION ONLY (n=36 in high-D is singular for distance claims)."""
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

    X = np.vstack([aminergic_white, probe_white])
    y = np.array([0] * len(aminergic_white) + [1] * len(probe_white))
    try:
        lda = LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")
        lda.fit(X, y)
        w = lda.coef_.ravel()
    except Exception:
        # Deliberate tolerance: the LDA axis is VISUALISATION ONLY (spec §4 — it
        # is singular at n=36 in high-D and explicitly never a distance claim).
        # An empty dict means "no illustration axis"; the only consumer guards
        # with `if q1["lda_axis"]:` before writing the optional viz TSV. It
        # never feeds Q1/Q2 or the verdict, so collapsing here is harmless.
        return {}
    axis = {}
    for name, M in project.items():
        M = np.atleast_2d(M)
        axis[name] = (M @ w).tolist()
    return axis


# --------------------------------------------------------------------------- #
# Q2 — add-impact                                                             #
# --------------------------------------------------------------------------- #
def novelty_leave_in_out(
    candidate_emb: Dict[str, np.ndarray],
    anchor_emb: Dict[str, np.ndarray],
    anchor_fam: Dict[str, str],
    taar_emb: Dict[str, np.ndarray],
    cfg: Config,
) -> Dict[str, object]:
    """Exact leave-in/leave-out candidate novelty via the production scorer.

    WITHOUT: prototypes from anchors only. WITH: the 36 TAAR added AS THEIR OWN
    'TAAR' family (not merged into aminergic). Also the RelativeMahalanobis
    variant (novelty_score background_mean) for reporting, without changing the
    scorer. Returns per-candidate rows + the paired arrays + the with-channel
    nearest-family assignment (for localization)."""
    with_emb = {**anchor_emb, **taar_emb}
    with_fam = {**anchor_fam, **{t: "TAAR" for t in taar_emb}}

    ch_wo = mahalanobis_channel(candidate_emb, anchor_emb, anchor_fam, cfg.k)
    ch_wi = mahalanobis_channel(candidate_emb, with_emb, with_fam, cfg.k)
    ch_wo_rel = mahalanobis_channel(candidate_emb, anchor_emb, anchor_fam, cfg.k,
                                    relative=True)
    ch_wi_rel = mahalanobis_channel(candidate_emb, with_emb, with_fam, cfg.k,
                                    relative=True)

    ids = list(candidate_emb)
    without = np.array([ch_wo[c]["emb_novelty"] for c in ids])
    withv = np.array([ch_wi[c]["emb_novelty"] for c in ids])
    rows = []
    for i, c in enumerate(ids):
        rows.append({
            "candidate_id": c,
            "without": float(without[i]),
            "with": float(withv[i]),
            "delta": float(withv[i] - without[i]),
            "without_rel": float(ch_wo_rel[c]["emb_novelty"]),
            "with_rel": float(ch_wi_rel[c]["emb_novelty"]),
            "with_nearest_family": ch_wi[c]["emb_nonchemo_family"],
        })
    return {"ids": ids, "without": without, "with": withv, "rows": rows}


def nearest_anchor_localization(rows: List[Dict[str, object]]) -> Dict[str, object]:
    """How many candidates now assign to the TAAR prototype (nearest family ==
    'TAAR' under the WITH scenario)? If none -> the block is inert."""
    flipped = [r["candidate_id"] for r in rows if r["with_nearest_family"] == "TAAR"]
    return {"count": len(flipped), "ids": flipped}


def rank_agreement(without: np.ndarray, withv: np.ndarray, topk: Sequence[int]) -> Dict[str, object]:
    """Wilcoxon signed-rank on the paired novelty deltas, Kendall tau-b on the
    two rankings, and top-k Jaccard of the highest-novelty candidate sets."""
    from scipy.stats import kendalltau, wilcoxon

    d = withv - without
    if np.allclose(d, 0.0):
        w_stat, w_p = 0.0, 1.0
    else:
        try:
            res = wilcoxon(without, withv, zero_method="wilcox")
            w_stat, w_p = float(res.statistic), float(res.pvalue)
        except ValueError:
            w_stat, w_p = 0.0, 1.0
    tau, tau_p = kendalltau(without, withv)
    order_wo = np.argsort(-without)
    order_wi = np.argsort(-withv)
    jac = {}
    for k in topk:
        kk = min(k, len(without))
        a, b = set(order_wo[:kk]), set(order_wi[:kk])
        jac[k] = len(a & b) / len(a | b) if (a | b) else 1.0
    return {
        "wilcoxon_stat": w_stat, "wilcoxon_p": w_p,
        "kendall_tau": float(tau) if tau == tau else 0.0,
        "kendall_p": float(tau_p) if tau_p == tau_p else 1.0,
        "jaccard": jac,
    }


def random36_null(
    candidate_emb: Dict[str, np.ndarray],
    anchor_emb: Dict[str, np.ndarray],
    anchor_fam: Dict[str, str],
    without: np.ndarray,
    n_add: int,
    cfg: Config,
) -> Dict[str, object]:
    """Monte-Carlo: does the TAAR block suppress candidate novelty MORE than an
    arbitrary n_add-sequence block added as its own prototype?

    Draw n_add ids from the anchor pool, add them as a NEW own-family 'RAND'
    (parallel to the TAAR own-prototype mechanic), re-score candidate novelty,
    and measure the per-candidate and aggregate suppression S = without - with.
    D_taar comes from the caller (mean of the real WITH scenario). Returns the
    aggregate p, per-candidate p (+ BH), and the null mean/sd.
    """
    ids = list(candidate_emb)
    anchor_ids = list(anchor_emb)
    n_add = min(n_add, len(anchor_ids))
    rng = np.random.RandomState(cfg.seed + 1)

    S_rand = np.zeros((cfg.n_random, len(ids)))
    D_rand = np.zeros(cfg.n_random)
    for d in range(cfg.n_random):
        pick = rng.choice(len(anchor_ids), size=n_add, replace=False)
        rand_emb = dict(anchor_emb)
        rand_fam = dict(anchor_fam)
        for j, p in enumerate(pick):
            rid = f"RAND_{j}"
            rand_emb[rid] = anchor_emb[anchor_ids[p]]
            rand_fam[rid] = "RAND"
        ch = mahalanobis_channel(candidate_emb, rand_emb, rand_fam, cfg.k)
        withr = np.array([ch[c]["emb_novelty"] for c in ids])
        S_rand[d] = without - withr
        D_rand[d] = S_rand[d].mean()
    return {"S_rand": S_rand, "D_rand": D_rand,
            "null_mean": float(D_rand.mean()), "null_sd": float(D_rand.std())}


# --------------------------------------------------------------------------- #
# robustness                                                                  #
# --------------------------------------------------------------------------- #
def _knn_novelty(query: np.ndarray, ref: np.ndarray, k: int) -> np.ndarray:
    from sklearn.neighbors import NearestNeighbors

    k = min(k, len(ref))
    nn = NearestNeighbors(n_neighbors=k).fit(ref)
    dist, _ = nn.kneighbors(np.atleast_2d(query))
    return dist.mean(axis=1)


def _lof_novelty(query: np.ndarray, ref: np.ndarray, k: int) -> np.ndarray:
    from sklearn.neighbors import LocalOutlierFactor

    k = min(k, max(1, len(ref) - 1))
    clf = LocalOutlierFactor(n_neighbors=k, novelty=True).fit(ref)
    return -clf.score_samples(np.atleast_2d(query))  # higher = more novel


def robustness_novelty(
    cand_white: np.ndarray,
    anchor_white: np.ndarray,
    taar_white: np.ndarray,
    maha_without: np.ndarray,
    maha_with: np.ndarray,
    rel_without: np.ndarray,
    rel_with: np.ndarray,
    cfg: Config,
) -> Dict[str, object]:
    """Recompute candidate novelty under (a) tied-cov Mahalanobis (from the prod
    scorer, passed in), (b) kNN mean-distance, (c) LOF. Report whether the Q2
    suppression conclusion is stable across all three + the RelativeMahalanobis
    read-out (scorer unchanged)."""
    from scipy.stats import kendalltau

    with_ref = np.vstack([anchor_white, taar_white])
    methods: Dict[str, Dict[str, object]] = {}

    def _entry(wo, wi):
        supp = float(np.mean(wo - wi))
        tau, _ = kendalltau(wo, wi)
        return {
            "mean_suppression": supp,
            "suppression_sign": int(np.sign(supp)),
            "kendall_tau": float(tau) if tau == tau else 0.0,
        }

    methods["mahalanobis"] = _entry(maha_without, maha_with)
    methods["knn"] = _entry(
        _knn_novelty(cand_white, anchor_white, cfg.knn),
        _knn_novelty(cand_white, with_ref, cfg.knn),
    )
    methods["lof"] = _entry(
        _lof_novelty(cand_white, anchor_white, cfg.knn),
        _lof_novelty(cand_white, with_ref, cfg.knn),
    )
    ref_sign = methods["mahalanobis"]["suppression_sign"]
    stable = all(m["suppression_sign"] == ref_sign for m in methods.values())
    return {
        "methods": methods,
        "stable": bool(stable),
        "relative_maha": {
            "mean_suppression": float(np.mean(rel_without - rel_with)),
            "note": "RelativeMahalanobis (Ren 2021); reported only, scorer unchanged",
        },
    }


# --------------------------------------------------------------------------- #
# orchestration: per model                                                    #
# --------------------------------------------------------------------------- #
def run_per_model(
    embeddings: Dict[str, np.ndarray],
    anchor_families: Dict[str, str],
    model_name: str,
    cfg: Optional[Config] = None,
) -> Dict[str, object]:
    """Full Q1 + Q2 + robustness analysis for ONE model. Pure compute (no IO)."""
    cfg = cfg or Config()
    if not embeddings:
        raise ValueError("empty embeddings")
    sets = split_by_prefix(embeddings)
    validate_sets(sets, anchor_families, cfg)
    coverage = _coverage(sets, anchor_families)

    anchor_emb = sets["anchor"]
    cand_emb = sets["candidate"]
    probe_emb = sets["probe"]
    taar_emb = sets["taar"]

    # production tied-cov metric (anchors only) -> whitening.
    precision = tied_precision(anchor_emb, anchor_families)
    L = whitening_matrix(precision)

    anchor_ids = list(anchor_emb)
    anchor_fams = [anchor_families[a] for a in anchor_ids]
    anchor_white = _stack(anchor_emb, anchor_ids) @ L
    amin_ids = [a for a in anchor_ids if anchor_families[a] == AMINERGIC]
    amin_white = _stack(anchor_emb, amin_ids) @ L
    probe_white = _stack(probe_emb, list(probe_emb)) @ L
    taar_white = _stack(taar_emb, list(taar_emb)) @ L
    cand_white = _stack(cand_emb, list(cand_emb)) @ L

    # ---- Q1 ----
    conf_p = conformal_pvalues(taar_white, amin_white, cfg)
    conf_bh = bh_fdr(conf_p)
    rejected = conf_bh < cfg.alpha
    conf_rows = [
        {"taar_id": t, "p_conformal": float(conf_p[i]),
         "p_bh": float(conf_bh[i]), "rejected": bool(rejected[i])}
        for i, t in enumerate(taar_emb)
    ]
    mmd_amin = two_sample_tests(taar_white, amin_white, cfg.n_perm, cfg.seed)
    mmd_probe = two_sample_tests(taar_white, probe_white, cfg.n_perm, cfg.seed)

    post, pred, calib_method, class_order = family_probe_posteriors(
        anchor_white, anchor_fams, taar_white, cfg)
    probe_rows = [
        {"taar_id": t, "p_aminergic": float(post[i].get(AMINERGIC, 0.0)),
         "predicted_family": pred[i], "calib_method": calib_method,
         "posterior": post[i]}
        for i, t in enumerate(taar_emb)
    ]
    purity = knn_class_purity(anchor_white, anchor_fams, taar_white, cfg.knn)
    purity_rows = [{"taar_id": t, **purity[i]} for i, t in enumerate(taar_emb)]
    lda_axis = lda_illustration_axis(
        amin_white, probe_white,
        {"taar": taar_white, "candidate": cand_white})

    indistinguishable = int(rejected.sum()) == 0
    q1 = {
        "conformal": conf_rows,
        "n_rejected": int(rejected.sum()),
        "indistinguishable": bool(indistinguishable),
        "mmd_taar_aminergic": {"stat": mmd_amin["mmd_stat"], "p": mmd_amin["mmd_p"]},
        "energy_taar_aminergic": {"stat": mmd_amin["energy_stat"], "p": mmd_amin["energy_p"]},
        "mmd_taar_probes": {"stat": mmd_probe["mmd_stat"], "p": mmd_probe["mmd_p"]},
        "energy_taar_probes": {"stat": mmd_probe["energy_stat"], "p": mmd_probe["energy_p"]},
        "probe": probe_rows,
        "probe_class_order": class_order,
        "knn_purity": purity_rows,
        "lda_axis": lda_axis,
    }

    # ---- Q2 ----
    lio = novelty_leave_in_out(cand_emb, anchor_emb, anchor_families, taar_emb, cfg)
    without, withv = lio["without"], lio["with"]
    ra = rank_agreement(without, withv, cfg.topk)
    loc = nearest_anchor_localization(lio["rows"])

    D_taar_per = without - withv
    D_taar_agg = float(D_taar_per.mean())
    null = random36_null(cand_emb, anchor_emb, anchor_families, without,
                         n_add=len(taar_emb), cfg=cfg)
    S_rand, D_rand = null["S_rand"], null["D_rand"]
    p_agg = (1 + int(np.sum(D_rand >= D_taar_agg))) / (cfg.n_random + 1)
    per_cand_p = np.array([
        (1 + int(np.sum(S_rand[:, j] >= D_taar_per[j]))) / (cfg.n_random + 1)
        for j in range(len(without))
    ])
    per_cand_bh = bh_fdr(per_cand_p)
    random36_rows = [
        {"candidate_id": lio["ids"][j], "S_taar": float(D_taar_per[j]),
         "p_random36": float(per_cand_p[j]), "p_bh": float(per_cand_bh[j])}
        for j in range(len(without))
    ]
    no_excess = p_agg >= cfg.alpha
    q2 = {
        "novelty": lio["rows"],
        "wilcoxon": {"stat": ra["wilcoxon_stat"], "p": ra["wilcoxon_p"]},
        "kendall": {"tau": ra["kendall_tau"], "p": ra["kendall_p"]},
        "jaccard": ra["jaccard"],
        "nearest_anchor_flips": loc,
        "random36": {
            "D_taar": D_taar_agg, "p_agg": float(p_agg),
            "null_mean": null["null_mean"], "null_sd": null["null_sd"],
            "n_bh_sig": int(np.sum(per_cand_bh < cfg.alpha)),
            "per_candidate": random36_rows,
        },
        "no_excess_suppression": bool(no_excess),
    }

    # ---- robustness ----
    rel_wo = np.array([r["without_rel"] for r in lio["rows"]])
    rel_wi = np.array([r["with_rel"] for r in lio["rows"]])
    robustness = robustness_novelty(
        cand_white, anchor_white, taar_white,
        without, withv, rel_wo, rel_wi, cfg)

    verdict = "ADD" if (indistinguishable and no_excess) else "HOLD"
    reason = (
        f"Q1 indistinguishable={indistinguishable} "
        f"(conformal rejects {int(rejected.sum())}/{len(taar_emb)} at FDR {cfg.alpha}); "
        f"Q2 no_excess_suppression={no_excess} "
        f"(random-36 aggregate p={p_agg:.3g}, nearest-anchor flips="
        f"{loc['count']}). Rule: ADD iff both true."
    )
    return {
        "model": model_name,
        "coverage": coverage,
        "q1": q1,
        "q2": q2,
        "robustness": robustness,
        "verdict": verdict,
        "verdict_reason": reason,
        "config": vars(cfg) | {"topk": list(cfg.topk)},
    }


# --------------------------------------------------------------------------- #
# consensus                                                                   #
# --------------------------------------------------------------------------- #
def run_consensus(results: Sequence[Dict[str, object]]) -> Dict[str, object]:
    """Combine per-model results the way the project's rank aggregation would:
    report each model's verdict, a simple agreement summary, and a Borda
    mean-rank aggregation of per-candidate suppression across models. A verdict
    disagreement is surfaced explicitly (it is itself a finding)."""
    per_model = {r["model"]: r["verdict"] for r in results}
    verdicts = set(per_model.values())
    agreement = len(verdicts) == 1
    consensus_verdict = next(iter(verdicts)) if agreement else "DISAGREEMENT"

    # rank-aggregate per-candidate suppression (S = without - with) across models.
    rankagg: List[Dict[str, object]] = []
    have_q2 = all("q2" in r and r["q2"].get("novelty") for r in results)
    if have_q2 and results:
        common = set.intersection(*[
            {row["candidate_id"] for row in r["q2"]["novelty"]} for r in results
        ])
        supp_by_model = {}
        for r in results:
            s = {row["candidate_id"]: (row["without"] - row["with"])
                 for row in r["q2"]["novelty"]}
            supp_by_model[r["model"]] = s
        ranks = {c: [] for c in common}
        for r in results:
            s = supp_by_model[r["model"]]
            ids = [c for c in common]
            order = sorted(ids, key=lambda c: -s[c])          # most suppressed first
            for rank, c in enumerate(order):
                ranks[c].append(rank)
        rankagg = sorted(
            [{"candidate_id": c, "mean_rank": float(np.mean(ranks[c])),
              **{f"suppression_{m}": float(supp_by_model[m][c]) for m in per_model}}
             for c in common],
            key=lambda d: d["mean_rank"],
        )
    return {
        "per_model": per_model,
        "agreement": bool(agreement),
        "consensus_verdict": consensus_verdict,
        "candidate_suppression_rankagg": rankagg,
    }


# --------------------------------------------------------------------------- #
# output writing                                                              #
# --------------------------------------------------------------------------- #
def _write_tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if c is None else str(c) for c in row) + "\n")


def _maybe_umap(out_dir: Path, points: Dict[str, np.ndarray]) -> None:
    """UMAP illustration; wrapped so its absence never fails the run."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import umap  # noqa: F401

        labels, mats = [], []
        for name, M in points.items():
            M = np.atleast_2d(M)
            labels += [name] * len(M)
            mats.append(M)
        X = np.vstack(mats)
        emb2d = umap.UMAP(random_state=0).fit_transform(X)
        fig, ax = plt.subplots(figsize=(6, 5))
        for name in points:
            mask = np.array(labels) == name
            ax.scatter(emb2d[mask, 0], emb2d[mask, 1], s=8, label=name, alpha=0.6)
        ax.legend()
        ax.set_title("UMAP (illustration only — not for distance claims)")
        fig.savefig(out_dir / "umap.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
    except Exception as exc:
        # Deliberate tolerance mandated by the spec (§7): the UMAP png is an
        # OPTIONAL illustration and its absence must never fail the run. We
        # leave a trace (a note file + stderr) so a skipped viz is legible and
        # not confused with a computed result — it never feeds the verdict.
        (out_dir / "umap.SKIPPED.txt").write_text(f"UMAP viz skipped: {exc!r}\n")
        print(f"[viz] UMAP skipped ({exc!r}); analysis outputs unaffected",
              file=sys.stderr, flush=True)


def write_outputs(results: Dict[str, object], out_dir: str,
                  whitened: Optional[Dict[str, np.ndarray]] = None) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    q1, q2 = results["q1"], results["q2"]

    _write_tsv(out / "q1_conformal.tsv",
               ["taar_id", "p_conformal", "p_bh", "rejected"],
               [[r["taar_id"], r["p_conformal"], r["p_bh"], r["rejected"]]
                for r in q1["conformal"]])

    ts_rows = []
    for comp, key in [("taar_vs_aminergic", "aminergic"), ("taar_vs_probes", "probes")]:
        ts_rows.append([comp, "mmd", q1[f"mmd_taar_{key}"]["stat"], q1[f"mmd_taar_{key}"]["p"]])
        ts_rows.append([comp, "energy", q1[f"energy_taar_{key}"]["stat"], q1[f"energy_taar_{key}"]["p"]])
    _write_tsv(out / "q1_twosample.tsv", ["comparison", "test", "statistic", "p_perm"], ts_rows)

    order = q1["probe_class_order"]
    _write_tsv(out / "q1_probe_posterior.tsv",
               ["taar_id", "p_aminergic", "predicted_family", "calib_method"] + order,
               [[r["taar_id"], r["p_aminergic"], r["predicted_family"], r["calib_method"]]
                + [r["posterior"].get(c, 0.0) for c in order] for r in q1["probe"]])

    _write_tsv(out / "q1_knn_purity.tsv",
               ["taar_id", "k", "aminergic_fraction", "majority_family", "composition"],
               [[r["taar_id"], r["k"], r["aminergic_fraction"], r["majority_family"],
                 json.dumps(r["composition"])] for r in q1["knn_purity"]])

    if q1["lda_axis"]:
        lda_rows = []
        for name, coords in q1["lda_axis"].items():
            for i, c in enumerate(coords):
                lda_rows.append([f"{name}_{i}", name, c])
        _write_tsv(out / "q1_lda_axis_illustration.tsv", ["point", "set", "lda_coord"], lda_rows)

    _write_tsv(out / "q2_novelty_deltas.tsv",
               ["candidate_id", "novelty_without", "novelty_with", "delta",
                "novelty_without_rel", "novelty_with_rel", "with_nearest_family"],
               [[r["candidate_id"], r["without"], r["with"], r["delta"],
                 r["without_rel"], r["with_rel"], r["with_nearest_family"]]
                for r in q2["novelty"]])

    _write_tsv(out / "q2_rank_agreement.tsv", ["metric", "value"],
               [["wilcoxon_stat", q2["wilcoxon"]["stat"]],
                ["wilcoxon_p", q2["wilcoxon"]["p"]],
                ["kendall_tau", q2["kendall"]["tau"]],
                ["kendall_p", q2["kendall"]["p"]]]
               + [[f"jaccard_top{k}", v] for k, v in q2["jaccard"].items()])

    loc = q2["nearest_anchor_flips"]
    _write_tsv(out / "q2_nearest_anchor.tsv", ["candidate_id", "nearest_family_with"],
               [[cid, "TAAR"] for cid in loc["ids"]]
               or [[f"# flips={loc['count']}", ""]])

    r36 = q2["random36"]
    _write_tsv(out / "q2_random36_null.tsv",
               ["candidate_id", "S_taar", "p_random36", "p_bh"],
               [[r["candidate_id"], r["S_taar"], r["p_random36"], r["p_bh"]]
                for r in r36["per_candidate"]]
               + [["__aggregate__", r36["D_taar"], r36["p_agg"], ""]])

    rob = results["robustness"]
    _write_tsv(out / "robustness.tsv",
               ["method", "mean_suppression", "suppression_sign", "kendall_tau"],
               [[m, v["mean_suppression"], v["suppression_sign"], v["kendall_tau"]]
                for m, v in rob["methods"].items()]
               + [["relative_maha", rob["relative_maha"]["mean_suppression"], "", ""],
                  ["__stable_across_methods__", rob["stable"], "", ""]])

    # machine-readable summary (drops the big per-item lists for compactness but
    # keeps per-candidate suppression so consensus rank-agg works from files).
    summary = {
        "model": results["model"],
        "verdict": results["verdict"],
        "verdict_reason": results["verdict_reason"],
        "coverage": results["coverage"],
        "q1": {k: v for k, v in q1.items()
               if k not in ("conformal", "probe", "knn_purity", "lda_axis")},
        "q1_n_rejected": q1["n_rejected"],
        "q1_indistinguishable": q1["indistinguishable"],
        "q2_no_excess_suppression": q2["no_excess_suppression"],
        "q2_summary": {"wilcoxon": q2["wilcoxon"], "kendall": q2["kendall"],
                       "jaccard": q2["jaccard"],
                       "nearest_anchor_flips": q2["nearest_anchor_flips"]["count"],
                       "random36_p_agg": r36["p_agg"], "random36_D_taar": r36["D_taar"]},
        "robustness": rob,
        "candidate_suppression": {r["candidate_id"]: r["without"] - r["with"]
                                  for r in q2["novelty"]},
        "config": results["config"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    _write_verdict_md(out / "VERDICT.md", results)
    if whitened:
        _maybe_umap(out, whitened)


def _write_verdict_md(path: Path, results: Dict[str, object]) -> None:
    q1, q2, cov = results["q1"], results["q2"], results["coverage"]
    r36 = q2["random36"]
    lines = [
        f"# TAAR embedding-placement verdict — {results['model']}",
        "",
        f"## VERDICT: **{results['verdict']}**",
        "",
        results["verdict_reason"],
        "",
        "## Coverage (n per set)",
        f"- anchors: {cov['anchor']}  (aminergic={cov['anchor_families'].get('aminergic')})",
        f"- candidates: {cov['candidate']}",
        f"- probes: {cov['probe']}",
        f"- TAARs: {cov['taar']}",
        "",
        "## Q1 — placement (tied-cov Mahalanobis geometry)",
        f"- Conformal (answer of record): {q1['n_rejected']}/{cov['taar']} TAAR "
        f"rejected as non-aminergic at FDR {results['config']['alpha']} → "
        f"indistinguishable={q1['indistinguishable']}",
        f"- MMD TAAR-vs-aminergic: stat={q1['mmd_taar_aminergic']['stat']:.4g}, "
        f"p={q1['mmd_taar_aminergic']['p']:.3g}",
        f"- Energy TAAR-vs-aminergic: stat={q1['energy_taar_aminergic']['stat']:.4g}, "
        f"p={q1['energy_taar_aminergic']['p']:.3g}",
        f"- MMD TAAR-vs-probes: p={q1['mmd_taar_probes']['p']:.3g}; "
        f"Energy TAAR-vs-probes: p={q1['energy_taar_probes']['p']:.3g}",
        "",
        "## Q2 — add-impact",
        f"- Wilcoxon signed-rank on paired novelty deltas: p={q2['wilcoxon']['p']:.3g}",
        f"- Kendall tau-b (without vs with): {q2['kendall']['tau']:.3g}",
        f"- top-k Jaccard: " + ", ".join(f"@{k}={v:.3g}" for k, v in q2["jaccard"].items()),
        f"- Nearest-anchor localization: {q2['nearest_anchor_flips']['count']} "
        f"candidate(s) now assign to the TAAR prototype",
        f"- Random-36 null: aggregate p={r36['p_agg']:.3g} "
        f"(D_taar={r36['D_taar']:.4g}, null_mean={r36['null_mean']:.4g}); "
        f"per-candidate BH-significant excess suppression: {r36['n_bh_sig']} → "
        f"no_excess_suppression={q2['no_excess_suppression']}",
        "",
        "## Robustness",
        f"- Q2 suppression conclusion stable across "
        f"Mahalanobis/kNN/LOF: {results['robustness']['stable']}",
        f"- RelativeMahalanobis mean suppression: "
        f"{results['robustness']['relative_maha']['mean_suppression']:.4g} "
        f"(reported only; scorer unchanged)",
        "",
        "## Decision rule",
        "ADD (as own prototype) iff Q1 indistinguishable AND Q2 no excess "
        "suppression beyond the random-36 null; else HOLD.",
    ]
    path.write_text("\n".join(lines) + "\n")


def _whitened_points(embeddings, anchor_families, cfg) -> Dict[str, np.ndarray]:
    sets = split_by_prefix(embeddings)
    precision = tied_precision(sets["anchor"], anchor_families)
    L = whitening_matrix(precision)
    return {
        "anchor": _stack(sets["anchor"], list(sets["anchor"])) @ L,
        "candidate": _stack(sets["candidate"], list(sets["candidate"])) @ L,
        "probe": _stack(sets["probe"], list(sets["probe"])) @ L,
        "taar": _stack(sets["taar"], list(sets["taar"])) @ L,
    }


def run_and_write(
    emb_npz: str, anchor_labels_tsv: str, out_dir: str,
    model_name: str, cfg: Optional[Config] = None,
) -> Dict[str, object]:
    """Load an npz, join anchor families from the TSV, analyse, write outputs."""
    cfg = cfg or Config()
    embeddings = load_embeddings(emb_npz)
    if not embeddings:
        raise ValueError(f"no embeddings loaded from {emb_npz}")
    sets = split_by_prefix(embeddings)
    anchor_families = map_anchor_families(list(sets["anchor"]), anchor_labels_tsv)
    print(f"[{model_name}] coverage: "
          f"anchors={len(sets['anchor'])} candidates={len(sets['candidate'])} "
          f"probes={len(sets['probe'])} TAARs={len(sets['taar'])}", flush=True)
    results = run_per_model(embeddings, anchor_families, model_name, cfg)
    whitened = _whitened_points(embeddings, anchor_families, cfg)
    write_outputs(results, out_dir, whitened=whitened)
    print(f"[{model_name}] VERDICT: {results['verdict']}", flush=True)
    return results


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #
def _cfg_from_args(a) -> Config:
    return Config(
        k=a.k, n_perm=a.n_perm, n_random=a.n_random, knn=a.knn, alpha=a.alpha,
        seed=a.seed, split_frac=a.split_frac, min_aminergic=a.min_aminergic,
        topk=tuple(a.topk), calib=a.calib,
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--emb-npz", help="one npz of {id: vector} over all ~2395 ids")
    ap.add_argument("--anchor-labels", help="anchor_set_PROD.tsv (accession/family)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--model", default="model", help="model tag (proteinclip3b/protrek)")
    ap.add_argument("--consensus-from", nargs="+", default=None,
                    help="per-model out-dirs (each with summary.json) to combine")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--n-random", type=int, default=500)
    ap.add_argument("--knn", type=int, default=15)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--split-frac", type=float, default=0.5)
    ap.add_argument("--min-aminergic", type=int, default=4)
    ap.add_argument("--topk", type=int, nargs="+", default=[20, 50])
    ap.add_argument("--calib", default="sigmoid", choices=["sigmoid", "isotonic"])
    a = ap.parse_args(argv)

    if a.consensus_from:
        summaries = []
        for d in a.consensus_from:
            s = json.loads((Path(d) / "summary.json").read_text())
            # reconstruct a minimal results-shaped dict for run_consensus.
            supp = s.get("candidate_suppression", {})
            summaries.append({
                "model": s["model"], "verdict": s["verdict"],
                "q2": {"novelty": [{"candidate_id": c, "without": v, "with": 0.0}
                                   for c, v in supp.items()]},
            })
        con = run_consensus(summaries)
        out = Path(a.out_dir)
        out.mkdir(parents=True, exist_ok=True)
        _write_tsv(out / "consensus_verdicts.tsv", ["model", "verdict"],
                   [[m, v] for m, v in con["per_model"].items()]
                   + [["__consensus__", con["consensus_verdict"]]])
        if con["candidate_suppression_rankagg"]:
            cols = ["candidate_id", "mean_rank"] + [
                k for k in con["candidate_suppression_rankagg"][0]
                if k.startswith("suppression_")]
            _write_tsv(out / "consensus_candidate_suppression.tsv", cols,
                       [[r[c] for c in cols] for r in con["candidate_suppression_rankagg"]])
        verdict_line = (f"CONSENSUS: {con['consensus_verdict']} "
                        f"(agreement={con['agreement']})")
        (out / "VERDICT.md").write_text(
            f"# Consensus\n\n{verdict_line}\n\n"
            + "\n".join(f"- {m}: {v}" for m, v in con["per_model"].items())
            + ("\n\n**Models DISAGREE — this disagreement is itself a finding.**\n"
               if not con["agreement"] else "\n"))
        print(verdict_line, flush=True)
        return 0

    if not a.emb_npz or not a.anchor_labels:
        ap.error("--emb-npz and --anchor-labels are required for a per-model run")
    run_and_write(a.emb_npz, a.anchor_labels, a.out_dir, a.model, _cfg_from_args(a))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
