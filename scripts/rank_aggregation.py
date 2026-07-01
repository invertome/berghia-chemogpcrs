"""Label-free rank aggregation over per-signal candidate rankings.

The council rejected trained/label-based candidate ranking: novel
lineage-specific-expansion (LSE) chemoreceptors are anti-correlated with
resemblance-to-known, and there are ~0 in-species positives to train or
validate against. The robust core is rank aggregation over the existing
per-signal rank-lists -- no weights, no labels.

Two methods, both label-free and parameter-light:
  - RRA (Robust Rank Aggregation, Kolde et al. 2012): treats each signal's
    normalized rank for a candidate as a draw from Uniform(0,1) under the
    null "this candidate is not special", and looks for the most extreme
    (smallest) order statistic across signals via the beta-distribution
    CDF. Primary method -- yields a significance-like score (lower =
    better).
  - RRF (Reciprocal Rank Fusion, Cormack et al. 2009): a simpler
    Sigma 1/(k+rank) sum. Cross-check (higher = better).

`scripts/audit_signal_ranking_independence.py` (Task 0) flags signal groups
that share a confound (e.g. phylo_score/og_confidence_score both derive
from the same OrthoFinder orthogroup/gene tree). Both methods here accept
those groups and FUSE each group into one list (mean normalized rank
across the group's present members) before scoring, so a shared confound
casts one vote, not one vote per redundant signal.
"""
from __future__ import annotations

import numpy as np
from scipy.special import betainc
from scipy.stats import rankdata


def _is_missing(value):
    if value is None:
        return True
    return value != value  # true only for NaN


def normalized_ranks(scores, higher_is_better=True):
    """Convert one signal's raw scores into normalized ranks in (0, 1].

    NaN/None entries are dropped -- a signal only "votes" where it has
    data. The best-scoring id gets the smallest normalized rank (1/n);
    the worst gets 1.0. Ties share the average rank.
    """
    present = {k: float(v) for k, v in scores.items() if not _is_missing(v)}
    n = len(present)
    if n == 0:
        return {}
    ids = list(present.keys())
    values = np.array([present[i] for i in ids], dtype=float)
    ranks = rankdata(-values if higher_is_better else values, method="average")
    return {i: r / n for i, r in zip(ids, ranks)}


def _fuse_group(member_lists):
    """Collapse several normalized-rank dicts into one by averaging each
    id's normalized rank across the members where it is present."""
    ids = set()
    for m in member_lists:
        ids.update(m.keys())
    return {i: sum(m[i] for m in member_lists if i in m) /
               len([m for m in member_lists if i in m])
            for i in ids}


def _effective_lists(per_signal_ranklists, groups=None):
    """Normalize every signal, then fuse per `groups` into effective lists.

    Each returned dict maps id -> normalized rank in (0, 1] (lower =
    better). A group's member signals collapse into ONE effective list.
    Signals not named in any group remain their own singleton list.
    """
    normalized = {name: normalized_ranks(scores)
                  for name, scores in per_signal_ranklists.items()}

    groups = [list(g) for g in groups] if groups else []
    grouped_names = {name for g in groups for name in g}
    for name in per_signal_ranklists:
        if name not in grouped_names:
            groups.append([name])

    effective = []
    for group in groups:
        members = [normalized[name] for name in group if name in normalized]
        if members:
            effective.append(_fuse_group(members))
    return effective


def rra_score(per_signal_ranklists, groups=None):
    """Robust Rank Aggregation (Kolde et al. 2012) score per id.

    For each id, gather its normalized rank from every effective list it
    appears in, sort ascending r_(1) <= ... <= r_(m), and take the Kolde
    beta-score of each order statistic: the regularized incomplete beta
    betainc(i, m-i+1, r_(i)) -- the probability the i-th of m Uniform(0,1)
    draws is <= r_(i) under the null. rho = min over i. The returned score
    applies the Bonferroni-style correction min(rho * m, 1.0). LOWER =
    more consistently top-ranked.
    """
    lists = _effective_lists(per_signal_ranklists, groups)
    ids = set()
    for lst in lists:
        ids.update(lst.keys())

    scores = {}
    for i in ids:
        r = sorted(lst[i] for lst in lists if i in lst)
        m = len(r)
        betas = [betainc(k, m - k + 1, r[k - 1]) for k in range(1, m + 1)]
        rho = min(betas)
        scores[i] = min(rho * m, 1.0)
    return scores


def rrf_score(per_signal_ranklists, k=60, groups=None):
    """Reciprocal Rank Fusion (Cormack et al. 2009) score per id.

    Sigma 1/(k + rank) over the effective lists an id appears in, where
    rank is the 1-based rank within that list (best = 1). HIGHER = better.
    """
    lists = _effective_lists(per_signal_ranklists, groups)
    scores = {}
    for lst in lists:
        ids = list(lst.keys())
        if not ids:
            continue
        values = np.array([lst[i] for i in ids], dtype=float)
        # smaller normalized rank -> better -> rank 1
        ranks = rankdata(values, method="average")
        for i, rnk in zip(ids, ranks):
            scores[i] = scores.get(i, 0.0) + 1.0 / (k + rnk)
    return scores


def aggregate(per_signal_ranklists, method="rra", groups=None):
    """Fuse per-signal rank-lists into one ordered id list, best first.

    method="rra" sorts ascending rho (lower = better); "rrf" sorts
    descending score (higher = better). Ties break by ascending id.
    """
    if method == "rra":
        scores = rra_score(per_signal_ranklists, groups=groups)
        return sorted(scores, key=lambda i: (scores[i], i))
    if method == "rrf":
        scores = rrf_score(per_signal_ranklists, k=60, groups=groups)
        return sorted(scores, key=lambda i: (-scores[i], i))
    raise ValueError(f"unknown method: {method!r} (expected 'rra' or 'rrf')")
