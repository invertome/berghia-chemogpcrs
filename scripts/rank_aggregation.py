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

import os

import numpy as np
import pandas as pd
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


# --------------------------------------------------------------------------- #
# Bridge to the candidate-ranking dataframe (shared by rank_candidates.py and
# compare_ranking_methods.py so the 12-signal spec lives in exactly one place)
# --------------------------------------------------------------------------- #
# The 12 core ranking signals as (signal_key, has_flag), PLUS (Task 6) five
# optional structural/embedding/microswitch evidence-channel signals appended
# below. has_flag=None marks the four base signals that are ALWAYS present;
# every other entry contributes a candidate's value only where its
# has_*_data boolean is True. The two non-obvious flag names (gprotein/ecl)
# mirror the production ranked CSV and
# scripts/audit_signal_ranking_independence.py's FLAG_OVERRIDES.
#
# Each entry is (key, flag) or (key, flag, column, invert):
#   - column: the exact df column to read. Omitted/None (the 12 core
#     signals) means "look up f'{key}_score_norm' then fall back to
#     f'{key}_score'" (rank_candidates.py's CSV convention). The Task-6
#     evidence channels don't follow that convention (they are not
#     normalized/weighted scores), so they name their own literal column.
#   - invert: True marks an EXCLUSION signal. Rank aggregation has no
#     weights to negate, so the ONLY way an exclusion signal can lower a
#     candidate's standing is to flip its polarity before it ever reaches
#     normalized_ranks(): the stored ranklist value becomes (1 - raw), so
#     "higher stored value = better" holds for every signal, exclusion or
#     not. See merge_evidence_channels() for how the Task-4/5/6 channel
#     outputs get joined onto the ranking df in the first place.
SIGNAL_SPEC = [
    ("phylo", None),
    ("purifying", None),
    ("positive", None),
    ("lse_depth", None),
    ("synteny", "has_synteny_data"),
    ("expression", "has_expression_data"),
    ("chemosensory_expr", "has_chemosensory_expr_data"),
    ("gprotein_coexpr", "has_gprotein_data"),
    ("ecl_divergence", "has_ecl_data"),
    ("expansion", "has_expansion_data"),
    ("og_confidence", "has_og_confidence_data"),
    ("tandem_cluster", "has_tandem_cluster_data"),
    # --- Task 6: optional structural/embedding/microswitch evidence
    # channels. Dormant by default (merge_evidence_channels() never ran /
    # channel columns absent / has_*_data all False) -- the 12 signals above
    # alone then reproduce the pre-Task-6 ranking exactly (see
    # tests/unit/test_channel_integration.py's regression tests).
    ("struct_novelty", "has_struct_data", "struct_novelty", False),
    ("struct_nonchemo_corrob", "has_struct_data", "struct_nonchemo_corrob", True),
    ("emb_classA_sim", "has_emb_data", "emb_classA_sim", False),
    ("emb_nonchemo_sim", "has_emb_data", "emb_nonchemo_sim", True),
    ("or_microswitch", "has_or_microswitch_data", "or_microswitch", False),
]


def build_ranklists_from_df(df, id_col="id"):
    """Build ``{signal: {id: score}}`` for every signal in SIGNAL_SPEC from a df.

    The 12 core ranking signals prefer the normalized ``<signal>_score_norm``
    column and fall back to the raw ``<signal>_score`` column (the written
    ranked CSV carries the raw columns, not the norm ones -- rank aggregation
    ranks within each signal, which is invariant to the monotonic min-max
    normalization, so either yields identical rankings). The Task-6
    evidence-channel signals instead read their own exact ``column`` named in
    SIGNAL_SPEC, since they are not normalized/weighted scores.

    A gated signal (``flag`` not None) contributes a candidate's value only
    where its ``has_*_data`` flag is True; the four base signals are always
    included. Higher stored value = better for every signal. An EXCLUSION
    signal (``invert=True`` in SIGNAL_SPEC) stores ``1 - raw_value`` instead
    of the raw value, so that invariant holds for it too -- rank aggregation
    has no weights, so inversion is the only way an exclusion signal can
    lower, never raise, a candidate's standing. Missing/NaN values and empty
    signals are dropped so a signal only "votes" where it has data; a signal
    whose column (or, for a gated signal, whose flag column) is entirely
    absent from ``df`` is skipped -- this is what keeps the Task-6 channels
    dormant until merge_evidence_channels() has actually joined them in.
    """
    ids = list(df[id_col])
    ranklists = {}
    for spec in SIGNAL_SPEC:
        key, flag = spec[0], spec[1]
        column = spec[2] if len(spec) > 2 else None
        invert = spec[3] if len(spec) > 3 else False
        if column is not None:
            col = column
        else:
            norm_col, raw_col = f"{key}_score_norm", f"{key}_score"
            col = norm_col if norm_col in df.columns else raw_col
        if col not in df.columns:
            continue
        values = list(df[col])
        if flag is None:
            keep = [True] * len(ids)
        elif flag in df.columns:
            keep = [bool(v) for v in df[flag]]
        else:
            # gated signal whose flag column is absent -> treat as no data
            continue
        signal_map = {}
        for idv, val, present in zip(ids, values, keep):
            if not present or _is_missing(val):
                continue
            val = float(val)
            if invert:
                val = 1.0 - val
            signal_map[idv] = val
        if signal_map:
            ranklists[key] = signal_map
    return ranklists


def merge_evidence_channels(df, struct_tsv=None, emb_tsv=None,
                             microswitch_tsv=None, id_col="id"):
    """Left-join the Task-4/5/6 evidence-channel outputs onto ``df`` by id.

    Each ``*_tsv`` path is a tab-separated table of an already-computed
    channel's per-candidate columns, keyed by ``id_col``:
        struct_tsv      -- id, struct_novelty, struct_nonchemo_corrob[, ...]
                           (structural_evidence.structural_channel() output)
        emb_tsv         -- id, emb_classA_sim, emb_nonchemo_sim[, ...]
                           (embedding_evidence.embedding_channel() output)
        microswitch_tsv -- id, or_microswitch
                           (or_microswitch.or_microswitch_flag() output)

    A ``None`` path, or one that doesn't exist on disk, leaves that channel
    entirely DORMANT: its ``has_*_data`` flag is False for every row and none
    of its score columns are added to the output -- the channel simply never
    ran (or hasn't been produced yet). This mirrors parse_foldseek() /
    load_embeddings() treating a missing input file as "no evidence", never
    an error.

    Sets ``has_struct_data`` / ``has_emb_data`` / ``has_or_microswitch_data``
    to True ONLY for the ``df`` rows whose id was found in the corresponding
    TSV. This is a LEFT join: every row of ``df`` is preserved regardless of
    channel coverage, and a row absent from a TSV keeps its score column(s)
    as NaN (dropped downstream by build_ranklists_from_df's missing-value
    handling) with that channel's flag False.
    """
    merged = df.copy()
    channels = (
        (struct_tsv, "has_struct_data",
         ["struct_novelty", "struct_nonchemo_corrob", "struct_state"]),
        (emb_tsv, "has_emb_data",
         ["emb_classA_sim", "emb_nonchemo_sim", "emb_nonchemo_family",
          "emb_leakage_flag"]),
        (microswitch_tsv, "has_or_microswitch_data", ["or_microswitch"]),
    )
    for tsv_path, flag_col, value_cols in channels:
        if tsv_path is None or not os.path.exists(tsv_path):
            merged[flag_col] = False
            continue
        chan = pd.read_csv(tsv_path, sep="\t").set_index(id_col)
        for col in value_cols:
            if col in chan.columns:
                merged[col] = merged[id_col].map(chan[col])
        merged[flag_col] = merged[id_col].isin(chan.index)
    return merged


def normalize_group_names(groups):
    """Translate audit group names (``<signal>_score``) to the signal keys
    used by :func:`build_ranklists_from_df` (``<signal>``).

    scripts/audit_signal_ranking_independence.py emits groups keyed by the raw
    ``*_score`` column names; the ranklists here are keyed by the bare signal
    key. Stripping the ``_score`` suffix bridges the two (all 12 map
    one-to-one). Idempotent: names without the suffix pass through unchanged.
    """
    def strip(name):
        return name[: -len("_score")] if name.endswith("_score") else name
    return [[strip(n) for n in group] for group in groups]


def rerank_output(df_sorted, method, groups=None):
    """Return ``df_sorted`` reordered by the chosen ranking method.

    ``method != "rankagg"`` (the default ``"weighted"``): returns ``df_sorted``
    UNCHANGED -- a strict identity so the weighted production ordering/output
    stays byte-identical. ``method == "rankagg"``: reorders rows by a Robust
    Rank Aggregation fusion of the SIGNAL_SPEC per-signal ranklists (the 12
    core signals, plus any Task-6 evidence channels present; ``groups`` fuses
    confounded signals into one vote); ids covered by no signal keep their
    incoming (weighted) order at the tail via a stable sort.
    """
    if method != "rankagg":
        return df_sorted
    order = aggregate(build_ranklists_from_df(df_sorted), method="rra", groups=groups)
    pos = {idv: i for i, idv in enumerate(order)}
    tail = len(pos)
    keyed = df_sorted.assign(_rankagg_pos=[pos.get(i, tail) for i in df_sorted["id"]])
    reordered = keyed.sort_values("_rankagg_pos", kind="stable")
    return reordered.drop(columns="_rankagg_pos")
