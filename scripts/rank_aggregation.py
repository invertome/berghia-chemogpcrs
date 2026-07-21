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
    CDF. That minimum, rho, is then referred to its EXACT null distribution
    (Kolde's own exact path, computed per distinct number-of-lists m -- see
    `rho_null_pvalue`), NOT to the Bonferroni bound min(rho*m, 1.0) this
    module used to return. Primary method -- yields a significance-like
    score (lower = better).
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

import math
import os
import sys
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.special import betainc, betaincinv
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


def _binomial_pmf(n, p):
    """P(Binomial(n, p) = d) for d = 0..n, as a plain list.

    scipy.stats.binom.pmf is ~15x slower per call than this at the sizes used
    here (n <= ~30) and rho_null_pvalue is called O(candidates x bootstrap
    replicates) times, so the dispatch overhead dominates.
    """
    if p <= 0.0:
        return [1.0] + [0.0] * n
    if p >= 1.0:
        return [0.0] * n + [1.0]
    q = 1.0 - p
    return [math.comb(n, d) * p ** d * q ** (n - d) for d in range(n + 1)]


@lru_cache(maxsize=65536)
def rho_null_pvalue(rho, m):
    """P(rho_null <= rho) -- the EXACT null distribution of Kolde's rho.

    rho_null is the minimum of the m order-statistic beta scores when the m
    normalized ranks really are i.i.d. Uniform(0,1), i.e. the statistic
    computed by :func:`rho_statistic` under the null "this candidate is not
    special". This is Kolde et al. (2012)'s own exact path, and it REPLACES
    the Bonferroni bound ``min(rho * m, 1.0)`` the module used to return
    (bead 8k8e).

    Why the bound had to go: it is a bound, so it CLIPS. On the real
    439-candidate cohort 204 candidates (46.5%) hit the 1.0 clip in ONE tie
    block that flattened 112 distinct underlying rho values into a single
    score, whereupon ``aggregate`` ordered them alphabetically by transcript
    id. Nearly half the "ranking" was the alphabet.

    THE NULL DEPENDS ON m, so this is computed per distinct m and memoized per
    (rho, m). m varies across candidates (measured: 7 for 431 candidates, 6 for
    8, because a signal only votes where it has data), and the same rho is less
    surprising with more lists -- scores are NOT comparable across candidates
    without conditioning on m.

    Method. rho <= x iff some order statistic U_(k) falls at or below
    a_k = the x-quantile of Beta(k, m-k+1), so the survival event is
    {U_(k) > a_k for all k}. The a_k are non-decreasing in k, so sweeping the
    thresholds in order and tracking N(a_j) = #{points <= a_j} gives an exact
    dynamic program. The probability is accumulated as the FIRST-CROSSING mass
    at each threshold -- a sum of strictly positive terms -- rather than as
    1 - survival, which would catastrophically cancel in the tail where the
    best candidates live (an id ranked top in all 7 lists has rho ~ 1e-19, and
    1 - (1 - 7e-19) is exactly 1.0 in double precision).

    Below rho ~ 1e-136 (m=6) / 1e-152 (m=7) scipy's ``betaincinv`` returns NaN.
    There the union bound m*rho IS the answer: the deficit of the exact value
    below m*rho scales as rho**(1/m), which at the crossover is already under
    1e-13 relative and shrinks from there. That branch is exact to double
    precision, not an approximation of convenience -- and it is unreachable by
    any real cohort (rho >= (1/n)**m, so m=7 would need n > 1e21 candidates).

    LOWER = more consistently top-ranked. Strictly increasing in rho on (0, 1),
    so it never ties two candidates whose evidence differs.
    """
    if m < 1:
        raise ValueError(f"rho_null_pvalue needs m >= 1, got {m!r}")
    if rho <= 0.0:
        return 0.0
    if rho >= 1.0:
        return 1.0

    ks = np.arange(1, m + 1)
    thresholds = betaincinv(ks, m - ks + 1, rho)
    if not np.all(np.isfinite(thresholds)):
        return min(m * rho, 1.0)
    a = [float(t) for t in thresholds]

    crossed = 0.0
    prev = 0.0
    # surviving[n] = P(N(a_j) = n and no threshold crossed at or before j)
    surviving = {0: 1.0}
    for j in range(1, m + 1):
        aj = max(a[j - 1], prev)          # a_k is non-decreasing; enforce it
        remaining_span = 1.0 - prev
        # each of the not-yet-placed points is uniform on (prev, 1]
        p = 1.0 if remaining_span <= 0.0 else min((aj - prev) / remaining_span, 1.0)
        nxt = {}
        for n, mass in surviving.items():
            for d, weight in enumerate(_binomial_pmf(m - n, p)):
                if weight == 0.0:
                    continue
                total = n + d
                if total >= j:            # U_(j) <= a_j: the null event occurred
                    crossed += mass * weight
                else:
                    nxt[total] = nxt.get(total, 0.0) + mass * weight
        surviving = nxt
        prev = aj
        if not surviving:
            break
    # crossed <= 1 mathematically; the min only absorbs float round-off at
    # rho -> 1, where the answer is 1 anyway. It is NOT the removed clip.
    return min(crossed, 1.0)


def rho_statistic(per_signal_ranklists, groups=None):
    """Kolde's raw rho statistic and list count per id, as ``{id: (rho, m)}``.

    For each id, gather its normalized rank from every effective list it
    appears in, sort ascending r_(1) <= ... <= r_(m), and take the Kolde
    beta-score of each order statistic: the regularized incomplete beta
    betainc(i, m-i+1, r_(i)) -- the probability the i-th of m Uniform(0,1)
    draws is <= r_(i) under the null. rho = min over i.

    Exposed separately from :func:`rra_score` because rho is not a p-value (it
    is the minimum of m dependent, marginally-uniform quantities) and because
    audits need the pre-null statistic and its m.
    """
    lists = _effective_lists(per_signal_ranklists, groups)
    ids = set()
    for lst in lists:
        ids.update(lst.keys())

    out = {}
    for i in ids:
        r = sorted(lst[i] for lst in lists if i in lst)
        m = len(r)
        out[i] = (min(betainc(k, m - k + 1, r[k - 1]) for k in range(1, m + 1)), m)
    return out


def rra_score(per_signal_ranklists, groups=None):
    """Robust Rank Aggregation (Kolde et al. 2012) score per id.

    The Kolde beta/rho step is :func:`rho_statistic`; the score is that rho
    mapped through its EXACT null distribution, :func:`rho_null_pvalue`,
    conditioned on the id's own m. LOWER = more consistently top-ranked.

    Bead 8k8e: this used to return the Bonferroni bound ``min(rho * m, 1.0)``,
    which clipped 46.5% of the real cohort into a single tie block broken by
    transcript id. The exact null is strictly monotone in rho and unclipped, so
    two candidates share a score only if their rank vectors genuinely agree --
    see :func:`rra_tied_block_size`, which reports those blocks rather than
    letting the alphabet resolve them silently.
    """
    return {i: rho_null_pvalue(rho, m)
            for i, (rho, m) in rho_statistic(per_signal_ranklists, groups).items()}


def rra_tied_block_size(per_signal_ranklists, groups=None, scores=None):
    """``{id: how many candidates share this id's exact RRA score}``.

    A block size > 1 means :func:`aggregate` cannot order those candidates on
    evidence -- it falls back to ascending id, which is arbitrary. Under the
    exact null that happens only for a GENUINE tie (identical sorted rank
    vectors and identical m), which is a real statement about the data and must
    stay a tie; what must not happen is for it to stay invisible. Stage 07
    writes this as the ``rra_tied_block_size`` column so a shortlist consumer
    can see that positions inside a block are interchangeable.

    ``scores`` reuses an already-computed :func:`rra_score` map (same
    ``groups``) instead of recomputing it.
    """
    if scores is None:
        scores = rra_score(per_signal_ranklists, groups=groups)
    counts = {}
    for value in scores.values():
        counts[value] = counts.get(value, 0) + 1
    return {i: counts[v] for i, v in scores.items()}


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
# The four formerly-ungated base signals now name their availability flag, so
# bead o98's phylo-absence rule (a candidate outside the class-A tree has no
# meaningful phylo/lse_divergence signal, and is judged on its other axes rather
# than scored a present 0.0) applies on BOTH ranking paths instead of only the
# weighted one, and dN/dS gates on whether aBSREL actually reported. Their
# flags are OPTIONAL: unlike the evidence channels, a df with no flag column
# (a legacy ranked CSV) keeps them voting rather than dropping them -- see
# _OPTIONAL_FLAG_SIGNALS in build_ranklists_from_df.
SIGNAL_SPEC = [
    ("phylo", "has_phylo_data"),
    ("purifying", "has_dnds_data"),
    ("positive", "has_dnds_data"),
    ("lse_divergence", "has_phylo_data"),
    # Bead hf3u: TWO depth axes, deliberately not fused. `lse_divergence` scores
    # cumulative branch length (divergence = rate x time); `lse_nesting_depth`
    # scores root-to-tip node count (duplication depth). Measured on this
    # repo's real trees they agree population-wide (spearman +0.897 / +0.604)
    # and barely at all inside the top quartile that alone affects ranking
    # (+0.078 / +0.034), so which one is used changes the scored set (62.7% /
    # 42.5% retained). There is no positive control to arbitrate them -- the
    # pipeline is subtractive -- so both vote and
    # audit_signal_ranking_independence.py decides empirically, via the group
    # map rank_candidates._load_signal_groups feeds back in here, whether they
    # are redundant enough to count as one. Unlike lse_divergence it gates on its
    # OWN flag: it is measured from its own population against its own
    # threshold, and is NOT in _OPTIONAL_FLAG_SIGNALS, so a CSV predating the
    # axis skips it rather than voting blind.
    ("lse_nesting_depth", "has_lse_nesting_depth_data"),
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
    # Phase-0 maha channel: S_novel as a POSITIVE voter (high = divergent from
    # known families = surface it). Dormant unless the maha channel populated
    # emb_novelty (cosine channel emits emb_nonchemo_sim/emb_classA_sim instead);
    # the two schemas are mutually dormant, so nothing changes unless maha is used.
    ("emb_novelty", "has_emb_data", "emb_novelty", False),
    ("or_microswitch", "has_or_microswitch_data", "or_microswitch", False),
]


# Signals whose has_*_data flag is honoured when the column is present but
# whose absence means "this CSV predates the flag", not "no data". Keeping the
# fallback is what stops a legacy ranked CSV from producing an EMPTY ranklist
# set (a gated signal with no flag column is skipped entirely).
_OPTIONAL_FLAG_SIGNALS = frozenset({"phylo", "purifying", "positive", "lse_divergence"})

# Explicit "exclude nothing" for RANKAGG_EXCLUDED_SIGNALS. A token is required
# because an empty string cannot be distinguished from an accidental empty
# export; no signal in SIGNAL_SPEC is named this, so it cannot collide.
_EXCLUDE_NOTHING_TOKEN = "none"


def excluded_signals_from_weights(weights, env_var="RANKAGG_EXCLUDED_SIGNALS"):
    """Signals that must NOT vote under rank aggregation.

    RRA is weight-free by design and stays so: this does NOT reinterpret
    weights as RRA inputs. It reads one specific, discrete piece of intent out
    of them -- a weight of exactly 0 is an EXCLUSION, not a small weight -- and
    turns it into an explicit named set, so the semantics are visible in the
    ranking's provenance rather than implied by a number the aggregator never
    sees.

    The motivating case is ``PURIFYING_WEIGHT=0`` (bead -ea9): whole-gene
    purifying selection is the wrong signal for chemoreceptor discovery, so the
    weighted scorer drops the axis entirely. Without this, `purifying` votes at
    full strength under the production default RANK_METHOD=rankagg.

    Policy is reversible without touching weights: setting
    ``RANKAGG_EXCLUDED_SIGNALS`` overrides the derivation verbatim, as a
    comma-separated signal list.

      unset               -> derive the exclusions from the weights
      ``""`` / whitespace -> no intent expressed; derive, and say so on stderr
      ``"none"``          -> explicitly exclude nothing
      ``"a,b"``           -> exclude exactly those

    An empty value used to mean "exclude nothing", which made this a trapdoor
    (bead wtwi): ``export RANKAGG_EXCLUDED_SIGNALS=`` yields ``""``, which is
    not None, so the override branch returned the empty set and the
    weight-derived exclusion was skipped -- re-enabling every zero-weighted
    signal at FULL strength, the exact failure this function exists to
    prevent. The trap is that the intended value and the accidental value (an
    empty export, an unset variable interpolated into a wrapper, a CI default
    that resolved to nothing) are the SAME STRING, so no care at the call site
    could distinguish them. An empty string is not a measurement of intent, so
    it no longer decides anything; "exclude nothing" now needs the explicit
    ``none`` token, which cannot be produced by accident.
    """
    override = os.getenv(env_var)
    if override is not None:
        names = {name.strip() for name in override.split(",") if name.strip()}
        if _EXCLUDE_NOTHING_TOKEN in {n.lower() for n in names}:
            if len(names) > 1:
                raise ValueError(
                    f"{env_var}={override!r} both excludes nothing and names "
                    f"signals to exclude. Use '{_EXCLUDE_NOTHING_TOKEN}' alone, "
                    f"or list only the signals to exclude."
                )
            return set()
        if names:
            return names
        # Empty / whitespace / bare separators: no intent expressed. Fall back
        # to the safe derivation rather than to the unsafe empty set, and make
        # the fallback visible so it cannot be mistaken for a deliberate policy.
        derived = {name for name, weight in weights.items() if float(weight) == 0.0}
        print(
            f"WARNING: {env_var} is set but empty, which states no policy. "
            f"Falling back to the weight-derived exclusion "
            f"{sorted(derived) or '(none)'}. Set {env_var}="
            f"'{_EXCLUDE_NOTHING_TOKEN}' if you really mean to exclude nothing.",
            file=sys.stderr,
        )
        return derived
    return {name for name, weight in weights.items() if float(weight) == 0.0}


def build_ranklists_from_df(df, id_col="id", excluded=None):
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
    dormant until merge_evidence_channels() has actually joined them in. The
    exception is _OPTIONAL_FLAG_SIGNALS, which fall back to voting when their
    flag column is missing so a legacy CSV still ranks.

    ``excluded`` names signals that must not vote at all (see
    :func:`excluded_signals_from_weights`); the default excludes nothing, so
    omitting it reproduces the previous behaviour exactly.
    """
    excluded = set(excluded or ())
    # Normalize the pre-rename `lse_depth_*` schema to `lse_divergence_*`
    # BEFORE any column lookup. This is the shared aggregation entry point
    # (rankagg, rank_confidence, compare_ranking_methods, shortlist_impact all
    # arrive here), so doing it once here is what stops a ranked CSV written
    # before the rename from quietly contributing NO lse_divergence ranklist --
    # which would read as "this candidate has no divergence signal" rather than
    # as the stale-schema problem it actually is. Announced, never silent; see
    # _rank_candidates_lib.apply_legacy_column_aliases.
    from _rank_candidates_lib import apply_legacy_column_aliases
    df = apply_legacy_column_aliases(df, source="rank-aggregation input")

    ids = list(df[id_col])
    ranklists = {}
    for spec in SIGNAL_SPEC:
        key, flag = spec[0], spec[1]
        if key in excluded:
            continue
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
        elif key in _OPTIONAL_FLAG_SIGNALS:
            # legacy CSV predating the flag column -> keep the signal voting
            keep = [True] * len(ids)
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
         ["emb_classA_sim", "emb_nonchemo_sim", "emb_novelty",
          "emb_nonchemo_family", "emb_leakage_flag"]),
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


def rerank_output(df_sorted, method, groups=None, excluded=None):
    """Return ``df_sorted`` reordered by the chosen ranking method.

    ``method != "rankagg"`` (the default ``"weighted"``): returns ``df_sorted``
    UNCHANGED -- a strict identity so the weighted production ordering/output
    stays byte-identical. ``method == "rankagg"``: reorders rows by a Robust
    Rank Aggregation fusion of the SIGNAL_SPEC per-signal ranklists (the 12
    core signals, plus any Task-6 evidence channels present; ``groups`` fuses
    confounded signals into one vote; ``excluded`` names signals barred from
    voting, e.g. a zero-weight axis -- see
    :func:`excluded_signals_from_weights`); ids covered by no signal keep their
    incoming (weighted) order at the tail via a stable sort.

    Under rankagg it also writes ``rra_tied_block_size`` (bead 8k8e): how many
    candidates share each row's exact RRA score. Anything above 1 is a genuine
    tie whose internal order is decided by ascending id, i.e. arbitrarily, and
    a shortlist consumer needs to see that rather than infer a precision the
    aggregation does not have. Ids with no voting signal get 0 -- they were
    never scored, so they are not tied with anything.
    """
    if method != "rankagg":
        return df_sorted
    ranklists = build_ranklists_from_df(df_sorted, excluded=excluded)
    scores = rra_score(ranklists, groups=groups)
    blocks = rra_tied_block_size(ranklists, groups=groups, scores=scores)
    order = sorted(scores, key=lambda i: (scores[i], i))
    pos = {idv: i for i, idv in enumerate(order)}
    tail = len(pos)
    keyed = df_sorted.assign(
        _rankagg_pos=[pos.get(i, tail) for i in df_sorted["id"]],
        rra_tied_block_size=[blocks.get(i, 0) for i in df_sorted["id"]],
    )
    reordered = keyed.sort_values("_rankagg_pos", kind="stable")
    return reordered.drop(columns="_rankagg_pos")
