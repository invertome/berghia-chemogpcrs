import math
import re

import pandas as pd
import pytest
from scipy.special import betainc

import rank_aggregation as ra


def test_normalized_ranks_drop_missing():
    nr = ra.normalized_ranks({"a": 0.9, "b": float("nan"), "c": 0.1})
    assert set(nr) == {"a", "c"}
    assert nr["a"] < nr["c"]  # a (highest raw score) -> smallest normalized rank


def test_normalized_ranks_lower_is_better():
    nr = ra.normalized_ranks({"a": 0.9, "b": 0.1}, higher_is_better=False)
    assert nr["b"] < nr["a"]  # smaller raw value wins when higher_is_better=False


def test_normalized_ranks_empty_and_singleton():
    assert ra.normalized_ranks({}) == {}
    assert ra.normalized_ranks({"a": 0.5}) == {"a": 1.0}


def test_rra_promotes_consistent_top():
    # 'a' is #1 in every signal; 'z' is last everywhere -> a must rank first.
    lists = {s: {"a": 0.99, "m": 0.5, "z": 0.01} for s in ["s1", "s2", "s3"]}
    order = ra.aggregate(lists, method="rra", groups=None)
    assert order[0] == "a" and order[-1] == "z"


def test_rra_matches_hand_computed_beta():
    # The Kolde et al. (2012) BETA step is unchanged by bead 8k8e: 'a' is
    # ranked #1-of-3 in three identical signals -> its normalized rank is 1/3
    # in all three; the minimum order-statistic beta occurs at the top order
    # statistic (i=m=3): betainc(3, 1, 1/3) = (1/3)**3 = 1/27.
    lists = {s: {"a": 0.99, "m": 0.5, "z": 0.01} for s in ["s1", "s2", "s3"]}
    m, r = 3, 1.0 / 3.0
    rho = min(betainc(k, m - k + 1, r) for k in range(1, m + 1))
    assert math.isclose(rho, 1.0 / 27.0, rel_tol=1e-9)

    stat = ra.rho_statistic(lists)
    assert stat["a"][1] == m
    assert math.isclose(stat["a"][0], 1.0 / 27.0, rel_tol=1e-9)

    # What CHANGED: the reported score is no longer the Bonferroni bound
    # min(rho*m, 1.0) = 1/9, it is the EXACT null P(rho_null <= rho) -- which
    # is strictly tighter than the bound for every 0 < rho < 1.
    scores = ra.rra_score(lists)
    assert math.isclose(scores["a"], ra.rho_null_pvalue(1.0 / 27.0, 3),
                        rel_tol=1e-12)
    assert scores["a"] < min(rho * m, 1.0)
    assert not math.isclose(scores["a"], 1.0 / 9.0, rel_tol=1e-3)


def test_rrf_matches_formula():
    lists = {"s1": {"a": 0.9, "b": 0.1}, "s2": {"a": 0.1, "b": 0.9}}
    sc = ra.rrf_score(lists, k=60)
    # a: rank1 in s1 (1/61) + rank2 in s2 (1/62); b symmetric -> tie
    assert math.isclose(sc["a"], 1.0 / 61.0 + 1.0 / 62.0, rel_tol=1e-9)
    assert math.isclose(sc["a"], sc["b"], rel_tol=1e-9)


def test_group_fusion_counts_confound_once():
    # 3 identical clone signals favor a>b>c; 1 independent signal favors
    # b>c>a. Ungrouped, the 3 clone votes outvote the 1 independent vote,
    # so 'a' wins. Grouped, the clones fuse into a SINGLE vote, making it
    # a fair 1-vs-1 fight -- 'b' wins because the independent signal ranks
    # 'a' dead last, while the clone group only ranks 'b' 2nd (not last).
    clones = {f"c{i}": {"a": 0.9, "b": 0.5, "c": 0.1} for i in range(3)}
    indep = {"ind": {"a": 0.1, "b": 0.9, "c": 0.5}}
    lists = {**clones, **indep}

    grouped = ra.aggregate(lists, method="rrf", groups=[["c0", "c1", "c2"], ["ind"]])
    ungrouped = ra.aggregate(lists, method="rrf", groups=None)

    assert grouped[0] == "b"     # 1 fused vote each -> independent's pick wins
    assert ungrouped[0] == "a"   # 3 raw clone votes outvote 1 independent vote


def test_partial_presence_rrf_matches_closed_form():
    # 'b' is absent from s2 -- it must not be treated as "worst" there,
    # it simply doesn't receive a vote from s2.
    lists = {"s1": {"a": 0.9, "b": 0.1}, "s2": {"a": 0.8}}
    rrf = ra.rrf_score(lists, k=60)
    assert math.isclose(rrf["a"], 2.0 / 61.0, rel_tol=1e-9)   # rank 1 in both signals
    assert math.isclose(rrf["b"], 1.0 / 62.0, rel_tol=1e-9)   # rank 2 in s1 only


def test_partial_presence_rra_invariant_to_unrelated_signal():
    # adding a signal that never mentions 'b' must not change b's score.
    only_a_and_b = {"s1": {"a": 0.9, "b": 0.1}}
    with_extra_signal = {"s1": {"a": 0.9, "b": 0.1}, "s2": {"a": 0.8, "x": 0.2}}
    base = ra.rra_score(only_a_and_b)
    extended = ra.rra_score(with_extra_signal)
    assert math.isclose(base["b"], extended["b"], rel_tol=1e-9)


def test_aggregate_dispatches_rra_ascending_and_rrf_descending():
    lists = {s: {"a": 0.99, "m": 0.5, "z": 0.01} for s in ["s1", "s2", "s3"]}
    assert ra.aggregate(lists, method="rra", groups=None) == ["a", "m", "z"]
    assert ra.aggregate(lists, method="rrf", groups=None) == ["a", "m", "z"]


def test_aggregate_tie_break_by_id():
    # Bead 8k8e requirement 3: a GENUINE tie (identical rank vectors) stays a
    # tie. aggregate() must still be DETERMINISTIC, so the arbitrary order
    # inside the block is fixed by id -- but the tie is no longer silent: the
    # block size is reported alongside (see the tie-block tests below).
    lists = {"s1": {"b": 0.5, "a": 0.5}}
    assert ra.aggregate(lists, method="rra", groups=None) == ["a", "b"]
    assert ra.aggregate(lists, method="rrf", groups=None) == ["a", "b"]
    assert ra.rra_tied_block_size(lists) == {"a": 2, "b": 2}


def test_aggregate_invalid_method_raises():
    with pytest.raises(ValueError):
        ra.aggregate({"s1": {"a": 0.5}}, method="bogus", groups=None)


# --------------------------------------------------------------------------- #
# build_ranklists_from_df: the single shared 12-signal spec
# --------------------------------------------------------------------------- #
def test_build_ranklists_base_always_present_gated_by_flag():
    df = pd.DataFrame(
        {
            "id": ["x", "y", "z"],
            "phylo_score_norm": [0.9, 0.5, 0.1],
            "purifying_score_norm": [0.1, 0.5, 0.9],
            "positive_score_norm": [0.2, 0.4, 0.6],
            "lse_divergence_score_norm": [0.3, 0.3, 0.3],
            "synteny_score_norm": [0.8, 0.2, 0.5],
            "has_synteny_data": [True, False, True],   # y gated OUT of synteny
            "gprotein_coexpr_score_norm": [0.7, 0.6, 0.5],
            "has_gprotein_data": [False, False, False],  # signal fully absent
        }
    )
    rl = ra.build_ranklists_from_df(df)
    # 4 base signals always present for every id
    for base in ("phylo", "purifying", "positive", "lse_divergence"):
        assert set(rl[base]) == {"x", "y", "z"}
    # gated signal: only ids whose has_*_data is True
    assert set(rl["synteny"]) == {"x", "z"}
    # a gated signal with NO True flags drops out entirely
    assert "gprotein_coexpr" not in rl
    # higher norm score = better is preserved as the stored value
    assert rl["phylo"]["x"] > rl["phylo"]["z"]


def test_build_ranklists_falls_back_to_raw_score_columns():
    # The written ranked CSV carries raw '<signal>_score' columns, not the
    # in-memory '<signal>_score_norm' ones; the builder must handle both.
    df = pd.DataFrame(
        {
            "id": ["x", "y"],
            "phylo_score": [0.9, 0.1],
            "purifying_score": [0.1, 0.9],
            "positive_score": [0.5, 0.5],
            "lse_divergence_score": [0.2, 0.8],
        }
    )
    rl = ra.build_ranklists_from_df(df)
    assert rl["phylo"] == {"x": 0.9, "y": 0.1}


def test_build_ranklists_uses_flag_overrides_for_gprotein_and_ecl():
    df = pd.DataFrame(
        {
            "id": ["x", "y"],
            "phylo_score_norm": [0.9, 0.1],
            "purifying_score_norm": [0.1, 0.9],
            "positive_score_norm": [0.5, 0.5],
            "lse_divergence_score_norm": [0.2, 0.8],
            "ecl_divergence_score_norm": [0.4, 0.6],
            "has_ecl_data": [True, True],            # NOT has_ecl_divergence_data
            "gprotein_coexpr_score_norm": [0.3, 0.7],
            "has_gprotein_data": [True, False],      # NOT has_gprotein_coexpr_data
        }
    )
    rl = ra.build_ranklists_from_df(df)
    assert set(rl["ecl_divergence"]) == {"x", "y"}
    assert set(rl["gprotein_coexpr"]) == {"x"}


def test_normalize_group_names_strips_score_suffix_idempotently():
    groups = [["phylo_score", "og_confidence_score"], ["synteny_score"]]
    assert ra.normalize_group_names(groups) == [["phylo", "og_confidence"], ["synteny"]]
    # already-stripped names pass through unchanged
    assert ra.normalize_group_names([["phylo"]]) == [["phylo"]]


# --------------------------------------------------------------------------- #
# Bead 8k8e: the exact null distribution of rho replaces min(rho*m, 1.0)
#
# The Bonferroni bound saturated 46.5% of the real 439-candidate cohort at
# exactly 1.0 -- ONE tie block holding 112 DISTINCT underlying rho values,
# ordered alphabetically by transcript id. The replacement is Kolde et al.
# (2012)'s own exact path: P(rho_null <= rho_obs), where rho_null is the
# minimum of the m order-statistic beta scores under m i.i.d. Uniform(0,1)
# draws. It is exact, strictly monotone in rho, never clipped, and -- because
# the null depends on m -- computed per distinct m.
# --------------------------------------------------------------------------- #
def _rho_of(sorted_ranks):
    """The Kolde rho for one candidate's sorted normalized ranks (reference
    implementation written straight from the paper, independent of the module)."""
    m = len(sorted_ranks)
    return min(betainc(k, m - k + 1, sorted_ranks[k - 1]) for k in range(1, m + 1))


def _exact_null_m2_closed_form(x):
    """P(rho_null <= x) for m = 2, derived by hand.

    For m = 2 the beta scores are 1-(1-U_(1))^2 and U_(2)^2, so rho <= x iff
    U_(1) <= a1 = 1-sqrt(1-x) OR U_(2) <= a2 = sqrt(x). The survival
    probability is P(both points > a1) minus P(both points in (a1, a2]) =
    (1-a1)^2 - (a2-a1)^2, hence

        P(rho <= x) = 1 - (1-x) + (sqrt(x) + sqrt(1-x) - 1)^2.
    """
    return x + (math.sqrt(x) + math.sqrt(1.0 - x) - 1.0) ** 2


def test_exact_null_matches_a_hand_derived_closed_form_for_m2():
    for x in (1e-8, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.9, 0.999):
        assert math.isclose(ra.rho_null_pvalue(x, 2),
                            _exact_null_m2_closed_form(x), rel_tol=1e-10), x


def test_exact_null_is_the_identity_for_a_single_list():
    # With m = 1 the only beta score is betainc(1, 1, U) = U ~ Uniform(0,1),
    # so rho IS a p-value already and the transform must be the identity.
    for x in (0.01, 0.1, 0.5, 0.9):
        assert math.isclose(ra.rho_null_pvalue(x, 1), x, rel_tol=1e-12)


def test_exact_null_is_uniform_under_the_null():
    """The defining property: feeding null-generated rho through the CDF must
    return Uniform(0,1). Checked as a KS test plus the two tail rates that
    actually matter for a ranking (5% and 1%)."""
    import numpy as np
    from scipy.stats import kstest

    rng = np.random.default_rng(20260721)
    for m in (2, 6, 7):
        u = np.sort(rng.random((40000, m)), axis=1)
        rho = np.column_stack(
            [betainc(k, m - k + 1, u[:, k - 1]) for k in range(1, m + 1)]
        ).min(axis=1)
        p = np.array([ra.rho_null_pvalue(float(v), m) for v in rho])
        assert kstest(p, "uniform").pvalue > 0.01, f"m={m} null is not uniform"
        assert abs((p <= 0.05).mean() - 0.05) < 0.005, m
        assert abs((p <= 0.01).mean() - 0.01) < 0.003, m


def test_exact_null_depends_on_m_so_scores_are_comparable_across_candidates():
    # The whole point of computing the null PER m: the same rho means
    # something different with 6 lists than with 7, and the real cohort mixes
    # both (m=7 for 431 candidates, m=6 for 8). More lists -> more chances for
    # some order statistic to be extreme -> the same rho is LESS surprising.
    rho = 0.02
    p6, p7 = ra.rho_null_pvalue(rho, 6), ra.rho_null_pvalue(rho, 7)
    assert p6 != p7
    assert p6 < p7
    assert [ra.rho_null_pvalue(rho, m) for m in range(1, 9)] == sorted(
        ra.rho_null_pvalue(rho, m) for m in range(1, 9)
    )


def test_exact_null_is_cached_per_distinct_m():
    ra.rho_null_pvalue.cache_clear()
    for _ in range(50):
        ra.rho_null_pvalue(0.03, 7)
        ra.rho_null_pvalue(0.03, 6)
    info = ra.rho_null_pvalue.cache_info()
    assert info.currsize == 2          # one cached null per (rho, m) seen
    assert info.misses == 2            # computed once per distinct m
    assert info.hits == 98


def test_exact_null_is_strictly_monotone_and_never_clipped():
    import numpy as np

    for m in (2, 6, 7):
        xs = np.unique(np.concatenate(
            [np.logspace(-18, -2, 500), np.linspace(0.01, 0.9999, 500)]))
        ps = np.array([ra.rho_null_pvalue(float(x), m) for x in xs])
        assert np.all(np.diff(ps) > 0), f"m={m} not strictly increasing"
        assert ps.max() < 1.0, f"m={m} clipped below rho=1"
    # The endpoints are the only saturated values, and they are exact.
    assert ra.rho_null_pvalue(0.0, 7) == 0.0
    assert ra.rho_null_pvalue(1.0, 7) == 1.0


def test_exact_null_is_bounded_above_by_the_bonferroni_bound_it_replaces():
    # min(rho*m, 1) is a valid but loose bound on the same probability; the
    # exact value must sit strictly below it everywhere in (0, 1).
    for m in (2, 6, 7, 13):
        for x in (1e-12, 1e-6, 1e-3, 0.01, 0.1, 0.4, 0.9):
            assert ra.rho_null_pvalue(x, m) < min(x * m, 1.0), (m, x)


def test_exact_null_survives_the_extreme_tail_without_going_non_finite():
    # scipy's betaincinv returns NaN below ~1e-136 (m=6) / ~1e-152 (m=7). In
    # that regime the union bound m*rho is the exact answer to better than
    # 1e-13 relative (the deficit scales as rho**(1/m)), so the tail branch is
    # exact-to-double-precision rather than a fallback guess. Verified here at
    # the crossover, where BOTH branches are computable.
    for m in (6, 7):
        near = 10.0 ** -(130 if m == 6 else 145)
        assert math.isclose(ra.rho_null_pvalue(near, m), m * near, rel_tol=1e-12)
        deep = 1e-200
        p = ra.rho_null_pvalue(deep, m)
        assert math.isfinite(p) and p > 0.0
        assert math.isclose(p, m * deep, rel_tol=1e-12)


def test_exact_null_resolves_scores_the_bonferroni_clip_flattened():
    """The regression that defines the bead.

    Ten candidates whose rho values are all DISTINCT but all large enough that
    rho*m >= 1 -- under the old formula every one of them scored exactly 1.0
    and the 'ranking' among them was alphabetical. The exact null must give
    ten distinct, correctly ordered scores.
    """
    # Twelve candidates ranked by two disagreeing signals. With m=2 the cap
    # binds for every rho >= 0.5, which here is the bottom half of the list.
    names = list("abcdefghijkl")
    n = len(names)
    shuffled = [1, 0, 8, 2, 10, 9, 7, 6, 4, 3, 5, 11]
    lists = {
        "s1": {name: float(n - i) for i, name in enumerate(names)},
        "s2": {name: float(n - shuffled[i]) for i, name in enumerate(names)},
    }
    stat = ra.rho_statistic(lists)
    capped = [i for i, (rho, m) in stat.items() if min(rho * m, 1.0) >= 1.0]
    assert len(capped) >= 4, "fixture does not reproduce the saturation"
    assert len({stat[i][0] for i in capped}) == len(capped), \
        "fixture must give the capped candidates DISTINCT underlying rho"

    scores = ra.rra_score(lists)
    assert len({scores[i] for i in capped}) == len(capped), \
        "exact null must not collapse distinct rho into one score"
    # and the resolved order must follow rho, not the alphabet
    by_score = sorted(capped, key=lambda i: scores[i])
    by_rho = sorted(capped, key=lambda i: stat[i][0])
    assert by_score == by_rho


# --------------------------------------------------------------------------- #
# Bead 8k8e requirement 3: genuine ties survive, and are SURFACED
# --------------------------------------------------------------------------- #
def test_genuine_ties_stay_tied_and_report_their_block_size():
    # x and y have identical rank vectors in every signal -> a genuine tie
    # that no amount of exactness can or should break. z is distinct.
    lists = {
        "s1": {"x": 0.5, "y": 0.5, "z": 0.9},
        "s2": {"x": 0.5, "y": 0.5, "z": 0.9},
    }
    scores = ra.rra_score(lists)
    assert scores["x"] == scores["y"]
    assert scores["z"] != scores["x"]
    assert ra.rra_tied_block_size(lists) == {"x": 2, "y": 2, "z": 1}


def test_ties_from_a_shared_rho_are_real_and_are_reported_not_broken():
    """rho is a MINIMUM over order statistics, so two candidates can share it
    while their full rank vectors differ -- they simply agree on the single
    most extreme coordinate, and RRA looks at nothing else.

    This is a property of Kolde's statistic, not of the null transform: no
    exactness can separate them, and inventing a separation would be inventing
    evidence. On the real 439-candidate cohort every residual tie block is of
    this kind (23 blocks, zero of which collapse distinct rho). So the
    contract is: stay tied, and SAY SO.
    """
    # Five candidates, two signals. 'p' is 1st of 5 in s1 and 4th in s2 ->
    # ranks (0.2, 0.8); 'q' is 1st of 5 in s2 and last in s1 -> (0.2, 1.0).
    # Different vectors, but both minimise at the FIRST order statistic with
    # the same value 0.2, so both get rho = 1-(1-0.2)^2 = 0.36.
    lists = {
        "s1": {"p": 5.0, "r": 4.0, "s": 3.0, "t": 2.0, "q": 1.0},
        "s2": {"q": 5.0, "r": 4.0, "s": 3.0, "p": 2.0, "t": 1.0},
    }
    stat = ra.rho_statistic(lists)
    ranks = ra._effective_lists(lists, None)
    vec = {i: tuple(sorted(lst[i] for lst in ranks if i in lst)) for i in stat}
    assert vec["p"] != vec["q"], "fixture must give the two DIFFERENT rank vectors"
    assert stat["p"][0] == stat["q"][0], "fixture must give them the SAME rho"

    scores = ra.rra_score(lists)
    assert scores["p"] == scores["q"]
    # 's' (ranks 0.6, 0.6) lands on the same rho from a third distinct vector,
    # so the reported block is all three -- exactly the fact a shortlist needs.
    blocks = ra.rra_tied_block_size(lists)
    assert blocks["p"] == blocks["q"] == blocks["s"] == 3
    assert scores["r"] != scores["p"] and blocks["r"] == 1


def test_tied_block_size_is_one_for_every_candidate_when_nothing_ties():
    lists = {"s1": {"a": 0.9, "b": 0.5, "c": 0.1}}
    assert set(ra.rra_tied_block_size(lists).values()) == {1}


def test_tied_block_size_honours_groups_like_the_scorer():
    lists = {
        "s1": {"a": 0.9, "b": 0.5},
        "s1_copy": {"a": 0.9, "b": 0.5},
        "s2": {"a": 0.5, "b": 0.9},
    }
    groups = [["s1", "s1_copy"], ["s2"]]
    scores = ra.rra_score(lists, groups=groups)
    blocks = ra.rra_tied_block_size(lists, groups=groups)
    assert set(blocks) == set(scores)
    # fusing the clones makes a and b symmetric -> a genuine 2-way tie
    assert scores["a"] == scores["b"]
    assert blocks == {"a": 2, "b": 2}


def test_rerank_output_surfaces_the_tied_block_size_column():
    df = pd.DataFrame(
        {
            "id": ["x", "y", "z"],
            # x and y are identical on every voting signal -> genuine tie
            "phylo_score": [0.5, 0.5, 0.9],
            "purifying_score": [0.5, 0.5, 0.9],
            "positive_score": [0.5, 0.5, 0.9],
            "lse_divergence_score": [0.5, 0.5, 0.9],
        }
    )
    out = ra.rerank_output(df, "rankagg")
    assert "rra_tied_block_size" in out.columns
    blocks = dict(zip(out["id"], out["rra_tied_block_size"]))
    assert blocks == {"x": 2, "y": 2, "z": 1}
    # the weighted path is untouched: strict identity, no new column
    weighted = ra.rerank_output(df, "weighted")
    assert "rra_tied_block_size" not in weighted.columns
    assert weighted.equals(df)


def test_rank_candidates_writes_the_tied_block_column_under_rankagg_only():
    """rank_candidates.py writes ``df_sorted[output_cols]``, so a column the
    reranker computes but output_cols omits never reaches the CSV. Pinned at
    source level (running main() needs the full pipeline inputs), mirroring
    test_rank_candidates_emits_final_rank.py.

    It must be appended CONDITIONALLY: the unconditional default-fill further
    down would otherwise write 0.0 on the weighted path, which reads as
    "nothing is tied" when the truth is "no rank aggregation ran".
    """
    from pathlib import Path

    src = (Path(__file__).resolve().parents[2]
           / "scripts" / "rank_candidates.py").read_text()
    assert re.search(
        r"if RANK_METHOD == 'rankagg' and 'rra_tied_block_size' in df_sorted\.columns:"
        r"\s*\n\s*output_cols\.append\('rra_tied_block_size'\)", src), \
        "rra_tied_block_size must be appended to output_cols under rankagg only"
    # and it must NOT be in the static list (which the weighted path also writes)
    static = re.search(r"output_cols\s*=\s*\[(.*?)\n\]", src, re.DOTALL)
    assert static and "rra_tied_block_size" not in static.group(1)
