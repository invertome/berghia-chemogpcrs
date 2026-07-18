"""Unit tests for the fusion-consensus harness (bead cw3.16).

Covers the pure statistical primitives (combiners, confound-independence gate,
family vote + LOO vote-famAcc, dual selection, bootstrap CI) and the
contract-emitting `build_consensus_channel`. No torch/esm/npz/GPU — every test
runs on toy in-memory inputs. `main` (IO glue) is exercised only for
importability; its correctness is IO-shaped and not unit-tested here.
"""
import sys

import numpy as np

sys.path.insert(0, "scripts")

# ---------------------------------------------------------------------------
# Task 1 — consensus combiners
# ---------------------------------------------------------------------------
from fusion_consensus import rank_average_consensus, robust_rank_aggregation


def _novelty(**models):  # {model: {cid: novelty}}
    return dict(models)


def test_rank_average_orders_by_cross_model_novelty():
    nov = _novelty(
        m1={"a": 9.0, "b": 5.0, "c": 1.0},
        m2={"a": 8.0, "b": 4.0, "c": 2.0},
    )
    s = rank_average_consensus(nov)
    assert s["a"] > s["b"] > s["c"]


def test_rra_consistent_top_gets_small_p():
    nov = _novelty(
        m1={"a": 9, "b": 8, "c": 7, "d": 1},
        m2={"a": 9, "b": 8, "c": 7, "d": 1},
        m3={"a": 9, "b": 8, "c": 7, "d": 1},
    )
    p = robust_rank_aggregation(nov)
    assert 0 < p["a"] <= 1.0
    assert p["a"] < p["d"]          # consistently-top < consistently-bottom


def test_rra_missing_candidate_uses_available_models():
    nov = _novelty(m1={"a": 9, "b": 1}, m2={"a": 8})  # b only in m1
    p = robust_rank_aggregation(nov)
    assert set(p) == {"a", "b"}


# ---------------------------------------------------------------------------
# Task 2 — confound-independence gate
# ---------------------------------------------------------------------------
from fusion_consensus import confound_residuals, residual_independence


def test_common_mode_confounds_trigger_stop():
    rng = np.random.default_rng(0)
    length = {f"c{i}": float(i) for i in range(60)}
    ident = {f"c{i}": float(i) for i in range(60)}          # identity == length here
    # both models' novelty are the SAME function of length -> shared confound
    m = {f"c{i}": length[f"c{i}"] + rng.normal(0, 0.01) for i in range(60)}
    nov = {"m1": dict(m), "m2": {k: v + rng.normal(0, 0.01) for k, v in m.items()}}
    res = confound_residuals(nov, {"length": length, "identity": ident})
    dec = residual_independence(res, threshold=0.6)
    assert dec["independent"] is False


def test_independent_residuals_trigger_go():
    rng = np.random.default_rng(1)
    length = {f"c{i}": float(i) for i in range(60)}
    ident = {f"c{i}": float(rng.normal()) for i in range(60)}
    nov = {  # residual signal is independent noise per model
        "m1": {f"c{i}": length[f"c{i}"] + rng.normal(0, 5) for i in range(60)},
        "m2": {f"c{i}": length[f"c{i}"] + rng.normal(0, 5) for i in range(60)},
    }
    res = confound_residuals(nov, {"length": length, "identity": ident})
    dec = residual_independence(res, threshold=0.6)
    assert dec["independent"] is True


# ---------------------------------------------------------------------------
# Task 3 — consensus family vote + reference LOO vote-famAcc
# ---------------------------------------------------------------------------
from fusion_consensus import plurality_family, loo_vote_famacc


def test_plurality_family_majority_and_tiebreak():
    pmf = {"m1": {"x": "pep"}, "m2": {"x": "pep"}, "m3": {"x": "amine"}}
    assert plurality_family(pmf)["x"] == "pep"
    tie = {"m1": {"y": "amine"}, "m2": {"y": "pep"}}      # tie -> alphabetical
    assert plurality_family(tie)["y"] == "amine"


def test_loo_vote_famacc_majority_corrects_single_error():
    truth = {"r1": "pep", "r2": "amine"}
    per_model = {
        "m1": {"r1": "pep",   "r2": "amine"},
        "m2": {"r1": "pep",   "r2": "amine"},
        "m3": {"r1": "amine", "r2": "opsin"},           # both wrong, out-voted
    }
    assert loo_vote_famacc(per_model, truth) == 1.0


# ---------------------------------------------------------------------------
# Task 4 — dual selection criterion with bootstrap CIs
# ---------------------------------------------------------------------------
from fusion_consensus import bootstrap_ci_diff, select_consensus


def test_bootstrap_ci_diff_brackets_the_point_difference():
    # a has a mean ~1 higher than b; the CI of mean(a)-mean(b) must bracket ~1
    rng = np.random.default_rng(7)
    a = rng.normal(1.0, 0.5, 200)
    b = rng.normal(0.0, 0.5, 200)
    lo, hi = bootstrap_ci_diff(a, b, np.mean, n=500, seed=7)
    assert lo < hi
    assert lo < 1.0 < hi


def test_select_ci_rejects_significantly_worse_and_accepts_noninferior():
    # CI-based dual criterion (cw3.16.1): reject a combiner only if it is
    # SIGNIFICANTLY worse than the best single model — famAcc CI entirely below 0,
    # or confound CI entirely above 0. "Don't chase noise": non-inferior combiners
    # are accepted; the noise combiner is rejected because its famAcc is robustly
    # worse, not because its point estimate merely differs.
    rng = np.random.default_rng(3)
    cids = [f"c{i}" for i in range(120)]
    ident = {c: float(i) for i, c in enumerate(cids)}
    length = {c: float(rng.normal()) for c in cids}
    confounds = {"identity": ident, "length": length}
    refs = [f"r{i}" for i in range(100)]

    # best single: moderately identity-confounded novelty, famAcc 0.80
    best_nov = {c: ident[c] * 0.3 + rng.normal(0, 5) for c in cids}
    best_correct = {r: int(i < 80) for i, r in enumerate(refs)}
    # signal combiner: LESS confounded (weaker identity coupling) + higher famAcc
    signal_nov = {c: ident[c] * 0.05 + rng.normal(0, 5) for c in cids}
    signal_correct = {r: int(i < 88) for i, r in enumerate(refs)}
    # noise combiner: near-zero confound but robustly worse famAcc (0.55) -> reject
    noise_nov = {c: rng.normal() for c in cids}
    noise_correct = {r: int(i < 55) for i, r in enumerate(refs)}

    result = select_consensus(
        combiner_novelty={"rra": signal_nov, "noise": noise_nov},
        combiner_ref_correct={"rra": signal_correct, "noise": noise_correct},
        confounds=confounds,
        best_single_novelty=best_nov,
        best_single_ref_correct=best_correct,
        n_boot=400, seed=3,
    )
    assert result["selected"] == "rra"
    assert result["rejected"]["noise"]["reason"] == "famacc"
    # the accepted signal combiner carries its CIs for auditability
    assert "famacc_ci" in result["accepted"]["rra"]
    assert "confound_ci" in result["accepted"]["rra"]


# ---------------------------------------------------------------------------
# Task 5 — contract-compatible consensus channel + CLI importability
# ---------------------------------------------------------------------------
from fusion_consensus import build_consensus_channel


def test_channel_emits_mahalanobis_contract():
    nov = {"m1": {"a": 9, "b": 1}, "m2": {"a": 8, "b": 2}}
    fam = {"m1": {"a": "pep", "b": "amine"}, "m2": {"a": "pep", "b": "opsin"}}
    ch = build_consensus_channel(nov, fam)
    assert set(ch["a"]) == {"emb_novelty", "emb_nonchemo_family",
                            "has_emb_data", "emb_leakage_flag"}
    assert ch["a"]["emb_nonchemo_family"] == "pep"
    assert ch["a"]["has_emb_data"] is True and ch["a"]["emb_leakage_flag"] is True
    assert ch["a"]["emb_novelty"] > ch["b"]["emb_novelty"]


def test_main_is_importable_and_self_test_runs():
    import fusion_consensus
    assert callable(fusion_consensus.main)
    # the synthetic --self-test path runs the full gate->combiners->selection chain
    fusion_consensus.main(["--self-test"])


# ---------------------------------------------------------------------------
# Task 6 (cw3.6) — length-deconfounded consensus channel
# ---------------------------------------------------------------------------
def test_deconfound_removes_length_only_novelty():
    # A pool where novelty tracks length exactly (c0..c9) plus one candidate
    # ("real") whose novelty far exceeds what its length predicts. RAW RRA ranks
    # the longest pure-length candidate ("c9") above "real"; after deconfounding
    # novelty on seq_len, the length-only signal is removed and "real" (novel
    # beyond its length) must rise above "c9". This is the locked cw3.6 fix:
    # the consensus channel must combine length-DECONFOUNDED novelties, not raw.
    length = {f"c{i}": float((i + 1) * 10) for i in range(10)}   # 10..100
    length["real"] = 25.0
    base_nov = dict(length)                                      # novelty == length
    base_nov["real"] = 95.0                                      # novel beyond length
    nov = {"m1": dict(base_nov), "m2": dict(base_nov)}
    fam = {m: {c: "amine" for c in base_nov} for m in nov}

    raw = build_consensus_channel(nov, fam)
    dec = build_consensus_channel(nov, fam, deconfound={"seq_len": length})

    # every candidate is still emitted with the full contract after deconfounding
    assert set(dec) == set(base_nov)
    assert set(dec["real"]) == {"emb_novelty", "emb_nonchemo_family",
                                "has_emb_data", "emb_leakage_flag"}
    # raw: the longest pure-length candidate outranks the real-novel one
    assert raw["c9"]["emb_novelty"] > raw["real"]["emb_novelty"]
    # deconfounded: length-only novelty removed -> real signal rises above c9
    assert dec["real"]["emb_novelty"] > dec["c9"]["emb_novelty"]


def test_deconfound_none_matches_raw_channel():
    # deconfound=None (default) must be identical to the raw combiner path.
    nov = {"m1": {"a": 9, "b": 5, "c": 1}, "m2": {"a": 8, "b": 4, "c": 2}}
    fam = {"m1": {c: "amine" for c in "abc"}, "m2": {c: "amine" for c in "abc"}}
    explicit_none = build_consensus_channel(nov, fam, deconfound=None)
    default = build_consensus_channel(nov, fam)
    for c in "abc":
        assert explicit_none[c]["emb_novelty"] == default[c]["emb_novelty"]


# ---------------------------------------------------------------------------
# A1 (v4bs.2) — phylogeny-residualized novelty as a SECOND (dormant) channel
# ---------------------------------------------------------------------------
def test_residual_deconfound_emits_excess_beyond_phylogeny_channel():
    # A pool where novelty tracks tree-distance exactly (t0..t9), plus 'beyond'
    # whose novelty far exceeds what its (small) tree-distance predicts. RAW
    # emb_novelty ranks the most phylo-distant t9 top; the residual channel
    # (novelty residualized on tree_distance) removes the phylogeny-explained
    # part, so 'beyond' (novel beyond phylogenetic expectation) rises above t9.
    tdist = {f"t{i}": float((i + 1) * 10) for i in range(10)}    # 10..100
    tdist["beyond"] = 25.0
    base = dict(tdist)
    base["beyond"] = 95.0                                       # novel beyond phylo
    nov = {"m1": dict(base), "m2": dict(base)}
    fam = {m: {c: "amine" for c in base} for m in nov}

    ch = build_consensus_channel(
        nov, fam, residual_deconfound={"tree_distance": tdist}
    )
    # both the raw AND the residual novelty are emitted per candidate
    assert "emb_novelty" in ch["beyond"] and "emb_novelty_residual" in ch["beyond"]
    # raw: the most phylo-distant pure-phylo candidate outranks 'beyond'
    assert ch["t9"]["emb_novelty"] > ch["beyond"]["emb_novelty"]
    # residual: phylogeny removed -> the beyond-phylo signal rises above t9
    assert ch["beyond"]["emb_novelty_residual"] > ch["t9"]["emb_novelty_residual"]


def test_no_residual_deconfound_keeps_original_four_key_contract():
    # Backward-compatible: without residual_deconfound, no extra key is added.
    nov = {"m1": {"a": 9, "b": 1}, "m2": {"a": 8, "b": 2}}
    fam = {"m1": {"a": "pep", "b": "amine"}, "m2": {"a": "pep", "b": "amine"}}
    ch = build_consensus_channel(nov, fam)
    assert "emb_novelty_residual" not in ch["a"]


def test_residual_is_none_for_candidates_missing_the_confound():
    # A candidate with no tree_distance gets a None residual but keeps its
    # primary emb_novelty (dormant partial coverage, like the seq_len path).
    tdist = {"a": 10.0, "b": 20.0}                              # 'nocover' absent
    nov = {"m1": {"a": 9, "b": 5, "nocover": 1},
           "m2": {"a": 8, "b": 4, "nocover": 2}}
    fam = {m: {c: "amine" for c in ("a", "b", "nocover")} for m in nov}
    ch = build_consensus_channel(
        nov, fam, residual_deconfound={"tree_distance": tdist}
    )
    assert ch["nocover"]["emb_novelty_residual"] is None
    assert isinstance(ch["nocover"]["emb_novelty"], float)
