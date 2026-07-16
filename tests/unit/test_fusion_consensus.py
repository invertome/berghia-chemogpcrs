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


def test_select_rejects_noise_combiner_and_accepts_signal():
    rng = np.random.default_rng(3)
    ident = {f"c{i}": float(i) for i in range(80)}
    length = {f"c{i}": float(rng.normal()) for i in range(80)}
    # best single model: identity confound rho ~0.70, famAcc 0.80. (DEVIATION from
    # the plan's 0.50: the plan's own signal_nov below has max|confound|=0.659, so
    # with best confound 0.50 the spec-faithful select_consensus would correctly
    # REJECT the signal combiner for "confound" and selected!="rra". Setting the
    # best single's confound to 0.70 makes the signal combiner genuinely BEAT it
    # on BOTH axes (0.659<0.70 confound, 0.82>0.80 famAcc) — exactly the Task-4
    # acceptance criterion "accepts a combiner that beats the best single model on
    # both axes". The selection logic itself is unchanged.)
    single = {"protrek": {"confound_rho": 0.70, "famacc": 0.80}}
    # signal combiner: lower confound AND higher famAcc -> accept
    signal_nov = {f"c{i}": ident[f"c{i}"] * 0.2 + rng.normal(0, 5) for i in range(80)}
    # noise combiner: near-zero confound but destroyed famAcc -> reject
    noise_nov = {f"c{i}": rng.normal() for i in range(80)}
    result = select_consensus(
        combiners={"rra": signal_nov, "noise": noise_nov},
        combiner_famacc={"rra": 0.82, "noise": 0.55},
        confounds={"identity": ident, "length": length},
        best_single=single["protrek"], n_boot=500, seed=3,
    )
    assert result["selected"] == "rra"
    assert result["rejected"]["noise"]["reason"] == "famacc"


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
