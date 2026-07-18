import numpy as np
from rank_confidence import bootstrap_rank_intervals, assign_tiers


def test_dominant_candidate_has_tight_top_rank():
    per_signal = {f"s{j}": {"win": 10.0, "mid": 5.0, "lose": 1.0} for j in range(6)}
    res = bootstrap_rank_intervals(per_signal, k=1, n_boot=200, seed=0)
    assert res["win"]["rank_ci_hi"] == 1
    assert res["win"]["p_top_k"] == 1.0
    assert res["lose"]["p_top_k"] == 0.0


def test_tiers_partition_by_p_top_k():
    res = {"a": {"p_top_k": 0.95}, "b": {"p_top_k": 0.4}, "c": {"p_top_k": 0.05}}
    assert assign_tiers(res) == {"a": "high", "b": "plausible", "c": "tail"}
