import numpy as np
from permutation_null import topk_separation, permutation_null_pvalue


def _signals_with_shared_truth(n, rng):
    truth = rng.normal(size=n)
    ids = [f"c{i}" for i in range(n)]
    return {f"s{j}": {ids[i]: truth[i] + rng.normal(0, 0.5) for i in range(n)}
            for j in range(5)}, ids


def test_separation_positive_when_top_is_separated():
    order = ["a", "b", "c", "d"]
    agg = {"a": 3.0, "b": 2.5, "c": 0.2, "d": 0.1}
    assert topk_separation(order, agg, k=2) > 0


def test_real_signal_beats_null():
    rng = np.random.default_rng(0)
    signals, _ = _signals_with_shared_truth(60, rng)
    assert permutation_null_pvalue(signals, k=10, n_perm=200, seed=1) < 0.05


def test_pure_noise_does_not_beat_null():
    rng = np.random.default_rng(2)
    ids = [f"c{i}" for i in range(60)]
    signals = {f"s{j}": {i: rng.normal() for i in ids} for j in range(5)}
    assert permutation_null_pvalue(signals, k=10, n_perm=200, seed=3) > 0.2
