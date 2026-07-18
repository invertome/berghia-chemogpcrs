"""Unit tests for the A3 model-role-split producer (epic v4bs, task v4bs.4).

Pure-function tests only (no torch/esm/GPU). The module is imported by bare
name because ``tests/unit/conftest.py`` puts ``scripts/`` on ``sys.path``.
"""
import numpy as np

from model_role_split import main, novelty_percentiles, role_gap


def test_novelty_percentiles_range_and_direction():
    # Highest input score must map to the highest percentile; all in (0, 1].
    scores = {"a": -3.0, "b": 0.0, "c": 10.0, "d": 2.5}
    pct = novelty_percentiles(scores)
    assert set(pct) == set(scores)
    assert all(0.0 < v <= 1.0 for v in pct.values())
    # Ordering of percentiles follows ordering of raw scores.
    assert pct["c"] > pct["d"] > pct["b"] > pct["a"]
    # The maximum-score candidate sits at the top of the (0, 1] range.
    assert pct["c"] == max(pct.values())
    assert np.isclose(pct["c"], 1.0)


def test_novelty_percentiles_ties_are_average_ranks():
    # Two equal scores must receive one identical, deterministic percentile
    # (the average of their tied ranks), not an order-dependent value.
    scores = {"a": 1.0, "b": 5.0, "c": 5.0, "d": 9.0}
    pct = novelty_percentiles(scores)
    assert pct["b"] == pct["c"]
    # ranks 1,(2,3)->2.5,(2,3)->2.5,4 over n=4  ->  0.625 for the tied pair.
    assert np.isclose(pct["b"], 2.5 / 4)
    assert np.isclose(pct["a"], 1.0 / 4)
    assert np.isclose(pct["d"], 4.0 / 4)


def test_role_gap_identical_maps_all_zero():
    # If phylo and function novelty are the same map, the two roles agree
    # perfectly and every gap must be exactly 0.
    scores = {"a": 1.0, "b": 4.0, "c": 2.0, "d": 9.0, "e": -1.0}
    gap = role_gap(scores, dict(scores))
    assert set(gap) == set(scores)
    assert all(g == 0.0 for g in gap.values())


def test_role_gap_function_top_phylo_bottom_is_largest_positive():
    # Candidate "x" is the MOST novel in function but the LEAST novel in phylo:
    # it must carry the single largest positive gap in the set.
    phylo = {"x": 0.0, "y": 5.0, "z": 8.0, "w": 3.0}
    function = {"x": 100.0, "y": 5.0, "z": 8.0, "w": 3.0}
    gap = role_gap(phylo, function)
    assert gap["x"] == max(gap.values())
    assert gap["x"] > 0.0
    assert all(gap["x"] >= g for g in gap.values())


def test_role_gap_scale_invariance_multiply_and_shift():
    # Percentiles are rank-based, so a positive-affine transform of ALL of one
    # model's scores must leave every gap unchanged.
    phylo = {"a": 1.0, "b": 3.0, "c": 7.0, "d": 2.0}
    function = {"a": 4.0, "b": 1.0, "c": 9.0, "d": 6.0}
    base = role_gap(phylo, function)
    scaled = role_gap({k: v * 13.0 for k, v in phylo.items()}, function)
    shifted = role_gap(phylo, {k: v + 100.0 for k, v in function.items()})
    for cid in base:
        assert np.isclose(base[cid], scaled[cid])
        assert np.isclose(base[cid], shifted[cid])


def test_role_gap_uses_intersection_only():
    # A candidate present in only one map must be absent from the output.
    phylo = {"a": 1.0, "b": 2.0, "c": 3.0, "only_phylo": 9.0}
    function = {"a": 2.0, "b": 1.0, "c": 3.0, "only_func": 9.0}
    gap = role_gap(phylo, function)
    assert set(gap) == {"a", "b", "c"}
    assert "only_phylo" not in gap
    assert "only_func" not in gap


def test_role_gap_definition_matches_percentile_difference():
    # Gap == function percentile minus phylo percentile, each computed WITHIN
    # its own model's full candidate set (rank-based, scale-free).
    phylo = {"a": 1.0, "b": 2.0, "c": 3.0}
    function = {"a": 3.0, "b": 2.0, "c": 1.0}
    p_pct = novelty_percentiles(phylo)
    f_pct = novelty_percentiles(function)
    gap = role_gap(phylo, function)
    for cid in gap:
        assert np.isclose(gap[cid], f_pct[cid] - p_pct[cid])


def _write_novelty_tsv(path, mapping):
    import pandas as pd
    pd.DataFrame(
        {"candidate_id": list(mapping), "novelty": list(mapping.values())}
    ).to_csv(path, sep="\t", index=False)


def test_main_writes_expected_schema_intersection_and_sort(tmp_path):
    import pandas as pd

    phylo = {"a": 0.0, "b": 5.0, "c": 8.0, "d": 3.0, "only_phylo": 9.0}
    function = {"a": 100.0, "b": 5.0, "c": 8.0, "d": 3.0, "only_func": 9.0}
    phylo_tsv = tmp_path / "novelty_phylo_PROD.tsv"
    function_tsv = tmp_path / "novelty_function_PROD.tsv"
    out_tsv = tmp_path / "role_split.tsv"
    _write_novelty_tsv(phylo_tsv, phylo)
    _write_novelty_tsv(function_tsv, function)

    main([
        "--phylo-tsv", str(phylo_tsv),
        "--function-tsv", str(function_tsv),
        "--out", str(out_tsv),
    ])

    df = pd.read_csv(out_tsv, sep="\t")
    assert list(df.columns) == [
        "id", "emb_novelty_phylo", "emb_novelty_function", "emb_role_gap"
    ]
    # Only the intersection is emitted.
    assert set(df["id"]) == {"a", "b", "c", "d"}
    # Sorted by emb_role_gap descending; "a" (top function / bottom phylo) leads.
    assert list(df["emb_role_gap"]) == sorted(df["emb_role_gap"], reverse=True)
    assert df.iloc[0]["id"] == "a"
    # Column contents match the pure functions.
    p_pct = novelty_percentiles(phylo)
    f_pct = novelty_percentiles(function)
    row_a = df[df["id"] == "a"].iloc[0]
    assert np.isclose(row_a["emb_novelty_phylo"], p_pct["a"])
    assert np.isclose(row_a["emb_novelty_function"], f_pct["a"])
    assert np.isclose(row_a["emb_role_gap"], f_pct["a"] - p_pct["a"])
