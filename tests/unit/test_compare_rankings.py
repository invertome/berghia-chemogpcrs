"""Tests for scripts/compare_rankings.py."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import compare_rankings as cr


def _write_ranked(path: Path, ids: list[str]) -> None:
    df = pd.DataFrame({
        "id": ids,
        "rank_score": [1.0 - 0.01 * i for i in range(len(ids))],
    })
    df.to_csv(path, index=False)


def test_jaccard_top_n_identical_lists() -> None:
    a = ["x", "y", "z", "w"]
    assert cr.jaccard_top_n(a, a, 4) == 1.0
    assert cr.jaccard_top_n(a, a, 2) == 1.0


def test_jaccard_top_n_disjoint() -> None:
    a = ["a", "b", "c"]
    b = ["x", "y", "z"]
    assert cr.jaccard_top_n(a, b, 3) == 0.0


def test_jaccard_top_n_partial_overlap() -> None:
    a = ["a", "b", "c", "d"]
    b = ["a", "b", "x", "y"]
    # top-2 of each: {a,b} ∩ {a,b} = 2; ∪ = 2; Jaccard = 1.0
    assert cr.jaccard_top_n(a, b, 2) == 1.0
    # top-4: {a,b,c,d} ∩ {a,b,x,y} = 2; ∪ = 6; Jaccard = 2/6
    assert cr.jaccard_top_n(a, b, 4) == pytest_approx(2 / 6)


def pytest_approx(value: float, tol: float = 1e-9) -> object:
    """Inline approx since pytest fixtures aren't imported here."""
    class _A:
        def __init__(self, v: float, t: float) -> None:
            self.v = v
            self.t = t
        def __eq__(self, other: object) -> bool:
            return abs(float(other) - self.v) < self.t
        def __repr__(self) -> str:
            return f"approx({self.v}±{self.t})"
    return _A(value, tol)


def test_median_rank_shift_identical() -> None:
    a = ["x", "y", "z", "w"]
    assert cr.median_rank_shift_top_n(a, a, 4) == 0.0


def test_median_rank_shift_reversed() -> None:
    a = ["x", "y", "z"]
    b = ["z", "y", "x"]
    # rank shifts in top-3: x: 0->2 (|2|), y: 1->1 (|0|), z: 2->0 (|2|)
    # sorted: [0, 2, 2]; median = 2
    assert cr.median_rank_shift_top_n(a, b, 3) == 2.0


def test_compare_two_identical_inputs() -> None:
    a = ["g1", "g2", "g3", "g4", "g5"]
    rows = cr.compare_two("A", a, "B", a, [3, 5])
    for row in rows:
        # scipy spearmanr/kendalltau may return 0.9999999... for identical
        # ranks due to float rounding; allow tolerance.
        assert abs(row["spearman_rho"] - 1.0) < 1e-9
        assert abs(row["kendall_tau"] - 1.0) < 1e-9
        assert row["jaccard_topN"] == 1.0
        assert row["median_rank_shift_topN"] == 0.0
        assert row["n_common"] == 5


def test_compare_two_disjoint_inputs() -> None:
    a = ["a1", "a2", "a3"]
    b = ["b1", "b2", "b3"]
    rows = cr.compare_two("A", a, "B", b, [3])
    assert rows[0]["n_common"] == 0
    assert rows[0]["jaccard_topN"] == 0.0
    assert rows[0]["spearman_rho"] is None


def test_parse_csv_args() -> None:
    out = cr.parse_csv_args(["default=path/a.csv", "no_dnds=path/b.csv"])
    assert out == {"default": "path/a.csv", "no_dnds": "path/b.csv"}


def test_parse_csv_args_handles_paths_with_equals() -> None:
    """If the path contains '=', only split on the first one."""
    out = cr.parse_csv_args(["cfg=path/with=equals.csv"])
    assert out == {"cfg": "path/with=equals.csv"}


def test_load_ranked_sorts_by_score(tmp_path: Path) -> None:
    """If rank_score is present, the loader returns IDs in score-descending order
    even if the CSV is not pre-sorted."""
    csv = tmp_path / "ranked.csv"
    pd.DataFrame({
        "id": ["g1", "g2", "g3"],
        "rank_score": [0.5, 0.9, 0.7],
    }).to_csv(csv, index=False)
    assert cr.load_ranked(str(csv)) == ["g2", "g3", "g1"]


def test_load_ranked_without_score_column(tmp_path: Path) -> None:
    """If rank_score is missing, return IDs in file order."""
    csv = tmp_path / "no_score.csv"
    pd.DataFrame({"id": ["g1", "g2", "g3"]}).to_csv(csv, index=False)
    assert cr.load_ranked(str(csv)) == ["g1", "g2", "g3"]


def test_load_ranked_missing_id_column_raises(tmp_path: Path) -> None:
    csv = tmp_path / "bad.csv"
    pd.DataFrame({"name": ["a", "b"]}).to_csv(csv, index=False)
    try:
        cr.load_ranked(str(csv))
    except ValueError as e:
        assert "id" in str(e)
    else:
        raise AssertionError("Expected ValueError for missing id column")
