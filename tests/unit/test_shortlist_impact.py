"""Tests for the shortlist-impact diagnostic (does emb_novelty move the top-k).

The ground-truth WITH/WITHOUT RRA orders are computed here directly from
``rank_aggregation`` (never hardcoded), and the tool's outputs are checked
against them. One synthetic df is engineered so that removing the novelty
voter demonstrably changes the top-3 (a candidate enters, another leaves);
another df simply lacks the novelty column, so the toggle is a no-op.
"""
import pandas as pd
import pytest

import rank_aggregation as ra
import shortlist_impact as si


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
def _pivotal_df():
    """10 candidates where emb_novelty is pivotal to the top-3.

    Per-signal ranks (1 = best) are laid out so that ``Y`` is strong on one
    base signal (positive) plus novelty, while ``B`` sits inside the top-3 on
    base alone. Removing novelty drops ``Y`` and lets ``B`` back in.
    """
    ids = ["A", "B", "Y", "Z", "F1", "F2", "F3", "F4", "F5", "F6"]
    fillers = ids[4:]
    n = len(ids)

    def base(top4):
        d = dict(top4)
        for i, f in enumerate(fillers):
            d[f] = 5 + i
        return d

    ranks = {
        "phylo":      base({"A": 1, "B": 2, "Z": 3, "Y": 4}),
        "purifying":  base({"B": 1, "A": 2, "Z": 3, "Y": 4}),
        "positive":   base({"Y": 1, "Z": 2, "A": 3, "B": 4}),
        "lse_divergence":  base({"A": 1, "B": 2, "Z": 3, "Y": 4}),
        # Y strong on novelty; fillers middling; A/B/Z weak on novelty
        "emb_novelty": {"Y": 1, "F1": 2, "F2": 3, "F3": 4, "F4": 5,
                        "F5": 6, "F6": 7, "Z": 8, "A": 9, "B": 10},
    }
    data = {"id": ids}
    for sig, rk in ranks.items():
        col = sig if sig == "emb_novelty" else f"{sig}_score"
        data[col] = [n - rk[i] + 1 for i in ids]  # higher score = better rank
    data["has_emb_data"] = [True] * n
    return pd.DataFrame(data)


def _base_only_df():
    """4 always-present base signals, no emb_novelty column at all."""
    ids = [f"c{i}" for i in range(6)]
    n = len(ids)
    lin = list(range(n, 0, -1))
    return pd.DataFrame({
        "id": ids,
        "phylo_score": lin,
        "purifying_score": lin,
        "positive_score": lin,
        "lse_divergence_score": lin,
    })


def _ground_truth_orders(df, toggle="emb_novelty"):
    rl = ra.build_ranklists_from_df(df)
    order_full = ra.aggregate(rl, method="rra")
    order_ablated = ra.aggregate(
        {k: v for k, v in rl.items() if k != toggle}, method="rra")
    return order_full, order_ablated


# --------------------------------------------------------------------------- #
# topk_jaccard
# --------------------------------------------------------------------------- #
def test_topk_jaccard_matches_set_definition():
    a = ["a", "b", "c", "d"]
    b = ["a", "c", "x", "y"]
    # top-2: {a,b} vs {a,c} -> inter {a}=1, union {a,b,c}=3
    assert si.topk_jaccard(a, b, 2) == pytest.approx(1 / 3)


def test_topk_jaccard_identical_is_one_and_disjoint_is_zero():
    a = ["a", "b", "c", "d"]
    assert si.topk_jaccard(a, a, 3) == pytest.approx(1.0)
    assert si.topk_jaccard(["a", "b"], ["x", "y"], 2) == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# topk_movement
# --------------------------------------------------------------------------- #
def test_topk_movement_entered_left_and_biggest_mover():
    order_full, order_ablated = _ground_truth_orders(_pivotal_df())
    # order_a = counterfactual (novelty removed), order_b = full (novelty in)
    mv = si.topk_movement(order_ablated, order_full, k=3)

    # novelty pulls Y INTO the top-3 and pushes B OUT
    assert mv["n_entered"] == 1 and mv["entered"] == ["Y"]
    assert mv["n_left"] == 1 and mv["left"] == ["B"]

    # Y is the biggest mover; positions are read from the ground-truth orders
    expected_top = ("Y",
                    order_ablated.index("Y") + 1,
                    order_full.index("Y") + 1)
    assert mv["movers"][0] == expected_top
    # movers are sorted by absolute movement, descending
    deltas = [abs(ra_ - rb_) for _, ra_, rb_ in mv["movers"]]
    assert deltas == sorted(deltas, reverse=True)
    # every union-of-top-k id is represented
    union = set(order_ablated[:3]) | set(order_full[:3])
    assert {m[0] for m in mv["movers"]} == union


# --------------------------------------------------------------------------- #
# shortlist_impact
# --------------------------------------------------------------------------- #
def test_shortlist_impact_novelty_changes_topk():
    df = _pivotal_df()
    order_full, order_ablated = _ground_truth_orders(df)
    a3, b3 = set(order_full[:3]), set(order_ablated[:3])
    expected_jac = len(a3 & b3) / len(a3 | b3)   # independent of si.topk_jaccard

    res = si.shortlist_impact(df, signal_to_toggle="emb_novelty", k=3)

    assert res["jaccard"] == pytest.approx(expected_jac)
    assert res["jaccard"] < 1.0
    assert res["signal_present"] is True
    assert res["movement"]["entered"] == ["Y"]
    assert res["movement"]["left"] == ["B"]
    assert res["movement"]["movers"][0][0] == "Y"


def test_shortlist_impact_absent_signal_is_a_noop():
    df = _base_only_df()
    res = si.shortlist_impact(df, signal_to_toggle="emb_novelty", k=3)

    assert res["jaccard"] == 1.0
    assert res["signal_present"] is False
    assert "note" in res
    assert res["movement"]["n_entered"] == 0
    assert res["movement"]["n_left"] == 0


# --------------------------------------------------------------------------- #
# main / CLI
# --------------------------------------------------------------------------- #
def test_main_prints_report(tmp_path, capsys):
    csv_path = tmp_path / "ranked.csv"
    _pivotal_df().to_csv(csv_path, index=False)

    res = si.main(["--ranked-csv", str(csv_path), "--signal", "emb_novelty",
                   "--k", "3"])
    out = capsys.readouterr().out

    assert res["jaccard"] == pytest.approx(0.5)
    assert "Jaccard" in out
    assert "0.5000" in out
    assert "Y" in out  # the candidate novelty pulls in
