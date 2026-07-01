"""Tests for the active-learning probe-batch selector (Task 7).

Covers: disagreement_score's dispersion behavior, select_batch's classifier
gate + mode ordering + diversity cap, write_batch's exact column contract,
and an end-to-end main() CLI smoke test.
"""
from collections import Counter

import pandas as pd
import pytest

import rank_aggregation as ra
import select_probe_batch as spb


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _ranks_from_df(df):
    """Mirror main()'s own composition: raw ranklists -> per-signal normalized ranks."""
    raw = ra.build_ranklists_from_df(df)
    return {sig: ra.normalized_ranks(scores) for sig, scores in raw.items()}


# --------------------------------------------------------------------------- #
# disagreement_score
# --------------------------------------------------------------------------- #
def test_disagreement_score_high_for_conflicting_signals():
    # 'a' is #1 in s1 (1/3) and last in s2 (3/3); 'b' is consistently middle.
    per_signal_ranks = {
        "s1": {"a": 1 / 3, "b": 2 / 3, "c": 3 / 3},
        "s2": {"a": 3 / 3, "b": 2 / 3, "c": 1 / 3},
    }
    assert spb.disagreement_score(per_signal_ranks, "a") > spb.disagreement_score(per_signal_ranks, "b")


def test_disagreement_score_zero_when_consistent():
    per_signal_ranks = {"s1": {"a": 0.5}, "s2": {"a": 0.5}, "s3": {"a": 0.5}}
    assert spb.disagreement_score(per_signal_ranks, "a") == pytest.approx(0.0)


def test_disagreement_score_needs_at_least_two_signals():
    # present in exactly one signal -> nothing to disagree with -> 0.0
    assert spb.disagreement_score({"s1": {"a": 0.2}}, "a") == 0.0
    # present in one signal, absent from the other (not "0 vs 0.2") -> 0.0
    assert spb.disagreement_score({"s1": {"a": 0.2}, "s2": {}}, "a") == 0.0
    # absent from every signal entirely -> 0.0
    assert spb.disagreement_score({"s1": {"b": 0.2}}, "missing") == 0.0
    assert spb.disagreement_score({}, "missing") == 0.0


def test_disagreement_score_matches_population_variance():
    # Locks the dispersion metric choice (population variance, ddof=0).
    import numpy as np

    per_signal_ranks = {"s1": {"a": 0.1}, "s2": {"a": 0.9}, "s3": {"a": 0.5}}
    expected = float(np.var([0.1, 0.9, 0.5]))
    assert spb.disagreement_score(per_signal_ranks, "a") == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# select_batch: classifier gate
# --------------------------------------------------------------------------- #
def test_select_batch_excludes_non_chemoreceptor_and_likely_non_chemoreceptor():
    df = pd.DataFrame(
        {
            "id": ["c0", "c1", "c2", "c3"],
            "phylo_score": [4, 3, 2, 1],
            "classification": [
                "non-chemoreceptor",
                "likely-non-chemoreceptor",
                "chemoreceptor-candidate",
                "chemoreceptor-candidate",
            ],
        }
    )
    per_signal_ranks = _ranks_from_df(df)
    selected = spb.select_batch(df, per_signal_ranks, batch_size=4, mode="topscore")
    assert "c0" not in selected
    assert "c1" not in selected
    assert set(selected) == {"c2", "c3"}


def test_select_batch_missing_classification_is_not_excluded():
    df = pd.DataFrame(
        {
            "id": ["c0", "c1"],
            "phylo_score": [2, 1],
            "classification": [None, "non-chemoreceptor"],
        }
    )
    per_signal_ranks = _ranks_from_df(df)
    selected = spb.select_batch(df, per_signal_ranks, batch_size=2, mode="topscore")
    assert selected == ["c0"]


# --------------------------------------------------------------------------- #
# select_batch: batch_size + mode ordering
# --------------------------------------------------------------------------- #
def test_select_batch_respects_batch_size_upper_bound():
    df = pd.DataFrame(
        {"id": [f"c{i}" for i in range(10)], "phylo_score": list(range(10, 0, -1))}
    )
    per_signal_ranks = _ranks_from_df(df)
    selected = spb.select_batch(df, per_signal_ranks, batch_size=3, mode="topscore")
    assert len(selected) == 3
    assert selected == ["c0", "c1", "c2"]


def test_select_batch_topscore_matches_pure_classifier_gated_aggregate():
    df = pd.DataFrame(
        {
            "id": ["c0", "c1", "c2", "c3"],
            "phylo_score": [4, 3, 2, 1],
            "purifying_score": [4, 3, 2, 1],
            "classification": [
                "chemoreceptor-candidate",
                "non-chemoreceptor",
                "chemoreceptor-candidate",
                "chemoreceptor-candidate",
            ],
        }
    )
    per_signal_ranks = _ranks_from_df(df)
    full_order = ra.aggregate(ra.build_ranklists_from_df(df), method="rra")
    expected = [i for i in full_order if i != "c1"][:2]

    got = spb.select_batch(df, per_signal_ranks, batch_size=2, mode="topscore")
    assert got == expected


def test_select_batch_active_mode_promotes_high_disagreement_candidate():
    # c2 is worst-of-5 under phylo but flips to best-of-5 under purifying --
    # a big swing that leaves its rough standing clearly behind c0 and c1,
    # yet gives it far higher disagreement than everyone else.
    df = pd.DataFrame(
        {
            "id": ["c0", "c1", "c2", "c3", "c4"],
            "phylo_score": [5, 4, 1, 3, 2],
            "purifying_score": [5, 4, 10, 3, 2],
        }
    )
    per_signal_ranks = _ranks_from_df(df)

    topscore = spb.select_batch(df, per_signal_ranks, batch_size=2, mode="topscore")
    assert "c2" not in topscore  # sanity precondition: not a top-2 pick on standing alone

    active = spb.select_batch(df, per_signal_ranks, batch_size=2, mode="active")
    assert "c2" in active  # promoted because of its far larger cross-signal disagreement


def test_select_batch_invalid_mode_raises():
    df = pd.DataFrame({"id": ["a"], "phylo_score": [1.0]})
    per_signal_ranks = _ranks_from_df(df)
    with pytest.raises(ValueError):
        spb.select_batch(df, per_signal_ranks, batch_size=1, mode="bogus")


# --------------------------------------------------------------------------- #
# select_batch: diversity cap
# --------------------------------------------------------------------------- #
def test_select_batch_diversity_cap_limits_per_cluster():
    df = pd.DataFrame(
        {
            "id": ["a0", "a1", "a2", "a3", "b0", "b1"],
            "phylo_score": [6, 5, 4, 3, 2, 1],  # cluster A ids all rank above cluster B
            "cluster": ["A", "A", "A", "A", "B", "B"],
        }
    )
    per_signal_ranks = _ranks_from_df(df)
    selected = spb.select_batch(
        df, per_signal_ranks, batch_size=4, cluster_col="cluster", mode="topscore"
    )
    assert len(selected) == 4
    cluster_of = dict(zip(df["id"], df["cluster"]))
    counts = Counter(cluster_of[i] for i in selected)
    assert counts["A"] <= 2  # ceil(4/2) = 2
    assert counts["B"] <= 2
    assert set(selected) == {"a0", "a1", "b0", "b1"}


def test_select_batch_active_mode_respects_diversity_cap_too():
    df = pd.DataFrame(
        {
            "id": ["a0", "a1", "a2", "a3", "b0", "b1"],
            "phylo_score": [6, 5, 4, 3, 2, 1],
            "purifying_score": [1, 2, 3, 4, 5, 6],  # disagreement for everyone
            "cluster": ["A", "A", "A", "A", "B", "B"],
        }
    )
    per_signal_ranks = _ranks_from_df(df)
    selected = spb.select_batch(
        df, per_signal_ranks, batch_size=4, cluster_col="cluster", mode="active"
    )
    assert len(selected) == 4
    cluster_of = dict(zip(df["id"], df["cluster"]))
    counts = Counter(cluster_of[i] for i in selected)
    assert counts["A"] <= 2
    assert counts["B"] <= 2


def test_select_batch_default_cluster_col_matches_pipeline_schema():
    # Regression: the real ranked CSV names the tandem-array column
    # 'tandem_cluster_id' (rank_candidates.py output_cols), NEVER 'cluster'.
    # The diversity cap -- the tool's headline feature -- must fire by DEFAULT
    # against that schema; a wrong default silently no-ops it and a batch could
    # end up 100% from one tandem cluster. This df has NO 'cluster' column.
    df = pd.DataFrame(
        {
            "id": ["a0", "a1", "a2", "a3", "b0", "b1"],
            "phylo_score": [6, 5, 4, 3, 2, 1],  # cluster A ids all rank above B
            "tandem_cluster_id": ["A", "A", "A", "A", "B", "B"],
        }
    )
    per_signal_ranks = _ranks_from_df(df)
    selected = spb.select_batch(df, per_signal_ranks, batch_size=4, mode="topscore")
    assert len(selected) == 4
    cluster_of = dict(zip(df["id"], df["tandem_cluster_id"]))
    counts = Counter(cluster_of[i] for i in selected)
    assert counts["A"] <= 2  # ceil(4/2) = 2 -- cap enforced BY DEFAULT
    assert counts["B"] <= 2
    assert set(selected) == {"a0", "a1", "b0", "b1"}


def test_select_batch_diversity_cap_skipped_when_cluster_col_absent():
    df = pd.DataFrame({"id": [f"c{i}" for i in range(4)], "phylo_score": [4, 3, 2, 1]})
    per_signal_ranks = _ranks_from_df(df)
    selected = spb.select_batch(
        df, per_signal_ranks, batch_size=4, cluster_col="nonexistent_col", mode="topscore"
    )
    assert selected == ["c0", "c1", "c2", "c3"]


# --------------------------------------------------------------------------- #
# write_batch
# --------------------------------------------------------------------------- #
def test_write_batch_emits_expected_columns_and_populated_rationale(tmp_path):
    # Input names the tandem-array column 'tandem_cluster_id' (the real
    # rank_candidates.py schema, and write_batch's default cluster_col); the
    # OUTPUT column is always 'cluster' regardless.
    df = pd.DataFrame(
        {
            "id": ["c0", "c1", "c2"],
            "phylo_score": [3, 2, 1],
            "purifying_score": [1, 2, 3],
            "tandem_cluster_id": ["g1", "g1", "g2"],
        }
    )
    per_signal_ranks = _ranks_from_df(df)
    selected = ["c0", "c2"]
    out = tmp_path / "batch.tsv"

    spb.write_batch(selected, df, per_signal_ranks, str(out))

    written = pd.read_csv(out, sep="\t")
    assert list(written.columns) == [
        "id",
        "rank_agg_position",
        "disagreement",
        "cluster",
        "top_supporting_signals",
        "conflicting_signals",
        "rationale",
    ]
    assert list(written["id"]) == selected
    assert written["rationale"].map(lambda s: isinstance(s, str) and len(s) > 0).all()

    by_id = written.set_index("id")
    for cand_id in selected:
        expected_dis = spb.disagreement_score(per_signal_ranks, cand_id)
        assert by_id.loc[cand_id, "disagreement"] == pytest.approx(expected_dis)

    assert by_id.loc["c0", "cluster"] == "g1"
    assert by_id.loc["c2", "cluster"] == "g2"
    # c0: phylo best (norm 1/3), purifying worst (norm 3/3) -> exact best/worst signal names
    assert by_id.loc["c0", "top_supporting_signals"] == "phylo"
    assert by_id.loc["c0", "conflicting_signals"] == "purifying"


def test_write_batch_creates_parent_directories(tmp_path):
    df = pd.DataFrame({"id": ["a", "b"], "phylo_score": [2, 1]})
    per_signal_ranks = _ranks_from_df(df)
    out = tmp_path / "nested" / "dir" / "batch.tsv"

    spb.write_batch(["a"], df, per_signal_ranks, str(out))

    assert out.exists()


# --------------------------------------------------------------------------- #
# main() CLI
# --------------------------------------------------------------------------- #
def test_main_cli_end_to_end(tmp_path):
    df = pd.DataFrame(
        {
            "id": [f"c{i}" for i in range(6)],
            "phylo_score": [6, 5, 4, 3, 2, 1],
            "purifying_score": [1, 2, 3, 4, 5, 6],
            "classification": ["chemoreceptor-candidate"] * 5 + ["non-chemoreceptor"],
            "tandem_cluster_id": ["A", "A", "B", "B", "C", "C"],
        }
    )
    csv_path = tmp_path / "ranked.csv"
    df.to_csv(csv_path, index=False)
    out_path = tmp_path / "out" / "batch.tsv"

    rc = spb.main(
        [
            "--ranked-csv",
            str(csv_path),
            "--batch-size",
            "3",
            "--mode",
            "active",
            "--out",
            str(out_path),
        ]
    )

    assert rc == 0
    written = pd.read_csv(out_path, sep="\t")
    assert len(written) == 3
    assert "c5" not in set(written["id"])  # non-chemoreceptor excluded
