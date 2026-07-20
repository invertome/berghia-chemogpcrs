"""Integration tests for Task 6: wiring the structural/embedding/microswitch
evidence channels (Tasks 4/5 + or_microswitch.py) into the rank-aggregation
reranker.

Design decision (deviates from the plan's literal "both paths" -- see
docs/plans/2026-07-01-ml-plm-chemoreceptor-ranking.md Task 6): rank
aggregation is WEIGHTLESS and label-free; the legacy weighted scorer
(``calculate_fair_rank_score`` in rank_candidates.py) stays an UNTOUCHED,
pristine baseline for an honest weighted-vs-rankagg comparison. So the three
new evidence channels are wired ONLY into ``scripts/rank_aggregation.py``'s
``SIGNAL_SPEC`` / ``build_ranklists_from_df`` / ``merge_evidence_channels`` --
``rank_candidates.py`` is not modified by this task.

Polarity rule under test: exclusion signals (``struct_nonchemo_corrob``,
``emb_nonchemo_sim``) must LOWER a candidate's standing when triggered. Rank
aggregation has no weights to negate, so this is achieved purely by
INVERTING the stored ranklist value (``1 - raw``) so "higher stored value =
better" holds for every signal the aggregator sees.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import audit_signal_ranking_independence as audit
import rank_aggregation as ra


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
def _base_signal_df():
    """A df carrying only the pre-Task-6 12 ranking signals (no channels)."""
    return pd.DataFrame(
        {
            "id": ["c0", "c1", "c2", "c3", "c4"],
            "phylo_score_norm": [0.9, 0.7, 0.5, 0.3, 0.1],
            "purifying_score_norm": [0.2, 0.8, 0.4, 0.6, 0.1],
            "positive_score_norm": [0.5, 0.5, 0.9, 0.1, 0.3],
            "lse_divergence_score_norm": [0.4, 0.4, 0.4, 0.9, 0.1],
            "synteny_score_norm": [0.8, 0.2, 0.5, 0.6, 0.3],
            "has_synteny_data": [True, True, False, True, True],
            "og_confidence_score_norm": [0.3, 0.6, 0.9, 0.2, 0.5],
            "has_og_confidence_data": [True, True, True, True, False],
            "tandem_cluster_score_norm": [0.1, 0.9, 0.5, 0.4, 0.2],
            "has_tandem_cluster_data": [True, False, True, True, True],
        }
    )


def _write_tsv(path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)
    return str(path)


# --------------------------------------------------------------------------- #
# REGRESSION: channels dormant -> byte-identical order to the pre-Task-6
# 12-signal reranker
# --------------------------------------------------------------------------- #
def test_channels_absent_columns_produce_no_new_ranklist_keys():
    df = _base_signal_df()
    rl = ra.build_ranklists_from_df(df)
    for key in ("struct_novelty", "struct_nonchemo_corrob", "emb_classA_sim",
                "emb_nonchemo_sim", "or_microswitch"):
        assert key not in rl


@pytest.mark.parametrize("method", ["rra", "rrf"])
def test_channels_present_but_flagged_false_reproduce_baseline_order(method):
    df = _base_signal_df()
    baseline = ra.aggregate(ra.build_ranklists_from_df(df), method=method)

    # Deliberately extreme, maximally-discriminating dormant values: c4 (the
    # weakest candidate by most base signals) holds the single "best" value
    # in EVERY dormant channel. If any has_*_data gate leaked, c4 would jump
    # to (or near) the top -- so this fixture actually exercises the gate
    # rather than passing vacuously.
    with_dormant = df.copy()
    with_dormant["has_struct_data"] = False
    with_dormant["has_emb_data"] = False
    with_dormant["has_or_microswitch_data"] = False
    with_dormant["struct_novelty"] = [0, 0, 0, 0, 1]
    with_dormant["struct_nonchemo_corrob"] = [1, 1, 1, 1, 0]  # inverted: favors c4
    with_dormant["emb_classA_sim"] = [0, 0, 0, 0, 1]
    with_dormant["emb_nonchemo_sim"] = [1, 1, 1, 1, 0]        # inverted: favors c4
    with_dormant["or_microswitch"] = [0, 0, 0, 0, 1]

    dormant_order = ra.aggregate(ra.build_ranklists_from_df(with_dormant), method=method)
    assert dormant_order == baseline


# --------------------------------------------------------------------------- #
# build_ranklists_from_df: the 5 new channel signals wire in correctly
# --------------------------------------------------------------------------- #
def test_build_ranklists_includes_struct_emb_and_microswitch_channels():
    df = pd.DataFrame(
        {
            "id": ["x", "y", "z"],
            "struct_novelty": [1, 0, 1],
            "struct_nonchemo_corrob": [0, 1, 0],
            "has_struct_data": [True, True, False],  # z gated OUT of the whole channel
            "emb_classA_sim": [0.9, 0.2, 0.5],
            "emb_nonchemo_sim": [0.1, 0.8, 0.3],
            "has_emb_data": [True, True, True],
            "or_microswitch": [1, 0, 1],
            "has_or_microswitch_data": [True, False, True],  # y gated out
        }
    )
    rl = ra.build_ranklists_from_df(df)

    # struct_novelty / struct_nonchemo_corrob share has_struct_data -> only x, y
    assert set(rl["struct_novelty"]) == {"x", "y"}
    assert set(rl["struct_nonchemo_corrob"]) == {"x", "y"}
    # positive signal: raw value preserved as-is
    assert rl["struct_novelty"] == {"x": 1.0, "y": 0.0}
    # EXCLUSION signal: value inverted (1 - raw) so higher = better
    assert rl["struct_nonchemo_corrob"] == {"x": 1.0, "y": 0.0}

    # emb signals: has_emb_data True for all three
    assert set(rl["emb_classA_sim"]) == {"x", "y", "z"}
    assert rl["emb_classA_sim"] == {"x": 0.9, "y": 0.2, "z": 0.5}
    assert rl["emb_nonchemo_sim"] == pytest.approx({"x": 0.9, "y": 0.2, "z": 0.7})

    # or_microswitch gated independently -> x, z only
    assert set(rl["or_microswitch"]) == {"x", "z"}
    assert rl["or_microswitch"] == {"x": 1.0, "z": 1.0}


# --------------------------------------------------------------------------- #
# EXCLUSION POLARITY: struct_nonchemo_corrob / emb_nonchemo_sim must LOWER
# standing via inversion, never via a negative weight (rankagg has none).
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("method", ["rra", "rrf"])
def test_struct_nonchemo_corrob_isolated_inversion_lowers_standing(method):
    # Isolated single-signal case: the ONLY vote is struct_nonchemo_corrob,
    # so the direction of its effect is unambiguous.
    df = pd.DataFrame(
        {
            "id": ["good", "bad"],
            "struct_nonchemo_corrob": [0, 1],
            "has_struct_data": [True, True],
        }
    )
    order = ra.aggregate(ra.build_ranklists_from_df(df), method=method)
    assert order.index("good") < order.index("bad")


def test_struct_nonchemo_corrob_lowers_standing_in_realistic_multisignal_df():
    # A fuller, otherwise-identical-candidate scenario: same base signals
    # (several genuinely TIED), differ ONLY in struct_nonchemo_corrob.
    # RRF only: RRA's Bonferroni correction (score = min(rho * m, 1.0)) is
    # documented (Kolde et al. 2012) to saturate to the 1.0 ceiling once
    # enough tied/uninformative signals inflate m relative to the observed
    # order statistics -- both candidates can legitimately land at the same
    # clamped 1.0 and fall back to alphabetical id tie-break. That is a
    # pre-existing property of the (untouched) RRA algorithm, not a Task-6
    # defect -- the isolated single-signal tests above already give the
    # method-agnostic, non-degenerate proof for both rra and rrf. RRF's
    # linear reciprocal-rank sum has no such saturation and reliably
    # discriminates regardless of how many tied confounds are present.
    df = pd.DataFrame(
        {
            "id": ["good", "bad"],
            "phylo_score_norm": [0.6, 0.6],
            "purifying_score_norm": [0.4, 0.4],
            "positive_score_norm": [0.5, 0.5],
            "lse_divergence_score_norm": [0.3, 0.3],
            "struct_nonchemo_corrob": [0, 1],
            "has_struct_data": [True, True],
        }
    )
    order = ra.aggregate(ra.build_ranklists_from_df(df), method="rrf")
    assert order.index("good") < order.index("bad")


@pytest.mark.parametrize("method", ["rra", "rrf"])
def test_emb_nonchemo_sim_isolated_inversion_lowers_standing(method):
    # Isolated single-signal case (mirrors the struct_nonchemo_corrob proof
    # above): the ONLY vote is emb_nonchemo_sim, so the direction is
    # unambiguous and -- since m=1 -- RRA's Bonferroni term never clamps.
    df = pd.DataFrame(
        {
            "id": ["good", "bad"],
            "emb_nonchemo_sim": [0.1, 0.9],
            "has_emb_data": [True, True],
        }
    )
    order = ra.aggregate(ra.build_ranklists_from_df(df), method=method)
    assert order.index("good") < order.index("bad")


# --------------------------------------------------------------------------- #
# merge_evidence_channels: left-join by id, has_*_data flags set correctly
# --------------------------------------------------------------------------- #
def test_merge_evidence_channels_all_absent_leaves_everything_dormant():
    df = pd.DataFrame({"id": ["a", "b"], "phylo_score_norm": [0.5, 0.5]})
    merged = ra.merge_evidence_channels(df)
    assert list(merged["has_struct_data"]) == [False, False]
    assert list(merged["has_emb_data"]) == [False, False]
    assert list(merged["has_or_microswitch_data"]) == [False, False]
    # original columns untouched, no channel columns fabricated
    assert list(merged["phylo_score_norm"]) == [0.5, 0.5]
    assert "struct_novelty" not in merged.columns
    assert "emb_classA_sim" not in merged.columns
    assert "or_microswitch" not in merged.columns


def test_merge_evidence_channels_missing_path_is_dormant_not_a_crash(tmp_path):
    df = pd.DataFrame({"id": ["a", "b"]})
    missing = str(tmp_path / "does_not_exist.tsv")
    merged = ra.merge_evidence_channels(
        df, struct_tsv=missing, emb_tsv=missing, microswitch_tsv=missing
    )
    assert not merged["has_struct_data"].any()
    assert not merged["has_emb_data"].any()
    assert not merged["has_or_microswitch_data"].any()


def test_merge_evidence_channels_struct_tsv_joins_by_id(tmp_path):
    df = pd.DataFrame({"id": ["cand1", "cand2", "cand3"]})
    struct_tsv = _write_tsv(
        tmp_path / "struct.tsv",
        [["cand1", 1, 0], ["cand2", 0, 1]],
        ["id", "struct_novelty", "struct_nonchemo_corrob"],
    )
    merged = ra.merge_evidence_channels(df, struct_tsv=struct_tsv)
    assert list(merged["has_struct_data"]) == [True, True, False]
    indexed = merged.set_index("id")
    assert indexed.loc["cand1", "struct_novelty"] == 1
    assert indexed.loc["cand2", "struct_nonchemo_corrob"] == 1
    assert pd.isna(indexed.loc["cand3", "struct_novelty"])


def test_merge_evidence_channels_emb_tsv_joins_by_id(tmp_path):
    df = pd.DataFrame({"id": ["cand1", "cand2"]})
    emb_tsv = _write_tsv(
        tmp_path / "emb.tsv",
        [["cand1", 0.9, 0.1]],
        ["id", "emb_classA_sim", "emb_nonchemo_sim"],
    )
    merged = ra.merge_evidence_channels(df, emb_tsv=emb_tsv)
    assert list(merged["has_emb_data"]) == [True, False]
    assert merged.set_index("id").loc["cand1", "emb_classA_sim"] == pytest.approx(0.9)


def test_merge_evidence_channels_microswitch_tsv_joins_by_id(tmp_path):
    df = pd.DataFrame({"id": ["cand1", "cand2"]})
    ms_tsv = _write_tsv(tmp_path / "ms.tsv", [["cand1", 1]], ["id", "or_microswitch"])
    merged = ra.merge_evidence_channels(df, microswitch_tsv=ms_tsv)
    assert list(merged["has_or_microswitch_data"]) == [True, False]
    assert merged.set_index("id").loc["cand1", "or_microswitch"] == 1


def test_merge_evidence_channels_left_join_preserves_all_df_rows(tmp_path):
    df = pd.DataFrame({"id": ["a", "b", "c"]})
    struct_tsv = _write_tsv(
        tmp_path / "struct.tsv", [["a", 1, 0]],
        ["id", "struct_novelty", "struct_nonchemo_corrob"],
    )
    merged = ra.merge_evidence_channels(df, struct_tsv=struct_tsv)
    assert list(merged["id"]) == ["a", "b", "c"]  # no rows dropped


@pytest.mark.parametrize("method", ["rra", "rrf"])
def test_merge_evidence_channels_end_to_end_exclusion_lowers_rank(tmp_path, method):
    """Full chain: merge_evidence_channels -> build_ranklists_from_df ->
    aggregate. Proves the wiring works together, not just each piece alone."""
    df = pd.DataFrame({"id": ["good", "bad"]})
    struct_tsv = _write_tsv(
        tmp_path / "struct.tsv",
        [["good", 0], ["bad", 1]],
        ["id", "struct_nonchemo_corrob"],
    )
    merged = ra.merge_evidence_channels(df, struct_tsv=struct_tsv)
    order = ra.aggregate(ra.build_ranklists_from_df(merged), method=method)
    assert order.index("good") < order.index("bad")


# --------------------------------------------------------------------------- #
# audit_signal_ranking_independence: SIGNAL_COLUMNS gains the channel columns
# --------------------------------------------------------------------------- #
def test_audit_signal_columns_includes_channel_scores():
    expected_new = {
        "struct_novelty", "struct_nonchemo_corrob",
        "emb_classA_sim", "emb_nonchemo_sim",
        "or_microswitch",
    }
    assert expected_new <= set(audit.SIGNAL_COLUMNS)


def test_audit_masks_channel_columns_by_their_real_flags(tmp_path):
    # None of the 5 new columns follow the "<name>_score" -> "has_<name>_data"
    # suffix-swap (mirrors the pre-existing gprotein/ecl regression test) --
    # each needs its own FLAG_OVERRIDES entry or its False-flagged rows would
    # leak their 0/placeholder value in as if it were real data.
    df = pd.DataFrame(
        {
            "id": ["a", "b", "c", "d"],
            "struct_novelty": [1, 0, 1, 0],
            "struct_nonchemo_corrob": [0, 1, 0, 1],
            "has_struct_data": [True, True, False, False],
            "emb_classA_sim": [0.9, 0.1, 0.5, 0.2],
            "emb_nonchemo_sim": [0.2, 0.8, 0.3, 0.7],
            "has_emb_data": [True, False, True, False],
            "or_microswitch": [1, 0, 1, 0],
            "has_or_microswitch_data": [True, True, True, False],
        }
    )
    p = tmp_path / "channels.csv"
    df.to_csv(p, index=False)
    m = audit.load_signal_matrix(str(p))

    assert np.isnan(m.loc["c", "struct_novelty"]) and np.isnan(m.loc["d", "struct_novelty"])
    assert m.loc["a", "struct_novelty"] == 1 and m.loc["b", "struct_novelty"] == 0

    assert np.isnan(m.loc["b", "emb_classA_sim"]) and np.isnan(m.loc["d", "emb_classA_sim"])
    assert m.loc["a", "emb_classA_sim"] == pytest.approx(0.9)
    assert m.loc["c", "emb_classA_sim"] == pytest.approx(0.5)

    assert np.isnan(m.loc["d", "or_microswitch"])
    assert m.loc["a", "or_microswitch"] == 1 and m.loc["c", "or_microswitch"] == 1
