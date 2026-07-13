"""Tests for wiring emb_novelty (S_novel, the Mahalanobis channel's positive
ranking axis) into rank_aggregation as a positive voter.

Phase-0 wire-up (bead cw3): the maha embedding channel emits emb_novelty (high =
divergent from known families = the candidates the user wants surfaced). It joins
as a POSITIVE weightless voter, and stays dormant when the cosine channel (no
emb_novelty column) ran instead — so nothing changes unless the maha channel is
used.
"""
from __future__ import annotations

import pandas as pd

from rank_aggregation import (
    SIGNAL_SPEC,
    build_ranklists_from_df,
    merge_evidence_channels,
)


def test_emb_novelty_is_a_positive_voter_in_signal_spec():
    entries = [s for s in SIGNAL_SPEC if s[0] == "emb_novelty"]
    assert entries, "emb_novelty missing from SIGNAL_SPEC"
    key, flag, column, invert = entries[0]
    assert flag == "has_emb_data"
    assert column == "emb_novelty"
    assert invert is False   # positive: high novelty ranks UP (not inverted)


def test_build_ranklists_ranks_higher_novelty_up_when_present():
    df = pd.DataFrame(
        {"id": ["a", "b"], "emb_novelty": [0.2, 0.9], "has_emb_data": [True, True]}
    )
    ranklists = build_ranklists_from_df(df)
    assert "emb_novelty" in ranklists
    assert ranklists["emb_novelty"]["b"] > ranklists["emb_novelty"]["a"]


def test_emb_novelty_dormant_when_column_absent():
    # cosine channel / no maha channel: no emb_novelty column -> signal dormant,
    # nothing added (existing rankings unchanged).
    df = pd.DataFrame({"id": ["a"], "phylo_score": [0.5]})
    assert "emb_novelty" not in build_ranklists_from_df(df)


def test_merge_joins_emb_novelty_from_a_maha_channel_tsv(tmp_path):
    df = pd.DataFrame({"id": ["a", "b"]})
    emb = tmp_path / "emb.tsv"
    pd.DataFrame({
        "id": ["a", "b"],
        "emb_nonchemo_family": ["opsin", "peptide"],
        "emb_novelty": [0.1, 0.8],
        "has_emb_data": [True, True],
        "emb_leakage_flag": [True, True],
    }).to_csv(emb, sep="\t", index=False)
    merged = merge_evidence_channels(df, emb_tsv=str(emb))
    assert "emb_novelty" in merged.columns
    assert merged.set_index("id").loc["b", "emb_novelty"] == 0.8
    assert bool(merged.set_index("id").loc["a", "has_emb_data"]) is True
