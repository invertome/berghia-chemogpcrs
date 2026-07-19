"""Tests for scripts/add_embedding_columns.py — always-present emb_novelty column.

cw3.6 item 3. The consensus embedding channel (fusion_consensus.py /
build_embedding_channel.py) writes a per-candidate TSV with an ``emb_novelty``
score. Under the default RANK_METHOD=weighted path that channel is never merged
into the ranked CSV (it only feeds the rankagg voter), so novelty is invisible
downstream. This augmenter left-joins ``emb_novelty`` (+ companion emb columns)
into the ranked CSV as an ALWAYS-PRESENT sortable column, regardless of
RANK_METHOD, so both ranked views (and the discovery score) can use it.

Values are read as strings in production (dtype=str, keep_default_na=False).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import add_embedding_columns as aec


def _ranked(tmp_path: Path, rows) -> str:
    p = tmp_path / "ranked.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return str(p)


def _channel(tmp_path: Path, rows) -> str:
    p = tmp_path / "embedding_channel.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return str(p)


def test_joins_emb_novelty_from_channel(tmp_path: Path):
    ranked = _ranked(tmp_path, [
        {"id": "c1", "rank_score": "0.9"},
        {"id": "c2", "rank_score": "0.5"},
    ])
    channel = _channel(tmp_path, [
        {"id": "c1", "emb_nonchemo_family": "amine", "emb_novelty": "7.5",
         "has_emb_data": "True", "emb_leakage_flag": "True"},
        {"id": "c2", "emb_nonchemo_family": "opsin", "emb_novelty": "1.2",
         "has_emb_data": "True", "emb_leakage_flag": "True"},
    ])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked, channel_tsv_path=channel,
                              out_path=str(out))
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert dict(zip(df["id"], df["emb_novelty"])) == {"c1": "7.5", "c2": "1.2"}
    assert dict(zip(df["id"], df["emb_nonchemo_family"])) == {"c1": "amine", "c2": "opsin"}


def test_joins_emb_novelty_residual_when_present(tmp_path: Path):
    # A1 (v4bs.2): the phylogeny-residualized novelty is a dormant descriptive
    # column — when the channel carries emb_novelty_residual it must be surfaced
    # in the ranked CSV alongside emb_novelty.
    ranked = _ranked(tmp_path, [{"id": "c1", "rank_score": "0.9"},
                                {"id": "c2", "rank_score": "0.5"}])
    channel = _channel(tmp_path, [
        {"id": "c1", "emb_nonchemo_family": "amine", "emb_novelty": "7.5",
         "emb_novelty_residual": "6.1", "has_emb_data": "True",
         "emb_leakage_flag": "True"},
        {"id": "c2", "emb_nonchemo_family": "opsin", "emb_novelty": "1.2",
         "emb_novelty_residual": "", "has_emb_data": "True",
         "emb_leakage_flag": "True"},
    ])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked, channel_tsv_path=channel,
                              out_path=str(out))
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert "emb_novelty_residual" in df.columns
    assert dict(zip(df["id"], df["emb_novelty_residual"])) == {"c1": "6.1", "c2": ""}


def _extra(tmp_path: Path, name: str, rows) -> str:
    p = tmp_path / name
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return str(p)


def test_extra_tsv_joins_its_columns_as_is(tmp_path: Path):
    # A3 (model_role_split) columns are already self-describing -> join unprefixed.
    ranked = _ranked(tmp_path, [{"id": "c1", "rank_score": "0.9"},
                                {"id": "c2", "rank_score": "0.5"}])
    extra = _extra(tmp_path, "role.tsv", [
        {"id": "c1", "emb_novelty_phylo": "0.2", "emb_role_gap": "0.7"},
        {"id": "c2", "emb_novelty_phylo": "0.8", "emb_role_gap": "-0.3"},
    ])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked, channel_tsv_path="",
                              out_path=str(out), extra_tsvs=[extra])
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert dict(zip(df["id"], df["emb_role_gap"])) == {"c1": "0.7", "c2": "-0.3"}
    assert dict(zip(df["id"], df["emb_novelty_phylo"])) == {"c1": "0.2", "c2": "0.8"}


def test_extra_tsv_prefix_namespaces_generic_column_names(tmp_path: Path):
    # A4 (embedding_nulls) emits generic names (p_value, collapsed) that would be
    # ambiguous in a ranked CSV -> PREFIX:PATH namespaces them.
    ranked = _ranked(tmp_path, [{"id": "c1", "rank_score": "0.9"}])
    extra = _extra(tmp_path, "nulls.tsv", [
        {"id": "c1", "p_value": "0.01", "collapsed": "False"},
    ])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked, channel_tsv_path="",
                              out_path=str(out), extra_tsvs=[f"emb_null:{extra}"])
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert df["emb_null_p_value"].iloc[0] == "0.01"
    assert df["emb_null_collapsed"].iloc[0] == "False"
    assert "p_value" not in df.columns          # un-namespaced name not leaked


def test_missing_extra_tsv_is_dormant_not_an_error(tmp_path: Path):
    ranked = _ranked(tmp_path, [{"id": "c1", "rank_score": "0.9"}])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked, channel_tsv_path="",
                              out_path=str(out),
                              extra_tsvs=[str(tmp_path / "nope.tsv")])
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert list(df["id"]) == ["c1"]             # rows preserved, no crash


def test_candidate_missing_from_channel_gets_blank(tmp_path: Path):
    ranked = _ranked(tmp_path, [{"id": "c1", "rank_score": "0.9"},
                                {"id": "c2", "rank_score": "0.5"}])
    channel = _channel(tmp_path, [
        {"id": "c1", "emb_nonchemo_family": "amine", "emb_novelty": "7.5",
         "has_emb_data": "True", "emb_leakage_flag": "True"},
    ])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked, channel_tsv_path=channel,
                              out_path=str(out))
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    row = df[df["id"] == "c2"].iloc[0]
    assert row["emb_novelty"] == ""          # no channel row -> blank, not a crash


def test_column_always_present_when_channel_absent(tmp_path: Path):
    # No channel TSV (embedding producer dormant): emb_novelty must STILL exist as
    # an (empty) column so downstream views/scores never KeyError on it.
    ranked = _ranked(tmp_path, [{"id": "c1", "rank_score": "0.9"}])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked,
                              channel_tsv_path=str(tmp_path / "nope.tsv"),
                              out_path=str(out))
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert "emb_novelty" in df.columns
    assert list(df["emb_novelty"]) == [""]


def test_channel_is_authoritative_no_duplicate_column(tmp_path: Path):
    # rankagg may already have merged an emb_novelty column; re-running the
    # augmenter must OVERWRITE from the channel (authoritative) and leave exactly
    # one emb_novelty column, never two.
    ranked = _ranked(tmp_path, [{"id": "c1", "rank_score": "0.9",
                                 "emb_novelty": "999"}])
    channel = _channel(tmp_path, [
        {"id": "c1", "emb_nonchemo_family": "amine", "emb_novelty": "7.5",
         "has_emb_data": "True", "emb_leakage_flag": "True"},
    ])
    out = tmp_path / "out.csv"
    aec.add_embedding_columns(ranked_csv_path=ranked, channel_tsv_path=channel,
                              out_path=str(out))
    df = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert list(df.columns).count("emb_novelty") == 1
    assert df.iloc[0]["emb_novelty"] == "7.5"
