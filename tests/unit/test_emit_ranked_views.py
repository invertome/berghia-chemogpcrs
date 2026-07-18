"""Tests for scripts/emit_ranked_views.py — the two ranked views (bead 1nr).

The pipeline emits one composite-sorted ranked CSV. A single sort conflates
two priorities, so this script splits the SAME scored data into two views:

  * CONFIDENCE view — the safe-bet HCR shortlist: chemoreceptor-candidate,
    not flagged for manual review, evidence-complete; sorted by rank_score.
  * DISCOVERY view — high-novelty divergent-LSE candidates a single composite
    would bury (reference-poor / manual-review / high tandem signal); sorted
    by a divergence-rewarding discovery score.

Column names mirror the ACTUAL augmented ranked CSV (rank_candidates.py
output_cols + add_classification_columns.py + add_og_coverage_columns.py):
id, rank_score, evidence_completeness, classification, needs_manual_review,
og_dnds_reliability, tandem_cluster_score_norm, positive_score_norm.
Values are read as strings in production (dtype=str, keep_default_na=False),
so the fixtures below use string values too.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import emit_ranked_views as erv


def _row(**kw):
    """A ranked-CSV row dict with sane defaults (all values as strings,
    matching the production read). Override any field via keyword."""
    base = {
        "id": "c",
        "rank_score": "0.5",
        "evidence_completeness": "0.9",
        "classification": "chemoreceptor-candidate",
        "needs_manual_review": "",
        "og_dnds_reliability": "high",
        "tandem_cluster_score_norm": "0.0",
        "positive_score_norm": "0.0",
    }
    base.update({k: str(v) for k, v in kw.items()})
    return base


def _df(*rows):
    """Build a string-typed DataFrame from _row() dicts."""
    return pd.DataFrame(list(rows))


# --------------------------------------------------------------------------
# CONFIDENCE view
# --------------------------------------------------------------------------

def test_confidence_keeps_complete_chemoreceptor_candidate():
    df = _df(_row(id="keep", evidence_completeness="0.8"))
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["keep"]


def test_confidence_excludes_below_completeness_threshold():
    df = _df(
        _row(id="hi", evidence_completeness="0.8"),
        _row(id="lo", evidence_completeness="0.5"),
    )
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["hi"]


def test_confidence_excludes_needs_manual_review():
    df = _df(
        _row(id="clean", needs_manual_review=""),
        _row(id="review", needs_manual_review="yes"),
    )
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["clean"]


def test_confidence_excludes_non_chemoreceptors():
    df = _df(
        _row(id="chemo", classification="chemoreceptor-candidate"),
        _row(id="nonchemo", classification="non-chemoreceptor",
             evidence_completeness="0.99"),
        _row(id="likely", classification="likely-non-chemoreceptor",
             evidence_completeness="0.99"),
    )
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["chemo"]


def test_confidence_sorted_by_rank_score_desc():
    df = _df(
        _row(id="mid", rank_score="0.5"),
        _row(id="hi", rank_score="0.9"),
        _row(id="lo", rank_score="0.1"),
    )
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["hi", "mid", "lo"]


def test_confidence_view_honors_final_rank_over_rank_score():
    # rank_candidates emits final_rank = the PRODUCTION order (weighted OR
    # rankagg). The confidence view must sort by it, not re-impose the weighted
    # rank_score order. final_rank here disagrees with rank_score to prove it.
    df = _df(
        _row(id="a", rank_score="0.9", final_rank="3"),
        _row(id="b", rank_score="0.5", final_rank="1"),
        _row(id="c", rank_score="0.7", final_rank="2"),
    )
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["b", "c", "a"]   # final_rank order, not rank_score


def test_confidence_view_falls_back_to_rank_score_without_final_rank():
    # No final_rank column (e.g. an older ranked CSV) -> unchanged behavior:
    # sort by rank_score descending.
    df = _df(
        _row(id="mid", rank_score="0.5"),
        _row(id="hi", rank_score="0.9"),
        _row(id="lo", rank_score="0.1"),
    )
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["hi", "mid", "lo"]


def test_confidence_threshold_is_configurable():
    df = _df(_row(id="c", evidence_completeness="0.6"))
    assert len(erv.build_confidence_view(df, min_completeness=0.7)) == 0
    assert list(erv.build_confidence_view(df, min_completeness=0.5)["id"]) == ["c"]


# --------------------------------------------------------------------------
# DISCOVERY view
# --------------------------------------------------------------------------

def test_discovery_includes_low_dnds_reliability():
    df = _df(_row(id="lowrel", og_dnds_reliability="low"))
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["lowrel"]


def test_discovery_includes_needs_manual_review():
    df = _df(_row(id="rev", og_dnds_reliability="high",
                  needs_manual_review="yes"))
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["rev"]


def test_discovery_includes_high_tandem_signal():
    df = _df(_row(id="tand", og_dnds_reliability="high",
                  needs_manual_review="", tandem_cluster_score_norm="0.9"))
    out = erv.build_discovery_view(df, tandem_high=0.5)
    assert list(out["id"]) == ["tand"]


def test_discovery_excludes_plain_high_confidence_candidate():
    # chemoreceptor-candidate, high reliability, no review, low tandem ->
    # belongs to the confidence shortlist, NOT the discovery view.
    df = _df(_row(id="plain", og_dnds_reliability="high",
                  needs_manual_review="", tandem_cluster_score_norm="0.1"))
    out = erv.build_discovery_view(df, tandem_high=0.5)
    assert list(out["id"]) == []


def test_discovery_excludes_non_chemoreceptors():
    # A non-chemoreceptor that trips every divergence disjunct is STILL
    # excluded — the classification gate is a hard AND.
    df = _df(_row(id="nonchemo", classification="non-chemoreceptor",
                  og_dnds_reliability="low", needs_manual_review="yes",
                  tandem_cluster_score_norm="0.9"))
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == []


def test_discovery_sorted_by_discovery_score_desc():
    df = _df(
        _row(id="low_disc", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.1", positive_score_norm="0.1"),
        _row(id="high_disc", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.9", positive_score_norm="0.9"),
        _row(id="mid_disc", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5"),
    )
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["high_disc", "mid_disc", "low_disc"]


def test_discovery_score_is_weighted_average_of_norm_signals():
    df = _df(_row(id="c", og_dnds_reliability="low",
                  tandem_cluster_score_norm="0.8", positive_score_norm="0.4"))
    out = erv.build_discovery_view(df, tandem_weight=2.0, positive_weight=1.0)
    # (2.0*0.8 + 1.0*0.4) / (2.0 + 1.0) = 2.0 / 3.0
    assert float(out.iloc[0]["discovery_score"]) == pytest.approx(2.0 / 3.0)


# --------------------------------------------------------------------------
# DISCOVERY score — consensus embedding-novelty term (cw3.6)
# --------------------------------------------------------------------------

def test_discovery_score_rewards_novelty_when_present():
    # emb_novelty (unbounded -log10 RRA p) is min-max normalized to [0,1] over the
    # input and added as a weighted discovery-score term. Two otherwise-identical
    # divergent candidates: the one with higher emb_novelty must score higher.
    df = _df(
        _row(id="hi_nov", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5",
             emb_novelty="10.0"),
        _row(id="lo_nov", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5",
             emb_novelty="1.0"),
    )
    out = erv.build_discovery_view(df, tandem_weight=2.0, positive_weight=1.0,
                                   novelty_weight=1.0)
    assert list(out["id"]) == ["hi_nov", "lo_nov"]
    s = dict(zip(out["id"], out["discovery_score"].astype(float)))
    # min-max: hi_nov->1.0, lo_nov->0.0; denom = 2+1+1 = 4
    assert s["hi_nov"] == pytest.approx((2 * .5 + 1 * .5 + 1 * 1.0) / 4)   # 0.625
    assert s["lo_nov"] == pytest.approx((2 * .5 + 1 * .5 + 1 * 0.0) / 4)   # 0.375


def test_discovery_score_ignores_absent_novelty_column():
    # No emb_novelty column -> novelty term dropped from BOTH numerator and
    # denominator; the score is the pre-existing tandem/positive weighted mean.
    df = _df(_row(id="c", og_dnds_reliability="low",
                  tandem_cluster_score_norm="0.8", positive_score_norm="0.4"))
    out = erv.build_discovery_view(df, tandem_weight=2.0, positive_weight=1.0,
                                   novelty_weight=1.0)
    assert float(out.iloc[0]["discovery_score"]) == pytest.approx(2.0 / 3.0)


def test_discovery_novelty_normalized_over_candidates_with_data_only():
    # A no-emb-data candidate (blank emb_novelty) must NOT define the min-max
    # floor: normalization is over candidates WITH data. Real novelties {5, 9};
    # the blank row contributes 0 and the min-data candidate maps to 0 (not a
    # positive value inflated by the blank sitting at 0).
    df = _df(
        _row(id="blank", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.0", positive_score_norm="0.0",
             emb_novelty=""),
        _row(id="lo", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.0", positive_score_norm="0.0",
             emb_novelty="5.0"),
        _row(id="hi", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.0", positive_score_norm="0.0",
             emb_novelty="9.0"),
    )
    out = erv.build_discovery_view(df, tandem_weight=0.0, positive_weight=0.0,
                                   novelty_weight=1.0)
    s = dict(zip(out["id"], out["discovery_score"].astype(float)))
    assert s["hi"] == pytest.approx(1.0)     # max of the with-data range
    assert s["lo"] == pytest.approx(0.0)     # min of the with-data range
    assert s["blank"] == pytest.approx(0.0)  # no data -> 0, doesn't set the floor


def test_discovery_score_ignores_all_blank_novelty_column():
    # Dormant embedding channel: emb_novelty column present but entirely blank.
    # It must be treated as ABSENT (dropped from BOTH numerator and denominator),
    # so the score equals the tandem/positive-only weighted mean — an empty
    # column must not dilute the other signals.
    df = _df(_row(id="c", og_dnds_reliability="low",
                  tandem_cluster_score_norm="0.8", positive_score_norm="0.4",
                  emb_novelty=""))
    out = erv.build_discovery_view(df, tandem_weight=2.0, positive_weight=1.0,
                                   novelty_weight=1.0)
    assert float(out.iloc[0]["discovery_score"]) == pytest.approx(2.0 / 3.0)


def test_main_reads_novelty_weight_from_env(tmp_path: Path, monkeypatch):
    ranked = tmp_path / "ranked.csv"
    # lo_nov listed FIRST: with equal tandem/positive, only the novelty term can
    # reorder it below hi_nov, so this asserts the env weight actually took effect.
    _df(
        _row(id="lo_nov", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5",
             emb_novelty="1.0"),
        _row(id="hi_nov", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5",
             emb_novelty="10.0"),
    ).to_csv(ranked, index=False)
    conf_out = tmp_path / "c.csv"
    disc_out = tmp_path / "d.csv"
    monkeypatch.setenv("DISCOVERY_NOVELTY_WEIGHT", "1.0")
    erv.main(["--ranked-csv", str(ranked),
              "--confidence-out", str(conf_out), "--discovery-out", str(disc_out)])
    disc = pd.read_csv(disc_out, keep_default_na=False, dtype=str)
    assert list(disc["id"]) == ["hi_nov", "lo_nov"]   # novelty broke the tie


# --------------------------------------------------------------------------
# Robustness: a missing OPTIONAL column must not crash
# --------------------------------------------------------------------------

def test_missing_classification_column_does_not_crash():
    # No 'classification' column -> the classification gate is skipped (warn),
    # rows are retained by both views rather than the run crashing.
    df = pd.DataFrame([{
        "id": "c", "rank_score": "0.5", "evidence_completeness": "0.9",
        "needs_manual_review": "", "og_dnds_reliability": "low",
        "tandem_cluster_score_norm": "0.0", "positive_score_norm": "0.0",
    }])
    conf = erv.build_confidence_view(df, min_completeness=0.7)
    disc = erv.build_discovery_view(df)
    assert list(conf["id"]) == ["c"]
    assert list(disc["id"]) == ["c"]


def test_missing_positive_norm_column_discovery_uses_tandem_only():
    # positive_score_norm absent -> that signal is dropped from the discovery
    # score; the denominator collapses to the tandem weight alone (no crash).
    df = pd.DataFrame([{
        "id": "c", "classification": "chemoreceptor-candidate",
        "og_dnds_reliability": "low", "needs_manual_review": "",
        "tandem_cluster_score_norm": "0.6",
    }])
    out = erv.build_discovery_view(df, tandem_weight=2.0, positive_weight=1.0)
    assert float(out.iloc[0]["discovery_score"]) == pytest.approx(0.6)


def test_missing_evidence_completeness_column_skips_that_filter():
    df = pd.DataFrame([{
        "id": "c", "classification": "chemoreceptor-candidate",
        "needs_manual_review": "", "rank_score": "0.5",
    }])
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert list(out["id"]) == ["c"]


# --------------------------------------------------------------------------
# Empty result -> header-only CSV (never crash on zero rows)
# --------------------------------------------------------------------------

def test_confidence_empty_result_writes_header_only(tmp_path: Path):
    df = _df(_row(id="nonchemo", classification="non-chemoreceptor"))
    out = erv.build_confidence_view(df, min_completeness=0.7)
    assert len(out) == 0
    p = tmp_path / "conf.csv"
    out.to_csv(p, index=False)
    lines = p.read_text().splitlines()
    assert len(lines) == 1                 # header only, no data rows
    assert "id" in lines[0]


def test_discovery_empty_result_writes_header_only(tmp_path: Path):
    df = _df(_row(id="plain", og_dnds_reliability="high",
                  needs_manual_review="", tandem_cluster_score_norm="0.0"))
    out = erv.build_discovery_view(df, tandem_high=0.5)
    assert len(out) == 0
    p = tmp_path / "disc.csv"
    out.to_csv(p, index=False)
    lines = p.read_text().splitlines()
    assert len(lines) == 1
    assert "discovery_score" in lines[0]   # score column present even when empty


# --------------------------------------------------------------------------
# main() end-to-end
# --------------------------------------------------------------------------

def test_main_writes_both_views(tmp_path: Path):
    ranked = tmp_path / "ranked.csv"
    _df(
        _row(id="safe", classification="chemoreceptor-candidate",
             evidence_completeness="0.9", needs_manual_review="",
             og_dnds_reliability="high", rank_score="0.9",
             tandem_cluster_score_norm="0.1", positive_score_norm="0.1"),
        _row(id="divergent", classification="chemoreceptor-candidate",
             evidence_completeness="0.3", needs_manual_review="yes",
             og_dnds_reliability="low", rank_score="0.4",
             tandem_cluster_score_norm="0.8", positive_score_norm="0.6"),
        _row(id="nonchemo", classification="non-chemoreceptor",
             evidence_completeness="0.95", needs_manual_review="",
             og_dnds_reliability="high", rank_score="0.99"),
    ).to_csv(ranked, index=False)

    conf_out = tmp_path / "confidence.csv"
    disc_out = tmp_path / "discovery.csv"
    rc = erv.main([
        "--ranked-csv", str(ranked),
        "--confidence-out", str(conf_out),
        "--discovery-out", str(disc_out),
    ])
    assert rc == 0
    conf = pd.read_csv(conf_out, keep_default_na=False, dtype=str)
    disc = pd.read_csv(disc_out, keep_default_na=False, dtype=str)
    assert list(conf["id"]) == ["safe"]        # safe-bet shortlist
    assert list(disc["id"]) == ["divergent"]   # divergent-LSE view


def test_main_reads_threshold_from_env(tmp_path: Path, monkeypatch):
    ranked = tmp_path / "ranked.csv"
    _df(_row(id="c", evidence_completeness="0.6", rank_score="0.5")).to_csv(
        ranked, index=False)
    conf_out = tmp_path / "c.csv"
    disc_out = tmp_path / "d.csv"
    monkeypatch.setenv("CONFIDENCE_MIN_COMPLETENESS", "0.5")
    erv.main([
        "--ranked-csv", str(ranked),
        "--confidence-out", str(conf_out),
        "--discovery-out", str(disc_out),
    ])
    conf = pd.read_csv(conf_out, keep_default_na=False, dtype=str)
    assert list(conf["id"]) == ["c"]   # 0.6 >= env threshold 0.5
