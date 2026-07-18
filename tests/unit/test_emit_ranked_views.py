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
        "emb_novelty": "0.0",
        "lse_depth_score": "0.0",
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
    # Monotone: high_disc dominates every signal, low_disc is worst -> the RRA
    # rho is strictly ordered, so discovery_score sorts high > mid > low.
    df = _df(
        _row(id="low_disc", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.1", positive_score_norm="0.1",
             emb_novelty="0.1", lse_depth_score="0.1"),
        _row(id="high_disc", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.9", positive_score_norm="0.9",
             emb_novelty="0.9", lse_depth_score="0.9"),
        _row(id="mid_disc", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5",
             emb_novelty="0.5", lse_depth_score="0.5"),
    )
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["high_disc", "mid_disc", "low_disc"]


# --------------------------------------------------------------------------
# DISCOVERY score — weight-free Robust Rank Aggregation (Task 2, cw-lkhu)
# --------------------------------------------------------------------------

def test_discovery_score_is_rra_over_ranklists():
    # The discovery view is ordered by a WEIGHT-FREE Robust Rank Aggregation over
    # {tandem, positive, novelty, lse_depth}, with discovery_score = -log10(rho).
    # Each candidate's score must equal -log10 of rank_aggregation.rra_score over
    # exactly those four signals -- proving all four participate and no weight
    # enters. rho is computed WITH lse_depth, so an impl that dropped lse_depth
    # (as the old weighted mean did) could not reproduce these values.
    import numpy as np
    from rank_aggregation import rra_score
    df = _df(
        _row(id="a", og_dnds_reliability="low", tandem_cluster_score_norm="0.9",
             positive_score_norm="0.8", emb_novelty="0.7", lse_depth_score="0.9"),
        _row(id="b", og_dnds_reliability="low", tandem_cluster_score_norm="0.5",
             positive_score_norm="0.5", emb_novelty="0.5", lse_depth_score="0.5"),
        _row(id="c", og_dnds_reliability="low", tandem_cluster_score_norm="0.1",
             positive_score_norm="0.2", emb_novelty="0.3", lse_depth_score="0.1"),
    )
    out = erv.build_discovery_view(df)
    per_signal = {
        "tandem":    {"a": 0.9, "b": 0.5, "c": 0.1},
        "positive":  {"a": 0.8, "b": 0.5, "c": 0.2},
        "novelty":   {"a": 0.7, "b": 0.5, "c": 0.3},
        "lse_depth": {"a": 0.9, "b": 0.5, "c": 0.1},
    }
    rho = rra_score(per_signal)
    got = dict(zip(out["id"], out["discovery_score"].astype(float)))
    for cid in ("a", "b", "c"):
        assert got[cid] == pytest.approx(-np.log10(max(rho[cid], 1e-300)))
    # ordering follows the RRA: lower rho == higher discovery_score, best first
    assert list(out["id"]) == sorted(rho, key=lambda i: rho[i])


def test_discovery_score_includes_novelty_signal():
    # emb_novelty participates in the RRA. With novelty the only varying signal
    # (the other signal columns dropped), the discovery order follows it. Input
    # order is scrambled so a no-op impl would NOT reproduce the expected order.
    df = _df(
        _row(id="lo_nov", og_dnds_reliability="low", emb_novelty="1.0"),
        _row(id="hi_nov", og_dnds_reliability="low", emb_novelty="9.0"),
        _row(id="mid_nov", og_dnds_reliability="low", emb_novelty="5.0"),
    ).drop(columns=["tandem_cluster_score_norm", "positive_score_norm",
                    "lse_depth_score"])
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["hi_nov", "mid_nov", "lo_nov"]


def test_discovery_score_includes_lse_depth_signal():
    # lse_depth participates in the RRA -- it was EXCLUDED from the old weighted
    # mean (the gap Task 2 closes). With lse_depth the only varying signal, the
    # discovery order follows it. Scrambled input order proves it is real sorting.
    df = _df(
        _row(id="shallow", og_dnds_reliability="low", lse_depth_score="0.1"),
        _row(id="deep", og_dnds_reliability="low", lse_depth_score="0.9"),
        _row(id="mid", og_dnds_reliability="low", lse_depth_score="0.5"),
    ).drop(columns=["tandem_cluster_score_norm", "positive_score_norm",
                    "emb_novelty"])
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["deep", "mid", "shallow"]


def test_discovery_score_handles_absent_novelty_column():
    # No emb_novelty column -> RRA runs over the remaining signals (tandem,
    # positive, lse_depth); graceful, ordered, finite scores.
    df = _df(
        _row(id="hi", og_dnds_reliability="low", tandem_cluster_score_norm="0.9",
             positive_score_norm="0.9", lse_depth_score="0.9"),
        _row(id="lo", og_dnds_reliability="low", tandem_cluster_score_norm="0.1",
             positive_score_norm="0.1", lse_depth_score="0.1"),
    ).drop(columns=["emb_novelty"])
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["hi", "lo"]
    assert out["discovery_score"].astype(float).notna().all()


def test_discovery_score_all_blank_novelty_column():
    # Dormant embedding channel: emb_novelty present but entirely blank. No id
    # votes in novelty, so RRA runs over the remaining signals -- no crash, no
    # dilution of the order.
    df = _df(
        _row(id="hi", og_dnds_reliability="low", tandem_cluster_score_norm="0.9",
             positive_score_norm="0.9", lse_depth_score="0.9", emb_novelty=""),
        _row(id="lo", og_dnds_reliability="low", tandem_cluster_score_norm="0.1",
             positive_score_norm="0.1", lse_depth_score="0.1", emb_novelty=""),
    )
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["hi", "lo"]


def test_discovery_blank_novelty_row_does_not_crash():
    # A blank emb_novelty on ONE candidate (mixed with data rows) simply drops
    # that id from the novelty vote; both rows are retained and ranked.
    df = _df(
        _row(id="hasnov", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5",
             lse_depth_score="0.5", emb_novelty="9.0"),
        _row(id="blank", og_dnds_reliability="low",
             tandem_cluster_score_norm="0.5", positive_score_norm="0.5",
             lse_depth_score="0.5", emb_novelty=""),
    )
    out = erv.build_discovery_view(df)
    assert set(out["id"]) == {"hasnov", "blank"}


def test_build_discovery_view_has_no_weight_params():
    # The hand-picked discovery weights are removed from the signature; only the
    # tandem membership cutoff remains.
    import inspect
    params = inspect.signature(erv.build_discovery_view).parameters
    assert "tandem_weight" not in params
    assert "positive_weight" not in params
    assert "novelty_weight" not in params
    assert "tandem_high" in params


def test_discovery_order_ignores_weight_env(tmp_path: Path, monkeypatch):
    # The DISCOVERY_*_WEIGHT env vars no longer feed the discovery order emitted
    # by main(): cranking DISCOVERY_TANDEM_WEIGHT (which under the old weighted
    # mean would surface the tandem-strong candidate 'a') must NOT reorder the
    # weight-free RRA output.
    ranked = tmp_path / "ranked.csv"
    _df(
        _row(id="a", og_dnds_reliability="low", tandem_cluster_score_norm="0.9",
             positive_score_norm="0.1", emb_novelty="0.1", lse_depth_score="0.1"),
        _row(id="b", og_dnds_reliability="low", tandem_cluster_score_norm="0.1",
             positive_score_norm="0.9", emb_novelty="0.9", lse_depth_score="0.9"),
    ).to_csv(ranked, index=False)

    def run():
        c = tmp_path / "c.csv"
        d = tmp_path / "d.csv"
        erv.main(["--ranked-csv", str(ranked), "--confidence-out", str(c),
                  "--discovery-out", str(d)])
        return list(pd.read_csv(d, keep_default_na=False, dtype=str)["id"])

    base = run()
    monkeypatch.setenv("DISCOVERY_TANDEM_WEIGHT", "99")
    assert run() == base


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


def test_missing_positive_norm_column_discovery_no_crash():
    # positive_score_norm / emb_novelty / lse_depth all absent -> RRA runs over
    # the single present signal (tandem); no crash, the row is retained with a
    # finite discovery_score.
    df = pd.DataFrame([{
        "id": "c", "classification": "chemoreceptor-candidate",
        "og_dnds_reliability": "low", "needs_manual_review": "",
        "tandem_cluster_score_norm": "0.6",
    }])
    out = erv.build_discovery_view(df)
    assert list(out["id"]) == ["c"]
    assert pd.notna(float(out.iloc[0]["discovery_score"]))


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
