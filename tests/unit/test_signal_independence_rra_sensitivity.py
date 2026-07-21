"""Tests for scripts/audit_rra_correlation_sensitivity.py.

Context: an expert council objected that the ~12 per-signal ranklists feeding
the production Robust Rank Aggregation (RRA) are NOT independent, and that
correlated ranklists "add weight to whatever the shared upstream artifact
says" rather than adding evidence.

`scripts/audit_signal_ranking_independence.py` already measures the empirical
correlation matrix, but it has two blind spots this module's subject closes:

  1. It reads the RAW ``<signal>_score`` columns, while rank aggregation votes
     on the ``<signal>_score_norm`` columns where those exist. The dN/dS
     reliability shrink (rank_candidates.py, `reliability_shrink`) writes ONLY
     the _norm columns, and that shrink is itself a shared computed quantity
     applied to BOTH selection signals -- so it is precisely a correlation
     source the raw-column audit cannot see.
  2. It never quantifies the CONSEQUENCE. A flagged pair at |rho|=0.75 may move
     nothing; an unflagged pair at |rho|=0.6 may move the shortlist. The 0.7
     grouping threshold is otherwise unfalsifiable.

So the subject measures the values RRA ACTUALLY votes on (by reusing
rank_aggregation.build_ranklists_from_df, which resolves _norm preference,
has_*_data gating and exclusion-signal inversion in one place -- no drift), and
reports what correlation does to THIS RRA implementation: duplication
inflation, tie saturation, and the fused-vs-unfused order delta.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import audit_rra_correlation_sensitivity as sens

SCRIPT = Path(__file__).resolve().parent.parent.parent / "scripts" / \
    "audit_rra_correlation_sensitivity.py"


# ---------------------------------------------------------------------------
# fact 1: the matrix reflects what RRA votes on, not the raw CSV columns
# ---------------------------------------------------------------------------

def test_signal_matrix_uses_the_norm_column_rra_actually_votes_on():
    # purifying_score_norm is the reliability-SHRUNK column; the raw
    # purifying_score is pre-shrink. RRA reads the _norm one, so the audit
    # matrix must too, or it measures a signal nobody votes with.
    df = pd.DataFrame({
        "id": ["a", "b", "c"],
        "purifying_score": [0.1, 0.2, 0.3],
        "purifying_score_norm": [0.9, 0.8, 0.7],
    })
    m = sens.signal_matrix(df)
    assert list(m["purifying"]) == [0.9, 0.8, 0.7], (
        "signal_matrix must prefer <signal>_score_norm, matching "
        "rank_aggregation.build_ranklists_from_df's column resolution"
    )


def test_signal_matrix_applies_exclusion_signal_inversion():
    # emb_nonchemo_sim is an EXCLUSION signal (invert=True in SIGNAL_SPEC):
    # RRA stores 1-raw so that "higher = better" holds for every signal.
    # Inversion flips the sign of every correlation this signal takes part in.
    df = pd.DataFrame({
        "id": ["a", "b"],
        "emb_nonchemo_sim": [0.1, 0.9],
        "has_emb_data": [True, True],
    })
    m = sens.signal_matrix(df)
    assert m["emb_nonchemo_sim"].tolist() == pytest.approx([0.9, 0.1])


def test_signal_matrix_masks_rows_whose_flag_is_false():
    # a gated signal only votes where has_*_data is True; elsewhere it must be
    # NaN, never a 0.0 no-data sentinel masquerading as a real low score.
    df = pd.DataFrame({
        "id": ["a", "b", "c"],
        "ecl_divergence_score": [0.0, 0.5, 0.7],
        "has_ecl_data": [False, True, True],
    })
    m = sens.signal_matrix(df)
    assert np.isnan(m.loc["a", "ecl_divergence"])
    assert m.loc["b", "ecl_divergence"] == 0.5


def test_signal_matrix_is_indexed_by_candidate_id():
    df = pd.DataFrame({"id": ["x", "y"], "phylo_score": [1.0, 2.0]})
    m = sens.signal_matrix(df)
    assert list(m.index) == ["x", "y"]


# ---------------------------------------------------------------------------
# fact 2: correlation measurement
# ---------------------------------------------------------------------------

def test_rank_correlation_matrix_detects_a_perfect_duplicate():
    m = pd.DataFrame({"a": [1.0, 2, 3, 4, 5], "b": [10.0, 20, 30, 40, 50],
                      "c": [5.0, 1, 4, 2, 3]}, index=list("vwxyz"))
    corr = sens.rank_correlation_matrix(m)
    assert corr.loc["a", "b"] == pytest.approx(1.0)
    assert abs(corr.loc["a", "c"]) < 0.9


def test_rank_correlation_matrix_uses_pairwise_complete_rows():
    # a NaN in one signal must not poison that signal's correlation with
    # others -- only the affected PAIR's overlapping rows are used.
    m = pd.DataFrame({"a": [1.0, 2, 3, 4, 5],
                      "b": [np.nan, 20, 30, 40, 50],
                      "c": [1.0, 2, 3, 4, 5]}, index=list("vwxyz"))
    corr = sens.rank_correlation_matrix(m)
    assert corr.loc["a", "b"] == pytest.approx(1.0)
    assert corr.loc["a", "c"] == pytest.approx(1.0)


def test_flag_correlated_pairs_reports_only_pairs_above_threshold():
    m = pd.DataFrame({"a": [1.0, 2, 3, 4, 5], "b": [1.0, 2, 3, 4, 5],
                      "c": [5.0, 1, 4, 2, 3]}, index=list("vwxyz"))
    pairs = sens.flag_correlated_pairs(sens.rank_correlation_matrix(m),
                                       threshold=0.7)
    assert [(p[0], p[1]) for p in pairs] == [("a", "b")]


def test_flag_correlated_pairs_catches_negative_correlation():
    # an EXCLUSION signal that mirrors another axis is just as redundant as a
    # positively correlated one; grouping keys on |rho|.
    m = pd.DataFrame({"a": [1.0, 2, 3, 4, 5], "b": [5.0, 4, 3, 2, 1]},
                     index=list("vwxyz"))
    pairs = sens.flag_correlated_pairs(sens.rank_correlation_matrix(m),
                                       threshold=0.7)
    assert len(pairs) == 1 and pairs[0][2] == pytest.approx(-1.0)


# ---------------------------------------------------------------------------
# fact 3: what correlation does to THIS RRA implementation
# ---------------------------------------------------------------------------

def test_duplicating_a_signal_inflates_its_favourites_rra_significance():
    # The council's mechanism, made concrete. RRA's beta order statistic
    # assumes INDEPENDENT Uniform(0,1) draws. Feeding the same ranklist twice
    # violates that, so a candidate the duplicated axis likes gets a smaller
    # (more "significant") rho than the evidence warrants.
    base = {
        "s1": {c: v for c, v in zip("abcdef", [6, 5, 4, 3, 2, 1])},
        "s2": {c: v for c, v in zip("abcdef", [1, 2, 3, 4, 5, 6])},
    }
    dup = dict(base, s1_copy=dict(base["s1"]))
    assert sens.rra_with(dup)["a"] < sens.rra_with(base)["a"], (
        "duplicating s1 must make s1's favourite 'a' score better (lower), "
        "which is exactly the double-counting the council warned about"
    )


def test_group_fusion_undoes_the_duplication_inflation():
    # rank_aggregation fuses a group into ONE effective list, so the duplicate
    # casts one vote. Fusing s1 with its copy must restore the ungrouped score.
    base = {
        "s1": {c: v for c, v in zip("abcdef", [6, 5, 4, 3, 2, 1])},
        "s2": {c: v for c, v in zip("abcdef", [1, 2, 3, 4, 5, 6])},
    }
    dup = dict(base, s1_copy=dict(base["s1"]))
    fused = sens.rra_with(dup, groups=[["s1", "s1_copy"]])
    assert fused["a"] == pytest.approx(sens.rra_with(base)["a"]), (
        "fusing an exact duplicate must reproduce the un-duplicated RRA score"
    )


def test_saturation_report_counts_ties_not_the_removed_bonferroni_cap():
    # Bead 8k8e replaced min(rho*m, 1.0) with the exact null, so "saturated at
    # 1.0" is now reachable ONLY by a candidate ranked dead last in every list
    # (rho == 1 exactly) -- here that is 'f'. What remains a first-class health
    # metric is TIE saturation: candidates aggregate() can only order by id.
    # Six candidates ranked identically by two identical signals -> no ties at
    # all, and exactly one genuinely-saturated candidate.
    lists = {
        "s1": {c: v for c, v in zip("abcdef", [6, 5, 4, 3, 2, 1])},
        "s2": {c: v for c, v in zip("abcdef", [6, 5, 4, 3, 2, 1])},
    }
    rep = sens.saturation_report(lists)
    assert rep["n_candidates"] == 6
    assert rep["n_saturated"] == 1          # only 'f', whose rho is exactly 1
    assert 0.0 <= rep["fraction_saturated"] <= 1.0
    assert rep["largest_tie_block"] == 1    # every score distinct
    assert rep["n_tied"] == 0
    assert rep["fraction_tied"] == 0.0


def test_saturation_report_counts_genuinely_tied_candidates():
    # Two candidates with identical rank vectors are a GENUINE tie: they stay
    # tied, and the report must say how many candidates sit in such blocks.
    lists = {
        "s1": {"a": 9.0, "b": 5.0, "c": 5.0, "d": 1.0},
        "s2": {"a": 9.0, "b": 5.0, "c": 5.0, "d": 1.0},
    }
    rep = sens.saturation_report(lists)
    assert rep["n_candidates"] == 4
    assert rep["largest_tie_block"] == 2
    assert rep["n_tied"] == 2
    assert rep["fraction_tied"] == pytest.approx(0.5)


def test_fusion_impact_is_a_no_op_for_singleton_groups():
    lists = {
        "s1": {c: v for c, v in zip("abcdef", [6, 5, 4, 3, 2, 1])},
        "s2": {c: v for c, v in zip("abcdef", [1, 3, 2, 6, 4, 5])},
    }
    imp = sens.fusion_impact(lists, groups=[["s1"], ["s2"]], top_k=3)
    assert imp["max_displacement"] == 0
    assert imp["top_k_jaccard"] == pytest.approx(1.0)
    assert imp["order_spearman"] == pytest.approx(1.0)


def test_fusion_impact_detects_reordering_when_a_real_group_fuses():
    # three lists, two of them near-duplicates: fusing them must change the
    # order, and the tool must report a non-zero displacement.
    lists = {
        "s1": {c: v for c, v in zip("abcdef", [6, 5, 4, 3, 2, 1])},
        "s1b": {c: v for c, v in zip("abcdef", [6, 5, 4, 3, 1, 2])},
        "s2": {c: v for c, v in zip("abcdef", [1, 2, 3, 4, 5, 6])},
    }
    imp = sens.fusion_impact(lists, groups=[["s1", "s1b"]], top_k=3)
    assert imp["max_displacement"] > 0
    assert imp["n_top_k_changed"] >= 0
    assert -1.0 <= imp["order_spearman"] <= 1.0


# ---------------------------------------------------------------------------
# fact 4: reporting + CLI
# ---------------------------------------------------------------------------

def _ranked_csv(tmp_path):
    rng = np.random.default_rng(0)
    n = 12
    ids = [f"c{i:02d}" for i in range(n)]
    phylo = np.linspace(1, 0, n)
    df = pd.DataFrame({
        "id": ids,
        "phylo_score": phylo,
        # og_confidence deliberately near-duplicates phylo -> must be flagged
        "og_confidence_score": phylo + rng.normal(0, 0.01, n),
        "has_og_confidence_data": [True] * n,
        "synteny_score": rng.permutation(n).astype(float),
        "has_synteny_data": [True] * n,
    })
    p = tmp_path / "ranked.csv"
    df.to_csv(p, index=False)
    return str(p)


def test_write_report_emits_tsv_json_and_markdown(tmp_path):
    df = pd.read_csv(_ranked_csv(tmp_path))
    res = sens.analyze(df, threshold=0.7, top_k=3)
    sens.write_report(res, str(tmp_path / "rra"))
    assert (tmp_path / "rra_correlation.tsv").exists()
    payload = json.loads((tmp_path / "rra_sensitivity.json").read_text())
    assert "flagged_pairs" in payload and "saturation" in payload
    assert "fusion_impact" in payload and payload["threshold"] == 0.7
    assert (tmp_path / "rra_sensitivity.md").exists()


def test_analyze_flags_the_planted_near_duplicate(tmp_path):
    df = pd.read_csv(_ranked_csv(tmp_path))
    res = sens.analyze(df, threshold=0.7, top_k=3)
    flagged = {frozenset((a, b)) for a, b, _ in res["flagged_pairs"]}
    assert frozenset(("phylo", "og_confidence")) in flagged


def test_cli_smoke_runs_end_to_end(tmp_path):
    csv = _ranked_csv(tmp_path)
    out = tmp_path / "cli"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--ranked-csv", csv,
         "--out-prefix", str(out), "--threshold", "0.7", "--top-k", "3"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "cli_sensitivity.json").exists()


# --------------------------------------------------------------------------- #
# Bead 8k8e requirement 4: detection that does not report is not detection.
#
# saturation_report() had ZERO callers -- stage 07 ran only
# audit_signal_ranking_independence.py, so the tie/saturation health of the
# shipped ranking was measured nowhere. Pin the wiring.
# --------------------------------------------------------------------------- #
STAGE_07 = Path(__file__).resolve().parent.parent.parent / "07_candidate_ranking.sh"


def test_stage07_invokes_the_rra_sensitivity_audit():
    text = STAGE_07.read_text()
    assert "audit_rra_correlation_sensitivity.py" in text, (
        "stage 07 must run the RRA sensitivity audit -- saturation_report() is "
        "the only thing that measures how much of the shipped order is decided "
        "by ties rather than by evidence"
    )


def test_stage07_audits_the_final_ranked_csv_with_the_derived_groups():
    text = STAGE_07.read_text()
    idx = text.index("audit_rra_correlation_sensitivity.py")
    # the invocation's argv is assembled just above the call, so window both ways
    block = text[max(0, idx - 900):idx + 400]
    assert "--ranked-csv" in block and "RANKED_CSV" in block, block
    assert "--out-prefix" in block, block
    # must reuse the independence audit's groups so the two audits agree
    assert "signal_independence_groups.json" in text
