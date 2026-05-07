"""Tests for scripts/validate_classification_against_aplysia.py.

Phase 6 Task 6.3 — retrospective validation harness. Reads a ranked
candidates CSV (with classification columns from stage 07) and the
curated Aplysia truth CSV; computes per-entry agreement and overall
accuracy / precision / recall metrics. Outputs a markdown-style
validation report.

Tests cover the local logic — no pipeline-output CSV needed.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import validate_classification_against_aplysia as vca


def _make_truth_csv(tmp_path: Path) -> str:
    """Compact truth CSV (3 chemoreceptors + 2 non-chemo)."""
    f = tmp_path / "truth.csv"
    f.write_text(
        "aplysia_protein_id,expected_classification,expected_family,"
        "expected_subfamily,source_citation,notes\n"
        "NP_001191513.1,chemoreceptor-candidate,,,Cummins 2009,A\n"
        "NP_001191494.1,chemoreceptor-candidate,,,Cummins 2009,B\n"
        "NP_001191495.1,chemoreceptor-candidate,,,Cummins 2009,C\n"
        "Q16950,non-chemoreceptor,aminergic,5HT,Swiss-Prot,5HT-1\n"
        "Q16951,non-chemoreceptor,aminergic,5HT,Swiss-Prot,5HT-2\n"
    )
    return str(f)


def _make_ranked_csv(tmp_path: Path, rows: list[dict]) -> str:
    f = tmp_path / "ranked.csv"
    pd.DataFrame(rows).to_csv(f, index=False)
    return str(f)


# ---- truth-set loading --------------------------------------------------

def test_load_truth_set(tmp_path: Path) -> None:
    truth = vca.load_truth(_make_truth_csv(tmp_path))
    assert len(truth) == 5
    assert truth["NP_001191513.1"]["expected_classification"] == "chemoreceptor-candidate"
    assert truth["Q16950"]["expected_family"] == "aminergic"
    assert truth["Q16950"]["expected_subfamily"] == "5HT"


# ---- per-entry verdict --------------------------------------------------

def test_verdict_chemoreceptor_correctly_unmarked() -> None:
    """Chemoreceptor truth + classifier kept it as chemoreceptor-candidate
    -> correct (true negative for non-chemo classification)."""
    truth = {"expected_classification": "chemoreceptor-candidate",
             "expected_family": "", "expected_subfamily": ""}
    obs = {"classification": "chemoreceptor-candidate",
           "classification_family": "", "classification_subfamily": ""}
    v = vca.verdict_for_entry(truth, obs)
    assert v["correct"] is True
    assert v["error_type"] == "TN"  # true negative (correctly NOT marked non-chemo)


def test_verdict_chemoreceptor_falsely_marked_as_non_chemo() -> None:
    """Chemoreceptor truth but classifier marked it non-chemo -> FALSE POSITIVE
    (classifier polluted the would-be wet-lab shortlist)."""
    truth = {"expected_classification": "chemoreceptor-candidate",
             "expected_family": "", "expected_subfamily": ""}
    obs = {"classification": "non-chemoreceptor",
           "classification_family": "aminergic",
           "classification_subfamily": "5HT"}
    v = vca.verdict_for_entry(truth, obs)
    assert v["correct"] is False
    assert v["error_type"] == "FP"


def test_verdict_non_chemo_correctly_marked() -> None:
    """Non-chemoreceptor truth + classifier marked it non-chemo with the
    right family -> correct (true positive)."""
    truth = {"expected_classification": "non-chemoreceptor",
             "expected_family": "aminergic",
             "expected_subfamily": "5HT"}
    obs = {"classification": "non-chemoreceptor",
           "classification_family": "aminergic",
           "classification_subfamily": "5HT"}
    v = vca.verdict_for_entry(truth, obs)
    assert v["correct"] is True
    assert v["error_type"] == "TP"
    assert v["family_correct"] is True
    assert v["subfamily_correct"] is True


def test_verdict_non_chemo_correct_family_wrong_subfamily() -> None:
    """Family right but subfamily wrong: counts as TP at the coarse level
    (classifier called non-chemo correctly) but flags the subfamily error."""
    truth = {"expected_classification": "non-chemoreceptor",
             "expected_family": "aminergic",
             "expected_subfamily": "5HT"}
    obs = {"classification": "non-chemoreceptor",
           "classification_family": "aminergic",
           "classification_subfamily": "dopamine"}
    v = vca.verdict_for_entry(truth, obs)
    assert v["correct"] is True   # family-level call was right
    assert v["error_type"] == "TP"
    assert v["family_correct"] is True
    assert v["subfamily_correct"] is False


def test_verdict_non_chemo_missed() -> None:
    """Non-chemoreceptor truth but classifier left it as candidate ->
    FALSE NEGATIVE (classifier missed a known non-chemo)."""
    truth = {"expected_classification": "non-chemoreceptor",
             "expected_family": "aminergic", "expected_subfamily": "5HT"}
    obs = {"classification": "chemoreceptor-candidate",
           "classification_family": "", "classification_subfamily": ""}
    v = vca.verdict_for_entry(truth, obs)
    assert v["correct"] is False
    assert v["error_type"] == "FN"


def test_verdict_likely_non_chemo_treated_as_correct_for_non_chemo_truth() -> None:
    """'likely-non-chemoreceptor' (medium confidence, 2-of-3 evidence)
    counts as correct for non-chemo truth — it's still 'flagged'."""
    truth = {"expected_classification": "non-chemoreceptor",
             "expected_family": "aminergic", "expected_subfamily": "5HT"}
    obs = {"classification": "likely-non-chemoreceptor",
           "classification_family": "aminergic",
           "classification_subfamily": "5HT"}
    v = vca.verdict_for_entry(truth, obs)
    assert v["correct"] is True
    assert v["error_type"] == "TP"


# ---- end-to-end summary -------------------------------------------------

def test_compute_metrics(tmp_path: Path) -> None:
    """End-to-end: truth + ranked CSV -> metrics dict."""
    truth_csv = _make_truth_csv(tmp_path)
    ranked = _make_ranked_csv(tmp_path, [
        # 2 chemoreceptors correctly unmarked
        {"id": "NP_001191513.1", "classification": "chemoreceptor-candidate",
         "classification_family": "", "classification_subfamily": ""},
        {"id": "NP_001191494.1", "classification": "chemoreceptor-candidate",
         "classification_family": "", "classification_subfamily": ""},
        # 1 chemoreceptor falsely marked non-chemo (FP)
        {"id": "NP_001191495.1", "classification": "non-chemoreceptor",
         "classification_family": "peptide", "classification_subfamily": ""},
        # 2 non-chemo correctly marked
        {"id": "Q16950", "classification": "non-chemoreceptor",
         "classification_family": "aminergic", "classification_subfamily": "5HT"},
        {"id": "Q16951", "classification": "likely-non-chemoreceptor",
         "classification_family": "aminergic", "classification_subfamily": "5HT"},
    ])
    out_md = tmp_path / "report.md"
    metrics = vca.run_validation(ranked, truth_csv, str(out_md))
    # 5 truth entries: 2 TN + 1 FP + 2 TP = 4/5 correct
    assert metrics["n_total"] == 5
    assert metrics["n_correct"] == 4
    assert metrics["accuracy"] == 4 / 5
    # FP rate (chemoreceptors falsely marked) = 1/3
    assert metrics["chemoreceptor_correctness"] == 2 / 3
    # Non-chemo recall = 2/2
    assert metrics["non_chemoreceptor_recall"] == 1.0
    # Markdown report should exist
    assert out_md.exists()
    text = out_md.read_text()
    assert "Aplysia" in text or "validation" in text.lower()


def test_truth_entry_missing_from_ranked(tmp_path: Path) -> None:
    """Truth entry not in ranked CSV (pipeline didn't process it) -> warn,
    don't fail. Counted as 'not_classified'."""
    truth_csv = _make_truth_csv(tmp_path)
    ranked = _make_ranked_csv(tmp_path, [
        {"id": "NP_001191513.1", "classification": "chemoreceptor-candidate",
         "classification_family": "", "classification_subfamily": ""},
    ])
    out_md = tmp_path / "report.md"
    metrics = vca.run_validation(ranked, truth_csv, str(out_md))
    assert metrics["n_not_classified"] == 4  # 4 truth entries missing
