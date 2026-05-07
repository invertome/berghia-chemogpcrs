"""Tests for scripts/add_classification_columns.py.

Phase 5 Task 5.1 — joins the 06c consensus TSV into the ranked
candidates CSV, adding columns:
    classification, classification_confidence, classification_family,
    classification_subfamily, classification_evidence

Mirrors the pattern of scripts/add_hcr_columns.py and
scripts/add_og_coverage_columns.py — reads ranked CSV, joins on
candidate_id, writes back. Default 'chemoreceptor-candidate' for any
candidate not in the consensus TSV.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import add_classification_columns as acc


def test_add_columns_joins_classifications(tmp_path: Path) -> None:
    """Standard end-to-end: ranked CSV gets the 5 classification columns."""
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text(
        "id,gene_name,rank_score\n"
        "bste_g1,Bste_g1,0.95\n"
        "bste_g2,Bste_g2,0.92\n"
        "bste_g3,Bste_g3,0.85\n"
    )
    consensus_tsv = tmp_path / "classifications.tsv"
    consensus_tsv.write_text(
        "candidate_id\tclassification\tclassification_confidence\t"
        "classification_family\tclassification_subfamily\t"
        "classification_evidence\n"
        "bste_g1\tnon-chemoreceptor\thigh\taminergic\t5HT\t"
        "hmm:aminergic;og:aminergic;placement:aminergic\n"
        "bste_g2\tlikely-non-chemoreceptor\tmedium\tpeptide\ttachykinin\t"
        "hmm:peptide;og:peptide;placement:opsin\n"
        # bste_g3 not in classifications -> default chemoreceptor-candidate
    )
    out_csv = tmp_path / "ranked_with_class.csv"
    acc.add_classification_columns(
        str(ranked_csv), str(consensus_tsv), str(out_csv))

    df = pd.read_csv(out_csv, keep_default_na=False, dtype=str)
    cols = list(df.columns)
    for c in ["classification", "classification_confidence",
              "classification_family", "classification_subfamily",
              "classification_evidence"]:
        assert c in cols
    row1 = df[df["id"] == "bste_g1"].iloc[0]
    assert row1["classification"] == "non-chemoreceptor"
    assert row1["classification_family"] == "aminergic"
    row3 = df[df["id"] == "bste_g3"].iloc[0]
    # Default for unclassified candidates
    assert row3["classification"] == "chemoreceptor-candidate"
    assert row3["classification_confidence"] == "NA"


def test_missing_consensus_tsv_falls_back_gracefully(tmp_path: Path) -> None:
    """If the consensus TSV doesn't exist (e.g. 06c hasn't run), all
    candidates default to chemoreceptor-candidate."""
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text("id,gene_name\nbste_g1,X\n")
    out_csv = tmp_path / "out.csv"
    acc.add_classification_columns(
        str(ranked_csv), str(tmp_path / "missing.tsv"), str(out_csv))
    df = pd.read_csv(out_csv, keep_default_na=False, dtype=str)
    assert df.iloc[0]["classification"] == "chemoreceptor-candidate"


def test_preserves_original_columns_and_order(tmp_path: Path) -> None:
    """The 5 new columns are appended; original columns stay in order."""
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text("id,col_a,col_b,rank_score\nbste_g1,A,B,0.5\n")
    consensus_tsv = tmp_path / "c.tsv"
    consensus_tsv.write_text(
        "candidate_id\tclassification\tclassification_confidence\t"
        "classification_family\tclassification_subfamily\t"
        "classification_evidence\n"
        "bste_g1\tnon-chemoreceptor\thigh\topsin\t\thmm:opsin;og:opsin\n"
    )
    out_csv = tmp_path / "out.csv"
    acc.add_classification_columns(
        str(ranked_csv), str(consensus_tsv), str(out_csv))
    df = pd.read_csv(out_csv, keep_default_na=False, dtype=str)
    cols = list(df.columns)
    # Original columns first, in order
    assert cols.index("id") < cols.index("col_a") < cols.index("col_b")
    assert cols.index("col_b") < cols.index("rank_score")
    # New columns appended after
    assert cols.index("rank_score") < cols.index("classification")


def test_subfamily_empty_when_no_subfamily_call(tmp_path: Path) -> None:
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text("id\nbste_g1\n")
    consensus_tsv = tmp_path / "c.tsv"
    consensus_tsv.write_text(
        "candidate_id\tclassification\tclassification_confidence\t"
        "classification_family\tclassification_subfamily\t"
        "classification_evidence\n"
        "bste_g1\tnon-chemoreceptor\thigh\topsin\t\thmm:opsin\n"
    )
    out_csv = tmp_path / "out.csv"
    acc.add_classification_columns(
        str(ranked_csv), str(consensus_tsv), str(out_csv))
    df = pd.read_csv(out_csv, keep_default_na=False, dtype=str)
    assert df.iloc[0]["classification_family"] == "opsin"
    assert df.iloc[0]["classification_subfamily"] == ""
