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
import pytest

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


# --- f2e: classification-based rank suppression -----------------------------
# 06c classification is a post-hoc annotation, not a hard filter. Without this,
# a non-chemoreceptor with strong non-phylo signals could surface in the
# shortlist (made more likely by o98 dropping the crude phylo penalty). Suppress
# (don't drop) 06c-classified non-chemoreceptors in rank_score, tunably.

def _consensus_tsv(tmp_path, rows):
    """rows: list of (candidate_id, classification). Other cols left blank."""
    lines = ["candidate_id\tclassification\tclassification_confidence\t"
             "classification_family\tclassification_subfamily\tclassification_evidence"]
    for cid, cls in rows:
        lines.append(f"{cid}\t{cls}\t\t\t\t")
    p = tmp_path / "classifications.tsv"
    p.write_text("\n".join(lines) + "\n")
    return p


def _run(tmp_path, ranked_text, consensus_rows, **kw):
    ranked = tmp_path / "ranked.csv"
    ranked.write_text(ranked_text)
    consensus = _consensus_tsv(tmp_path, consensus_rows)
    out = tmp_path / "out.csv"
    acc.add_classification_columns(str(ranked), str(consensus), str(out), **kw)
    return pd.read_csv(out)


def test_nonchemoreceptor_rank_score_suppressed(tmp_path):
    df = _run(tmp_path, "id,rank_score\nbste_g1,0.9\n",
              [("bste_g1", "non-chemoreceptor")])
    row = df[df["id"] == "bste_g1"].iloc[0]
    assert float(row["rank_score"]) == pytest.approx(0.9 * 0.1)
    assert float(row["rank_score_prefilter"]) == pytest.approx(0.9)
    assert float(row["classification_rank_factor"]) == pytest.approx(0.1)


def test_likely_nonchemoreceptor_moderately_suppressed(tmp_path):
    df = _run(tmp_path, "id,rank_score\nbste_g1,0.8\n",
              [("bste_g1", "likely-non-chemoreceptor")])
    row = df[df["id"] == "bste_g1"].iloc[0]
    assert float(row["rank_score"]) == pytest.approx(0.8 * 0.5)


def test_chemoreceptor_candidate_unchanged(tmp_path):
    # Default (not in consensus) and explicit chemoreceptor-candidate are x1.0.
    df = _run(tmp_path, "id,rank_score\nbste_g1,0.7\nbste_g2,0.6\n",
              [("bste_g2", "chemoreceptor-candidate")])
    assert float(df[df["id"] == "bste_g1"].iloc[0]["rank_score"]) == pytest.approx(0.7)
    assert float(df[df["id"] == "bste_g2"].iloc[0]["rank_score"]) == pytest.approx(0.6)


def test_suppression_reorders_shortlist(tmp_path):
    # A high-composite non-chemoreceptor must fall below a lower-composite
    # chemoreceptor-candidate after suppression.
    df = _run(tmp_path, "id,rank_score\nbste_hi,0.9\nbste_lo,0.5\n",
              [("bste_hi", "non-chemoreceptor")])
    assert df.iloc[0]["id"] == "bste_lo"   # 0.5 now outranks 0.9*0.1=0.09


def test_factors_configurable_can_disable(tmp_path):
    df = _run(tmp_path, "id,rank_score\nbste_g1,0.9\n",
              [("bste_g1", "non-chemoreceptor")],
              nonchemo_rank_factor=1.0, likely_nonchemo_rank_factor=1.0)
    assert float(df[df["id"] == "bste_g1"].iloc[0]["rank_score"]) == pytest.approx(0.9)


def test_no_rank_score_column_is_noop(tmp_path):
    # A ranked CSV without rank_score must not error (annotation still applied).
    df = _run(tmp_path, "id,gene_name\nbste_g1,Bste_g1\n",
              [("bste_g1", "non-chemoreceptor")])
    assert df[df["id"] == "bste_g1"].iloc[0]["classification"] == "non-chemoreceptor"


# --- cl2: needs_manual_review flag ------------------------------------------
# Flag candidates that defaulted to 'chemoreceptor-candidate' only because BOTH
# annotation-based sources (HMM + OG) went dark — the divergent/unclassifiable
# case where automated classification gives no signal, so a human should eyeball
# them before a real divergent chemoreceptor is silently deprioritized. Placement
# alone does NOT rescue (least reliable source).

def _consensus_evidence_tsv(tmp_path, rows):
    """rows: list of (candidate_id, classification, evidence)."""
    header = ("candidate_id\tclassification\tclassification_confidence\t"
              "classification_family\tclassification_subfamily\tclassification_evidence")
    lines = [header]
    for cid, cls, ev in rows:
        lines.append(f"{cid}\t{cls}\t\t\t\t{ev}")
    p = tmp_path / "classifications.tsv"
    p.write_text("\n".join(lines) + "\n")
    return p


def _run_ev(tmp_path, ranked_text, consensus_rows):
    ranked = tmp_path / "ranked.csv"
    ranked.write_text(ranked_text)
    consensus = _consensus_evidence_tsv(tmp_path, consensus_rows)
    out = tmp_path / "out.csv"
    acc.add_classification_columns(str(ranked), str(consensus), str(out))
    return pd.read_csv(out, keep_default_na=False, dtype=str)


def _flag(df, cid):
    return df[df["id"] == cid].iloc[0]["needs_manual_review"]


def test_flagged_when_both_hmm_and_og_dark(tmp_path):
    df = _run_ev(tmp_path, "id,rank_score\nc1,0.5\n",
                 [("c1", "chemoreceptor-candidate",
                   "hmm:unclassified-hmm;og:unclassified-og;placement:unclassified-placement")])
    assert _flag(df, "c1") == "yes"


def test_not_flagged_when_hmm_makes_a_call(tmp_path):
    df = _run_ev(tmp_path, "id,rank_score\nc1,0.5\n",
                 [("c1", "chemoreceptor-candidate",
                   "hmm:aminergic;og:unclassified-og;placement:NA")])
    assert _flag(df, "c1") == ""


def test_not_flagged_when_og_makes_a_call(tmp_path):
    df = _run_ev(tmp_path, "id,rank_score\nc1,0.5\n",
                 [("c1", "chemoreceptor-candidate",
                   "hmm:unclassified-hmm;og:peptide;placement:NA")])
    assert _flag(df, "c1") == ""


def test_placement_alone_does_not_rescue(tmp_path):
    # HMM+OG dark, placement classifies -> still 'chemoreceptor-candidate' and
    # still flagged (placement alone is the least-reliable source).
    df = _run_ev(tmp_path, "id,rank_score\nc1,0.5\n",
                 [("c1", "chemoreceptor-candidate",
                   "hmm:unclassified-hmm;og:unclassified-og;placement:aminergic")])
    assert _flag(df, "c1") == "yes"


def test_not_flagged_when_classified_non_chemoreceptor(tmp_path):
    df = _run_ev(tmp_path, "id,rank_score\nc1,0.5\n",
                 [("c1", "non-chemoreceptor", "hmm:aminergic;og:aminergic;placement:aminergic")])
    assert _flag(df, "c1") == ""


def test_flagged_when_no_consensus_entry(tmp_path):
    # Candidate absent from the consensus TSV -> default chemoreceptor-candidate,
    # empty evidence -> both sources dark -> flag for review.
    df = _run_ev(tmp_path, "id,rank_score\nc1,0.5\n", [])
    assert _flag(df, "c1") == "yes"
