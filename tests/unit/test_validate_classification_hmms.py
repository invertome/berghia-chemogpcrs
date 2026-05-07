"""Tests for scripts/validate_classification_hmms.py.

Phase 2 Task 2.2 — leave-one-out (LOO) cross-validation of the per-family
HMM library built in Task 2.1.

LOO procedure: for each HMM, hold out one training sequence at a time,
rebuild the HMM from remaining members, hmmscan the held-out sequence
against the FULL HMM library (with the rebuilt HMM), and record whether
the top hit's family matches the held-out sequence's true family.

Tests cover the local (no-subprocess) parts:
  - parsing manifest TSV
  - parsing HMM alignment FASTA
  - confusion-matrix accumulation
  - per-family recall / precision calculation
  - threshold-selection logic
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

import validate_classification_hmms as vch


# ---- parse_manifest -----------------------------------------------------

def test_parse_manifest_basic(tmp_path: Path) -> None:
    """Manifest TSV from Task 2.1 has columns family/subfamily/n_train/
    hmm_path/alignment_path. Coarse rows have empty subfamily."""
    m = tmp_path / "manifest.tsv"
    m.write_text(
        "family\tsubfamily\tn_train\thmm_path\talignment_path\n"
        "aminergic\t\t165\tresults/aminergic.hmm\tresults/aminergic.aln\n"
        "aminergic\t5HT\t39\tresults/aminergic_5HT.hmm\tresults/aminergic_5HT.aln\n"
        "peptide\t\t189\tresults/peptide.hmm\tresults/peptide.aln\n"
    )
    rows = vch.parse_manifest(str(m))
    assert len(rows) == 3
    assert rows[0]["family"] == "aminergic"
    assert rows[0]["subfamily"] == ""
    assert rows[0]["n_train"] == 165
    assert rows[1]["subfamily"] == "5HT"


# ---- parse_alignment_fasta ----------------------------------------------

def test_parse_alignment_fasta(tmp_path: Path) -> None:
    """Alignment FASTA is the unaligned-style FASTA hmmbuild reads (sequences
    may have gaps in real alignments — preserved for structural reuse but
    stripped before LOO rebuild)."""
    fa = tmp_path / "aln.fasta"
    fa.write_text(
        ">P28335\nMVNL---RNAVHSF\n"
        ">P50406\nMSEKLREALV----LV\n"
    )
    seqs = vch.parse_alignment_fasta(str(fa))
    assert seqs["P28335"] == "MVNL---RNAVHSF"
    assert seqs["P50406"] == "MSEKLREALV----LV"


def test_strip_gaps() -> None:
    """LOO rebuild needs ungapped sequences (re-aligned)."""
    assert vch.strip_gaps("MVNL---RNAVHSF") == "MVNLRNAVHSF"
    assert vch.strip_gaps("M.AB-CD..") == "MABCD"


# ---- combined hit label normalization -----------------------------------

def test_label_for_hmm_coarse() -> None:
    """For coarse-only HMMs (no subfamily), the family label is used directly."""
    assert vch.label_for_hmm("aminergic") == ("aminergic", "")


def test_label_for_hmm_medium() -> None:
    """Medium-granularity HMMs are named family_subfamily."""
    assert vch.label_for_hmm("aminergic_5HT") == ("aminergic", "5HT")


def test_label_for_hmm_with_hyphen_in_family() -> None:
    """Families with hyphens (class-B-secretin, class-F-frizzled) must
    parse correctly even though `_` is the family/subfamily separator."""
    assert vch.label_for_hmm("class-B-secretin") == ("class-B-secretin", "")
    assert vch.label_for_hmm("class-F-frizzled") == ("class-F-frizzled", "")


def test_label_for_hmm_peptide_subfamily_with_hyphen() -> None:
    """peptide_NPY-NPF, peptide_vasopressin-oxytocin contain hyphens in the
    subfamily portion — split on FIRST underscore only."""
    assert vch.label_for_hmm("peptide_NPY-NPF") == ("peptide", "NPY-NPF")
    assert vch.label_for_hmm("peptide_vasopressin-oxytocin") == (
        "peptide", "vasopressin-oxytocin")


# ---- confusion matrix + metrics -----------------------------------------

def test_confusion_matrix_basic() -> None:
    """LOO results: list of (true_family, predicted_family) per held-out seq."""
    results = [
        ("aminergic", "aminergic"),
        ("aminergic", "aminergic"),
        ("aminergic", "peptide"),  # misclassified
        ("peptide", "peptide"),
        ("peptide", "peptide"),
        ("opsin", "opsin"),
        ("opsin", "unclassified"),  # missed
    ]
    cm = vch.confusion_matrix(results)
    assert cm[("aminergic", "aminergic")] == 2
    assert cm[("aminergic", "peptide")] == 1
    assert cm[("peptide", "peptide")] == 2
    assert cm[("opsin", "opsin")] == 1
    assert cm[("opsin", "unclassified")] == 1


def test_recall_per_family() -> None:
    """Recall = TP / (TP + FN) per family."""
    results = [
        ("aminergic", "aminergic"),  # TP
        ("aminergic", "aminergic"),  # TP
        ("aminergic", "peptide"),    # FN (missed)
        ("aminergic", "unclassified"),  # FN
        ("peptide", "peptide"),       # TP
    ]
    metrics = vch.per_family_metrics(results)
    assert metrics["aminergic"]["recall"] == 2 / 4  # 2 of 4 aminergic hits
    assert metrics["peptide"]["recall"] == 1.0


def test_precision_per_family() -> None:
    """Precision = TP / (TP + FP) per family."""
    results = [
        ("aminergic", "aminergic"),  # TP for aminergic
        ("aminergic", "aminergic"),  # TP for aminergic
        ("peptide", "aminergic"),    # FP for aminergic (real peptide, called aminergic)
        ("peptide", "peptide"),      # TP for peptide
        ("opsin", "peptide"),        # FP for peptide
    ]
    metrics = vch.per_family_metrics(results)
    # aminergic: 2 TP + 1 FP = precision 2/3
    assert metrics["aminergic"]["precision"] == 2 / 3
    # peptide: 1 TP + 1 FP = precision 1/2
    assert metrics["peptide"]["precision"] == 1 / 2


def test_metrics_handle_zero_division() -> None:
    """Family with no positive predictions = precision NaN (not crash)."""
    results = [
        ("aminergic", "aminergic"),
        ("aminergic", "aminergic"),
    ]
    metrics = vch.per_family_metrics(results)
    # peptide has no TP and no FP; precision undefined
    assert "peptide" not in metrics or metrics.get("peptide", {}).get(
        "precision") is None or metrics["peptide"].get("n_total", 0) == 0


# ---- threshold selection ------------------------------------------------

def test_select_threshold_picks_loosest_passing() -> None:
    """Given LOO scan results (true_family, predicted_family, evalue) and
    a target recall, pick the loosest E-value threshold that still keeps
    recall above target."""
    # aminergic: 5 held-out; sorted by evalue (ascending = best to worst)
    # E-vals where the top hit is aminergic:
    #   1e-50 → aminergic (TP)
    #   1e-30 → aminergic (TP)
    #   1e-20 → peptide (FN)  <- the threshold matters here
    #   1e-10 → aminergic (TP)
    #   1e-5  → aminergic (TP)
    # If target recall = 0.8 (4/5), threshold needs to keep 4 TPs;
    # easiest interpretation: pick the loosest E-val cutoff such that
    # AT LEAST 4 of the 5 held-out get a TP at evalue <= cutoff.
    rows = [
        ("aminergic", "aminergic", 1e-50),
        ("aminergic", "aminergic", 1e-30),
        ("aminergic", "peptide",   1e-20),
        ("aminergic", "aminergic", 1e-10),
        ("aminergic", "aminergic", 1e-5),
    ]
    # 4 TPs total (recall = 4/5 = 0.8), all at evalues <= 1e-5
    threshold = vch.select_evalue_threshold(
        rows, family="aminergic", target_recall=0.8)
    # Should be >= 1e-5 (the loosest TP) and <= some reasonable cap
    assert threshold >= 1e-5
    assert threshold <= 1.0  # not arbitrarily loose


def test_select_threshold_default_when_no_data() -> None:
    """If no held-out sequences for a family, return a conservative default."""
    threshold = vch.select_evalue_threshold([], family="missing", target_recall=0.9)
    assert threshold == 1e-10  # reasonable default


# ---- ID handling --------------------------------------------------------

def test_held_out_id_resolves_to_known_family() -> None:
    """Map seq_id -> (family, subfamily) for the held-out check."""
    manifest = [
        {"family": "aminergic", "subfamily": "5HT", "alignment_path": "/p/x.aln",
         "n_train": 39, "hmm_path": "/p/x.hmm"},
    ]
    seqs_by_alignment = {"/p/x.aln": {"P28335": "MVNL", "P50406": "MSEK"}}
    id_to_family = vch.build_id_to_family_map(manifest, seqs_by_alignment)
    assert id_to_family["P28335"] == ("aminergic", "5HT")
    assert id_to_family["P50406"] == ("aminergic", "5HT")
