"""Tests for per-candidate embedding diagnostics + the scoring-confound report.

The confound report is the scientifically load-bearing part (bead cw3.11): it
must be able to DETECT when a model's novelty score is really tracking a nuisance
variable (sequence identity to nearest reference, or length) rather than biology.
So the decisive tests synthesize novelty AS a function of a confound and assert
the report recovers the (anti)correlation — a detector that can't detect a planted
confound is worthless.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from embedding_candidate_diagnostics import candidate_diagnostics, confound_report


def _reference():
    """Three characterized class-A families as tight, well-separated clusters,
    plus an out-of-class (class-B) and an orphan family that MUST be excluded
    from the novelty prototypes by the class-A restriction."""
    rng = np.random.RandomState(0)
    emb, labels = {}, {}
    plan = [("peptide", "A", np.eye(8)[0] * 6),
            ("aminergic", "A", np.eye(8)[1] * 6),
            ("opsin", "A", np.eye(8)[2] * 6),
            ("class-B-secretin", "B", np.eye(8)[3] * 6),
            ("orphan", "A", np.eye(8)[4] * 6)]
    i = 0
    for fam, cls, center in plan:
        for _ in range(8):
            k = f"ANCHOR_{cls}_1_P{i:03d}"
            emb[k] = rng.randn(8) * 0.15 + center
            labels[k] = fam
            i += 1
    return emb, labels


def test_candidate_near_a_family_scores_low_novelty_and_names_it():
    ref, labels = _reference()
    cand = {"on_peptide": np.eye(8)[0] * 6, "far": np.eye(8)[7] * 30}
    df = candidate_diagnostics(cand, ref, labels,
                               seq_len={"on_peptide": 400, "far": 400},
                               identity_to_nearest={"on_peptide": 0.9, "far": 0.2})
    r = df.set_index("candidate_id")
    assert r.loc["on_peptide", "nearest_family"] == "peptide"
    assert r.loc["far", "novelty"] > r.loc["on_peptide", "novelty"]


def test_nearest_family_never_out_of_class_or_orphan():
    ref, labels = _reference()
    # candidates sitting exactly on the class-B and orphan centroids
    cand = {"on_B": np.eye(8)[3] * 6, "on_orphan": np.eye(8)[4] * 6}
    df = candidate_diagnostics(cand, ref, labels,
                               seq_len={"on_B": 400, "on_orphan": 400},
                               identity_to_nearest={"on_B": 0.5, "on_orphan": 0.5})
    fams = set(df["nearest_family"])
    assert "class-B-secretin" not in fams and "orphan" not in fams
    assert fams <= {"peptide", "aminergic", "opsin"}


def test_margin_is_larger_for_an_unambiguous_candidate():
    ref, labels = _reference()
    cand = {"clear": np.eye(8)[0] * 6,                     # squarely peptide
            "between": (np.eye(8)[0] + np.eye(8)[1]) * 3}  # midway peptide/aminergic
    df = candidate_diagnostics(cand, ref, labels,
                               seq_len={"clear": 400, "between": 400},
                               identity_to_nearest={"clear": 0.9, "between": 0.5})
    r = df.set_index("candidate_id")
    assert r.loc["clear", "margin"] > r.loc["between", "margin"]


def test_confound_report_detects_a_planted_identity_confound():
    # novelty synthesized as a strict decreasing function of identity -> the
    # report MUST recover a strong negative Spearman (the ESM-C failure mode).
    df = pd.DataFrame({
        "candidate_id": [f"c{i}" for i in range(20)],
        "novelty": np.arange(20, 0, -1, dtype=float),
        "identity_to_nearest": np.linspace(0.1, 0.9, 20),
        "seq_len": np.random.RandomState(1).permutation(np.arange(300, 500, 10)),
    })
    rep = confound_report(df)
    assert rep["identity_to_nearest"]["spearman"] < -0.95
    # length was shuffled independently -> no strong correlation
    assert abs(rep["seq_len"]["spearman"]) < 0.6


def test_confound_report_flags_a_clean_scorer_as_uncorrelated():
    rng = np.random.RandomState(2)
    df = pd.DataFrame({
        "candidate_id": [f"c{i}" for i in range(50)],
        "novelty": rng.randn(50),
        "identity_to_nearest": rng.rand(50),
        "seq_len": rng.randint(300, 600, 50),
    })
    rep = confound_report(df)
    assert abs(rep["identity_to_nearest"]["spearman"]) < 0.4
    assert abs(rep["seq_len"]["spearman"]) < 0.4


def test_diagnostics_schema():
    ref, labels = _reference()
    cand = {"x": np.eye(8)[0] * 6}
    df = candidate_diagnostics(cand, ref, labels, seq_len={"x": 400},
                               identity_to_nearest={"x": 0.9})
    assert list(df.columns) == ["candidate_id", "novelty", "nearest_family",
                                "margin", "seq_len", "identity_to_nearest"]
