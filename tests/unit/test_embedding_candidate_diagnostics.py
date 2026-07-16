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


def test_confound_report_role_labels_and_low_identity_block():
    # length is a pure ARTIFACT axis (want ~0); identity is a REDUNDANCY axis
    # (moderate correlation is expected/healthy — the check is added-value, and it
    # carries a low-identity-subset sub-analysis).
    rng = np.random.RandomState(0)
    df = pd.DataFrame({
        "candidate_id": [f"c{i}" for i in range(30)],
        "novelty": rng.randn(30),
        "identity_to_nearest": np.linspace(0, 90, 30),
        "seq_len": rng.randint(300, 600, 30),
    })
    rep = confound_report(df)
    assert rep["seq_len"]["role"] == "artifact"
    assert rep["identity_to_nearest"]["role"] == "redundancy"
    assert "low_identity" in rep["identity_to_nearest"]
    assert "low_identity" not in rep["seq_len"]


def test_length_confound_flagged_as_artifact():
    n = 60
    seq_len = np.linspace(200, 800, n)
    df = pd.DataFrame({
        "candidate_id": [f"c{i}" for i in range(n)],
        "novelty": -seq_len + np.random.RandomState(3).normal(0, 5, n),  # tracks length
        "identity_to_nearest": np.random.RandomState(4).uniform(0, 90, n),
        "seq_len": seq_len,
    })
    assert confound_report(df)["seq_len"]["verdict"] == "artifact-confound"


def test_identity_redundant_when_novelty_is_pure_identity_proxy():
    n = 80
    ident = np.linspace(0, 90, n)
    df = pd.DataFrame({
        "candidate_id": [f"c{i}" for i in range(n)],
        "novelty": 90 - ident,                       # strictly monotone in identity
        "identity_to_nearest": ident,
        "seq_len": np.random.RandomState(5).randint(300, 600, n),
    })
    # near-perfect global correlation => the embedding just re-derives identity
    assert confound_report(df)["identity_to_nearest"]["verdict"] == "redundant-with-identity"


def test_identity_adds_value_when_decorrelated_in_low_regime():
    # moderate global correlation (divergent IS more novel — expected) but within
    # the low-identity tail novelty is INDEPENDENT of residual identity => the
    # embedding adds orthogonal signal exactly where alignment can't.
    n = 300
    rng = np.random.default_rng(0)
    ident = rng.uniform(0, 90, n)
    q = np.quantile(ident, 0.33)
    novelty = np.empty(n)
    hi = ident > q
    novelty[hi] = -0.5 * ident[hi] + rng.normal(0, 20, hi.sum())
    novelty[~hi] = rng.normal(50, 20, (~hi).sum())
    df = pd.DataFrame({
        "candidate_id": [f"c{i}" for i in range(n)],
        "novelty": novelty, "identity_to_nearest": ident,
        "seq_len": rng.integers(300, 600, n),
    })
    rep = confound_report(df)["identity_to_nearest"]
    assert rep["verdict"] == "adds-value"
    assert abs(rep["low_identity"]["spearman"]) < 0.4


def test_identity_low_regime_artifact_when_novelty_tracks_identity_among_orphans():
    # global correlation is only moderate (not flagged redundant), but WITHIN the
    # low-identity tail novelty still strongly tracks residual identity => it is an
    # identity proxy in the regime that matters (the Dijkhof failure).
    n = 300
    rng = np.random.default_rng(1)
    ident = rng.uniform(0, 90, n)
    q = np.quantile(ident, 0.33)
    novelty = np.empty(n)
    lo = ident <= q
    novelty[lo] = 90 - ident[lo] + rng.normal(0, 1, lo.sum())
    novelty[~lo] = rng.normal(45, 20, (~lo).sum())
    df = pd.DataFrame({
        "candidate_id": [f"c{i}" for i in range(n)],
        "novelty": novelty, "identity_to_nearest": ident,
        "seq_len": rng.integers(300, 600, n),
    })
    rep = confound_report(df)["identity_to_nearest"]
    assert rep["verdict"] == "low-identity-artifact"
    assert abs(rep["low_identity"]["spearman"]) > 0.4


def test_diagnostics_schema():
    ref, labels = _reference()
    cand = {"x": np.eye(8)[0] * 6}
    df = candidate_diagnostics(cand, ref, labels, seq_len={"x": 400},
                               identity_to_nearest={"x": 0.9})
    assert list(df.columns) == ["candidate_id", "novelty", "nearest_family",
                                "margin", "seq_len", "identity_to_nearest"]
