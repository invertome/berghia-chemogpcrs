"""Pin the validated-vs-unvalidated threshold boundary in the HMM classifier.

Audit finding #12. `assign_classification` used
`thresholds.get(family, DEFAULT_THRESHOLD)` with DEFAULT_THRESHOLD = 1e-10,
while the real leave-one-out thresholds in
results/classification/loo/loo_metrics.tsv span 2.4e-105 .. 2.5e-198.

Only the 8 custom families have a LOO row. Every Pfam-fallback label
therefore got a cutoff 95-190 orders of magnitude looser than any
validated family. Because the same Pfam library also annotates the
reference set that the orthogroup vote reads, the HMM and OG votes agree
BY CONSTRUCTION -- yielding 'likely-non-chemoreceptor' and a 0.5x score
multiplier on the strength of a path that was never benchmarked. In a
pipeline whose entire purpose is finding chemoreceptors, silently
demoting a real one is the expensive error.

Decision pinned here: an unvalidated call is still EMITTED (the
annotation is useful) but it is TAGGED, and a tagged call may not carry
the consensus into a demotion tier on its own. Fail-open on annotation,
fail-closed on demotion.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import classify_consensus as cc
import classify_via_hmm as cvh


# Real LOO thresholds, read from results/classification/loo/loo_metrics.tsv
# on 2026-07-20. The loosest validated threshold is aminergic at 2.4e-105.
LOOSEST_VALIDATED_THRESHOLD = 2.4e-105


# ---------------------------------------------------------------------------
# 1. The gap the audit describes is real and pinned
# ---------------------------------------------------------------------------

def test_unvalidated_fallback_is_far_looser_than_any_validated_threshold():
    """Documents the gap rather than hiding it: this is exactly why an
    unvalidated call may not be treated as equivalent evidence."""
    assert cvh.UNVALIDATED_FALLBACK_THRESHOLD > LOOSEST_VALIDATED_THRESHOLD * 1e50


# ---------------------------------------------------------------------------
# 2. Threshold provenance is reported per call
# ---------------------------------------------------------------------------

def test_validated_family_is_tagged_loo():
    hits = [("aminergic_5HT", 1e-120, 400.0)]
    result = cvh.assign_classification(hits, {"aminergic": 2.4e-105})
    assert result["family"] == "aminergic"
    assert result["threshold_source"] == "loo"


def test_family_without_loo_row_is_tagged_unvalidated():
    """A Pfam-fallback label has no LOO row — the call is emitted but tagged."""
    hits = [("class-A-7tm", 1e-30, 120.0)]
    result = cvh.assign_classification(hits, {"aminergic": 2.4e-105})
    assert result["family"] == "class-A-7tm"
    assert result["threshold_source"] == "unvalidated-default"


def test_unclassified_call_reports_no_threshold_source():
    hits = [("aminergic_5HT", 1e-3, 20.0)]
    result = cvh.assign_classification(hits, {"aminergic": 2.4e-105})
    assert result["family"] == "unclassified-hmm"
    assert result["threshold_source"] == ""


def test_validated_hit_is_preferred_over_an_unvalidated_better_evalue():
    """A validated family that passes must win over an unvalidated label,
    even when the unvalidated hit has a better E-value — otherwise the
    loose default lets the unbenchmarked path outrank the benchmarked one."""
    hits = [
        ("class-A-7tm", 1e-200, 700.0),     # unvalidated, best E-value
        ("aminergic_5HT", 1e-120, 400.0),   # validated, passes its threshold
    ]
    result = cvh.assign_classification(hits, {"aminergic": 2.4e-105})
    assert result["family"] == "aminergic"
    assert result["threshold_source"] == "loo"


def test_unvalidated_call_survives_when_no_validated_family_passes():
    """Fail-open on annotation: we still report what was found."""
    hits = [("class-A-7tm", 1e-30, 120.0), ("aminergic_5HT", 1e-3, 20.0)]
    result = cvh.assign_classification(hits, {"aminergic": 2.4e-105})
    assert result["family"] == "class-A-7tm"
    assert result["threshold_source"] == "unvalidated-default"


def test_empty_threshold_table_tags_everything_unvalidated():
    """With no LOO metrics at all, nothing is validated — the whole run
    must be visibly unvalidated rather than silently authoritative."""
    hits = [("aminergic_5HT", 1e-120, 400.0)]
    result = cvh.assign_classification(hits, {})
    assert result["threshold_source"] == "unvalidated-default"


# ---------------------------------------------------------------------------
# 3. Consensus must not demote on unvalidated evidence alone
# ---------------------------------------------------------------------------

def _write_tsv(path: Path, rows: list[dict]) -> str:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return str(path)


def test_unvalidated_hmm_call_cannot_produce_a_medium_demotion(tmp_path):
    """The core bug: HMM + OG 'agree' because the same unbenchmarked Pfam
    library annotated both, yielding likely-non-chemoreceptor + 0.5x."""
    h = _write_tsv(tmp_path / "h.tsv", [{
        "candidate_id": "g1", "hmm_family": "aminergic", "hmm_subfamily": "",
        "evalue": "1e-30", "score": "120", "hmm_target": "class-A-7tm",
        "threshold_source": "unvalidated-default",
    }])
    o = _write_tsv(tmp_path / "o.tsv", [{
        "candidate_id": "g1", "og_vote_family": "aminergic",
        "og_vote_subfamily": "",
    }])
    out = tmp_path / "out.tsv"
    cc.run_consensus(h, o, None, str(out))
    row = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert row["classification"] == "chemoreceptor-candidate"
    assert row["classification_confidence"] == "NA"


def test_validated_hmm_call_still_produces_a_medium_demotion(tmp_path):
    """Guard against over-correction: the validated path must be unchanged."""
    h = _write_tsv(tmp_path / "h.tsv", [{
        "candidate_id": "g1", "hmm_family": "aminergic", "hmm_subfamily": "",
        "evalue": "1e-120", "score": "400", "hmm_target": "aminergic_5HT",
        "threshold_source": "loo",
    }])
    o = _write_tsv(tmp_path / "o.tsv", [{
        "candidate_id": "g1", "og_vote_family": "aminergic",
        "og_vote_subfamily": "",
    }])
    out = tmp_path / "out.tsv"
    cc.run_consensus(h, o, None, str(out))
    row = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert row["classification"] == "likely-non-chemoreceptor"
    assert row["classification_confidence"] == "medium"


def test_unvalidated_hmm_call_cannot_produce_a_high_demotion(tmp_path):
    h = _write_tsv(tmp_path / "h.tsv", [{
        "candidate_id": "g1", "hmm_family": "opsin", "hmm_subfamily": "",
        "evalue": "1e-30", "score": "120", "hmm_target": "class-A-7tm",
        "threshold_source": "unvalidated-default",
    }])
    o = _write_tsv(tmp_path / "o.tsv", [{
        "candidate_id": "g1", "og_vote_family": "opsin", "og_vote_subfamily": "",
    }])
    p = _write_tsv(tmp_path / "p.tsv", [{
        "candidate_id": "g1", "placement_family": "opsin",
        "placement_subfamily": "",
    }])
    out = tmp_path / "out.tsv"
    cc.run_consensus(h, o, p, str(out))
    row = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert row["classification"] == "chemoreceptor-candidate"


def test_unvalidated_evidence_is_still_visible_in_the_evidence_string(tmp_path):
    """The annotation must not vanish — a reviewer has to be able to see
    that an unvalidated HMM call happened and was discounted."""
    h = _write_tsv(tmp_path / "h.tsv", [{
        "candidate_id": "g1", "hmm_family": "aminergic", "hmm_subfamily": "",
        "evalue": "1e-30", "score": "120", "hmm_target": "class-A-7tm",
        "threshold_source": "unvalidated-default",
    }])
    o = _write_tsv(tmp_path / "o.tsv", [{
        "candidate_id": "g1", "og_vote_family": "aminergic",
        "og_vote_subfamily": "",
    }])
    out = tmp_path / "out.tsv"
    cc.run_consensus(h, o, None, str(out))
    row = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert "aminergic" in row["classification_evidence"]
    assert "unvalidated" in row["classification_evidence"]


def test_missing_threshold_source_column_is_treated_as_validated(tmp_path):
    """Backwards compatibility: a TSV written before this column existed
    must keep its old meaning rather than silently losing all demotions."""
    h = _write_tsv(tmp_path / "h.tsv", [{
        "candidate_id": "g1", "hmm_family": "aminergic", "hmm_subfamily": "",
        "evalue": "1e-120", "score": "400", "hmm_target": "aminergic_5HT",
    }])
    o = _write_tsv(tmp_path / "o.tsv", [{
        "candidate_id": "g1", "og_vote_family": "aminergic",
        "og_vote_subfamily": "",
    }])
    out = tmp_path / "out.tsv"
    cc.run_consensus(h, o, None, str(out))
    row = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert row["classification"] == "likely-non-chemoreceptor"
