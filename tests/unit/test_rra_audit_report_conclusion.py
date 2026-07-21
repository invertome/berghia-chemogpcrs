"""The audit report's conclusion must follow from its own numbers (bead pu86).

``write_report`` appended this sentence unconditionally, four lines below the
computed effect block::

    An order Spearman of 1.0 with 0 displacement means the measured correlation
    does NOT change the ranking, and the council's double-counting objection,
    while structurally real, is inert here.

On the real 439-candidate cohort the numbers printed directly above it are
Spearman 0.6122, top-20 Jaccard 0.333, all 20 top-k ids changed, max
displacement 380. The report asserted a conclusion its own numbers refute.

The sentence predates this work but was dead code -- the script had no callers
until stage 07 was wired to it, at which point it began shipping every run.

A conclusion that cannot be false is not a conclusion, so this pins BOTH
directions: it must say "inert" only when the measurement actually supports it,
and must say the opposite when it does not.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from audit_rra_correlation_sensitivity import write_report


def _result(spearman: float, displacement: int, jaccard: float, changed: int) -> dict:
    return {
        "signals": ["a", "b"],
        "threshold": 0.7,
        "flagged_pairs": [],
        "correlation": pd.DataFrame({"a": [1.0]}, index=["a"]),
        "saturation": {"n_candidates": 439, "n_distinct_scores": 219,
                       "n_tied": 243, "fraction_tied": 0.554,
                       "largest_tie_block": 97, "n_saturated": 0,
                       "fraction_saturated": 0.0},
        "fusion_impact": {"order_spearman": spearman, "top_k": 20,
                          "top_k_jaccard": jaccard, "n_top_k_changed": changed,
                          "max_displacement": displacement},
    }


def _md(tmp_path: Path, **kw) -> str:
    write_report(_result(**kw), str(tmp_path / "rra"))
    return (tmp_path / "rra_sensitivity.md").read_text()


def test_does_not_claim_inert_when_the_ranking_actually_moves(tmp_path: Path):
    """The real-cohort case that exposed this."""
    md = _md(tmp_path, spearman=0.6122, displacement=380, jaccard=0.333, changed=20)
    assert "is inert here" not in md, (
        "the report claims the correlation is inert while its own numbers show "
        "Spearman 0.6122 and 380 positions of displacement")
    assert "0.6122" in md


def test_says_the_ranking_moved_when_it_moved(tmp_path: Path):
    """Must-accept the opposite direction: silence is not enough, the report
    has to state the finding."""
    md = _md(tmp_path, spearman=0.6122, displacement=380, jaccard=0.333, changed=20)
    assert "does change the ranking" in md.lower() or "changes the ranking" in md.lower(), (
        "the report drops the conclusion entirely instead of reporting that the "
        "correlation moved the order")


def test_claims_inert_only_when_the_measurement_supports_it(tmp_path: Path):
    """Must-reject control: a guard tuned only one way is useless."""
    md = _md(tmp_path, spearman=1.0, displacement=0, jaccard=1.0, changed=0)
    assert "inert" in md, (
        "with Spearman 1.0 and zero displacement the correlation genuinely is "
        "inert, and the report should say so")


def test_the_tie_fraction_reaches_the_summary_not_only_the_json(tmp_path: Path):
    """fraction_saturated is ~0 by construction under the exact null, so a
    summary quoting only that reads healthy while 55% of the order is tied."""
    md = _md(tmp_path, spearman=1.0, displacement=0, jaccard=1.0, changed=0)
    assert "55.4%" in md or "0.554" in md, "the tie fraction is not in the report body"


def test_json_and_markdown_do_not_disagree(tmp_path: Path):
    write_report(_result(0.6122, 380, 0.333, 20), str(tmp_path / "rra"))
    payload = json.loads((tmp_path / "rra_sensitivity.json").read_text())
    md = (tmp_path / "rra_sensitivity.md").read_text()
    assert f"{payload['fusion_impact']['order_spearman']:.4f}" in md
