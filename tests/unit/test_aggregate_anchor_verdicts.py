"""Unit tests for aggregate_anchor_verdicts.py — combine per-class C3 verdicts
into a single report-only summary (no config is mutated)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import aggregate_anchor_verdicts as agg  # noqa: E402


def _verdict(klass, verdict, **kw):
    d = {"class": klass, "verdict": verdict, "n_infiltrations": 0, "rf": 0.0,
         "support_drop": 0.0, "reasons": []}
    d.update(kw)
    return d


def test_aggregate_summarizes_excluded_classes():
    verdicts = [
        _verdict("A", "include"),
        _verdict("B", "exclude", n_infiltrations=2, reasons=["2 anchors nest in a supported in-group clade"]),
        _verdict("C", "include"),
        _verdict("F", "exclude", rf=0.2, reasons=["in-group RF 0.200 > 0.1"]),
    ]
    report = agg.aggregate_verdicts(verdicts)
    assert report["keep_outgroup_anchors"] == ["A", "C"]
    assert report["exclude_outgroup_anchors"] == ["B", "F"]
    assert report["n_classes"] == 4
    assert {r["class"] for r in report["per_class"]} == {"A", "B", "C", "F"}


def test_render_markdown_has_table_and_recommendation():
    report = agg.aggregate_verdicts([
        _verdict("A", "include"),
        _verdict("B", "exclude", n_infiltrations=1, reasons=["infiltration"]),
    ])
    md = agg.render_markdown(report)
    assert "| class |" in md.lower() or "| Class |" in md
    assert "A" in md and "B" in md
    assert "include" in md and "exclude" in md
    # report-only: surfaces a recommendation, does not claim to have applied it
    assert "ANCHOR_OUTGROUP_CLASSES" in md or "keep" in md.lower()


def test_load_and_aggregate_dir(tmp_path):
    for k, v in [("A", "include"), ("B", "exclude")]:
        (tmp_path / f"verdict_class_{k}.json").write_text(json.dumps(_verdict(k, v)))
    report = agg.aggregate_dir(str(tmp_path))
    assert report["keep_outgroup_anchors"] == ["A"]
    assert report["exclude_outgroup_anchors"] == ["B"]
