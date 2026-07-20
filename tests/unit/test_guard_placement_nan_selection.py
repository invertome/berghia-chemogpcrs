"""Guard: a non-finite LWR must never shadow, nor masquerade as, a measurement.

``test_classification_vocab_placement.py`` already pins the GATE: a NaN LWR
reaching ``classify_placement_for_id`` no longer falls through
``NaN < 0.80 == False`` into a confident call. This file pins the two holes
upstream and downstream of that gate.

1. SELECTION (``parse_jplace``). The best placement was chosen with
   ``max(..., key=lambda p: p[lwr_idx])``. ``max`` keeps its running winner
   unless a later item compares strictly greater, and every comparison with NaN
   is False -- so a NaN placement appearing FIRST is never displaced, and it
   shadows a perfectly good measured placement behind it. The candidate then
   hits the (correct) NaN gate and is reported ``unclassified-placement``: a
   real, confident measurement silently downgraded to "could not measure".

2. REPORTING. The gate reports a degenerate LWR as ``0.0``, which the output
   TSV wrote as ``0.000`` -- indistinguishable from a genuine, measured LWR of
   zero. "No data" and "measured a low value" must be distinguishable.

Why this matters: placement is one of three consensus sources in the
non-chemoreceptor classifier. Losing a real placement moves a candidate from
3-of-3 to 2-of-3 agreement, changing its confidence tier -- and exclusion is
the pipeline's load-bearing output.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import classify_via_placement as cvp  # noqa: E402


EDGE_MAP = {3: ("aminergic", "5HT"), 7: ("opsin", "")}


def _jplace(tmp_path: Path, placements: list[dict]) -> str:
    p = tmp_path / "epa_result.jplace"
    p.write_text(json.dumps({
        "fields": ["edge_num", "like_weight_ratio"],
        "placements": placements,
        "tree": "();",
    }))
    return str(p)


# ---------------------------------------------------------------------------
# 1. A non-finite LWR must not shadow a real one
# ---------------------------------------------------------------------------

def test_nan_listed_first_does_not_shadow_a_confident_placement(tmp_path: Path):
    """The bug: max() returns the NaN entry (edge 7), losing the real 0.97
    placement on edge 3, and the candidate is reported unclassified."""
    path = _jplace(tmp_path, [
        {"n": ["cand1"], "p": [[7, float("nan")], [3, 0.97]]},
    ])
    got = cvp.parse_jplace(path)
    assert got["cand1"]["edge_num"] == 3
    assert got["cand1"]["lwr"] == pytest.approx(0.97)

    call = cvp.classify_placement_for_id("cand1", got, EDGE_MAP)
    assert call["family"] == "aminergic", "a measured placement was lost to NaN"


def test_infinities_do_not_shadow_either(tmp_path: Path):
    path = _jplace(tmp_path, [
        {"n": ["cand1"], "p": [[7, float("inf")], [3, 0.91]]},
    ])
    got = cvp.parse_jplace(path)
    assert got["cand1"]["edge_num"] == 3


def test_best_finite_placement_wins_among_several(tmp_path: Path):
    path = _jplace(tmp_path, [
        {"n": ["cand1"], "p": [[7, 0.20], [3, 0.55], [7, float("nan")]]},
    ])
    got = cvp.parse_jplace(path)
    assert got["cand1"]["edge_num"] == 3
    assert got["cand1"]["lwr"] == pytest.approx(0.55)


def test_json_nan_literal_is_the_real_ingress(tmp_path: Path):
    """EPA-ng writing a bare NaN token is accepted by json.load, which is how a
    NaN gets in at all. Pinned so nobody 'simplifies' this away as impossible."""
    p = tmp_path / "epa_result.jplace"
    p.write_text('{"fields": ["edge_num", "like_weight_ratio"], '
                 '"placements": [{"n": ["c"], "p": [[7, NaN], [3, 0.99]]}], '
                 '"tree": "();"}')
    got = cvp.parse_jplace(str(p))
    assert got["c"]["edge_num"] == 3


# ---------------------------------------------------------------------------
# 2. Undetermined must be an explicit state, not a plausible number
# ---------------------------------------------------------------------------

def test_all_placements_non_finite_is_explicitly_undetermined(tmp_path: Path):
    path = _jplace(tmp_path, [
        {"n": ["cand1"], "p": [[7, float("nan")], [3, float("nan")]]},
    ])
    got = cvp.parse_jplace(path)
    assert got["cand1"]["lwr_determined"] is False
    call = cvp.classify_placement_for_id("cand1", got, EDGE_MAP)
    assert call["family"] == "unclassified-placement"
    assert call["lwr_determined"] is False


def test_measured_placement_is_marked_determined(tmp_path: Path):
    path = _jplace(tmp_path, [{"n": ["cand1"], "p": [[3, 0.97]]}])
    got = cvp.parse_jplace(path)
    assert got["cand1"]["lwr_determined"] is True
    assert cvp.classify_placement_for_id(
        "cand1", got, EDGE_MAP)["lwr_determined"] is True


def test_a_genuinely_measured_zero_is_not_undetermined():
    """The distinction the TSV has to carry: 0.0 measured != 0.0 unmeasured."""
    call = cvp.classify_placement_for_id(
        "q", {"q": {"edge_num": 3, "lwr": 0.0, "lwr_determined": True}}, EDGE_MAP)
    assert call["family"] == "unclassified-placement"   # below threshold
    assert call["lwr_determined"] is True               # but it WAS measured


def test_candidate_absent_from_placements_is_undetermined():
    """EPA-ng returning nothing for a query is 'no data', not 'LWR = 0'."""
    call = cvp.classify_placement_for_id("nope", {}, EDGE_MAP)
    assert call["family"] == "unclassified-placement"
    assert call["lwr_determined"] is False


def test_nan_still_reports_a_non_nan_score():
    """Pre-existing contract (test_classification_vocab_placement.py): the NaN
    must not leak out as a score, or it compares False downstream too."""
    call = cvp.classify_placement_for_id(
        "q", {"q": {"edge_num": 3, "lwr": float("nan")}}, EDGE_MAP)
    assert not math.isnan(call["lwr"])
    assert call["lwr"] == 0.0
    assert call["lwr_determined"] is False


# ---------------------------------------------------------------------------
# 3. Malformed rows must not crash the parser
# ---------------------------------------------------------------------------

def test_short_placement_row_is_skipped_not_fatal(tmp_path: Path):
    path = _jplace(tmp_path, [{"n": ["cand1"], "p": [[7], [3, 0.93]]}])
    got = cvp.parse_jplace(path)
    assert got["cand1"]["edge_num"] == 3


def test_empty_placement_list_yields_no_entry(tmp_path: Path):
    path = _jplace(tmp_path, [{"n": ["cand1"], "p": []}])
    assert cvp.parse_jplace(path) == {}


# ---------------------------------------------------------------------------
# 4. The TSV must distinguish the two
# ---------------------------------------------------------------------------

def test_format_lwr_marks_undetermined_distinctly():
    assert cvp.format_lwr({"lwr": 0.0, "lwr_determined": False}) == "NA"
    assert cvp.format_lwr({"lwr": 0.0, "lwr_determined": True}) == "0.000"
    assert cvp.format_lwr({"lwr": 0.955, "lwr_determined": True}) == "0.955"


def test_format_lwr_defaults_to_determined_for_legacy_rows():
    """Rows written before the flag existed keep their historical meaning."""
    assert cvp.format_lwr({"lwr": 0.42}) == "0.420"
