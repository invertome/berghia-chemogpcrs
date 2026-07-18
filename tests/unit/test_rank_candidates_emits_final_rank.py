"""Regression: rank_candidates.py must emit a final_rank column carrying the
production order (weighted OR rankagg), so emit_ranked_views' confidence view
sorts by it rather than re-deriving order from the weighted rank_score.

Running main() needs the full pipeline inputs, so (mirroring the lightweight
source-level check in test_rank_candidates_emits_orthogroup.py) this asserts the
two lines of the contract are present in the script source:
  1. df_sorted['final_rank'] is assigned as a 1-based range, AND
  2. 'final_rank' is listed in output_cols (so it reaches the CSV).
The behavioral half (confidence view honoring final_rank) is covered in
test_emit_ranked_views.py::test_confidence_view_honors_final_rank_over_rank_score.
"""
from __future__ import annotations

import os
import re


def _src() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.normpath(os.path.join(here, "..", ".."))
    with open(os.path.join(repo, "scripts", "rank_candidates.py")) as f:
        return f.read()


def test_final_rank_assigned_as_1_based_range():
    assert re.search(r"df_sorted\[['\"]final_rank['\"]\]\s*=\s*range\(\s*1,", _src()), \
        "rank_candidates.py must assign a 1-based final_rank on the production order"


def test_final_rank_in_output_cols():
    m = re.search(r"output_cols\s*=\s*\[(.*?)\]", _src(), re.DOTALL)
    assert m, "could not find output_cols list in rank_candidates.py"
    assert "'final_rank'" in m.group(1) or '"final_rank"' in m.group(1), \
        "final_rank must be in output_cols so it reaches the ranked CSV"
