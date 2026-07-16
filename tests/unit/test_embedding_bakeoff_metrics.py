"""Tests for the harness-log -> tidy bake-off metrics parser.

The proven harness (scratch_lofo_bakeoff.py) PRINTS a fixed-width table per
(model, refset); a cross-model view needs those numbers as a tidy table. Parsing
the stdout logs (rather than patching the harness) works on the 16 runs already
completed AND future ones, and never touches the validated metric code.
"""
from __future__ import annotations

import textwrap

from embedding_bakeoff_metrics import BAKEOFF_COLUMNS, parse_harness_log

# Two refsets of one model; the round-2-style `nan` invariance (no candidates)
# must parse, and both maha and cosine config rows must be captured.
SAMPLE_LOG = textwrap.dedent("""\
    === model=protrek ref=reference_protrek_PROD.npz cand=candidates_protrek_classA.npz ===

    ########## protrek :: FULL-clean ##########
    # reference_protrek_PROD.npz: 1094 refs, 11 families, 790 candidates, dim 1024

    score     calib   cent     LOFO-AUROC  invariance  LOO-famAcc
    --------------------------------------------------------------
    cos       raw     single        0.900       0.780       0.657
    maha      raw     single        0.994       1.000       0.760
    maha      raw     multi         0.994       1.000       0.794

    ########## protrek :: VERIFIED-clean ##########
    # reference_protrek_PROD.npz: 723 refs, 11 families, 790 candidates, dim 1024

    score     calib   cent     LOFO-AUROC  invariance  LOO-famAcc
    --------------------------------------------------------------
    maha      raw     multi         0.997       nan         0.791

    === harness_clean protrek DONE ===
    """)


def test_parse_returns_one_row_per_model_refset_config():
    rows = parse_harness_log(SAMPLE_LOG)
    # 3 config rows in FULL + 1 in VERIFIED
    assert len(rows) == 4


def test_parse_extracts_model_refset_and_metadata():
    rows = parse_harness_log(SAMPLE_LOG)
    full = [r for r in rows if r["refset"] == "FULL"]
    assert {r["model"] for r in full} == {"protrek"}
    r = full[0]
    assert r["n_refs"] == 1094 and r["n_families"] == 11 and r["dim"] == 1024


def test_parse_reads_the_maha_multi_row_values():
    rows = parse_harness_log(SAMPLE_LOG)
    mm = next(r for r in rows if r["refset"] == "FULL" and r["score"] == "maha" and r["cent"] == "multi")
    assert mm["lofo_auroc"] == 0.994
    assert mm["loo_famacc"] == 0.794
    assert mm["invariance"] == 1.000


def test_parse_accepts_nan_invariance_from_candidateless_runs():
    rows = parse_harness_log(SAMPLE_LOG)
    ver = next(r for r in rows if r["refset"] == "VERIFIED")
    assert ver["invariance"] != ver["invariance"]  # NaN
    assert ver["n_refs"] == 723 and ver["lofo_auroc"] == 0.997


def test_rows_match_the_declared_schema():
    rows = parse_harness_log(SAMPLE_LOG)
    for r in rows:
        assert list(r.keys()) == list(BAKEOFF_COLUMNS)


def test_empty_or_headerless_log_yields_no_rows():
    assert parse_harness_log("") == []
    assert parse_harness_log("nothing to see here\n") == []
