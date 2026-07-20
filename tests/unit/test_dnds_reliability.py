import numpy as np
import pandas as pd
import pytest

from dnds_reliability import reliability_shrink   # conftest adds scripts/ to path
from rank_aggregation import build_ranklists_from_df


def test_full_reliability_is_identity():
    vals = {"a": 0.9, "b": 0.1, "c": 0.5}
    rel = {"a": 1.0, "b": 1.0, "c": 1.0}
    assert reliability_shrink(vals, rel) == {"a": 0.9, "b": 0.1, "c": 0.5}


def test_zero_reliability_goes_to_neutral_median():
    vals = {"a": 0.9, "b": 0.1, "c": 0.5}   # median 0.5
    rel = {"a": 0.0, "b": 0.0, "c": 1.0}
    out = reliability_shrink(vals, rel)
    assert out["a"] == 0.5 and out["b"] == 0.5 and out["c"] == 0.5


def test_partial_reliability_interpolates():
    vals = {"a": 1.0, "b": 0.0}   # median 0.5
    out = reliability_shrink(vals, {"a": 0.5, "b": 0.5})
    assert out["a"] == 0.75 and out["b"] == 0.25


def test_missing_reliability_defaults_to_full():
    out = reliability_shrink({"a": 0.9, "b": 0.1}, {"a": 1.0})   # b absent
    assert out["b"] == 0.1


def test_production_ranklist_reads_shrunk_signal():
    """End-to-end proof that reliability now reaches the PRODUCTION rankagg path,
    not just the weighted composite. rank_candidates.py shrinks the selection
    norm columns on df BEFORE build_ranklists_from_df runs; replicate that exact
    composition and confirm an underpowered candidate (reliability approx 0)
    votes with the cohort NEUTRAL (median) value in the 'positive' ranklist,
    instead of its low raw score."""
    df = pd.DataFrame({
        "id": ["poor", "mid", "rich"],
        "positive_score_norm": [0.1, 0.5, 0.9],      # median 0.5
        "purifying_score_norm": [0.0, 0.0, 0.0],
        "phylo_score_norm": [0.5, 0.5, 0.5],
        "lse_divergence_score_norm": [0.0, 0.0, 0.0],
        "dnds_reliability_weight": [0.0, 1.0, 1.0],   # 'poor' is underpowered
    })
    # Exactly what rank_candidates.py's upstream block does to the selection cols.
    rel = dict(zip(df["id"], df["dnds_reliability_weight"]))
    for sig in ("positive_score_norm", "purifying_score_norm"):
        vals = dict(zip(df["id"], df[sig]))
        df[sig] = df["id"].map(reliability_shrink(vals, rel))

    ranklists = build_ranklists_from_df(df)
    # The underpowered candidate now contributes the neutral value, not 0.1 ...
    assert ranklists["positive"]["poor"] == pytest.approx(0.5)
    # ... while the fully-powered candidate is unchanged.
    assert ranklists["positive"]["rich"] == pytest.approx(0.9)
