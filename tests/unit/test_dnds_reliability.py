import numpy as np
from dnds_reliability import reliability_shrink   # conftest adds scripts/ to path


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
