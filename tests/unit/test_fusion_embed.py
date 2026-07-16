import sys, numpy as np
sys.path.insert(0, "scripts")
from fusion_embed import concat_fuse

def test_concat_fuse_intersects_and_reports_dropped():
    a = {"x": np.array([3.0, 0.0]), "y": np.array([1.0, 1.0]), "z": np.array([1.0, 0.0])}
    b = {"x": np.array([0.0, 4.0]), "y": np.array([2.0, 0.0])}   # no z
    fused, dropped = concat_fuse(a, b)
    assert set(fused) == {"x", "y"} and dropped == ["z"]
    assert fused["x"].shape == (4,)
    assert np.isclose(np.linalg.norm(fused["x"][:2]), 1.0)      # each half unit-norm

def test_concat_fuse_pca_dim_and_no_overlap():
    a = {"x": np.array([1.0, 0.0]), "y": np.array([0.0, 1.0])}
    b = {"p": np.array([1.0, 0.0])}
    fused, dropped = concat_fuse(a, b)
    assert fused == {} and sorted(dropped) == ["p", "x", "y"]
