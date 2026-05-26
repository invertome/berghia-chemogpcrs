"""Tests for scripts/montage_structure_figures.py.

Assembles per-structure PNG renders into a single labelled grid montage for
slides/manuscripts. The grid geometry is deterministic so we can assert exact
canvas dimensions.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import montage_structure_figures as msf


def _dummy_pngs(tmp_path: Path, n: int) -> list[str]:
    from PIL import Image

    paths = []
    for i in range(n):
        p = tmp_path / f"struct_{i}.png"
        Image.new("RGB", (100, 100), (10 * i, 10 * i, 10 * i)).save(p)
        paths.append(str(p))
    return paths


def test_montage_grid_dimensions(tmp_path: Path) -> None:
    """3 images at cols=2 -> 2x2 grid; canvas size follows the layout formula."""
    paths = _dummy_pngs(tmp_path, 3)
    out = tmp_path / "montage.png"

    msf.build_montage(paths, ["a", "b", "c"], str(out), cols=2,
                      cell=(50, 50), pad=10, label_h=20)

    assert out.exists()
    from PIL import Image
    im = Image.open(out)
    # width  = pad + cols*(cell_w + pad)            = 10 + 2*(50+10) = 130
    # height = pad + rows*(cell_h + label_h + pad)  = 10 + 2*(50+20+10) = 170
    assert im.size == (130, 170)


def test_montage_empty_raises(tmp_path: Path) -> None:
    """No input images is an error, not a blank canvas."""
    with pytest.raises(ValueError):
        msf.build_montage([], [], str(tmp_path / "m.png"))
