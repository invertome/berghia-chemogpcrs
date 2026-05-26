#!/usr/bin/env python3
"""Assemble per-structure PNG renders into one labelled grid montage.

Stage 08 renders each AlphaFold 3 prediction to a PNG (via PyMOL); this pools
them into a single figure suitable for slides / supplementary panels. Labels
default to each file's stem (the candidate id).
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Sequence


def build_montage(
    image_paths: Sequence[str],
    labels: Sequence[str],
    output_path: str,
    cols: int | None = None,
    cell: tuple[int, int] = (400, 400),
    pad: int = 10,
    label_h: int = 20,
    bg: tuple[int, int, int] = (255, 255, 255),
) -> str:
    """Tile ``image_paths`` into a grid PNG with a label under each cell.

    Canvas geometry (deterministic):
        width  = pad + cols * (cell_w + pad)
        height = pad + rows * (cell_h + label_h + pad)
    """
    from PIL import Image, ImageDraw

    n = len(image_paths)
    if n == 0:
        raise ValueError("build_montage requires at least one image")
    if cols is None:
        cols = math.ceil(math.sqrt(n))
    cols = max(1, cols)
    rows = math.ceil(n / cols)
    cw, ch = cell

    width = pad + cols * (cw + pad)
    height = pad + rows * (ch + label_h + pad)
    canvas = Image.new("RGB", (width, height), bg)
    draw = ImageDraw.Draw(canvas)

    for idx, path in enumerate(image_paths):
        r, c = divmod(idx, cols)
        x = pad + c * (cw + pad)
        y = pad + r * (ch + label_h + pad)
        img = Image.open(path).convert("RGB").resize((cw, ch))
        canvas.paste(img, (x, y))
        label = labels[idx] if idx < len(labels) else ""
        if label:
            draw.text((x, y + ch + 4), label, fill=(0, 0, 0))

    canvas.save(output_path)
    return output_path


def main(argv: Sequence[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("images", nargs="+", help="per-structure PNG files")
    ap.add_argument("--out", required=True, help="output montage PNG")
    ap.add_argument("--cols", type=int, default=None, help="columns (default: ~sqrt(n))")
    args = ap.parse_args(argv)
    labels = [Path(p).stem for p in args.images]
    build_montage(args.images, labels, args.out, cols=args.cols)


if __name__ == "__main__":
    main()
