# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///
"""Render a predicted protein structure to a publication PNG (+ PyMOL session).

Default coloring is the canonical AlphaFold pLDDT confidence palette — pLDDT is
stored per-atom in the mmCIF B-factor column. Use ``--color-mode ss`` for
secondary-structure coloring (e.g. experimental reference structures that have
no pLDDT). Headless OSMesa software rendering: no GPU, display, or X server
needed (suitable for an HPC compute node). Generated with the `pymol` skill.

Run standalone via uv (auto-installs pymol-open-source-whl from the header):
    uv run scripts/render_structure_figure.py model.cif out.png --session out.pse
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Sequence


def plddt_color_rules() -> list[tuple[str, tuple[float, float, float], str]]:
    """Canonical AlphaFold pLDDT confidence bands (very-low -> very-high).

    Each rule is (color_name, (r,g,b) in 0-1, PyMOL B-factor selection).
    Thresholds match the AlphaFold DB legend: <50 orange, 50-70 yellow,
    70-90 light-blue, >=90 dark-blue. Selections are CUMULATIVE and MUST be
    applied in list order (low->high) so each higher band overwrites the
    previous fill. PyMOL's selection lexer rejects ">="/"<=", so we use the
    "color all, then b > X upwards" idiom rather than half-open ranges.
    """
    return [
        ("af_vlow", (0xFF / 255, 0x7D / 255, 0x45 / 255), "all"),
        ("af_low", (0xFF / 255, 0xDB / 255, 0x13 / 255), "b > 50"),
        ("af_conf", (0x65 / 255, 0xCB / 255, 0xF3 / 255), "b > 70"),
        ("af_high", (0x00 / 255, 0x53 / 255, 0xD6 / 255), "b > 90"),
    ]


def render(
    structure_path: str,
    png_path: str,
    session_path: str | None = None,
    color_mode: str = "plddt",
    width: int = 1600,
    height: int = 1200,
) -> int:
    """Render one structure; returns 0 on success, non-zero on failure."""
    os.environ["PYOPENGL_PLATFORM"] = "osmesa"
    import pymol  # deferred so unit tests can import this module without PyMOL

    pymol.pymol_argv = ["pymol", "-cq"]
    pymol.finish_launching()
    from pymol import cmd

    if not os.path.exists(structure_path):
        print(f"[render] ERROR: structure not found: {structure_path}", file=sys.stderr)
        cmd.quit()
        return 2
    cmd.load(structure_path, "structure")
    if cmd.count_atoms("all") == 0:
        print(f"[render] ERROR: no atoms loaded from {structure_path}", file=sys.stderr)
        cmd.quit()
        return 3

    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)

    if color_mode == "plddt":
        for name, rgb, sel in plddt_color_rules():
            cmd.set_color(name, list(rgb))
            cmd.color(name, sel)
    else:  # secondary-structure coloring
        cmd.color("marine", "ss h")
        cmd.color("orange", "ss s")
        cmd.color("grey80", "ss l+''")

    cmd.orient()
    if os.path.dirname(png_path):
        os.makedirs(os.path.dirname(png_path), exist_ok=True)
    cmd.png(png_path, width=width, height=height, dpi=300, ray=1)
    if session_path:
        cmd.save(session_path)
    cmd.quit()
    return 0


def main(argv: Sequence[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("structure", help="input structure (mmCIF or PDB)")
    ap.add_argument("png", help="output PNG path")
    ap.add_argument("--session", default=None, help="also save a .pse session")
    ap.add_argument("--color-mode", choices=["plddt", "ss"], default="plddt")
    ap.add_argument("--width", type=int, default=1600)
    ap.add_argument("--height", type=int, default=1200)
    a = ap.parse_args(argv)
    rc = render(a.structure, a.png, a.session, a.color_mode, a.width, a.height)
    sys.exit(rc)


if __name__ == "__main__":
    main()
