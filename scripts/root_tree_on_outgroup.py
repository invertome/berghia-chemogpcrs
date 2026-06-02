#!/usr/bin/env python3
"""root_tree_on_outgroup.py — root a per-class tree on its swap-map outgroup.

Bead vo8.3 (per-class refactor follow-up). IQ-TREE leaves the per-class trees
unrooted. The per-class architecture adds a sister-class outgroup to each tree
(swap-map A<-C, B/C/F<-A), so outgroup rooting is valid here (unlike the old
combined multi-class tree, where midpoint rooting was used). This roots the
tree on the outgroup taxa, falling back to midpoint when the outgroup is absent
or non-monophyletic.

Produces a ``*.rooted.treefile`` sibling; the raw IQ-TREE ``*.treefile`` is left
untouched for consumers that want the unrooted topology.

Usage:
    python3 root_tree_on_outgroup.py \\
        --tree   class_A.treefile \\
        --outgroup-ids outgroup_class_A.fa \\
        --out    class_A.rooted.treefile \\
        [--fallback midpoint|none]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional


def parse_outgroup_ids(path: Path) -> list[str]:
    """Read outgroup leaf names from a plain id file or a FASTA.

    FASTA: each ``>`` header's first whitespace-delimited token is the id.
    Plain: one id per non-empty line.
    """
    text = path.read_text()
    if any(line.startswith(">") for line in text.splitlines()):
        ids = [line[1:].split()[0] for line in text.splitlines()
               if line.startswith(">") and line[1:].split()]
    else:
        ids = [line.strip() for line in text.splitlines() if line.strip()]
    # de-dupe preserving order
    seen: set[str] = set()
    out: list[str] = []
    for i in ids:
        if i not in seen:
            seen.add(i)
            out.append(i)
    return out


def load_tree(path: Path):
    """Load a Newick tree, trying ete3 formats tolerant of IQ-TREE support labels."""
    from ete3 import Tree
    last_exc: Optional[Exception] = None
    for fmt in (1, 0, 5, 2, 3):
        try:
            return Tree(str(path), format=fmt)
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
    raise ValueError(f"could not parse tree {path}: {last_exc}")


def root_tree(tree, outgroup_ids: list[str], fallback: str = "midpoint") -> str:
    """Root ``tree`` in place on the outgroup; return the mode used.

    Modes: outgroup_single, outgroup_clade, midpoint, unrooted.
    """
    leaves = set(tree.get_leaf_names())
    og = [n for n in outgroup_ids if n in leaves]

    if len(og) == 1:
        tree.set_outgroup(og[0])
        return "outgroup_single"

    if len(og) >= 2 and len(og) < len(leaves):
        is_mono = tree.check_monophyly(values=og, target_attr="name")[0]
        if is_mono:
            tree.set_outgroup(tree.get_common_ancestor(og))
            return "outgroup_clade"

    # fallback
    if fallback == "midpoint":
        try:
            mid = tree.get_midpoint_outgroup()
        except Exception:  # noqa: BLE001
            mid = None
        if mid is not None:
            tree.set_outgroup(mid)
            return "midpoint"
    return "unrooted"


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--tree", type=Path, required=True, help="Input IQ-TREE .treefile")
    p.add_argument("--outgroup-ids", type=Path, required=True,
                   help="Outgroup leaf ids — plain id file or FASTA of outgroup seqs")
    p.add_argument("--out", type=Path, required=True, help="Output rooted .treefile")
    p.add_argument("--fallback", choices=["midpoint", "none"], default="midpoint",
                   help="What to do when outgroup is absent/non-monophyletic")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    if not args.tree.exists():
        print(f"[root_tree_on_outgroup] ERROR: tree not found: {args.tree}",
              file=sys.stderr)
        return 1

    outgroup_ids = (parse_outgroup_ids(args.outgroup_ids)
                    if args.outgroup_ids.exists() else [])
    if not outgroup_ids:
        print(f"[root_tree_on_outgroup] WARN: no outgroup ids in "
              f"{args.outgroup_ids}", file=sys.stderr)

    tree = load_tree(args.tree)
    mode = root_tree(tree, outgroup_ids, fallback=args.fallback)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    tree.write(outfile=str(args.out), format=1)
    print(f"[root_tree_on_outgroup] rooted ({mode}) → {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
