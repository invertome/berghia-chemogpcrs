#!/usr/bin/env python3
"""build_og_class_majority.py — assign each orthogroup a majority GPCR class.

Bead vo8.2 (per-class refactor follow-up). Stage 04's per-OG routing reads
``results/classification/og_class_majority.tsv`` (orthogroup<TAB>class) to put
each orthogroup tree under its class subdirectory. Nothing produced that file,
so every OG fell back to class_A. This producer derives it.

Inputs:
  - OrthoFinder ``Orthogroups.tsv`` (OG_id<TAB>sp1 members<TAB>sp2 ...; members
    comma-space separated within a cell).
  - One or more classifier TSVs from classify_gpcr_by_class.py
    (seq_id<TAB>class<TAB>...; class in A/B/C/F). Typically candidate_classes.tsv
    plus class_berghia.tsv.

Output (only written when there is real data — see ``main``):
  - og_class_majority.tsv: ``orthogroup<TAB>class`` (header + one row per OG
    that has at least one classified member). LF line endings (a trailing CR
    would break stage 04's ``awk ... $1==og {print $2}`` match).

Usage:
    python3 build_og_class_majority.py \\
        --orthogroups results/orthogroups/OrthoFinder/Results*/Orthogroups/Orthogroups.tsv \\
        --classes results/.../classify/candidate_classes.tsv \\
        --classes results/.../classify/class_berghia.tsv \\
        --out results/classification/og_class_majority.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

# Tie-break order: deepest-diverging class first so a split OG roots sensibly.
CLASS_ORDER = ["A", "B", "C", "F"]
_MEMBER_SPLIT = re.compile(r",\s*")


def parse_class_tsv(path: Path) -> dict[str, str]:
    """Read a classifier TSV → {seq_id: class}. Skips header and empty classes."""
    out: dict[str, str] = {}
    with path.open() as fh:
        next(fh, None)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            seq_id, cls = parts[0].strip(), parts[1].strip()
            if seq_id and cls:
                out[seq_id] = cls
    return out


def parse_orthogroups_tsv(path: Path) -> dict[str, list[str]]:
    """Read OrthoFinder Orthogroups.tsv → {og_id: [member seq_ids]}."""
    out: dict[str, list[str]] = {}
    with path.open() as fh:
        next(fh, None)  # header
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            if not cells or not cells[0].strip():
                continue
            og = cells[0].strip()
            members: list[str] = []
            for cell in cells[1:]:
                cell = cell.strip()
                if not cell:
                    continue
                members.extend(m for m in _MEMBER_SPLIT.split(cell) if m)
            out[og] = members
    return out


def majority_class(classes: list[str]) -> Optional[str]:
    """Plurality class; ties broken by CLASS_ORDER. None if empty."""
    if not classes:
        return None
    counts = Counter(classes)
    top = max(counts.values())

    def order_key(c: str) -> int:
        return CLASS_ORDER.index(c) if c in CLASS_ORDER else len(CLASS_ORDER)

    return min((c for c, n in counts.items() if n == top), key=order_key)


def build(
    orthogroups: dict[str, list[str]],
    class_map: dict[str, str],
    *,
    min_fraction: float = 0.0,
) -> list[tuple[str, str]]:
    """For each OG, majority class of its classified members.

    OGs with no classified member are omitted (stage 04 then routes them to
    'unclassified'). If ``min_fraction`` > 0, the winning class must cover at
    least that fraction of the OG's classified members, else the OG is omitted.
    """
    rows: list[tuple[str, str]] = []
    for og in sorted(orthogroups):
        member_classes = [class_map[m] for m in orthogroups[og] if m in class_map]
        cls = majority_class(member_classes)
        if cls is None:
            continue
        if min_fraction > 0.0:
            frac = member_classes.count(cls) / len(member_classes)
            if frac < min_fraction:
                continue
        rows.append((og, cls))
    return rows


def write_tsv(rows: list[tuple[str, str]], out_path: Path) -> None:
    """Write ``orthogroup<TAB>class`` with a header. LF line endings."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        fh.write("orthogroup\tclass\n")
        for og, cls in rows:
            fh.write(f"{og}\t{cls}\n")


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--orthogroups", type=Path, required=True,
                   help="OrthoFinder Orthogroups.tsv")
    p.add_argument("--classes", type=Path, action="append", required=True,
                   help="Classifier seq_id/class TSV (repeatable)")
    p.add_argument("--out", type=Path, required=True,
                   help="Output og_class_majority.tsv")
    p.add_argument("--min-fraction", type=float, default=0.0,
                   help="Min fraction of classified members the winning class "
                        "must cover (default 0 = plurality)")
    p.add_argument("--force", action="store_true",
                   help="Overwrite --out if it already exists")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    if args.out.exists() and not args.force:
        print(f"[build_og_class_majority] {args.out} exists; skipping (use --force).",
              file=sys.stderr)
        return 0

    if not args.orthogroups.exists():
        print(f"[build_og_class_majority] Orthogroups not found ({args.orthogroups}); "
              "writing nothing so stage 04 keeps its class_A default.", file=sys.stderr)
        return 0

    class_map: dict[str, str] = {}
    for cf in args.classes:
        if cf.exists():
            class_map.update(parse_class_tsv(cf))
        else:
            print(f"[build_og_class_majority] WARN: class file missing: {cf}",
                  file=sys.stderr)

    if not class_map:
        print("[build_og_class_majority] no sequence classes available; writing "
              "nothing so stage 04 keeps its class_A default.", file=sys.stderr)
        return 0

    orthogroups = parse_orthogroups_tsv(args.orthogroups)
    rows = build(orthogroups, class_map, min_fraction=args.min_fraction)
    write_tsv(rows, args.out)
    print(f"[build_og_class_majority] {len(rows)}/{len(orthogroups)} orthogroups "
          f"classified → {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
