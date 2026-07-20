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

# classify_gpcr_by_class.py writes class in {A, B, C, F, unclassified}.
# 'unclassified' is a SENTINEL meaning "no class evidence", not a class: it must
# never cast a vote, never enter the min-fraction denominator, and never win a
# majority (which would emit a confident assignment whose content is "we do not
# know"). Matched case-insensitively.
#
# This is also the exact token stage 04 and stage 05 route on
# (results/phylogenies/protein/class_${OG_CLASS}/), so it doubles as the value
# written for an orthogroup with no real votes.
UNCLASSIFIED = "unclassified"


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


def real_votes(classes: list[str]) -> list[str]:
    """Drop the 'unclassified' sentinel — it is an absence of evidence."""
    return [c for c in classes if c.strip().lower() != UNCLASSIFIED]


def majority_class(classes: list[str]) -> Optional[str]:
    """Plurality class; ties broken by CLASS_ORDER. None if no real votes.

    The 'unclassified' sentinel is excluded before counting, so it can neither
    win nor dilute. ``None`` means "no class evidence" — the caller states that
    explicitly rather than letting an absent row imply it.
    """
    classes = real_votes(classes)
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
    """For each OG, majority class of its members that carry real class evidence.

    Every orthogroup gets exactly one row. An OG with no real votes — no
    classified member, or only 'unclassified' ones — is written as
    ``unclassified`` EXPLICITLY rather than omitted, so a failed lookup
    downstream means a genuinely missing or truncated table instead of
    doubling as "no evidence".

    If ``min_fraction`` > 0, the winning class must cover at least that
    fraction of the OG's REAL votes (the sentinel is not in the denominator);
    otherwise the OG is likewise reported ``unclassified``.
    """
    rows: list[tuple[str, str]] = []
    for og in sorted(orthogroups):
        member_classes = real_votes(
            [class_map[m] for m in orthogroups[og] if m in class_map])
        cls = majority_class(member_classes)
        if cls is not None and min_fraction > 0.0:
            if member_classes.count(cls) / len(member_classes) < min_fraction:
                cls = None
        rows.append((og, cls if cls is not None else UNCLASSIFIED))
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
        if not cf.exists():
            print(f"[build_og_class_majority] WARN: class file missing: {cf}",
                  file=sys.stderr)
            continue
        incoming = parse_class_tsv(cf)
        # Stage 04 globs class_*.tsv, so several tables can classify the same
        # seq_id. A bare .update() lets the last file silently win a genuine
        # disagreement; report it instead of resolving it invisibly.
        conflicts = sorted(s for s, c in incoming.items()
                           if s in class_map and class_map[s] != c)
        if conflicts:
            print(f"[build_og_class_majority] WARN: {cf.name} reclassifies "
                  f"{len(conflicts)} seq_id(s) already classified differently "
                  f"(last file wins): "
                  f"{', '.join(f'{s} {class_map[s]}->{incoming[s]}' for s in conflicts[:5])}"
                  f"{' ...' if len(conflicts) > 5 else ''}", file=sys.stderr)
        class_map.update(incoming)

    if not class_map:
        print("[build_og_class_majority] no sequence classes available; writing "
              "nothing so stage 04 keeps its class_A default.", file=sys.stderr)
        return 0

    orthogroups = parse_orthogroups_tsv(args.orthogroups)
    rows = build(orthogroups, class_map, min_fraction=args.min_fraction)
    write_tsv(rows, args.out)
    # Count real assignments, not rows: every OG now has a row, and an
    # 'unclassified' one is an explicit non-answer.
    n_classified = sum(1 for _og, cls in rows if cls != UNCLASSIFIED)
    print(f"[build_og_class_majority] {n_classified}/{len(orthogroups)} orthogroups "
          f"classified ({len(rows) - n_classified} explicitly unclassified) "
          f"→ {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
