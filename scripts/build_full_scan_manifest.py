#!/usr/bin/env python3
"""build_full_scan_manifest.py — build the full-557 (P6) scan manifest.

P6 prep (bead kib). The P5 manifest builder reads the Phase-1a species list;
P6 scans every consolidated proteome. consolidate_proteomes_for_genome_wide_og.py
writes one file per species-tree leaf to
``${GENOME_WIDE_ORTHOFINDER_DIR}/input/<taxid>_<sanitized_binomial>.fa``. This
enumerates that directory into the same 3-column manifest the scan array
consumes (via sbatch_run_*_scan.sh, ``cut -f3`` → proteome path):

    taxid  binomial  proteome_path

LF line endings (a trailing CR would survive the wrapper's cut -f3 and break
the per-row ``[ -f <path> ]`` check — the P5 bug).

Usage:
    python3 build_full_scan_manifest.py \\
        --proteomes-dir results/orthogroups_genome_wide/input \\
        --out p6_full_scan_manifest.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

_LEAF_RE = re.compile(r"^(\d+)_(.+)$")
_EXTS = (".aa.fna", ".fasta", ".faa", ".fna", ".fa")


@dataclass(frozen=True)
class ManifestRow:
    taxid: int
    binomial: str
    proteome_path: Path


def parse_leaf_filename(name: str) -> Optional[tuple[int, str]]:
    """Parse ``<taxid>_<sanitized_binomial><ext>`` → (taxid, 'Genus species').

    Returns None when the name has no leading numeric taxid + binomial.
    """
    stem = name
    for ext in _EXTS:
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    else:
        stem = Path(name).stem
    m = _LEAF_RE.match(stem)
    if not m:
        return None
    binomial = m.group(2).replace("_", " ").strip()
    if not binomial:
        return None
    return int(m.group(1)), binomial


def enumerate_proteomes(proteomes_dir: Path, pattern: str = "*.fa") -> list[ManifestRow]:
    """Enumerate proteome files → manifest rows (absolute paths)."""
    rows: list[ManifestRow] = []
    for p in sorted(proteomes_dir.glob(pattern)):
        if not p.is_file():
            continue
        parsed = parse_leaf_filename(p.name)
        if parsed is None:
            print(f"[build_full_scan_manifest] SKIP (no taxid_binomial): {p.name}",
                  file=sys.stderr)
            continue
        taxid, binomial = parsed
        rows.append(ManifestRow(taxid, binomial, p.resolve()))
    return rows


def write_manifest_tsv(rows: list[ManifestRow], out_path: Path) -> None:
    """Write taxid/binomial/proteome_path, sorted by taxid, LF line endings."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        fh.write("taxid\tbinomial\tproteome_path\n")
        for r in sorted(rows, key=lambda x: x.taxid):
            fh.write(f"{r.taxid}\t{r.binomial}\t{r.proteome_path}\n")


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--proteomes-dir", type=Path, required=True,
                   help="Directory of <taxid>_<binomial>.fa proteomes "
                        "(consolidation input dir).")
    p.add_argument("--out", type=Path, required=True, help="Output manifest TSV.")
    p.add_argument("--pattern", default="*.fa",
                   help="Glob for proteome files (default *.fa).")
    p.add_argument("--force", action="store_true",
                   help="Overwrite --out if it already exists.")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    if args.out.exists() and not args.force:
        print(f"[build_full_scan_manifest] {args.out} exists; skipping (use --force).",
              file=sys.stderr)
        return 0

    rows = enumerate_proteomes(args.proteomes_dir, pattern=args.pattern)
    write_manifest_tsv(rows, args.out)
    print(f"[build_full_scan_manifest] {len(rows)} proteomes → {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
