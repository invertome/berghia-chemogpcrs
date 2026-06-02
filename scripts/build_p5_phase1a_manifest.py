"""Build the P5 Phase 1a validation manifest.

Reads ``references/species_tree/proteome_manifest.tsv`` (Phase 1a, 91
annotated species), resolves proteome cache paths via the canonical
``sanitize_sample_name`` convention, and emits a three-column TSV:

    taxid  binomial  proteome_path

Used by the SLURM array wrapper
``sbatch_run_p5_phase1a_scan.sh`` to map array-task index → proteome path.

Path convention (must match consolidate_proteomes_for_genome_wide_og.py):
    <proteomes_dir>/<taxid>_<sanitized_binomial>.faa

Idempotent: skips output if it already exists unless --force.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Sanitization — MUST match build_braker4_samples_csv.sanitize_sample_name
# and consolidate_proteomes_for_genome_wide_og.sanitize_sample_name
# ---------------------------------------------------------------------------

_NON_ID_CHARS = re.compile(r"[^A-Za-z0-9_]+")
_MULTI_UNDERSCORE = re.compile(r"_+")


def sanitize_sample_name(name: str) -> str:
    """Collapse non-[A-Za-z0-9_] runs → '_', dedupe, strip leading/trailing."""
    s = _NON_ID_CHARS.sub("_", name.strip())
    s = _MULTI_UNDERSCORE.sub("_", s)
    return s.strip("_")


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ManifestRow:
    taxid: int
    binomial: str
    proteome_path: Path


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def read_phase1a_manifest(manifest_path: Path) -> list[ManifestRow]:
    """Parse the Phase 1a manifest TSV.

    Returns one ManifestRow per row that has a non-empty ``accession``
    and a parseable ``taxid``.  Does NOT check whether the proteome file
    exists — that is handled by :func:`build_manifest`.
    """
    rows: list[ManifestRow] = []
    with manifest_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            accession = (row.get("accession") or "").strip()
            if not accession:
                continue
            try:
                taxid = int(row["taxid"])
            except (KeyError, ValueError):
                continue
            binomial = (row.get("binomial") or "").strip()
            rows.append(ManifestRow(
                taxid=taxid,
                binomial=binomial,
                proteome_path=Path(),  # placeholder; resolved in build_manifest
            ))
    return rows


def resolve_proteome_path(
    taxid: int,
    binomial: str,
    proteomes_dir: Path,
) -> Path:
    """Return the expected .faa path for a species under ``proteomes_dir``.

    Uses the canonical ``<taxid>_<sanitized_binomial>.faa`` naming scheme
    that matches how Phase 1a proteomes are cached by the NCBI download
    pipeline and referenced by consolidate_proteomes_for_genome_wide_og.py.

    Always returns an absolute path.
    """
    leaf = f"{taxid}_{sanitize_sample_name(binomial)}"
    return (proteomes_dir / f"{leaf}.faa").resolve()


def build_manifest(
    manifest_path: Path,
    *,
    proteomes_dir: Path,
) -> list[ManifestRow]:
    """Read Phase 1a manifest, resolve paths, skip missing files.

    Emits a warning to stderr for each species whose proteome file is not
    found in ``proteomes_dir``.  Skipping is not fatal — validation simply
    won't run for that species.
    """
    parsed = read_phase1a_manifest(manifest_path)
    result: list[ManifestRow] = []
    for row in parsed:
        path = resolve_proteome_path(row.taxid, row.binomial, proteomes_dir)
        if not path.exists():
            print(
                f"[build_p5_phase1a_manifest] SKIP taxid={row.taxid} "
                f"({row.binomial}): proteome not found at {path}",
                file=sys.stderr,
            )
            continue
        result.append(ManifestRow(
            taxid=row.taxid,
            binomial=row.binomial,
            proteome_path=path,
        ))
    return result


def write_manifest_tsv(rows: list[ManifestRow], out_path: Path) -> None:
    """Write the three-column manifest TSV to ``out_path``.

    Always writes a header row; body may be empty if no proteomes were
    found. All proteome paths are written as absolute paths.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        # lineterminator="\n": csv.writer defaults to CRLF, which leaves a
        # trailing '\r' on col3 when the SLURM scan wrapper does `cut -f3`,
        # breaking `[ -f "<path>\r" ]` for every row.
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "binomial", "proteome_path"])
        for r in rows:
            writer.writerow([r.taxid, r.binomial, str(r.proteome_path.resolve())])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Build the P5 Phase 1a validation manifest TSV. "
            "Reads Phase 1a proteome_manifest.tsv, resolves proteome cache "
            "paths, emits taxid/binomial/proteome_path for each species that "
            "has a proteome file on disk."
        ),
    )
    p.add_argument(
        "--manifest", type=Path, required=True,
        help="Phase 1a manifest TSV (references/species_tree/proteome_manifest.tsv).",
    )
    p.add_argument(
        "--proteomes-dir", type=Path,
        default=Path("references/species_tree/cache/proteomes"),
        help=(
            "Directory containing <taxid>_<sanitized_binomial>.faa files. "
            "Default: references/species_tree/cache/proteomes"
        ),
    )
    p.add_argument(
        "--out", type=Path, required=True,
        help="Output manifest TSV path.",
    )
    p.add_argument(
        "--force", action="store_true",
        help="Re-run even if --out already exists.",
    )
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = _build_argparser().parse_args(argv)

    # Idempotency: skip if output exists and --force not given
    if args.out.exists() and not args.force:
        print(
            f"[build_p5_phase1a_manifest] Output already exists at {args.out}; "
            "skipping (use --force to overwrite).",
            file=sys.stderr,
        )
        return 0

    rows = build_manifest(args.manifest, proteomes_dir=args.proteomes_dir)
    write_manifest_tsv(rows, args.out)

    print(
        f"[build_p5_phase1a_manifest] {len(rows)} rows written to {args.out}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
