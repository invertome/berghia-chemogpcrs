"""Build the P5 Phase 1a validation manifest.

Reads ``references/species_tree/proteome_manifest.tsv`` (Phase 1a, 91
annotated species), resolves proteome cache paths via the canonical
``sanitize_sample_name`` convention, and emits a three-column TSV:

    taxid  binomial  proteome_path

Used by the SLURM array wrapper
``sbatch_run_p5_phase1a_scan.sh`` to map array-task index → proteome path.

Path convention: ``<proteomes_dir>/<taxid>_<sanitized_binomial><suffix>``.
The directories and suffixes are NOT re-declared here -- they are imported
from consolidate_proteomes_for_genome_wide_og, which is the single source of
truth for where each phase's proteomes live. Two copies of a path table
drift, and the drift is silent.

This module used to hardcode ``references/species_tree/cache/proteomes`` as
the default, a directory no producer writes and which does not exist on the
cluster at all (measured 2026-07-20). Phase 1a proteomes are written by
download_species_tree_phase1a.py into ``species_tree_data/ncbi_proteomes``
as ``<leaf>.faa`` (133 files). The wrong default did not crash: every
species reported individually "not found" and the run exited 0 with a
header-only manifest, so the downstream scan array had nothing to do and
also succeeded. main() now returns 2 when a non-empty manifest resolves
NOTHING -- see the guard there.

Idempotent: skips output if it already exists unless --force.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

sys.path.insert(0, str(Path(__file__).resolve().parent))

# Single source of truth for phase → directory and the probed protein-FASTA
# suffixes. _probe_proteome is imported (rather than reimplemented) for the
# same reason: the "which extension does the producer actually write" question
# has exactly one correct answer and it is already encoded there.
from consolidate_proteomes_for_genome_wide_og import (  # noqa: E402
    PHASE_PROTEOME_DIRS,
    PROTEOME_SUFFIXES,
    _probe_proteome,
)

PHASE = "1a"

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


def phase1a_proteome_dirs(
    base_dir: Path = Path("."),
    override: Optional[Path] = None,
) -> list[Path]:
    """Directories to probe for Phase 1a proteomes, canonical first.

    ``override`` (the --proteomes-dir flag) wins outright so an operator can
    point at an arbitrary location. Otherwise the phase → directory table is
    consulted, resolved against ``base_dir``.
    """
    if override is not None:
        return [Path(override)]
    return [Path(base_dir) / d for d in PHASE_PROTEOME_DIRS[PHASE]]


def resolve_proteome_path(
    taxid: int,
    binomial: str,
    proteomes_dir: Union[Path, Sequence[Path]],
) -> Path:
    """Return the proteome path for a species, probing known suffixes.

    ``proteomes_dir`` is a single directory or an ordered sequence of them
    (canonical producer location first). Each is probed over
    PROTEOME_SUFFIXES, matching how the consolidator locates the very same
    files -- the producer writes ``<leaf>.faa``, but probing means a future
    rename does not silently zero this script out.

    When nothing is found, returns the CANONICAL directory's ``<leaf>.faa``
    so the caller's "not found" message points an operator at the location
    they should go look in -- never at a directory no producer writes.

    Always returns an absolute path.
    """
    leaf = f"{taxid}_{sanitize_sample_name(binomial)}"
    dirs: list[Path] = (
        [Path(proteomes_dir)]
        if isinstance(proteomes_dir, (str, Path))
        else [Path(d) for d in proteomes_dir]
    )
    for d in dirs:
        found = _probe_proteome(d, leaf)
        if found is not None:
            return found.resolve()
    return (dirs[0] / f"{leaf}.faa").resolve()


def build_manifest(
    manifest_path: Path,
    *,
    proteomes_dir: Optional[Union[Path, Sequence[Path]]] = None,
    base_dir: Path = Path("."),
) -> list[ManifestRow]:
    """Read Phase 1a manifest, resolve paths, skip missing files.

    ``proteomes_dir`` may be a single directory, an ordered sequence of
    directories, or None — in which case the phase → directory table is used,
    resolved against ``base_dir``.

    Emits a warning to stderr for each species whose proteome file is not
    found.  Skipping an individual species is not fatal — validation simply
    won't run for it. A run in which NOTHING resolves is a different animal
    and is caught by main()'s guard, not here.
    """
    dirs = (
        phase1a_proteome_dirs(base_dir)
        if proteomes_dir is None
        else ([Path(proteomes_dir)]
              if isinstance(proteomes_dir, (str, Path))
              else [Path(d) for d in proteomes_dir])
    )
    parsed = read_phase1a_manifest(manifest_path)
    result: list[ManifestRow] = []
    seen: dict[int, ManifestRow] = {}
    for row in parsed:
        path = resolve_proteome_path(row.taxid, row.binomial, dirs)
        if not path.exists():
            print(
                f"[build_p5_phase1a_manifest] SKIP taxid={row.taxid} "
                f"({row.binomial}): proteome not found at {path}",
                file=sys.stderr,
            )
            continue
        # NCBI organism-name synonyms put one taxid on two manifest rows
        # (e.g. 6573 Mizuhopecten / Patinopecten yessoensis — same assembly,
        # byte-identical proteome). Keep the first row whose file exists so the
        # taxon isn't double-counted; name the dropped synonym loudly.
        kept = seen.get(row.taxid)
        if kept is not None:
            print(
                f"[build_p5_phase1a_manifest] DUPLICATE taxid {row.taxid}: "
                f"keeping '{kept.binomial}', dropping '{row.binomial}'",
                file=sys.stderr,
            )
            continue
        new_row = ManifestRow(
            taxid=row.taxid,
            binomial=row.binomial,
            proteome_path=path,
        )
        seen[row.taxid] = new_row
        result.append(new_row)
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
        "--proteomes-dir", type=Path, default=None,
        help=(
            "Override the directory containing "
            "<taxid>_<sanitized_binomial>.faa files. Default: the phase 1a "
            "entries of PHASE_PROTEOME_DIRS resolved against --base-dir "
            f"({', '.join(PHASE_PROTEOME_DIRS[PHASE])})."
        ),
    )
    p.add_argument(
        "--base-dir", type=Path, default=Path("."),
        help="Project root for resolving the phase 1a proteome directories.",
    )
    p.add_argument(
        "--out", type=Path, required=True,
        help="Output manifest TSV path.",
    )
    p.add_argument(
        "--force", action="store_true",
        help="Re-run even if --out already exists.",
    )
    p.add_argument(
        "--allow-empty", action="store_true",
        help=(
            "Permit a run in which the manifest lists species but NONE of "
            "their proteomes resolve (e.g. the Phase 1a download has not run "
            "yet). Without this the run fails loudly rather than emitting a "
            "header-only manifest that silently gives the scan array no work."
        ),
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

    dirs = phase1a_proteome_dirs(args.base_dir, args.proteomes_dir)
    n_parsed = len(read_phase1a_manifest(args.manifest))
    rows = build_manifest(args.manifest, proteomes_dir=dirs)
    write_manifest_tsv(rows, args.out)

    print(
        f"[build_p5_phase1a_manifest] {len(rows)}/{n_parsed} rows written "
        f"to {args.out}",
        file=sys.stderr,
    )

    # Zero-resolution guard. A manifest that lists species but resolves NONE
    # of them is a path/convention failure, not N independently missing
    # species — and it is invisible downstream, because a header-only manifest
    # just gives the SLURM scan array zero tasks, which then also succeeds.
    # Partial gaps stay a warning: Phase 1a downloads land incrementally, and
    # an individually missing species is a real (reported) measurement.
    if n_parsed and not rows and not args.allow_empty:
        print(
            f"[build_p5_phase1a_manifest] ERROR: resolved 0/{n_parsed} "
            f"proteomes. That is a path/convention failure, not {n_parsed} "
            f"independently missing species — check --base-dir, the manifest, "
            f"and that download_species_tree_phase1a.py really writes into "
            f"[{', '.join(str(d) for d in dirs)}] "
            f"(probed suffixes: {', '.join(PROTEOME_SUFFIXES)}). "
            f"Pass --allow-empty if the Phase 1a download genuinely has not "
            f"run yet.",
            file=sys.stderr,
        )
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
