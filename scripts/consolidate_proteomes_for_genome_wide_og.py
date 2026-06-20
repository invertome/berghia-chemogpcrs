"""Gather all species-tree proteomes into one OrthoFinder input directory.

Reads Phase 1a (annotated), Phase 1d/1e (BRAKER4 de-novo), Phase 1g
(denovotranscript, optional), and the Berghia RefSeq proteome; writes each
to ${GENOME_WIDE_ORTHOFINDER_DIR}/input/<taxid>_<sanitized_binomial>.fa with
sequence headers rewritten to >{leaf_name}|{original_id}.

Filename and sequence-prefix format matches the species-tree leaf names by
construction (same taxid_sanitized_binomial scheme), so OrthoFinder species
column names equal tree leaves with no post-hoc remapping needed.

By design this consolidation should produce exactly one output file per
species-tree leaf. Any `missing_proteome` / `empty_fasta` / `duplicate_taxid`
entries in the consolidation_report.tsv are project-level gaps to be fixed
upstream — 03c will refuse to run CAFE5 until the report is clean.

Path conventions (see docs/plans/2026-05-28-orthology-557-expansion.md):
  Phase 1a  : <base_dir>/references/species_tree/cache/proteomes/<leaf>.faa
  Phase 1f  : <base_dir>/references/species_tree/cache/proteomes_braker4/<leaf>.faa
              OR <braker4_output_dir>/<leaf>/results/braker.aa (fallback)
  Phase 1g  : <base_dir>/references/species_tree/cache/proteomes_denovotranscript/<leaf>.faa
  Berghia   : supplied via --berghia-proteome CLI flag (caller resolves env var)
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Sanitization — MUST match build_braker4_samples_csv.sanitize_sample_name
# ---------------------------------------------------------------------------

_NON_ID = re.compile(r"[^A-Za-z0-9_]+")
_MULTI_UNDER = re.compile(r"_+")

BERGHIA_TAXID = 1287507


def sanitize_sample_name(name: str) -> str:
    """Collapse non-[A-Za-z0-9_] runs → '_', dedupe, strip leading/trailing '_'."""
    s = _NON_ID.sub("_", name.strip())
    s = _MULTI_UNDER.sub("_", s)
    return s.strip("_")


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ProteomeSource:
    taxid: int
    binomial: str
    fasta_path: Path
    phase: str  # "1a" | "1f" | "1g" | "berghia" | "recovery"


@dataclass(frozen=True)
class ConsolidationStatus:
    taxid: int
    binomial: str
    phase: str
    status: str  # "ok" | "missing_proteome" | "empty_fasta" | "duplicate_taxid"
    n_seqs: int = 0
    target: Optional[Path] = None
    message: str = ""


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def rewrite_proteome(src: Path, dst: Path, leaf_name: str) -> int:
    """Stream src FASTA → dst with header rewritten to >{leaf_name}|{original_id}.

    Only the first whitespace-delimited token of the original header becomes
    `original_id`; the rest of the description line is dropped.
    Returns the number of sequences written.
    """
    n = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                original_id = line[1:].split()[0]
                fout.write(f">{leaf_name}|{original_id}\n")
                n += 1
            else:
                fout.write(line)
    return n


def consolidate_one(
    src: ProteomeSource,
    *,
    out_dir: Path,
    force: bool = False,
) -> ConsolidationStatus:
    """Write one proteome to out_dir with leaf-name filename + header prefix.

    Returns a ConsolidationStatus describing the outcome.  Idempotent: if the
    target already exists and has ≥1 sequence header, skips unless force=True.
    """
    leaf = f"{src.taxid}_{sanitize_sample_name(src.binomial)}"
    target = out_dir / f"{leaf}.fa"

    if not src.fasta_path.exists():
        return ConsolidationStatus(
            src.taxid, src.binomial, src.phase, "missing_proteome",
            message=f"not found: {src.fasta_path}",
        )

    # Idempotency check: skip if target already has at least one header line
    if target.exists() and not force:
        if target.stat().st_size > 0:
            with target.open() as fh:
                has_header = any(ln.startswith(">") for ln in fh)
            if has_header:
                return ConsolidationStatus(
                    src.taxid, src.binomial, src.phase, "ok",
                    target=target,
                    message="skipped (already present)",
                )

    out_dir.mkdir(parents=True, exist_ok=True)
    n = rewrite_proteome(src.fasta_path, target, leaf)

    if n == 0:
        target.unlink(missing_ok=True)
        return ConsolidationStatus(
            src.taxid, src.binomial, src.phase, "empty_fasta",
            message="0 sequences in source FASTA",
        )

    return ConsolidationStatus(
        src.taxid, src.binomial, src.phase, "ok",
        n_seqs=n, target=target,
    )


# ---------------------------------------------------------------------------
# Manifest loading + path resolution
# ---------------------------------------------------------------------------

def _leaf_name(taxid: int, binomial: str) -> str:
    return f"{taxid}_{sanitize_sample_name(binomial)}"


def _read_manifest(path: Path, phase: str, base_dir: Path, braker4_output_dir: Path) -> list[ProteomeSource]:
    """Parse a species-tree manifest TSV and return ProteomeSources.

    Handles three column schemas:
      - Phase 1a (proteome_manifest.tsv):        clade column = 'clade'
      - Phase 1d (extension_inventory.tsv):      clade column = 'clade_name'
      - Phase 1e (genome_inventory_unannotated.tsv): clade column = 'clade'

    Skips rows with empty accession.
    """
    sources: list[ProteomeSource] = []
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            accession = (row.get("accession") or "").strip()
            if not accession:
                continue
            try:
                taxid = int(row["taxid"])
            except (KeyError, ValueError):
                continue
            binomial = (row.get("binomial") or "").strip()
            if not binomial:
                continue

            leaf = _leaf_name(taxid, binomial)
            fasta_path = _locate_proteome(taxid, binomial, leaf, phase, base_dir, braker4_output_dir)
            sources.append(ProteomeSource(
                taxid=taxid,
                binomial=binomial,
                fasta_path=fasta_path,
                phase=phase,
            ))
    return sources


def _locate_proteome(
    taxid: int,
    binomial: str,
    leaf: str,
    phase: str,
    base_dir: Path,
    braker4_output_dir: Path,
) -> Path:
    """Return the expected proteome path for a species, given its phase.

    For Phase 1a: canonical cache at references/species_tree/cache/proteomes/
    For Phase 1f (Phase 1d/1e BRAKER4): canonical post-process cache at
        references/species_tree/cache/proteomes_braker4/; fall back to the
        BRAKER4 raw output braker.aa file.
    For Phase 1g: references/species_tree/cache/proteomes_denovotranscript/
    """
    fname = f"{leaf}.faa"

    if phase == "1a":
        return base_dir / "references" / "species_tree" / "cache" / "proteomes" / fname

    if phase == "1f":
        canonical = base_dir / "references" / "species_tree" / "cache" / "proteomes_braker4" / fname
        if canonical.exists():
            return canonical
        # Fallback: raw BRAKER4 output directory
        raw = braker4_output_dir / leaf / "results" / "braker.aa"
        return canonical if not raw.exists() else raw

    if phase == "1g":
        return base_dir / "references" / "species_tree" / "cache" / "proteomes_denovotranscript" / fname

    # berghia / recovery — caller provides fasta_path directly
    return base_dir / fname


def load_sources(args) -> tuple[list[ProteomeSource], list[ConsolidationStatus]]:
    """Read all manifests + locate proteome paths per phase.

    Returns (sources_to_process, dropped_duplicates) where:
      - sources_to_process: list of ProteomeSource for unique taxids (first-seen wins)
      - dropped_duplicates: list of ConsolidationStatus with status='duplicate_taxid'

    Args has attributes:
        phase1a_manifest  (Path or None)
        phase1d_manifest  (Path or None)
        phase1e_manifest  (Path or None)
        phase1g_manifest  (Path or None)
        berghia_proteome  (Path or None)
        base_dir          (Path) — project root for relative path resolution
        braker4_output_dir (Path) — BRAKER4 per-species output root
    """
    base_dir = Path(args.base_dir)
    braker4_output_dir = Path(args.braker4_output_dir)
    seen_taxids: dict[int, str] = {}  # taxid -> phase of first-seen
    sources: list[ProteomeSource] = []
    duplicates: list[ConsolidationStatus] = []

    def _add(srcs: list[ProteomeSource]) -> None:
        for s in srcs:
            if s.taxid not in seen_taxids:
                seen_taxids[s.taxid] = s.phase
                sources.append(s)
            else:
                # Record as duplicate_taxid status
                duplicates.append(ConsolidationStatus(
                    taxid=s.taxid,
                    binomial=s.binomial,
                    phase=s.phase,
                    status="duplicate_taxid",
                    message=f"already seen in phase {seen_taxids[s.taxid]}",
                ))

    # Add Berghia first (canonical source) so it wins any duplicate-taxid conflicts
    if getattr(args, "berghia_proteome", None):
        berghia_path = Path(args.berghia_proteome)
        _add([ProteomeSource(
            taxid=BERGHIA_TAXID,
            binomial="Berghia stephanieae",
            fasta_path=berghia_path,
            phase="berghia",
        )])

    if getattr(args, "phase1a_manifest", None):
        _add(_read_manifest(
            Path(args.phase1a_manifest), "1a", base_dir, braker4_output_dir,
        ))

    # Phase 1d/1e: prefer legacy per-phase manifests when explicitly supplied;
    # fall back to the unified genome_inventory.tsv when neither is given.
    _phase1d = getattr(args, "phase1d_manifest", None)
    _phase1e = getattr(args, "phase1e_manifest", None)
    _unified = getattr(args, "manifest", None)

    if _phase1d or _phase1e:
        # Legacy explicit paths — honour them as before.
        if _phase1d:
            _add(_read_manifest(
                Path(_phase1d), "1f", base_dir, braker4_output_dir,
            ))
        if _phase1e:
            _add(_read_manifest(
                Path(_phase1e), "1f", base_dir, braker4_output_dir,
            ))
    elif _unified:
        # Unified manifest covers both 1d and 1e species (BRAKER4-annotated).
        _add(_read_manifest(
            Path(_unified), "1f", base_dir, braker4_output_dir,
        ))

    if getattr(args, "phase1g_manifest", None):
        _add(_read_manifest(
            Path(args.phase1g_manifest), "1g", base_dir, braker4_output_dir,
        ))

    return sources, duplicates


# ---------------------------------------------------------------------------
# Report writing
# ---------------------------------------------------------------------------

def write_consolidation_report(
    statuses: list[ConsolidationStatus],
    path: Path,
) -> None:
    """Write a TSV report with one row per species."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "taxid", "binomial", "phase", "status",
            "n_seqs", "target", "message",
        ])
        for s in statuses:
            w.writerow([
                s.taxid, s.binomial, s.phase, s.status,
                s.n_seqs,
                str(s.target) if s.target else "",
                s.message,
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Consolidate all species-tree proteomes into one OrthoFinder "
            "input dir with leaf-name-matching filenames + sequence prefixes. "
            "Feeds the genome-wide OrthoFinder run (see plan 2026-05-28)."
        ),
    )
    p.add_argument(
        "--out-dir", required=True, type=Path,
        help="OrthoFinder input dir (GENOME_WIDE_ORTHOFINDER_DIR/input).",
    )
    p.add_argument(
        "--base-dir", type=Path, default=Path("."),
        help="Project root for resolving relative proteome cache paths.",
    )
    p.add_argument(
        "--braker4-output-dir", type=Path,
        default=Path("species_tree_data/braker4_run/output"),
        help="BRAKER4 per-species output root (fallback path for raw braker.aa).",
    )
    p.add_argument(
        "--phase1a-manifest", type=Path, default=None,
        help="references/species_tree/proteome_manifest.tsv",
    )
    p.add_argument(
        "--manifest", type=Path,
        default=Path("references/species_tree/genome_inventory.tsv"),
        help=(
            "Unified genome inventory TSV (default: references/species_tree/"
            "genome_inventory.tsv). Covers species previously split across "
            "extension_inventory.tsv (Phase 1d) and genome_inventory_unannotated.tsv "
            "(Phase 1e). Ignored when --phase1d-manifest or --phase1e-manifest "
            "are supplied explicitly."
        ),
    )
    p.add_argument(
        "--phase1d-manifest", type=Path, default=None,
        help=(
            "Legacy Phase 1d manifest (extension_inventory.tsv). "
            "When supplied, overrides --manifest for Phase 1d species."
        ),
    )
    p.add_argument(
        "--phase1e-manifest", type=Path, default=None,
        help=(
            "Legacy Phase 1e manifest (genome_inventory_unannotated.tsv). "
            "When supplied, overrides --manifest for Phase 1e species."
        ),
    )
    p.add_argument(
        "--phase1g-manifest", type=Path, default=None,
        help="Optional Phase 1g (denovotranscript) manifest.",
    )
    p.add_argument(
        "--berghia-proteome", type=Path, default=None,
        help="Berghia RefSeq proteome FASTA (genomes/${BERGHIA_FILE_PREFIX}.proteins.fa).",
    )
    p.add_argument(
        "--force", action="store_true",
        help="Overwrite existing target files (default: skip if already present).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    sources, duplicates = load_sources(args)
    if not sources and not duplicates:
        print("[consolidate] WARNING: no proteome sources found — check manifests/args.",
              file=sys.stderr)

    statuses: list[ConsolidationStatus] = []
    for src in sources:
        status = consolidate_one(src, out_dir=out_dir, force=args.force)
        statuses.append(status)

    # Append any duplicate_taxid statuses
    statuses.extend(duplicates)

    report_path = out_dir.parent / "consolidation_report.tsv"
    write_consolidation_report(statuses, report_path)

    n_ok = sum(1 for s in statuses if s.status == "ok")
    n_total = len(statuses)
    n_bad = n_total - n_ok

    print(f"[consolidate] {n_ok}/{n_total} ok — report at {report_path}")

    # Hard-warn (not fail) so the caller sees the gap; 03c has the hard assertion.
    if n_bad:
        bad_statuses = [s for s in statuses if s.status != "ok"]
        print(
            f"[consolidate] WARN: {n_bad} species not consolidated cleanly:\n"
            + "\n".join(
                f"  {s.status} taxid={s.taxid} {s.binomial}: {s.message}"
                for s in bad_statuses
            ),
            file=sys.stderr,
        )
        print(
            "[consolidate] WARN: 03c CAFE5 will refuse to run until the "
            "consolidation report is clean (see consolidation_report.tsv).",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
