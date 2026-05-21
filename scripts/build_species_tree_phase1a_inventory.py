"""Phase 1a of species-tree pipeline: proteome inventory (no downloads).

Bead -c5d (under -dnk epic). For each of 239 reference taxa under
references/nath_et_al/one_to_one_ortholog/<clade>/, query NCBI Datasets
CLI for the best annotated assembly per species. Output a download
manifest TSV so we can see expected coverage before committing to the
bulk-download step (Phase 1b).

Selection priority:
  RefSeq + annotated  >  GenBank + annotated  >  none (drop)

Runs on Unity (sbatch on cpu/cpu-preempt, ~1h modest resources).
Downloads + intermediates land in
  /scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05/species_tree_data/
The final manifest TSV is small and gets scp'd to the repo at
  references/species_tree/proteome_manifest.tsv
for git-tracking.

See docs/plans/2026-05-21-species-tree-design.md (gitignored) for the
7-phase design context.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Union


# ----------------------------------------------------------------------
# Filename parsing
# ----------------------------------------------------------------------

_FILENAME_RE = re.compile(r"^(\d+)_(.+)$")


def parse_reference_filename(filename: Union[str, Path]) -> tuple[int, str]:
    """Parse `<taxid>_<binomial>.faa` filename into (taxid, binomial).

    Accepts a bare filename or any path; only the basename is consulted.
    Binomial is returned with underscores in the filename converted to
    spaces (so `Aplysia_californica` -> `Aplysia californica`).

    Raises ValueError on:
      - missing .faa extension
      - missing or non-integer taxid prefix
    """
    name = Path(str(filename)).name
    if not name.endswith(".faa"):
        raise ValueError(f"missing .faa extension: {name}")
    stem = name[: -len(".faa")]
    m = _FILENAME_RE.match(stem)
    if not m:
        raise ValueError(f"missing integer taxid prefix: {name}")
    taxid = int(m.group(1))
    binomial = m.group(2).replace("_", " ")
    return taxid, binomial


# ----------------------------------------------------------------------
# Reference taxa collection
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class ReferenceTaxon:
    taxid: int
    binomial: str
    clade: str


def collect_reference_taxa(refs_root: Union[str, Path]) -> list[ReferenceTaxon]:
    """Scan `<refs_root>/<clade>/*.faa` files and return ReferenceTaxon
    list sorted by taxid (deterministic ordering for downstream use).

    Filenames that don't match the `<taxid>_<binomial>.faa` schema are
    skipped with a stderr warning rather than crashing — real refs dirs
    can accumulate stray files (README, fixtures, etc.).

    Raises FileNotFoundError if `refs_root` doesn't exist.
    """
    root = Path(refs_root)
    if not root.exists():
        raise FileNotFoundError(f"refs_root does not exist: {root}")

    taxa: list[ReferenceTaxon] = []
    for clade_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        clade = clade_dir.name
        for faa in sorted(clade_dir.glob("*.faa")):
            try:
                taxid, binomial = parse_reference_filename(faa.name)
            except ValueError as e:
                print(f"warning: skipping {faa.name}: {e}", file=sys.stderr)
                continue
            taxa.append(ReferenceTaxon(taxid=taxid, binomial=binomial, clade=clade))

    taxa.sort(key=lambda t: t.taxid)
    return taxa


# ----------------------------------------------------------------------
# NCBI Datasets summary record -> best annotated assembly
# ----------------------------------------------------------------------

_ASSEMBLY_LEVEL_RANK = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1,
    "": 0,
}


@dataclass(frozen=True)
class AssemblyChoice:
    """Best annotated assembly for one taxid, distilled from a list
    of `datasets summary genome taxon` records.
    """
    accession: str
    source: str               # "RefSeq" | "GenBank"
    assembly_level: str       # "Chromosome" | "Scaffold" | "Contig" | "" (missing)
    annotation_status: str    # "Current" | "" (missing)
    est_protein_count: int    # 0 when unknown
    submission_date: str      # "" when unknown


def _is_annotated(record: dict) -> bool:
    """A datasets record counts as annotated iff it carries an
    annotation_info object. NCBI omits the field entirely when no
    annotation exists.
    """
    return isinstance(record.get("annotation_info"), dict)


def _source_label(record: dict) -> str:
    src = record.get("source_database", "")
    if "REFSEQ" in src.upper():
        return "RefSeq"
    if "GENBANK" in src.upper():
        return "GenBank"
    return src or "Unknown"


def pick_best_assembly(records: list[dict]) -> AssemblyChoice | None:
    """Choose the best annotated assembly from a list of NCBI Datasets
    summary records (the per-line JSON objects returned by
    `datasets summary genome taxon <taxid> --as-json-lines`).

    Priority:
      1) source == RefSeq AND annotated
      2) source == GenBank AND annotated
      (records that aren't annotated are excluded)

    Tie-breaks within priority class:
      a) higher assembly_level (Complete Genome > Chromosome > Scaffold > Contig)
      b) more recent submission_date (lex sort on ISO-ish YYYY-MM-DD)

    Returns None if no annotated record exists.
    """
    annotated = [r for r in records if _is_annotated(r)]
    if not annotated:
        return None

    def sort_key(r: dict) -> tuple:
        src = _source_label(r)
        # Priority: RefSeq (rank 2) > GenBank (rank 1) > Unknown (0)
        src_rank = {"RefSeq": 2, "GenBank": 1}.get(src, 0)
        level = (r.get("assembly_info") or {}).get("assembly_level", "") or ""
        level_rank = _ASSEMBLY_LEVEL_RANK.get(level, 0)
        date = (r.get("assembly_info") or {}).get("submission_date", "") or ""
        # Negate for descending sort on first three; date is descending too
        return (-src_rank, -level_rank, _negate_date(date))

    annotated.sort(key=sort_key)
    best = annotated[0]
    ann = best.get("annotation_info") or {}
    ann_stats = (ann.get("stats") or {}).get("gene_counts") or {}
    return AssemblyChoice(
        accession=best.get("accession", ""),
        source=_source_label(best),
        assembly_level=(best.get("assembly_info") or {}).get("assembly_level", "") or "",
        annotation_status=ann.get("status", "") or "",
        est_protein_count=int(ann_stats.get("protein_coding", 0) or 0),
        submission_date=(best.get("assembly_info") or {}).get("submission_date", "") or "",
    )


def _negate_date(date_str: str) -> str:
    """For descending date sort via ascending lex order. The trick:
    later dates should sort first, so we invert by subtracting from
    a sentinel string. For ISO YYYY-MM-DD format, simple negation:
    return a string that sorts in reverse of `date_str`.
    """
    # The cheap approach: prepend a NEGATIVE inverter. Since dates are
    # ISO-formatted, computing complement-per-digit is fiddly. Use a
    # simple inversion: pad to length 10, then map each digit to (9 - d).
    if not date_str:
        # Empty dates sort last in descending order. Return a tag that
        # sorts AFTER any real date when used as ascending key.
        return "9" * 10
    padded = date_str.ljust(10, "0")[:10]
    inverted = "".join(
        str(9 - int(c)) if c.isdigit() else c
        for c in padded
    )
    return inverted


# ----------------------------------------------------------------------
# Manifest entry + TSV writer
# ----------------------------------------------------------------------

MANIFEST_COLUMNS = (
    "taxid",
    "binomial",
    "clade",
    "source",
    "accession",
    "assembly_level",
    "annotation_status",
    "est_protein_count",
    "submission_date",
    "drop_reason",
)


@dataclass(frozen=True)
class ManifestEntry:
    """One row of the Phase 1a output manifest.

    For successful taxa: `choice` is the picked AssemblyChoice and
    `drop_reason` is "".
    For dropped taxa: `choice` is None and `drop_reason` explains
    why (e.g. "no_proteome_in_ncbi").
    """
    taxon: ReferenceTaxon
    choice: AssemblyChoice | None
    drop_reason: str


def write_manifest_tsv(path: Union[str, Path], entries: list[ManifestEntry]) -> None:
    """Write a TSV with MANIFEST_COLUMNS header, then one row per entry,
    sorted ascending by taxid. Dropped entries have empty assembly fields.
    """
    rows = sorted(entries, key=lambda e: e.taxon.taxid)
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for e in rows:
            t = e.taxon
            c = e.choice
            if c is None:
                row = [
                    str(t.taxid), t.binomial, t.clade,
                    "", "", "", "", "", "",
                    e.drop_reason,
                ]
            else:
                row = [
                    str(t.taxid), t.binomial, t.clade,
                    c.source, c.accession, c.assembly_level,
                    c.annotation_status, str(c.est_protein_count),
                    c.submission_date,
                    e.drop_reason,
                ]
            f.write("\t".join(row) + "\n")


# ----------------------------------------------------------------------
# NCBI Datasets CLI wrapper
# ----------------------------------------------------------------------

def query_datasets_for_taxon(
    taxid: int,
    datasets_bin: str = "datasets",
    timeout: int = 60,
) -> list[dict]:
    """Call `datasets summary genome taxon <taxid> --as-json-lines` and
    return the parsed list of JSON records (one per line).

    - Empty stdout -> empty list.
    - Non-zero returncode with a "no assemblies" stderr -> empty list
      (NCBI's way of saying "no hits"; not an error).
    - Other non-zero returncode -> RuntimeError with stderr text.
    - A malformed JSON line is skipped with a stderr warning; other
      lines in the same response are still returned.
    """
    cmd = [
        datasets_bin, "summary", "genome", "taxon", str(taxid),
        "--as-json-lines",
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout,
    )
    if result.returncode != 0:
        # Treat "no assemblies / no results" as empty, not error
        stderr_lower = (result.stderr or "").lower()
        if "no assemblies" in stderr_lower or "no results" in stderr_lower:
            return []
        raise RuntimeError(
            f"datasets failed (returncode={result.returncode}) "
            f"for taxid={taxid}: {result.stderr.strip()}"
        )

    # Some NCBI versions emit "no assemblies" via stderr while returning 0.
    stderr_lower = (result.stderr or "").lower()
    if not result.stdout.strip() and (
        "no assemblies" in stderr_lower or "no results" in stderr_lower
    ):
        return []

    records: list[dict] = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError as e:
            print(
                f"warning: skipping malformed JSON line for taxid={taxid}: {e}",
                file=sys.stderr,
            )
            continue
    return records


# ----------------------------------------------------------------------
# Orchestrator
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class InventorySummary:
    n_total: int
    n_with_proteome: int
    n_dropped: int
    n_refseq: int
    n_genbank: int
    n_query_errors: int


def build_inventory_manifest(
    refs_root: Union[str, Path],
    out_path: Union[str, Path],
    query_fn: Callable[[int], list[dict]] = query_datasets_for_taxon,
    progress_every: int = 10,
) -> InventorySummary:
    """Walk `refs_root`, query NCBI Datasets per taxon (via `query_fn`,
    injectable for testing), pick best assembly, write the manifest TSV
    to `out_path`. Returns an InventorySummary.

    Query failures (RuntimeError from query_fn) are recorded as
    drop_reason="query_error" so the whole batch doesn't abort on a
    transient NCBI hiccup — those taxa can be retried later.
    """
    taxa = collect_reference_taxa(refs_root)
    entries: list[ManifestEntry] = []
    n_refseq = 0
    n_genbank = 0
    n_query_errors = 0

    for i, taxon in enumerate(taxa, start=1):
        try:
            records = query_fn(taxon.taxid)
        except RuntimeError as e:
            print(
                f"  [{i}/{len(taxa)}] taxid={taxon.taxid} {taxon.binomial}: "
                f"query failed ({e}) — recording as drop_reason=query_error",
                file=sys.stderr,
            )
            entries.append(ManifestEntry(taxon=taxon, choice=None, drop_reason="query_error"))
            n_query_errors += 1
            continue

        choice = pick_best_assembly(records)
        if choice is None:
            entries.append(ManifestEntry(
                taxon=taxon, choice=None, drop_reason="no_proteome_in_ncbi",
            ))
        else:
            entries.append(ManifestEntry(taxon=taxon, choice=choice, drop_reason=""))
            if choice.source == "RefSeq":
                n_refseq += 1
            elif choice.source == "GenBank":
                n_genbank += 1

        if i % progress_every == 0 or i == len(taxa):
            print(
                f"  [{i}/{len(taxa)}] processed",
                file=sys.stderr,
            )

    write_manifest_tsv(out_path, entries)

    n_total = len(taxa)
    n_dropped = sum(1 for e in entries if e.drop_reason)
    n_with_proteome = n_total - n_dropped
    summary = InventorySummary(
        n_total=n_total,
        n_with_proteome=n_with_proteome,
        n_dropped=n_dropped,
        n_refseq=n_refseq,
        n_genbank=n_genbank,
        n_query_errors=n_query_errors,
    )
    print(
        f"\nInventory summary: {n_total} total | {n_with_proteome} with proteome "
        f"({n_refseq} RefSeq + {n_genbank} GenBank) | {n_dropped} dropped "
        f"({n_query_errors} query errors)",
        file=sys.stderr,
    )
    return summary


# ----------------------------------------------------------------------
# CLI entrypoint
# ----------------------------------------------------------------------

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Phase 1a: NCBI Datasets proteome inventory for 239 reference "
            "taxa. Writes a download manifest TSV (no big downloads)."
        )
    )
    p.add_argument(
        "--refs-root",
        type=Path,
        default=Path("references/nath_et_al/one_to_one_ortholog"),
        help="Directory holding <clade>/<taxid>_<binomial>.faa files",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("references/species_tree/proteome_manifest.tsv"),
        help="Output TSV path",
    )
    p.add_argument(
        "--datasets-bin",
        default="datasets",
        help="Path to the NCBI Datasets CLI binary (default: 'datasets' on PATH)",
    )
    p.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="Per-query subprocess timeout in seconds",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)

    def query_fn(taxid: int) -> list[dict]:
        return query_datasets_for_taxon(
            taxid, datasets_bin=args.datasets_bin, timeout=args.timeout,
        )

    summary = build_inventory_manifest(
        refs_root=args.refs_root,
        out_path=args.out,
        query_fn=query_fn,
    )
    return 0 if summary.n_total > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
