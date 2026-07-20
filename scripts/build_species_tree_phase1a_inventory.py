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

    Phase 1e extends with `contig_n50` + `total_length_bp` from
    assembly_stats; these help Phase 1f prioritize species with
    higher-quality assemblies for BRAKER3 annotation.
    """
    accession: str
    source: str               # "RefSeq" | "GenBank"
    assembly_level: str       # "Chromosome" | "Scaffold" | "Contig" | "" (missing)
    annotation_status: str    # "Current" | "" (missing == unannotated)
    est_protein_count: int    # 0 when unknown
    submission_date: str      # "" when unknown
    contig_n50: int = 0       # 0 when unknown
    total_length_bp: int = 0  # 0 when unknown


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


def _assembly_stats(record: dict) -> tuple[int, int]:
    """Return (contig_n50, total_length_bp) from a datasets record.
    Both default to 0 when missing.
    """
    stats = record.get("assembly_stats") or {}
    try:
        n50 = int(stats.get("contig_n50", 0) or 0)
    except (TypeError, ValueError):
        n50 = 0
    try:
        total = int(stats.get("total_sequence_length", 0) or 0)
    except (TypeError, ValueError):
        total = 0
    return n50, total


def pick_best_assembly(
    records: list[dict],
    require_annotation: bool = True,
) -> AssemblyChoice | None:
    """Choose the best assembly from a list of NCBI Datasets summary
    records (the per-line JSON objects returned by
    `datasets summary genome taxon <taxid> --as-json-lines`).

    Default (Phase 1a, require_annotation=True) priority:
      1) source == RefSeq AND annotated
      2) source == GenBank AND annotated
      (records that aren't annotated are excluded)

    Phase 1e mode (require_annotation=False) extends priority to:
      1) RefSeq AND annotated
      2) GenBank AND annotated
      3) GenBank AND unannotated  <-- BRAKER3 candidates
      (RefSeq+unannotated effectively doesn't exist on NCBI)

    Tie-breaks within priority class:
      a) higher assembly_level (Complete Genome > Chromosome > Scaffold > Contig)
      b) more recent release_date (NCBI's own field name; the column it
         lands in is still called submission_date for schema compatibility)

    Returns None when no candidate passes the priority filter.
    """
    candidates = []
    for r in records:
        if _is_annotated(r):
            candidates.append((r, True))   # (record, is_annotated)
        elif not require_annotation:
            candidates.append((r, False))
    if not candidates:
        return None

    def sort_key(item) -> tuple:
        r, annotated = item
        src = _source_label(r)
        src_rank = {"RefSeq": 2, "GenBank": 1}.get(src, 0)
        # Annotation rank: annotated (rank 1) beats unannotated (rank 0)
        # within the same source. Combined with source rank this gives:
        # RefSeq+annot=22, GenBank+annot=12, GenBank+unannot=10, etc.
        ann_rank = 1 if annotated else 0
        level = (r.get("assembly_info") or {}).get("assembly_level", "") or ""
        level_rank = _ASSEMBLY_LEVEL_RANK.get(level, 0)
        # NCBI datasets v2alpha exposes 'release_date'; there is no
        # 'submission_date' key. Reading the wrong name made this tie-break a
        # permanent no-op, so which of two equal-ranked assemblies won could
        # flip between re-queries.
        date = (r.get("assembly_info") or {}).get("release_date", "") or ""
        # Encode as a single tuple where lower sorts first; negate ranks
        # for descending. Use (src*10+ann) to keep "annotated source X"
        # beating "unannotated source X" before assembly level matters.
        combined_rank = src_rank * 10 + ann_rank
        return (-combined_rank, -level_rank, _negate_date(date))

    candidates.sort(key=sort_key)
    best, best_annotated = candidates[0]

    ann = best.get("annotation_info") or {} if best_annotated else {}
    ann_stats = (ann.get("stats") or {}).get("gene_counts") or {}
    n50, total_len = _assembly_stats(best)
    return AssemblyChoice(
        accession=best.get("accession", ""),
        source=_source_label(best),
        assembly_level=(best.get("assembly_info") or {}).get("assembly_level", "") or "",
        annotation_status=(ann.get("status", "") or "") if best_annotated else "",
        est_protein_count=int(ann_stats.get("protein_coding", 0) or 0),
        submission_date=(best.get("assembly_info") or {}).get("release_date", "") or "",
        contig_n50=n50,
        total_length_bp=total_len,
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


def write_manifest_tsv(
    path: Union[str, Path],
    entries: list[ManifestEntry],
    include_assembly_stats: bool = False,
) -> None:
    """Write a TSV with MANIFEST_COLUMNS header, then one row per entry,
    sorted ascending by taxid. Dropped entries have empty assembly fields.

    When `include_assembly_stats=True` (Phase 1e), two extra columns
    `contig_n50` and `total_length_bp` are appended.
    """
    rows = sorted(entries, key=lambda e: e.taxon.taxid)
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    columns = list(MANIFEST_COLUMNS)
    if include_assembly_stats:
        columns += ["contig_n50", "total_length_bp"]
    with out.open("w") as f:
        f.write("\t".join(columns) + "\n")
        for e in rows:
            t = e.taxon
            c = e.choice
            if c is None:
                row = [
                    str(t.taxid), t.binomial, t.clade,
                    "", "", "", "", "", "",
                    e.drop_reason,
                ]
                if include_assembly_stats:
                    row += ["", ""]
            else:
                row = [
                    str(t.taxid), t.binomial, t.clade,
                    c.source, c.accession, c.assembly_level,
                    c.annotation_status, str(c.est_protein_count),
                    c.submission_date,
                    e.drop_reason,
                ]
                if include_assembly_stats:
                    row += [str(c.contig_n50), str(c.total_length_bp)]
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
    # NCBI Datasets uses non-zero rc for several "valid taxid, but no
    # genome / annotation / proteome" outcomes — these are NOT errors,
    # they just mean "no proteome for this taxon". Recognize them.
    NO_HIT_PHRASES = (
        "no assemblies",
        "no results",
        "no genome data",            # "...is valid for 'X', but no genome data..."
        "no genomes",                # plural form some versions emit
        "no annotation",             # rarer; assemblies present but no proteome
        "no proteins",
        "did not match any genomes",
    )

    stderr_lower = (result.stderr or "").lower()

    if result.returncode != 0:
        if any(p in stderr_lower for p in NO_HIT_PHRASES):
            return []
        raise RuntimeError(
            f"datasets failed (returncode={result.returncode}) "
            f"for taxid={taxid}: {result.stderr.strip()}"
        )

    # Some NCBI versions emit a no-hit phrase via stderr while returning 0.
    if not result.stdout.strip() and any(p in stderr_lower for p in NO_HIT_PHRASES):
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


def read_taxids_to_skip(
    manifest_path: Union[str, Path],
    mode: str = "dropped_only",
) -> set[int]:
    """Read a prior manifest TSV and return the set of taxids to SKIP
    when re-running the inventory.

    `mode="dropped_only"` (default for Phase 1e): skip taxa that had
    a successful Phase 1a entry (drop_reason == ""), so we only
    re-query the taxa Phase 1a dropped.

    `mode="successful_only"`: skip taxa that were dropped — useful if
    you want to re-run only the previously-successful subset (e.g.,
    to refresh assembly stats).
    """
    path = Path(manifest_path)
    if not path.exists():
        raise FileNotFoundError(f"manifest not found: {path}")
    skip: set[int] = set()
    import csv as _csv
    with path.open() as f:
        reader = _csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                taxid = int(row["taxid"])
            except (KeyError, ValueError):
                continue
            drop_reason = (row.get("drop_reason") or "").strip()
            if mode == "dropped_only" and not drop_reason:
                skip.add(taxid)
            elif mode == "successful_only" and drop_reason:
                skip.add(taxid)
    return skip


def build_inventory_manifest(
    refs_root: Union[str, Path],
    out_path: Union[str, Path],
    query_fn: Callable[[int], list[dict]] = query_datasets_for_taxon,
    progress_every: int = 10,
    require_annotation: bool = True,
    taxids_to_skip: set[int] | None = None,
    include_assembly_stats: bool = False,
) -> InventorySummary:
    """Walk `refs_root`, query NCBI Datasets per taxon (via `query_fn`,
    injectable for testing), pick best assembly, write the manifest TSV
    to `out_path`. Returns an InventorySummary.

    Phase 1e args:
    - `require_annotation=False`: also accept unannotated GenBank
      assemblies as fallback candidates (3rd-priority).
    - `taxids_to_skip`: set of taxids to omit entirely from this run
      (e.g., skip Phase 1a's successful 91 species when re-querying
      only the dropped 148).
    - `include_assembly_stats=True`: emit contig_n50 + total_length_bp
      columns.

    Query failures (RuntimeError from query_fn) are recorded as
    drop_reason="query_error" so the whole batch doesn't abort on a
    transient NCBI hiccup — those taxa can be retried later.
    """
    taxa_all = collect_reference_taxa(refs_root)
    if taxids_to_skip:
        taxa = [t for t in taxa_all if t.taxid not in taxids_to_skip]
        skipped_n = len(taxa_all) - len(taxa)
        if skipped_n:
            print(
                f"Skipping {skipped_n} taxa per --limit-taxids-from filter",
                file=sys.stderr,
            )
    else:
        taxa = taxa_all
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

        choice = pick_best_assembly(records, require_annotation=require_annotation)
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

    write_manifest_tsv(out_path, entries, include_assembly_stats=include_assembly_stats)

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
    p.add_argument(
        "--allow-unannotated",
        action="store_true",
        help=(
            "Phase 1e mode: also accept unannotated GenBank assemblies "
            "as fallback candidates (priority: RefSeq+annot > GenBank+annot > "
            "GenBank+unannotated). When set, output TSV gains contig_n50 + "
            "total_length_bp columns for Phase 1f BRAKER3 prioritization."
        ),
    )
    p.add_argument(
        "--limit-taxids-from",
        type=Path,
        default=None,
        help=(
            "Path to a prior manifest TSV. When provided, restrict this run "
            "to taxa NOT yet covered there. Default mode 'dropped_only' "
            "skips successful Phase 1a entries."
        ),
    )
    p.add_argument(
        "--limit-mode",
        choices=("dropped_only", "successful_only"),
        default="dropped_only",
        help="Which subset to keep when --limit-taxids-from is given",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)

    def query_fn(taxid: int) -> list[dict]:
        return query_datasets_for_taxon(
            taxid, datasets_bin=args.datasets_bin, timeout=args.timeout,
        )

    skip_set: set[int] | None = None
    if args.limit_taxids_from:
        skip_set = read_taxids_to_skip(args.limit_taxids_from, mode=args.limit_mode)
        print(
            f"--limit-taxids-from {args.limit_taxids_from} (mode={args.limit_mode}): "
            f"{len(skip_set)} taxa will be skipped",
            file=sys.stderr,
        )

    summary = build_inventory_manifest(
        refs_root=args.refs_root,
        out_path=args.out,
        query_fn=query_fn,
        require_annotation=not args.allow_unannotated,
        taxids_to_skip=skip_set,
        include_assembly_stats=args.allow_unannotated,
    )
    return 0 if summary.n_total > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
