"""Generate a BRAKER4-format samples.csv for Phase 1f (bead -p49).

Reads `references/species_tree/genome_inventory_unannotated.tsv` (134
unannotated-genome species from Phase 1e) + a cache dir of pre-downloaded
genome FASTAs (Phase 1f genome downloader output) + an OrthoDB metazoa
protein DB; emits the BRAKER4 samples.csv for protein-only EP mode.

BRAKER4 schema reference: https://github.com/Gaius-Augustus/BRAKER4 README.
The 14 columns must appear in the exact order below; for EP mode only
four are populated:

    sample_name      = <taxid>_<binomial> sanitized to [a-zA-Z0-9_]
    genome           = absolute path to <cache>/<taxid>_<binomial>.fasta
    protein_fasta    = absolute path to OrthoDB metazoa FASTA
    busco_lineage    = metazoa_odb12  (matches BRAKER4 default for metazoans)

`genome_masked` left empty triggers BRAKER4's built-in RepeatMasker step
(adds 4-12h per species; per user "do things properly" preference + the
Nath 2025 methodology this is the right default).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union


BRAKER4_COLUMNS = (
    "sample_name", "genome", "genome_masked", "protein_fasta",
    "bam_files", "fastq_r1", "fastq_r2", "sra_ids",
    "varus_genus", "varus_species", "isoseq_bam", "isoseq_fastq",
    "busco_lineage", "reference_gtf",
)


@dataclass(frozen=True)
class SpeciesTarget:
    """One row from the Phase 1e unannotated-genome inventory."""
    taxid: int
    binomial: str
    clade: str
    accession: str


_NON_ID_CHARS = re.compile(r"[^A-Za-z0-9_]+")
_MULTI_UNDERSCORE = re.compile(r"_+")


def sanitize_sample_name(name: str) -> str:
    """Normalize to [a-zA-Z0-9_] only — collapse runs of underscores,
    strip leading/trailing underscores. BRAKER4 uses sample_name as a
    directory + filename component so any character outside that set
    would break filesystem layout.
    """
    s = _NON_ID_CHARS.sub("_", name.strip())
    s = _MULTI_UNDERSCORE.sub("_", s)
    return s.strip("_")


def read_targets(manifest_path: Union[str, Path]) -> list[SpeciesTarget]:
    """Read Phase 1e manifest, return rows with non-empty accession.

    Only the 4 standard columns are read; Phase 1e's extra contig_n50 +
    total_length_bp are ignored.
    """
    targets: list[SpeciesTarget] = []
    with Path(manifest_path).open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            accession = (row.get("accession") or "").strip()
            if not accession:
                continue
            try:
                taxid = int(row["taxid"])
            except (KeyError, ValueError):
                continue
            targets.append(SpeciesTarget(
                taxid=taxid,
                binomial=(row.get("binomial") or "").strip(),
                clade=(row.get("clade") or "").strip(),
                accession=accession,
            ))
    targets.sort(key=lambda t: t.taxid)
    return targets


def build_samples_row(
    target: SpeciesTarget,
    genome_cache: Union[str, Path],
    protein_db: Union[str, Path],
    *,
    busco_lineage: str = "metazoa_odb12",
) -> Optional[dict[str, str]]:
    """Build one BRAKER4 samples.csv row for `target`.

    Returns None when the expected genome FASTA isn't present in
    `genome_cache` — caller logs a warning and skips that species.
    """
    cache = Path(genome_cache)
    sample_name = sanitize_sample_name(f"{target.taxid}_{target.binomial}")
    genome_fa = cache / f"{sample_name}.fasta"
    if not genome_fa.exists():
        return None
    protein_fa = Path(protein_db)
    row = {col: "" for col in BRAKER4_COLUMNS}
    row["sample_name"] = sample_name
    row["genome"] = str(genome_fa.resolve())
    row["protein_fasta"] = str(protein_fa.resolve())
    row["busco_lineage"] = busco_lineage
    return row


def write_samples_csv(
    out_path: Union[str, Path],
    rows: list[dict[str, str]],
) -> None:
    """Write a CSV with BRAKER4_COLUMNS header + one row per dict.

    Uses standard csv.DictWriter so any embedded commas in paths are
    quoted properly — though our paths should never contain commas.
    """
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(BRAKER4_COLUMNS))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Generate BRAKER4 samples.csv for Phase 1f EP-mode runs over "
            "the unannotated-genome species. Reads the Phase 1e manifest "
            "+ a cache dir of pre-downloaded genome FASTAs + an OrthoDB "
            "protein DB; emits one row per species that has a cached genome."
        ),
    )
    p.add_argument(
        "--manifest", type=Path, required=True,
        help="Phase 1e manifest TSV.",
    )
    p.add_argument(
        "--genome-cache", type=Path, required=True,
        help="Dir containing <taxid>_<binomial>.fasta files (Phase 1f download output).",
    )
    p.add_argument(
        "--protein-db", type=Path, required=True,
        help="OrthoDB metazoa protein FASTA path (single file).",
    )
    p.add_argument(
        "--busco-lineage", default="metazoa_odb12",
        help="BUSCO lineage for QC (default: metazoa_odb12).",
    )
    p.add_argument(
        "--out", type=Path, required=True,
        help="Output samples.csv path.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)

    if not args.protein_db.exists():
        print(
            f"[build_braker4_samples_csv] WARNING: protein DB not present at "
            f"{args.protein_db}. The samples.csv will still reference this path "
            f"so BRAKER4 will fail at runtime if it's not staged. Proceeding.",
            file=sys.stderr,
        )

    targets = read_targets(args.manifest)
    rows: list[dict[str, str]] = []
    n_missing = 0
    for target in targets:
        row = build_samples_row(
            target, args.genome_cache, args.protein_db,
            busco_lineage=args.busco_lineage,
        )
        if row is None:
            n_missing += 1
            print(
                f"  skip taxid={target.taxid} {target.binomial}: "
                f"genome not in cache",
                file=sys.stderr,
            )
            continue
        rows.append(row)

    write_samples_csv(args.out, rows)

    print(
        f"[build_braker4_samples_csv] {len(rows)} rows written to {args.out} "
        f"({n_missing} skipped for missing genome)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
