"""Post-process BRAKER4 outputs into the canonical Phase 1f cache layout
(bead -p49).

Walks BRAKER4's `output/{sample_name}/results/` per species in samples.csv,
gunzips `braker.aa.gz` + `braker.codingseq.gz`, and writes them at
canonical paths:

    references/species_tree/cache/proteomes_braker4/<sample_name>.aa.fna
    references/species_tree/cache/proteomes_braker4/<sample_name>.cds.fna

A per-species report TSV records status:
  ok        -> both .aa.fna + .cds.fna produced
  ok_no_cds -> .aa.fna only (BRAKER4 ran but didn't ship codingseq;
               shouldn't happen but defensive)
  failed    -> braker.aa.gz missing (BRAKER4 didn't produce output)
  skipped   -> canonical files already present (idempotent re-runs)

The script exits 0 if everything is ok/skipped; exits 2 if any species
failed (so a wrapping sbatch can flag the run for attention).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union


@dataclass(frozen=True)
class PostprocessResult:
    sample_name: str
    status: str
    aa_path: Optional[Path]
    cds_path: Optional[Path]
    error: str = ""
    n_proteins: int = 0
    n_cds: int = 0


REPORT_COLUMNS = ("sample_name", "status", "n_proteins", "n_cds", "error")


def read_sample_names(samples_csv: Union[str, Path]) -> list[str]:
    """Read sample_name column from BRAKER4 samples.csv. Empty / whitespace
    sample_names are dropped (defensive — shouldn't occur in practice).
    """
    out: list[str] = []
    with Path(samples_csv).open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = (row.get("sample_name") or "").strip()
            if name:
                out.append(name)
    return out


def _count_fasta_records(path: Path) -> int:
    if not path or not path.exists():
        return 0
    n = 0
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def _gunzip_to(src_gz: Path, dest: Path) -> None:
    """Gunzip src_gz to dest. dest's parent is created if missing."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src_gz, "rb") as fin, dest.open("wb") as fout:
        shutil.copyfileobj(fin, fout)


def postprocess_one(
    sample_name: str,
    braker4_output: Union[str, Path],
    cache_dir: Union[str, Path],
    *,
    force: bool = False,
) -> PostprocessResult:
    """Post-process one species' BRAKER4 output.

    - If `force=False` and both canonical files already exist with content,
      return status='skipped'.
    - Else: locate braker.aa.gz (required) + braker.codingseq.gz (optional)
      under `<braker4_output>/<sample_name>/results/`, gunzip both into
      canonical names under `<cache_dir>/`.
    """
    cache = Path(cache_dir)
    output = Path(braker4_output)
    cache.mkdir(parents=True, exist_ok=True)
    aa_dest = cache / f"{sample_name}.aa.fna"
    cds_dest = cache / f"{sample_name}.cds.fna"

    if not force:
        aa_present = aa_dest.exists() and aa_dest.stat().st_size > 0
        cds_present = cds_dest.exists() and cds_dest.stat().st_size > 0
        if aa_present and cds_present:
            return PostprocessResult(
                sample_name=sample_name, status="skipped",
                aa_path=aa_dest, cds_path=cds_dest,
                n_proteins=_count_fasta_records(aa_dest),
                n_cds=_count_fasta_records(cds_dest),
            )

    results_dir = output / sample_name / "results"
    aa_gz = results_dir / "braker.aa.gz"
    cds_gz = results_dir / "braker.codingseq.gz"

    if not aa_gz.exists() or aa_gz.stat().st_size == 0:
        return PostprocessResult(
            sample_name=sample_name, status="failed",
            aa_path=None, cds_path=None,
            error=f"braker.aa.gz missing or empty at {aa_gz}",
        )

    _gunzip_to(aa_gz, aa_dest)
    n_proteins = _count_fasta_records(aa_dest)

    if cds_gz.exists() and cds_gz.stat().st_size > 0:
        _gunzip_to(cds_gz, cds_dest)
        return PostprocessResult(
            sample_name=sample_name, status="ok",
            aa_path=aa_dest, cds_path=cds_dest,
            n_proteins=n_proteins,
            n_cds=_count_fasta_records(cds_dest),
        )

    return PostprocessResult(
        sample_name=sample_name, status="ok_no_cds",
        aa_path=aa_dest, cds_path=None,
        n_proteins=n_proteins, n_cds=0,
    )


def postprocess_all(
    sample_names: list[str],
    braker4_output: Union[str, Path],
    cache_dir: Union[str, Path],
    *,
    force: bool = False,
    progress_every: int = 10,
) -> list[PostprocessResult]:
    """Loop over samples, return PostprocessResults in input order."""
    results: list[PostprocessResult] = []
    for i, name in enumerate(sample_names, start=1):
        result = postprocess_one(
            name, braker4_output, cache_dir, force=force,
        )
        results.append(result)
        if i % progress_every == 0 or i == len(sample_names):
            n_ok = sum(1 for r in results if r.status == "ok")
            n_no_cds = sum(1 for r in results if r.status == "ok_no_cds")
            n_skip = sum(1 for r in results if r.status == "skipped")
            n_fail = sum(1 for r in results if r.status == "failed")
            print(
                f"  [{i}/{len(sample_names)}] ok={n_ok} ok_no_cds={n_no_cds} "
                f"skipped={n_skip} failed={n_fail}",
                file=sys.stderr,
            )
    return results


def write_report_tsv(
    report_path: Union[str, Path],
    results: list[PostprocessResult],
) -> None:
    path = Path(report_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(REPORT_COLUMNS) + "\n")
        for r in results:
            f.write("\t".join([
                r.sample_name,
                r.status,
                str(r.n_proteins),
                str(r.n_cds),
                r.error.replace("\t", " ").replace("\n", " "),
            ]) + "\n")


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Post-process BRAKER4 outputs: gunzip + rename per-species "
            "braker.aa.gz + braker.codingseq.gz into a canonical cache dir."
        ),
    )
    p.add_argument(
        "--braker4-output", type=Path, required=True,
        help="BRAKER4 output dir (contains per-sample subdirs).",
    )
    p.add_argument(
        "--samples-csv", type=Path, required=True,
        help="BRAKER4 samples.csv (used to enumerate species).",
    )
    p.add_argument(
        "--cache-dir", type=Path, required=True,
        help="Canonical cache dir (proteomes_braker4/).",
    )
    p.add_argument(
        "--report", type=Path, required=True,
        help="Per-species post-process report TSV.",
    )
    p.add_argument(
        "--force", action="store_true",
        help="Overwrite existing canonical files (default: skip if present).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    sample_names = read_sample_names(args.samples_csv)
    print(
        f"[postprocess_braker4_outputs] {len(sample_names)} samples from "
        f"{args.samples_csv}; cache={args.cache_dir}",
        file=sys.stderr,
    )

    results = postprocess_all(
        sample_names, args.braker4_output, args.cache_dir, force=args.force,
    )
    write_report_tsv(args.report, results)

    n_ok = sum(1 for r in results if r.status == "ok")
    n_no_cds = sum(1 for r in results if r.status == "ok_no_cds")
    n_skipped = sum(1 for r in results if r.status == "skipped")
    n_failed = sum(1 for r in results if r.status == "failed")
    print(
        f"[postprocess_braker4_outputs] done: ok={n_ok} ok_no_cds={n_no_cds} "
        f"skipped={n_skipped} failed={n_failed} (report: {args.report})",
        file=sys.stderr,
    )
    return 0 if n_failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
