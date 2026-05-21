"""Phase 1a-CDS of species-tree pipeline: download paired protein + CDS
FASTAs for the species that Phase 1a inventory identified as having a
usable NCBI annotated assembly.

Bead -9dn under -dnk species-tree epic. Replaces the miniprot-based
recover_cds workflow as stage 05's primary nucleotide-alignment source
(recover_cds stays as a resilience fallback per the 2026-05-21 design
decision).

Each species gets a paired output:

    <cache_dir>/<taxid>_<binomial>.faa          # protein FASTA
    <cache_dir>/<taxid>_<binomial>.cds.fna      # CDS FASTA

Some GenBank-only assemblies ship protein.faa but no cds_from_genomic.fna
— those get status="ok_no_cds" so stage 05 can fall back to recover_cds
just for them, instead of failing the whole batch.

CLI runs on Unity via an sbatch wrapper. Downloads land in scratch
(`/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05/species_tree_data/ncbi_proteomes/`).
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Union


@dataclass(frozen=True)
class DownloadTarget:
    """One species that has a usable NCBI assembly to download."""
    taxid: int
    binomial: str    # genus + species, separated by a single space
    clade: str
    accession: str   # e.g. GCF_000237925.1 (RefSeq) or GCA_021461655.1 (GenBank)


def target_output_paths(target: DownloadTarget, cache_dir: Union[str, Path]) -> tuple[Path, Path]:
    """Return the canonical (protein_path, cds_path) pair for a target.

    Spaces in the binomial become underscores so the filename matches
    the input schema used elsewhere in the species-tree pipeline.
    """
    cache = Path(cache_dir)
    stem = f"{target.taxid}_{target.binomial.replace(' ', '_')}"
    return cache / f"{stem}.faa", cache / f"{stem}.cds.fna"


@dataclass(frozen=True)
class FetchResult:
    """The two FASTA files produced by `datasets download genome accession`
    + unzip, as paths in some staging dir. The downloader copies them
    out to canonical names. CDS may be None when the NCBI archive ships
    only `protein.faa` (common on GenBank-only assemblies).
    """
    protein_faa: Optional[Path]
    cds_fna: Optional[Path]
    error: str = ""
    ok: bool = True


@dataclass(frozen=True)
class DownloadResult:
    """One species' download outcome. `status` is the source of truth
    for the per-species report TSV; downstream code (stage 05 refactor,
    is_already_downloaded) reads `status` to decide retry/skip behavior.

    Statuses:
      ok                 -> paired protein + CDS at canonical paths
      ok_no_cds          -> protein only; NCBI archive lacked CDS
      skipped            -> already in valid end state, no-op this run
      download_failed    -> NCBI download or unzip failed
    """
    target: DownloadTarget
    status: str
    faa_path: Optional[Path]
    cds_path: Optional[Path]
    error: str = ""
    n_proteins: int = 0
    n_cds: int = 0


def _count_fasta_records(path: Path) -> int:
    """Count '>'-prefixed records. Returns 0 when path is None or empty."""
    if not path or not path.exists():
        return 0
    n = 0
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def download_one(
    target: DownloadTarget,
    cache_dir: Union[str, Path],
    *,
    fetcher: Callable[[str, Path], FetchResult],
    prior_status: Optional[str] = None,
    work_root: Optional[Union[str, Path]] = None,
) -> DownloadResult:
    """Orchestrate a single species download: idempotent skip → fetch →
    place at canonical paths → report.

    `fetcher(accession, work_dir) -> FetchResult` is the swappable
    download backend. The default (built at CLI time) calls NCBI
    Datasets and unzips. Tests inject a fake that drops fixture files.

    `prior_status` lets re-runs honor a previously-recorded
    `ok_no_cds` end state (protein-only on a GenBank assembly).
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    faa_dest, cds_dest = target_output_paths(target, cache)

    if is_already_downloaded(target, cache, prior_status=prior_status):
        return DownloadResult(
            target=target, status="skipped",
            faa_path=faa_dest if faa_dest.exists() else None,
            cds_path=cds_dest if cds_dest.exists() else None,
            n_proteins=_count_fasta_records(faa_dest) if faa_dest.exists() else 0,
            n_cds=_count_fasta_records(cds_dest) if cds_dest.exists() else 0,
        )

    work_dir = Path(work_root) if work_root else cache / "_work" / target.accession
    fetch = fetcher(target.accession, work_dir)
    if not fetch.ok:
        return DownloadResult(
            target=target, status="download_failed",
            faa_path=None, cds_path=None,
            error=fetch.error,
        )

    if fetch.protein_faa is None or not fetch.protein_faa.exists():
        # Should not happen for a valid annotated assembly, but treat as
        # a failure rather than silently producing nothing.
        return DownloadResult(
            target=target, status="download_failed",
            faa_path=None, cds_path=None,
            error="protein.faa missing in NCBI download archive",
        )

    shutil.copyfile(fetch.protein_faa, faa_dest)
    n_proteins = _count_fasta_records(faa_dest)

    if fetch.cds_fna is not None and fetch.cds_fna.exists():
        shutil.copyfile(fetch.cds_fna, cds_dest)
        return DownloadResult(
            target=target, status="ok",
            faa_path=faa_dest, cds_path=cds_dest,
            n_proteins=n_proteins, n_cds=_count_fasta_records(cds_dest),
        )

    return DownloadResult(
        target=target, status="ok_no_cds",
        faa_path=faa_dest, cds_path=None,
        n_proteins=n_proteins, n_cds=0,
    )


def is_already_downloaded(
    target: DownloadTarget,
    cache_dir: Union[str, Path],
    prior_status: str | None = None,
) -> bool:
    """True iff this target's paired output is in its valid end state.

    Two distinct "done" states:

    - Both .faa and .cds.fna present and non-empty (the normal happy path).
    - Only .faa present, and a prior run recorded `prior_status='ok_no_cds'`
      for this accession (GenBank-only assemblies that ship protein but
      not CDS — re-running would still produce the same result).

    Anything else (no files, half-finished download, empty file) returns
    False so the downloader retries.
    """
    faa, cds = target_output_paths(target, cache_dir)
    faa_present = faa.exists() and faa.stat().st_size > 0
    cds_present = cds.exists() and cds.stat().st_size > 0
    if faa_present and cds_present:
        return True
    if faa_present and prior_status == "ok_no_cds":
        return True
    return False


REPORT_COLUMNS = (
    "taxid", "binomial", "clade", "accession",
    "status", "n_proteins", "n_cds", "error",
)


# ----------------------------------------------------------------------
# NCBI Datasets CLI primitives
# ----------------------------------------------------------------------

def build_datasets_download_argv(
    accession: str,
    output_zip: Path,
    datasets_bin: str = "datasets",
) -> list[str]:
    """Build the `datasets download genome accession ...` argv for one
    accession. The `--include protein,cds` flag is the central
    behavior change for bead -9dn vs. proteome-only Phase 1a.
    """
    return [
        datasets_bin, "download", "genome", "accession", accession,
        "--include", "protein,cds",
        "--filename", str(output_zip),
    ]


def find_extracted_files(
    work_dir: Path,
    accession: str,
) -> tuple[Optional[Path], Optional[Path]]:
    """Locate protein.faa + cds_from_genomic.fna in the NCBI Datasets
    post-unzip layout: `<work_dir>/ncbi_dataset/data/<accession>/...`.

    Returns (protein_path, cds_path). Either may be None: missing
    protein.faa means the archive was bad; missing cds_from_genomic.fna
    is the common GenBank-only case (handled by stage 05's recover_cds
    fallback).

    Looking only inside the named accession's subdir is defensive: it
    prevents accidentally returning some other assembly's files if the
    work_dir was reused.
    """
    data_dir = Path(work_dir) / "ncbi_dataset" / "data" / accession
    if not data_dir.is_dir():
        return None, None
    faa = data_dir / "protein.faa"
    cds = data_dir / "cds_from_genomic.fna"
    return (faa if faa.exists() else None,
            cds if cds.exists() else None)


def make_ncbi_fetcher(
    datasets_bin: str = "datasets",
    *,
    runner: Optional[Callable] = None,
    timeout: int = 600,
) -> Callable[[str, Path], FetchResult]:
    """Build the real fetcher: call NCBI Datasets, unzip, locate files.

    `runner(argv, ...) -> CompletedProcess`-compatible callable is
    injectable so tests can simulate the subprocess call without
    hitting the network.
    """
    if runner is None:
        runner = subprocess.run

    def fetch(accession: str, work_dir: Path) -> FetchResult:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        zip_path = work_dir / f"{accession}.zip"
        argv = build_datasets_download_argv(
            accession=accession, output_zip=zip_path, datasets_bin=datasets_bin,
        )
        try:
            result = runner(
                argv, capture_output=True, text=True, timeout=timeout,
            )
        except Exception as e:  # subprocess.TimeoutExpired etc.
            return FetchResult(
                protein_faa=None, cds_fna=None,
                error=f"datasets invocation raised {type(e).__name__}: {e}",
                ok=False,
            )

        if result.returncode != 0 or not zip_path.exists():
            stderr = (result.stderr or "").strip()
            return FetchResult(
                protein_faa=None, cds_fna=None,
                error=(f"datasets returncode={result.returncode}: {stderr}"
                       if stderr else
                       f"datasets returncode={result.returncode}, no zip produced"),
                ok=False,
            )

        try:
            with zipfile.ZipFile(zip_path) as z:
                z.extractall(work_dir)
        except zipfile.BadZipFile as e:
            return FetchResult(
                protein_faa=None, cds_fna=None,
                error=f"corrupt zip archive: {e}",
                ok=False,
            )

        faa, cds = find_extracted_files(work_dir, accession)
        if faa is None:
            return FetchResult(
                protein_faa=None, cds_fna=None,
                error="protein.faa missing from NCBI archive",
                ok=False,
            )
        return FetchResult(protein_faa=faa, cds_fna=cds, error="", ok=True)

    return fetch


# ----------------------------------------------------------------------
# CLI entrypoint
# ----------------------------------------------------------------------

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Phase 1a-CDS: download paired protein + CDS FASTAs for the "
            "species identified by the Phase 1a inventory. One zip per "
            "accession via NCBI Datasets CLI; outputs land at "
            "<cache_dir>/<taxid>_<binomial>.{faa,cds.fna}."
        ),
    )
    p.add_argument(
        "--manifest", type=Path, required=True,
        help=("Phase 1a inventory TSV (references/species_tree/proteome_manifest.tsv)."),
    )
    p.add_argument(
        "--cache-dir", type=Path, required=True,
        help=("Output dir. On Unity this lands at "
              "/scratch3/.../species_tree_data/ncbi_proteomes/."),
    )
    p.add_argument(
        "--datasets-bin", default="datasets",
        help="Path to NCBI Datasets CLI binary (default: 'datasets' on PATH).",
    )
    p.add_argument(
        "--report", type=Path, default=None,
        help=("Per-species report TSV. Defaults to "
              "<cache_dir>/download_report.tsv. Re-runs read this to "
              "honor 'ok_no_cds' as a legitimate end state."),
    )
    p.add_argument(
        "--timeout", type=int, default=600,
        help="Per-accession datasets-CLI timeout in seconds.",
    )
    p.add_argument(
        "--work-root", type=Path, default=None,
        help=("Optional shared scratch dir for per-accession unzip work. "
              "Defaults to <cache_dir>/_work/<accession>/."),
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    report_path = args.report or (args.cache_dir / "download_report.tsv")

    targets = read_download_targets(args.manifest)
    print(
        f"[download_species_tree_phase1a] {len(targets)} targets read from "
        f"{args.manifest}; cache={args.cache_dir}",
        file=sys.stderr,
    )

    prior = read_prior_status(report_path)
    if prior:
        print(
            f"[download_species_tree_phase1a] prior report at {report_path}: "
            f"{len(prior)} previously-recorded statuses",
            file=sys.stderr,
        )

    fetcher = make_ncbi_fetcher(
        datasets_bin=args.datasets_bin, timeout=args.timeout,
    )
    results = download_all(
        targets, args.cache_dir,
        fetcher=fetcher,
        prior_statuses=prior,
        work_root=args.work_root,
    )
    write_report_tsv(report_path, results)

    n_ok = sum(1 for r in results if r.status == "ok")
    n_no_cds = sum(1 for r in results if r.status == "ok_no_cds")
    n_skipped = sum(1 for r in results if r.status == "skipped")
    n_failed = sum(1 for r in results if r.status == "download_failed")
    print(
        f"[download_species_tree_phase1a] done: ok={n_ok} ok_no_cds={n_no_cds} "
        f"skipped={n_skipped} failed={n_failed} (report: {report_path})",
        file=sys.stderr,
    )
    return 0 if n_failed == 0 else 2


def write_report_tsv(
    report_path: Union[str, Path],
    results: list[DownloadResult],
) -> None:
    """Persist the per-species download report TSV. The `status` column
    is the source of truth for the next run's prior-status lookup.
    """
    path = Path(report_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(REPORT_COLUMNS) + "\n")
        for r in results:
            f.write("\t".join([
                str(r.target.taxid),
                r.target.binomial,
                r.target.clade,
                r.target.accession,
                r.status,
                str(r.n_proteins),
                str(r.n_cds),
                r.error.replace("\t", " ").replace("\n", " "),
            ]) + "\n")


def read_prior_status(report_path: Union[str, Path]) -> dict[str, str]:
    """Return {accession: status} from a prior report TSV.

    Missing/empty file returns an empty dict (first-run case).
    Used to honor 'ok_no_cds' as a legitimate end state across runs.
    """
    path = Path(report_path)
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = (row.get("accession") or "").strip()
            status = (row.get("status") or "").strip()
            if acc and status:
                out[acc] = status
    return out


def download_all(
    targets: list[DownloadTarget],
    cache_dir: Union[str, Path],
    *,
    fetcher: Callable[[str, Path], FetchResult],
    prior_statuses: Optional[dict[str, str]] = None,
    progress_every: int = 5,
    work_root: Optional[Union[str, Path]] = None,
) -> list[DownloadResult]:
    """Loop over targets, return one DownloadResult per target in input
    order. A single failure does NOT abort the batch — failed accessions
    are recorded with status='download_failed' so the next run can retry
    just those.
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    prior = prior_statuses or {}
    results: list[DownloadResult] = []
    for i, target in enumerate(targets, start=1):
        prior_status = prior.get(target.accession)
        result = download_one(
            target, cache,
            fetcher=fetcher,
            prior_status=prior_status,
            work_root=work_root,
        )
        results.append(result)
        if i % progress_every == 0 or i == len(targets):
            n_ok = sum(1 for r in results if r.status == "ok")
            n_no_cds = sum(1 for r in results if r.status == "ok_no_cds")
            n_skip = sum(1 for r in results if r.status == "skipped")
            n_fail = sum(1 for r in results if r.status == "download_failed")
            print(
                f"  [{i}/{len(targets)}] ok={n_ok} ok_no_cds={n_no_cds} "
                f"skipped={n_skip} failed={n_fail}",
                file=sys.stderr,
            )
    return results


def read_download_targets(manifest_path: Union[str, Path]) -> list[DownloadTarget]:
    """Parse the Phase 1a inventory TSV and return only the rows with a
    non-empty `accession` (the species we can actually download).
    Dropped rows (no_proteome_in_ncbi, query_error, ...) are skipped.

    Result is sorted by taxid for deterministic downstream ordering.
    """
    path = Path(manifest_path)
    targets: list[DownloadTarget] = []
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
            targets.append(DownloadTarget(
                taxid=taxid,
                binomial=(row.get("binomial") or "").strip(),
                clade=(row.get("clade") or "").strip(),
                accession=accession,
            ))
    targets.sort(key=lambda t: t.taxid)
    return targets


if __name__ == "__main__":
    sys.exit(main())
