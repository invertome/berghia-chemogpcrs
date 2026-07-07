"""Phase 1f of species-tree pipeline: download genome FASTAs for the 134
species that Phase 1e identified as having a public GenBank assembly but
no public annotation. BRAKER4 (bead -p49) re-annotates these from scratch.

Companion to scripts/download_species_tree_phase1a.py (bead -9dn):
- Phase 1a-CDS:  --include protein,cds  -> <taxid>_<binomial>.{faa,cds.fna}
- Phase 1f:      --include genome       -> <taxid>_<binomial>.fasta

Output dir on Unity:
    /scratch3/.../species_tree_data/braker4_genomes/<taxid>_<binomial>.fasta

Most behaviors mirror -9dn:
  - Idempotent re-runs via download_report.tsv prior-status lookup.
  - Failed accessions recorded as 'download_failed', batch continues.
  - Defensive: glob-matches `<accession>_*_genomic.fna` only (NOT
    `cds_from_genomic.fna`, which would otherwise match the suffix).
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

# Share the samples-builder's canonical stem so the genome filename this writes
# is byte-identical to the name build_braker4_samples_csv later looks up. Using
# a bare `binomial.replace(' ', '_')` here diverged for punctuated names
# (e.g. "Ctena cf. galapagana" / "Chaetoderma sp. LZ-2023a"): the builder strips
# '.'/'-' and collapses underscores, so those genomes were downloaded but never
# found -> silently dropped from the BRAKER set.
from build_braker4_samples_csv import sanitize_sample_name


@dataclass(frozen=True)
class DownloadTarget:
    """One species with a usable GenBank genome to download (unannotated)."""
    taxid: int
    binomial: str
    clade: str
    accession: str


def target_output_path(target: DownloadTarget, cache_dir: Union[str, Path]) -> Path:
    """Canonical genome FASTA path for one species: `<cache_dir>/<taxid>_<binomial>.fasta`.

    The stem is sanitized with the SAME function the samples-builder uses
    (sanitize_sample_name: non-[A-Za-z0-9_] -> '_', collapsed) so the two agree
    on the filename for every species, including punctuated names.
    """
    cache = Path(cache_dir)
    stem = sanitize_sample_name(f"{target.taxid}_{target.binomial}")
    return cache / f"{stem}.fasta"


@dataclass(frozen=True)
class FetchResult:
    """Result of one `datasets download genome accession ...` + unzip."""
    genome_fa: Optional[Path]
    error: str = ""
    ok: bool = True


@dataclass(frozen=True)
class DownloadResult:
    """Per-species outcome. `status` drives report TSV + idempotent re-runs.

    Statuses:
      ok               -> genome FASTA at canonical path, non-empty
      skipped          -> already in valid end state, no-op this run
      download_failed  -> NCBI download or unzip failed, or archive missing genome
    """
    target: DownloadTarget
    status: str
    genome_path: Optional[Path]
    error: str = ""
    n_seqs: int = 0
    size_bytes: int = 0


def _count_fasta_records(path: Path) -> int:
    if not path or not path.exists():
        return 0
    n = 0
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def is_already_downloaded(
    target: DownloadTarget,
    cache_dir: Union[str, Path],
) -> bool:
    """True iff the canonical genome FASTA exists and is non-empty."""
    path = target_output_path(target, cache_dir)
    return path.exists() and path.stat().st_size > 0


def build_datasets_download_argv(
    accession: str,
    output_zip: Path,
    datasets_bin: str = "datasets",
) -> list[str]:
    """`datasets download genome accession <ACC> --include genome --filename <ZIP>`.

    Only genome (not protein, not cds) — for BRAKER4 to re-annotate.
    """
    return [
        datasets_bin, "download", "genome", "accession", accession,
        "--include", "genome",
        "--filename", str(output_zip),
    ]


def find_extracted_genome(
    work_dir: Path,
    accession: str,
) -> Optional[Path]:
    """Locate the genome FASTA in the NCBI Datasets post-unzip layout:
    `<work_dir>/ncbi_dataset/data/<accession>/<accession>_<asm>_genomic.fna`.

    Defensive: matches ONLY files whose name starts with `<accession>_`
    AND ends with `_genomic.fna`. This excludes `cds_from_genomic.fna`
    (which doesn't start with the accession) so a polluted workdir
    can't return the wrong file.
    """
    data_dir = Path(work_dir) / "ncbi_dataset" / "data" / accession
    if not data_dir.is_dir():
        return None
    for f in sorted(data_dir.iterdir()):
        name = f.name
        if name.startswith(f"{accession}_") and name.endswith("_genomic.fna"):
            return f
    return None


def make_ncbi_fetcher(
    datasets_bin: str = "datasets",
    *,
    runner: Optional[Callable] = None,
    timeout: int = 900,
) -> Callable[[str, Path], FetchResult]:
    """Build the real fetcher: call NCBI Datasets, unzip, locate genome FASTA.

    Genomes are big (50 MB - 5 GB compressed); default timeout is 15 min
    per accession. The Phase 1e manifest's largest genome (Conus tribblei,
    ~2.2 GB) takes ~5 min to download on a fast link.
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
        except Exception as e:
            return FetchResult(
                genome_fa=None,
                error=f"datasets invocation raised {type(e).__name__}: {e}",
                ok=False,
            )

        if result.returncode != 0 or not zip_path.exists():
            stderr = (result.stderr or "").strip()
            return FetchResult(
                genome_fa=None,
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
                genome_fa=None,
                error=f"corrupt zip archive: {e}",
                ok=False,
            )

        genome_fa = find_extracted_genome(work_dir, accession)
        if genome_fa is None:
            return FetchResult(
                genome_fa=None,
                error="genome FASTA missing from NCBI archive",
                ok=False,
            )
        return FetchResult(genome_fa=genome_fa, error="", ok=True)

    return fetch


def download_one(
    target: DownloadTarget,
    cache_dir: Union[str, Path],
    *,
    fetcher: Callable[[str, Path], FetchResult],
    work_root: Optional[Union[str, Path]] = None,
) -> DownloadResult:
    """Orchestrate one species: idempotent skip → fetch → place at canonical
    path → report.
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    dest = target_output_path(target, cache)

    if is_already_downloaded(target, cache):
        return DownloadResult(
            target=target, status="skipped",
            genome_path=dest,
            n_seqs=_count_fasta_records(dest),
            size_bytes=dest.stat().st_size,
        )

    work_dir = Path(work_root) if work_root else cache / "_work" / target.accession
    fetch = fetcher(target.accession, work_dir)
    if not fetch.ok:
        return DownloadResult(
            target=target, status="download_failed",
            genome_path=None, error=fetch.error,
        )

    if fetch.genome_fa is None or not fetch.genome_fa.exists():
        return DownloadResult(
            target=target, status="download_failed",
            genome_path=None,
            error="genome FASTA missing in NCBI download archive",
        )

    shutil.copyfile(fetch.genome_fa, dest)
    return DownloadResult(
        target=target, status="ok",
        genome_path=dest,
        n_seqs=_count_fasta_records(dest),
        size_bytes=dest.stat().st_size,
    )


REPORT_COLUMNS = (
    "taxid", "binomial", "clade", "accession",
    "status", "n_seqs", "size_bytes", "error",
)


def write_report_tsv(
    report_path: Union[str, Path],
    results: list[DownloadResult],
) -> None:
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
                str(r.n_seqs),
                str(r.size_bytes),
                r.error.replace("\t", " ").replace("\n", " "),
            ]) + "\n")


def read_prior_status(report_path: Union[str, Path]) -> dict[str, str]:
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
    progress_every: int = 5,
    work_root: Optional[Union[str, Path]] = None,
) -> list[DownloadResult]:
    """Loop over targets, one DownloadResult per target in input order.
    Failures don't abort the batch.
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    results: list[DownloadResult] = []
    for i, target in enumerate(targets, start=1):
        result = download_one(
            target, cache, fetcher=fetcher, work_root=work_root,
        )
        results.append(result)
        if i % progress_every == 0 or i == len(targets):
            n_ok = sum(1 for r in results if r.status == "ok")
            n_skip = sum(1 for r in results if r.status == "skipped")
            n_fail = sum(1 for r in results if r.status == "download_failed")
            total_bytes = sum(r.size_bytes for r in results)
            print(
                f"  [{i}/{len(targets)}] ok={n_ok} skipped={n_skip} failed={n_fail} "
                f"(total {total_bytes/1e9:.1f} GB)",
                file=sys.stderr,
            )
    return results


def read_download_targets(*manifest_paths: Union[str, Path]) -> list[DownloadTarget]:
    """Read one or more manifest TSVs; return DownloadTargets only for
    rows with a non-empty accession.

    Supports both schemas:
      - Phase 1e (genome_inventory_unannotated.tsv): clade column
        named `clade`.
      - Phase 1d (extension_inventory.tsv): clade column named
        `clade_name`; its `policy_class` column is ignored.

    Dedup: when the same taxid appears in multiple manifests, the
    FIRST manifest in `manifest_paths` order wins. Pass Phase 1e
    first to keep its annotation lineage as the canonical source.
    """
    targets: list[DownloadTarget] = []
    seen: set[int] = set()
    for manifest_path in manifest_paths:
        path = Path(manifest_path)
        with path.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                accession = (row.get("accession") or "").strip()
                if not accession:
                    continue
                if (row.get("drop_reason") or "").strip():
                    continue
                try:
                    taxid = int(row["taxid"])
                except (KeyError, ValueError):
                    continue
                if taxid in seen:
                    continue
                seen.add(taxid)
                clade = (row.get("clade") or row.get("clade_name") or "").strip()
                targets.append(DownloadTarget(
                    taxid=taxid,
                    binomial=(row.get("binomial") or "").strip(),
                    clade=clade,
                    accession=accession,
                ))
    targets.sort(key=lambda t: t.taxid)
    return targets


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Phase 1f: download genome FASTAs for the 134 unannotated-genome "
            "species identified by Phase 1e inventory. Outputs land at "
            "<cache_dir>/<taxid>_<binomial>.fasta for BRAKER4 to consume."
        ),
    )
    p.add_argument(
        "--manifest", type=Path, nargs="+",
        default=[Path("references/species_tree/genome_inventory.tsv")],
        help=(
            "One or more inventory TSVs (default: references/species_tree/"
            "genome_inventory.tsv). Accepts the unified schema or the older "
            "Phase 1e (genome_inventory_unannotated.tsv) and/or Phase 1d "
            "(extension_inventory.tsv) schemas. When passing multiple files, "
            "list the higher-priority one first so its entries win per-taxid dedup."
        ),
    )
    p.add_argument(
        "--cache-dir", type=Path, required=True,
        help="Output dir (species_tree_data/braker4_genomes/ on Unity).",
    )
    p.add_argument(
        "--datasets-bin", default="datasets",
        help="Path to NCBI Datasets CLI binary.",
    )
    p.add_argument(
        "--report", type=Path, default=None,
        help="Report TSV (default <cache_dir>/download_report.tsv).",
    )
    p.add_argument(
        "--timeout", type=int, default=900,
        help="Per-accession datasets-CLI timeout in seconds (default 15 min).",
    )
    p.add_argument(
        "--work-root", type=Path, default=None,
        help="Optional scratch dir for unzip work (default <cache_dir>/_work/).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    report_path = args.report or (args.cache_dir / "download_report.tsv")

    targets = read_download_targets(*args.manifest)
    print(
        f"[download_species_tree_phase1f] {len(targets)} targets read from "
        f"{args.manifest}; cache={args.cache_dir}",
        file=sys.stderr,
    )

    fetcher = make_ncbi_fetcher(
        datasets_bin=args.datasets_bin, timeout=args.timeout,
    )
    results = download_all(
        targets, args.cache_dir, fetcher=fetcher, work_root=args.work_root,
    )
    write_report_tsv(report_path, results)

    n_ok = sum(1 for r in results if r.status == "ok")
    n_skipped = sum(1 for r in results if r.status == "skipped")
    n_failed = sum(1 for r in results if r.status == "download_failed")
    total_bytes = sum(r.size_bytes for r in results)
    print(
        f"[download_species_tree_phase1f] done: ok={n_ok} skipped={n_skipped} "
        f"failed={n_failed} total={total_bytes/1e9:.1f} GB "
        f"(report: {report_path})",
        file=sys.stderr,
    )
    return 0 if n_failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
