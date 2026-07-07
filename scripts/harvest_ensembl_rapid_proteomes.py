#!/usr/bin/env python3
"""harvest_ensembl_rapid_proteomes.py — harvest Ensembl rapid-release proteomes.

Closes the Ensembl-annotation gap in the Phase-1f BRAKER set (bead 95b):
some genome-only species in genome_inventory.tsv have their EXACT assembly
already annotated by Ensembl rapid-release, so de-novo BRAKER annotation is
redundant. This driver downloads those Ensembl protein sets (with full
provenance), stages them under the standard proteome cache naming
(<cache>/<taxid>_<binomial>.faa[+ .cds.fna]), and emits the rows to append to
proteome_manifest.tsv + a detailed provenance sidecar. A separate `mark` step
then writes drop_reason='harvested_ensembl' on the genome_inventory rows so both
BRAKER consumers skip them.

Input work table (--work), one row per species, TAB-separated with header:
    taxid  binomial  clade  assembly_level  accession  ens_production_name  ens_species

Provenance recorded per species (user standing rule): source, assembly
accession, ensembl_production_name, geneset id + release, release_date, exact
download URL, access date (UTC, passed in via --access-date), sha256, n_proteins.

Ensembl rapid-release FTP layout:
  /pub/rapid-release/species/<Ens_Species>/<accession>/ensembl/geneset/<rel>/
      <Ens_Species>-<accession>-<rel>-pep.fa.gz   (proteins)
      <Ens_Species>-<accession>-<rel>-cds.fa.gz   (CDS)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import subprocess
import sys
import time
from pathlib import Path

FTP_ROOT = "https://ftp.ensembl.org/pub/rapid-release/species"
REL_RE = re.compile(r'href="(\d{4}_\d{2})/"')


def sanitize(name: str) -> str:
    """<taxid>_<binomial> stem convention: spaces -> underscores."""
    return name.strip().replace(" ", "_")


def curl(url: str, dest: Path | None = None, retries: int = 4) -> bytes | None:
    """GET url. Returns body bytes (dest=None) or writes to dest (returns b'')."""
    for attempt in range(retries):
        cmd = ["curl", "-sSfL", "--max-time", "180"]
        if dest is not None:
            cmd += ["-o", str(dest)]
        cmd.append(url)
        r = subprocess.run(cmd, capture_output=(dest is None))
        if r.returncode == 0:
            return r.stdout if dest is None else b""
        time.sleep(2 * (attempt + 1))
    return None


PROVIDERS = ("ensembl", "braker")  # Ensembl rapid-release annotates via either pipeline


def resolve_provider_release(species_dir: str, accession: str) -> tuple[str, str] | None:
    """Find the (provider, latest YYYY_MM release) for a species.

    Ensembl rapid-release serves genesets under either `ensembl/` (its own
    gene-annotation pipeline) or `braker/` (BRAKER-annotated). Try each.
    """
    for provider in PROVIDERS:
        listing = curl(f"{FTP_ROOT}/{species_dir}/{accession}/{provider}/geneset/")
        if listing is None:
            continue
        rels = REL_RE.findall(listing.decode("utf-8", "replace"))
        if rels:
            return provider, sorted(rels)[-1]
    return None


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def count_fasta(path: Path) -> int:
    n = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
    return n


def gunzip_to(src_gz: Path, dest: Path) -> None:
    with gzip.open(src_gz, "rb") as fi, open(dest, "wb") as fo:
        for chunk in iter(lambda: fi.read(1 << 20), b""):
            fo.write(chunk)


def harvest_one(row: dict, cache: Path, work: Path, access_date: str,
                force: bool = False) -> dict:
    taxid, binom = row["taxid"], row["binomial"]
    acc, ens_sp, prod = row["accession"], row["ens_species"], row["ens_production_name"]
    species_dir = sanitize(ens_sp)
    stem = f"{taxid}_{sanitize(binom)}"
    faa = cache / f"{stem}.faa"
    result = {**row, "status": "", "n_proteins": "", "url": "", "release": "",
              "sha256": "", "error": ""}

    if not force and faa.exists() and faa.stat().st_size > 0:
        result.update(status="skipped", n_proteins=count_fasta(faa))
        return result

    pr = resolve_provider_release(species_dir, acc)
    if not pr:
        result.update(status="error", error="no geneset release dir found")
        return result
    provider, rel = pr

    base = f"{FTP_ROOT}/{species_dir}/{acc}/{provider}/geneset/{rel}"
    pep_url = f"{base}/{species_dir}-{acc}-{rel}-pep.fa.gz"
    cds_url = f"{base}/{species_dir}-{acc}-{rel}-cds.fa.gz"
    work.mkdir(parents=True, exist_ok=True)
    pep_gz = work / f"{stem}.pep.fa.gz"
    if curl(pep_url, pep_gz) is None:
        result.update(status="error", error=f"pep download failed: {pep_url}")
        return result
    cache.mkdir(parents=True, exist_ok=True)
    gunzip_to(pep_gz, faa)
    npep = count_fasta(faa)
    if npep == 0:
        faa.unlink(missing_ok=True)
        result.update(status="error", error="0 proteins after gunzip")
        return result

    # CDS is best-effort (mirrors NCBI ok_no_cds semantics).
    cds_gz = work / f"{stem}.cds.fa.gz"
    if curl(cds_url, cds_gz) is not None:
        gunzip_to(cds_gz, cache / f"{stem}.cds.fna")

    result.update(status="ok", n_proteins=npep, url=pep_url, release=rel,
                  provider=provider, sha256=sha256_of(faa), access_date=access_date)
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, type=Path, help="29-species work table TSV")
    ap.add_argument("--cache", required=True, type=Path, help="proteome cache dir")
    ap.add_argument("--tmp", type=Path, default=Path("/tmp/ens_harvest"),
                    help="scratch dir for .gz downloads")
    ap.add_argument("--access-date", required=True,
                    help="UTC access date (provenance); pass date -u +%%Y-%%m-%%d")
    ap.add_argument("--out-report", required=True, type=Path)
    ap.add_argument("--out-manifest-add", required=True, type=Path,
                    help="proteome_manifest rows to append (10-col)")
    ap.add_argument("--out-provenance", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=0, help="process only first N (testing)")
    ap.add_argument("--force", action="store_true",
                    help="re-download even if the cached .faa exists")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.work), delimiter="\t"))
    if args.limit:
        rows = rows[: args.limit]

    results = []
    for i, row in enumerate(rows, 1):
        res = harvest_one(row, args.cache, args.tmp, args.access_date, force=args.force)
        results.append(res)
        print(f"[{i}/{len(rows)}] {res['status']:>7}  {row['taxid']}_{row['binomial']}"
              f"  n={res['n_proteins']}  {res.get('error','')}", flush=True)
        time.sleep(0.3)

    ok = [r for r in results if r["status"] in ("ok", "skipped")]
    # Report
    with open(args.out_report, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["taxid", "binomial", "accession", "status", "n_proteins", "error"])
        for r in results:
            w.writerow([r["taxid"], r["binomial"], r["accession"], r["status"],
                        r["n_proteins"], r["error"]])
    # proteome_manifest additions (10-col schema, drop_reason empty = usable)
    with open(args.out_manifest_add, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        for r in results:
            if r["status"] != "ok":
                continue
            w.writerow([r["taxid"], r["binomial"], r["clade"], "Ensembl_rapid_release",
                        r["accession"], r["assembly_level"], "Ensembl rapid-release geneset",
                        r["n_proteins"], "", ""])
    # provenance sidecar (full detail)
    with open(args.out_provenance, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["taxid", "binomial", "accession", "source", "ens_production_name",
                    "ens_species", "annotation_provider", "geneset_release", "download_url",
                    "access_date", "sha256", "n_proteins"])
        for r in results:
            if r["status"] != "ok":
                continue
            w.writerow([r["taxid"], r["binomial"], r["accession"], "Ensembl_rapid_release",
                        r["ens_production_name"], r["ens_species"], r.get("provider", ""),
                        r["release"], r["url"], r.get("access_date", ""), r["sha256"],
                        r["n_proteins"]])

    n_ok = sum(1 for r in results if r["status"] == "ok")
    n_skip = sum(1 for r in results if r["status"] == "skipped")
    n_err = sum(1 for r in results if r["status"] == "error")
    print(f"\nDONE: ok={n_ok} skipped={n_skip} error={n_err} / {len(results)}")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
