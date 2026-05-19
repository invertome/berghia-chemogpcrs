#!/usr/bin/env python3
"""fetch_curated_chemoreceptors.py — Fetch curated invertebrate
chemoreceptor protein sequences from NCBI by accession.

Reads references/curated_chemoreceptors_lophotrochozoa/manifest.tsv
and writes per-family FASTA files. Caches: skips entries whose
sequence is already present in the output file (idempotent).

Manifest TSV columns:
    family  subfamily  gene_name  organism  taxid  accession_type
    accession  citation_doi  validation_method  notes

Per-family output:
    references/curated_chemoreceptors_lophotrochozoa/<family>/<family>_seqs.fa

The FASTA header for each fetched sequence is:
    >{family}__{gene_name}__{taxid}__{accession}

Slashes / spaces / brackets stripped. Family is the prefix to make
later per-family extraction grep-able.

Bead -k0g (curated invertebrate chemoreceptor HMM library scaffolding,
2026-05-19). Stage 04 of the scaffolding effort.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

from Bio import Entrez, SeqIO


def _safe_id_part(s: str) -> str:
    """Make a FASTA-header-safe token (no spaces, no special chars)."""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in s.strip())


def parse_manifest(manifest_path: str) -> list[dict]:
    """Read the curation TSV. Returns list of dicts (one per row).
    Skips comment lines starting with #."""
    rows = []
    with open(manifest_path) as f:
        # Filter # comment lines from the DictReader stream
        lines = (line for line in f if not line.lstrip().startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            if not row.get("accession", "").strip():
                continue
            rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def existing_accessions(fasta_path: str) -> set[str]:
    """Set of accession tokens already in a FASTA file (parsed from
    header pattern family__gene__taxid__accession)."""
    seen: set[str] = set()
    if not os.path.exists(fasta_path):
        return seen
    for rec in SeqIO.parse(fasta_path, "fasta"):
        parts = rec.id.split("__")
        if parts:
            seen.add(parts[-1])
    return seen


def fetch_one(accession: str, accession_type: str, retries: int = 3) -> str:
    """Fetch a single protein sequence from NCBI. Returns the raw
    sequence string (no header). Raises on failure after retries."""
    if accession_type in ("refseq", "genbank"):
        db = "protein"
    elif accession_type == "uniprot":
        # UniProt accessions can be fetched via NCBI protein DB too,
        # but the canonical route is uniprot.org. Stay on NCBI for
        # one-stop-shop; falls back to error if UniProt-only accession.
        db = "protein"
    elif accession_type == "ensembl":
        db = "protein"
    else:
        db = "protein"

    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(db=db, id=accession, rettype="fasta",
                                   retmode="text")
            content = handle.read()
            handle.close()
            if not content.startswith(">"):
                raise RuntimeError(
                    f"Unexpected response for {accession} (db={db}): "
                    f"{content[:80]!r}"
                )
            # Strip header, return only sequence
            lines = content.splitlines()
            seq = "".join(line.strip() for line in lines[1:] if line.strip())
            if len(seq) < 50:
                raise RuntimeError(
                    f"Sequence too short for {accession}: {len(seq)} aa"
                )
            return seq
        except Exception as e:  # noqa: BLE001
            last_err = e
            if attempt < retries - 1:
                wait = 2 ** attempt
                print(
                    f"  retry {attempt + 1}/{retries} for {accession} "
                    f"in {wait}s ({e})",
                    file=sys.stderr,
                )
                time.sleep(wait)
    raise RuntimeError(f"Failed to fetch {accession}: {last_err}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--manifest",
        default="references/curated_chemoreceptors_lophotrochozoa/manifest.tsv",
    )
    ap.add_argument(
        "--output-dir",
        default="references/curated_chemoreceptors_lophotrochozoa",
        help="Per-family output: <output-dir>/<family>/<family>_seqs.fa",
    )
    ap.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL", "xibalbanus@gmail.com"),
        help="NCBI Entrez requires an email contact",
    )
    ap.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY"),
        help="Optional NCBI API key for higher rate limits",
    )
    ap.add_argument(
        "--throttle", type=float, default=0.4,
        help="Seconds to sleep between requests (NCBI free-tier limit: "
             "3 req/sec without API key, 10 with)",
    )
    ap.add_argument("--dry-run", action="store_true",
                    help="Parse manifest + check existing cache, don't fetch")
    args = ap.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    rows = parse_manifest(args.manifest)
    if not rows:
        print(f"No entries in manifest {args.manifest}", file=sys.stderr)
        return 1

    by_family: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_family[r["family"]].append(r)

    total_fetched = 0
    total_skipped = 0
    total_failed = 0
    insufficient_curation: list[str] = []

    for family, entries in sorted(by_family.items()):
        family_dir = Path(args.output_dir) / family
        family_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = family_dir / f"{family}_seqs.fa"

        if len(entries) < 3:
            insufficient_curation.append(
                f"{family} ({len(entries)} entries; hmmbuild needs >=3)"
            )

        cached = existing_accessions(str(fasta_path))
        new_records: list[tuple[str, str]] = []
        for entry in entries:
            acc = entry["accession"]
            if acc in cached:
                total_skipped += 1
                continue
            if args.dry_run:
                print(f"  would fetch [{family}] {acc} "
                      f"({entry.get('gene_name', '?')})", file=sys.stderr)
                continue
            try:
                seq = fetch_one(acc, entry.get("accession_type", "refseq"))
                header_token = "__".join([
                    _safe_id_part(entry["family"]),
                    _safe_id_part(entry.get("gene_name", "")),
                    _safe_id_part(entry.get("taxid", "")),
                    _safe_id_part(acc),
                ])
                new_records.append((header_token, seq))
                total_fetched += 1
                print(f"  fetched [{family}] {acc} "
                      f"({entry.get('gene_name', '?')}, {len(seq)} aa)",
                      file=sys.stderr)
                time.sleep(args.throttle)
            except Exception as e:  # noqa: BLE001
                print(f"  ERROR [{family}] {acc}: {e}", file=sys.stderr)
                total_failed += 1

        # Append (don't truncate) — preserves cache from prior runs
        if new_records:
            with open(fasta_path, "a") as out:
                for header, seq in new_records:
                    out.write(f">{header}\n{seq}\n")

    print("", file=sys.stderr)
    print(f"Summary: fetched={total_fetched} skipped_cached={total_skipped} "
          f"failed={total_failed}", file=sys.stderr)
    if insufficient_curation:
        print(f"Insufficient curation (will be skipped by builder): "
              f"{', '.join(insufficient_curation)}", file=sys.stderr)
    return 0 if total_failed == 0 else 3


if __name__ == "__main__":
    sys.exit(main())
