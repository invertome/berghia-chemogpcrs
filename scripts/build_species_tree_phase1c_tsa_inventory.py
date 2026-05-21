"""Phase 1c of species-tree pipeline: TSA inventory for dropped species.

Bead -2x3 (under -dnk epic). Phase 1a (-c5d) found only 91/239 species
have public proteomes in NCBI Datasets (38.1% coverage). This phase
queries NCBI `nuccore` for TSA (Transcriptome Shotgun Assembly) entries
on the 148 dropped species. Result decides whether to commit to a
Phase 1d EvidentialGene processing step that converts TSAs into
proteome proxies for BUSCO.

Reads Phase 1a's proteome_manifest.tsv, filters to drop_reason !=""
rows, queries Bio.Entrez per taxid, writes a TSA inventory TSV.

Query format (manually confirmed 2026-05-21):
  esearch -db nuccore -query "txid<taxid>[Organism]
                              AND srcdb_genbank[Properties]
                              AND TSA[Keyword]"
On Platynereis dumerilii (taxid 6359): 10,928 records, master GBZT.
On Riftia pachyptila (taxid 6426): 0 records.

Runs on Unity (sbatch on cpu/cpu-preempt, modest resources, ~30 min
for ~148 serial Entrez queries with API key + 0.34s/call rate cap).
Downloads + manifest land on scratch; final TSV scp'd to repo.

See docs/plans/2026-05-21-species-tree-design.md.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Union

from Bio import Entrez


# ----------------------------------------------------------------------
# Read dropped taxa from Phase 1a manifest
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class DroppedTaxon:
    taxid: int
    binomial: str
    clade: str


def read_dropped_taxa_from_manifest(
    manifest_path: Union[str, Path],
) -> list[DroppedTaxon]:
    """Read Phase 1a proteome_manifest.tsv and return only the rows
    where drop_reason is non-empty (i.e. the species we want to recover
    via TSA in this phase). Sorted by taxid for determinism.
    """
    path = Path(manifest_path)
    if not path.exists():
        raise FileNotFoundError(f"manifest not found: {path}")

    dropped: list[DroppedTaxon] = []
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if not row.get("drop_reason"):
                continue
            try:
                taxid = int(row["taxid"])
            except (KeyError, ValueError):
                continue
            dropped.append(DroppedTaxon(
                taxid=taxid,
                binomial=row.get("binomial", ""),
                clade=row.get("clade", ""),
            ))
    dropped.sort(key=lambda t: t.taxid)
    return dropped


# ----------------------------------------------------------------------
# Master-prefix derivation
# ----------------------------------------------------------------------

_WGS_RE = re.compile(r"^([A-Z]{4})\d{2,}")


def derive_master_prefix(accession: str) -> str:
    """TSA contig accessions follow WGS format: 4-char-prefix + 2 digits
    + 6 digits (e.g. GBZT01000001). Master record uses 4-char prefix +
    8 zeros. Return the 4-char prefix for use as a stable group ID.

    Falls back to returning the input verbatim if the accession doesn't
    look like WGS — we want the inventory to surface something even on
    unusual accession formats.
    """
    if not accession:
        return ""
    m = _WGS_RE.match(accession)
    if m:
        return m.group(1)
    return accession


# ----------------------------------------------------------------------
# Bio.Entrez TSA query
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class TsaQueryResult:
    has_tsa: bool
    n_records: int
    master_prefix: str
    first_accession: str


def query_tsa_for_taxon(
    taxid: int,
    email: str | None = None,
    api_key: str | None = None,
) -> TsaQueryResult:
    """Query NCBI nuccore for TSA entries belonging to `taxid`. Returns
    a TsaQueryResult — None of has_tsa, n_records, master_prefix is
    populated when there's no hit.

    Sets Entrez.email + Entrez.api_key (when supplied) so NCBI lets us
    run >3 req/s. The orchestrator should still respect a small sleep
    between calls to be polite.

    Raises whatever Bio.Entrez raises on network failure (caller's job
    to record as a query_error drop reason without aborting the batch).
    """
    if email:
        Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    term = (
        f"txid{taxid}[Organism] AND "
        "srcdb_genbank[Properties] AND "
        "TSA[Keyword]"
    )
    handle = Entrez.esearch(db="nuccore", term=term, retmax=1)
    record = Entrez.read(handle)
    try:
        # esearch handle can be closed; some Bio.Entrez versions are lazy
        handle.close()  # type: ignore[union-attr]
    except Exception:
        pass

    n = int(record.get("Count", 0) or 0)
    idlist = record.get("IdList", []) or []

    if n == 0 or not idlist:
        return TsaQueryResult(
            has_tsa=False, n_records=0,
            master_prefix="", first_accession="",
        )

    # esummary on first hit to extract a real accession (we need the
    # 4-char prefix to identify the TSA project for Phase 1d download).
    sum_handle = Entrez.esummary(db="nuccore", id=idlist[0])
    summaries = list(Entrez.parse(sum_handle))
    try:
        sum_handle.close()  # type: ignore[union-attr]
    except Exception:
        pass

    first_acc = ""
    if summaries:
        s = summaries[0]
        # Bio.Entrez summary dicts use 'AccessionVersion' (newer XML)
        # or 'Caption' (legacy). Accept either.
        first_acc = str(s.get("AccessionVersion") or s.get("Caption") or "")

    return TsaQueryResult(
        has_tsa=True,
        n_records=n,
        master_prefix=derive_master_prefix(first_acc),
        first_accession=first_acc,
    )


# ----------------------------------------------------------------------
# Inventory TSV writer
# ----------------------------------------------------------------------

INVENTORY_COLUMNS = (
    "taxid",
    "binomial",
    "clade",
    "has_tsa",
    "n_tsa_records",
    "tsa_master_prefix",
    "tsa_first_accession",
    "query_error",
)


@dataclass(frozen=True)
class TsaInventoryEntry:
    taxon: DroppedTaxon
    result: TsaQueryResult | None
    query_error: str


def write_tsa_inventory_tsv(
    path: Union[str, Path],
    entries: list[TsaInventoryEntry],
) -> None:
    """Write tsa_inventory.tsv with INVENTORY_COLUMNS header. Rows
    sorted by taxid. has_tsa is 'yes' / 'no' / '' (empty when query
    failed).
    """
    rows = sorted(entries, key=lambda e: e.taxon.taxid)
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("\t".join(INVENTORY_COLUMNS) + "\n")
        for e in rows:
            t = e.taxon
            if e.result is None:
                row = [
                    str(t.taxid), t.binomial, t.clade,
                    "", "", "", "",
                    e.query_error,
                ]
            else:
                has = "yes" if e.result.has_tsa else "no"
                row = [
                    str(t.taxid), t.binomial, t.clade,
                    has, str(e.result.n_records),
                    e.result.master_prefix, e.result.first_accession,
                    e.query_error,
                ]
            f.write("\t".join(row) + "\n")


# ----------------------------------------------------------------------
# Orchestrator
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class TsaInventorySummary:
    n_dropped_total: int
    n_with_tsa: int
    n_without_tsa: int
    n_query_errors: int


def build_tsa_inventory(
    manifest_path: Union[str, Path],
    out_path: Union[str, Path],
    query_fn: Callable[[int], TsaQueryResult] = query_tsa_for_taxon,
    progress_every: int = 10,
    rate_limit_sleep: float = 0.0,
) -> TsaInventorySummary:
    """Walk Phase 1a's dropped taxa, query NCBI nuccore for TSA presence
    per taxid (via injectable `query_fn`), write tsa_inventory.tsv.

    `rate_limit_sleep` adds a per-call sleep — NCBI tolerates 10 req/s
    with an API key, so default 0 is fine for the live wrapper that
    overrides via the orchestrator caller.

    Query failures are recorded as `query_error` strings; the batch
    continues. Same pattern as Phase 1a.
    """
    dropped = read_dropped_taxa_from_manifest(manifest_path)
    entries: list[TsaInventoryEntry] = []
    n_query_errors = 0
    n_with_tsa = 0

    for i, taxon in enumerate(dropped, start=1):
        try:
            result = query_fn(taxon.taxid)
        except Exception as e:  # broad: any Entrez/network failure
            print(
                f"  [{i}/{len(dropped)}] taxid={taxon.taxid} {taxon.binomial}: "
                f"query failed ({e})",
                file=sys.stderr,
            )
            entries.append(TsaInventoryEntry(
                taxon=taxon, result=None, query_error=str(e),
            ))
            n_query_errors += 1
        else:
            entries.append(TsaInventoryEntry(
                taxon=taxon, result=result, query_error="",
            ))
            if result.has_tsa:
                n_with_tsa += 1

        if rate_limit_sleep:
            time.sleep(rate_limit_sleep)

        if i % progress_every == 0 or i == len(dropped):
            print(f"  [{i}/{len(dropped)}] processed", file=sys.stderr)

    write_tsa_inventory_tsv(out_path, entries)

    n_total = len(dropped)
    n_without_tsa = n_total - n_with_tsa - n_query_errors
    summary = TsaInventorySummary(
        n_dropped_total=n_total,
        n_with_tsa=n_with_tsa,
        n_without_tsa=n_without_tsa,
        n_query_errors=n_query_errors,
    )
    print(
        f"\nTSA inventory: {n_total} dropped species queried | "
        f"{n_with_tsa} with TSA | {n_without_tsa} without TSA | "
        f"{n_query_errors} query errors",
        file=sys.stderr,
    )
    return summary


# ----------------------------------------------------------------------
# CLI entrypoint
# ----------------------------------------------------------------------

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Phase 1c: NCBI TSA (Transcriptome Shotgun Assembly) inventory "
            "for the species dropped by Phase 1a (no public proteome)."
        )
    )
    p.add_argument(
        "--manifest",
        type=Path,
        default=Path("references/species_tree/proteome_manifest.tsv"),
        help="Phase 1a proteome_manifest.tsv input",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("references/species_tree/tsa_inventory.tsv"),
        help="Output TSA inventory TSV",
    )
    p.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL", ""),
        help="NCBI contact email (or env NCBI_EMAIL)",
    )
    p.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help="NCBI API key (or env NCBI_API_KEY)",
    )
    p.add_argument(
        "--rate-limit-sleep",
        type=float,
        default=0.11,
        help=(
            "Per-call sleep in seconds. Default 0.11 keeps us <10 req/s "
            "even with the API key, so we're nice to NCBI."
        ),
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    if not args.email:
        print(
            "warning: --email/NCBI_EMAIL not set; NCBI policy requires "
            "an email contact for E-utilities use",
            file=sys.stderr,
        )

    def query_fn(taxid: int) -> TsaQueryResult:
        return query_tsa_for_taxon(
            taxid, email=args.email or None, api_key=args.api_key or None,
        )

    summary = build_tsa_inventory(
        manifest_path=args.manifest,
        out_path=args.out,
        query_fn=query_fn,
        rate_limit_sleep=args.rate_limit_sleep,
    )
    return 0 if summary.n_dropped_total > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
