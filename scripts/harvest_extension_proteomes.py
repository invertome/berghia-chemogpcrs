#!/usr/bin/env python3
"""harvest_extension_proteomes.py — close the Phase-1d extension annotation gap.

Phase 1d (build_species_tree_phase1d_extension_inventory.py) ran with
require_annotation=False and skipped Phase 1a's identify-annotated ->
download-proteome -> proteome_manifest step. This driver reinstates it for the
extension set: select the NCBI-annotated extension species, hand them to
download_species_tree_phase1a.py, and append the successfully-downloaded ones to
proteome_manifest.tsv. The `mark` subcommand then writes a `drop_reason` on
those rows in genome_inventory.tsv so the two BRAKER consumers (build_braker4_
samples_csv.py, download_species_tree_phase1f_genomes.py) skip them on
re-derive. Bead berghia-chemogpcrs-w2x.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

PROTEOME_MANIFEST_COLUMNS = (
    "taxid", "binomial", "clade", "source", "accession", "assembly_level",
    "annotation_status", "est_protein_count", "submission_date", "drop_reason",
)
DOWNLOAD_MANIFEST_COLUMNS = ("taxid", "binomial", "clade", "accession")


def annotated_accessions(datasets_jsonl_path: str) -> set:
    """Accessions whose NCBI datasets summary record carries an annotation_info object.

    Matches build_species_tree_phase1a_inventory._is_annotated: NCBI omits the
    field entirely when there is no annotation, so an annotation_info dict (even
    empty) counts as annotated.
    """
    out = set()
    with open(datasets_jsonl_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if isinstance(rec.get("annotation_info"), dict):
                acc = rec.get("accession", "")
                if acc:
                    out.add(acc)
    return out


def select_annotated_extension(extension_tsv: str, datasets_jsonl: str) -> list:
    """Extension rows whose accession is NCBI-annotated. clade <- clade_name."""
    ann = annotated_accessions(datasets_jsonl)
    out = []
    with open(extension_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if acc not in ann:
                continue
            out.append({
                "taxid": (row.get("taxid") or "").strip(),
                "binomial": (row.get("binomial") or "").strip(),
                "clade": (row.get("clade_name") or row.get("clade") or "").strip(),
                "accession": acc,
                "source": (row.get("source") or "").strip(),
                "assembly_level": (row.get("assembly_level") or "").strip(),
                "annotation_status": (row.get("annotation_status") or "").strip(),
                "est_protein_count": (row.get("est_protein_count") or "0").strip(),
                "submission_date": (row.get("submission_date") or "").strip(),
            })
    return out


def write_download_manifest(targets: list, out_tsv: str) -> None:
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(DOWNLOAD_MANIFEST_COLUMNS)
        for t in targets:
            w.writerow([t["taxid"], t["binomial"], t["clade"], t["accession"]])


def write_staged_manifest(targets: list, out_tsv: str) -> None:
    """The 10-col proteome_manifest rows for the harvested subset (drop_reason='')."""
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(PROTEOME_MANIFEST_COLUMNS),
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for t in targets:
            w.writerow({c: t.get(c, "") for c in PROTEOME_MANIFEST_COLUMNS})


def _read_status(download_results_tsv: str) -> dict:
    if not Path(download_results_tsv).exists():
        raise FileNotFoundError(f"download_results_tsv not found: {download_results_tsv!r}")
    out = {}
    with open(download_results_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row.get("taxid") or "").strip()] = (row.get("status") or "").strip()
    return out


def append_to_proteome_manifest(staged_tsv: str, download_results_tsv: str,
                                proteome_manifest_tsv: str) -> int:
    """Append staged species with a present proteome (ok/ok_no_cds/skipped — the
    last meaning already-cached) to proteome_manifest.tsv (idempotent)."""
    if not Path(proteome_manifest_tsv).exists():
        raise FileNotFoundError(
            f"proteome_manifest_tsv not found: {proteome_manifest_tsv!r} (run Phase 1a before harvest)")
    status = _read_status(download_results_tsv)
    existing = set()
    with open(proteome_manifest_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            existing.add((row.get("taxid") or "").strip())
    to_add = []
    with open(staged_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tx = (row.get("taxid") or "").strip()
            if tx in existing:
                continue
            # "skipped" = proteome already cached (present) from a prior run, so it
            # is a usable proteome and must be recorded too — keeps the harvest
            # idempotent when the cache is populated but the manifest was rebuilt.
            if status.get(tx) not in ("ok", "ok_no_cds", "skipped"):
                continue
            to_add.append({c: row.get(c, "") for c in PROTEOME_MANIFEST_COLUMNS})
    if to_add:
        with open(proteome_manifest_tsv, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(PROTEOME_MANIFEST_COLUMNS),
                               delimiter="\t", lineterminator="\n")
            for r in to_add:
                w.writerow(r)
    return len(to_add)


def mark_harvested_in_genome_inventory(
    genome_inventory_tsv: str,
    proteome_manifest_tsv: str,
    reason: str = "harvested_annotated",
) -> int:
    """Mark genome_inventory rows for species that already have a usable proteome.

    Reads the taxids of species WITH a usable proteome (empty drop_reason) from
    proteome_manifest_tsv — NOT the confirmed-no-proteome entries
    (drop_reason="no_proteome_in_ncbi"), which still need BRAKER annotation — then
    rewrites genome_inventory_tsv in place: for each row whose taxid is in that set
    and whose current drop_reason is empty, sets drop_reason to `reason`. Every
    other row and every other column is left unchanged. The header, exact column
    order, and all columns are preserved (fieldnames are read from the file).

    Idempotent: a row that already carries any non-empty drop_reason is not
    overwritten. Returns the count of rows newly marked.

    The file's existing line terminator (CRLF or LF) is detected and preserved
    so a mark run never silently rewrites every line's ending.
    """
    # Collect taxids of species with a USABLE proteome (empty drop_reason).
    # proteome_manifest also tracks species it confirmed have NO proteome
    # (drop_reason="no_proteome_in_ncbi"); those still need BRAKER de-novo
    # annotation, so they must NOT be marked/excluded from the target set.
    manifest_taxids: set[str] = set()
    with open(proteome_manifest_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tx = (row.get("taxid") or "").strip()
            drop = (row.get("drop_reason") or "").strip()
            if tx and not drop:
                manifest_taxids.add(tx)

    # Detect the existing line terminator so the rewrite preserves it.
    with open(genome_inventory_tsv, "rb") as fh:
        linesep = "\r\n" if b"\r\n" in fh.read() else "\n"

    # Read genome_inventory (newline="" lets csv strip the \r from field values).
    with open(genome_inventory_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    n_marked = 0
    for row in rows:
        tx = (row.get("taxid") or "").strip()
        current_reason = (row.get("drop_reason") or "").strip()
        # For each matching row with no existing drop_reason, set it to `reason`.
        if tx in manifest_taxids and not current_reason:
            row["drop_reason"] = reason
            n_marked += 1

    # Rewrite in place, preserving the original line terminator.
    with open(genome_inventory_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames,
                           delimiter="\t", lineterminator=linesep)
        w.writeheader()
        w.writerows(rows)

    return n_marked


def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "harvest_extension_proteomes — close the Phase-1d extension annotation gap. "
            "Select NCBI-annotated extension species and append downloaded proteomes "
            "to proteome_manifest.tsv."
        )
    )
    sub = p.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    sel = sub.add_parser(
        "select",
        help="Select annotated extension species and write download + staged manifests.",
    )
    sel.add_argument("--extension-tsv", required=True,
                     help="Phase 1d extension inventory TSV.")
    sel.add_argument("--datasets-jsonl", required=True,
                     help="NCBI datasets summary JSONL for extension accessions.")
    sel.add_argument("--download-manifest-out", required=True,
                     help="Output path for the 4-col download manifest TSV.")
    sel.add_argument("--staged-out", required=True,
                     help="Output path for the 10-col staged proteome manifest TSV.")

    app = sub.add_parser(
        "append",
        help="Append successfully-downloaded extension species to proteome_manifest.tsv.",
    )
    app.add_argument("--staged", required=True,
                     help="Staged manifest TSV produced by the select subcommand.")
    app.add_argument("--download-results", required=True,
                     help="Per-species download report TSV (taxid, status, …).")
    app.add_argument("--proteome-manifest", required=True,
                     help="proteome_manifest.tsv to append to.")

    mrk = sub.add_parser(
        "mark",
        help=(
            "Mark harvested species in genome_inventory.tsv with a drop_reason so "
            "BRAKER consumers skip them on re-derive."
        ),
    )
    mrk.add_argument("--genome-inventory", required=True,
                     help="genome_inventory.tsv to update in place.")
    mrk.add_argument("--proteome-manifest", required=True,
                     help="proteome_manifest.tsv whose taxids will be marked.")

    return p


def main(argv=None) -> None:
    args = _build_argparser().parse_args(argv)

    if args.command == "select":
        for flag, path in (
            ("--extension-tsv", args.extension_tsv),
            ("--datasets-jsonl", args.datasets_jsonl),
        ):
            if not Path(path).exists():
                print(f"ERROR: {flag} path does not exist: {path!r}", file=sys.stderr)
                sys.exit(2)
        targets = select_annotated_extension(args.extension_tsv, args.datasets_jsonl)
        write_download_manifest(targets, args.download_manifest_out)
        write_staged_manifest(targets, args.staged_out)
        print(f"[harvest] selected {len(targets)} annotated extension species",
              file=sys.stderr)

    elif args.command == "append":
        for flag, path in (
            ("--staged", args.staged),
            ("--download-results", args.download_results),
            ("--proteome-manifest", args.proteome_manifest),
        ):
            if not Path(path).exists():
                print(f"ERROR: {flag} path does not exist: {path!r}", file=sys.stderr)
                sys.exit(2)
        n = append_to_proteome_manifest(
            args.staged, args.download_results, args.proteome_manifest
        )
        print(f"[harvest] appended {n} species to {args.proteome_manifest}",
              file=sys.stderr)

    elif args.command == "mark":
        for flag, path in (
            ("--genome-inventory", args.genome_inventory),
            ("--proteome-manifest", args.proteome_manifest),
        ):
            if not Path(path).exists():
                print(f"ERROR: {flag} path does not exist: {path!r}", file=sys.stderr)
                sys.exit(2)
        n = mark_harvested_in_genome_inventory(
            args.genome_inventory, args.proteome_manifest
        )
        print(f"[harvest] marked {n} rows in {args.genome_inventory}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
