#!/usr/bin/env python3
"""harvest_extension_proteomes.py — close the Phase-1d extension annotation gap.

Phase 1d (build_species_tree_phase1d_extension_inventory.py) ran with
require_annotation=False and skipped Phase 1a's identify-annotated ->
download-proteome -> proteome_manifest step. This driver reinstates it for the
extension set: select the NCBI-annotated extension species, hand them to
download_species_tree_phase1a.py, and append the successfully-downloaded ones to
proteome_manifest.tsv. The existing extension-inventory disjointness logic then
drops them from the BRAKER target set on re-derive. Bead berghia-chemogpcrs-w2x.

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations

import csv
import json
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
        w = csv.writer(fh, delimiter="\t")
        w.writerow(DOWNLOAD_MANIFEST_COLUMNS)
        for t in targets:
            w.writerow([t["taxid"], t["binomial"], t["clade"], t["accession"]])


def write_staged_manifest(targets: list, out_tsv: str) -> None:
    """The 10-col proteome_manifest rows for the harvested subset (drop_reason='')."""
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(PROTEOME_MANIFEST_COLUMNS), delimiter="\t")
        w.writeheader()
        for t in targets:
            w.writerow({c: t.get(c, "") for c in PROTEOME_MANIFEST_COLUMNS})
