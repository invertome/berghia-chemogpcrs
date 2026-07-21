#!/usr/bin/env python3
"""Pull alternate join keys from UniProt for anchors that failed the OrthoDB join.

The identifier-gap anchors are ones whose species IS in OrthoDB but whose
UniProt accession carries no UniProt cross-reference on the OrthoDB side. The
accession is not the only bridge: OrthoDB's gene_xrefs table also carries
RefSeq, Ensembl, NCBI gene ids and gene names, and its genes table carries the
original sequence id. So we ask UniProt what else each accession is known as,
and try those.

Also emits the anchor SEQUENCE and its md5, because exact sequence match is the
last-resort key that depends on no identifier agreeing at all.

Everything here is obtained by querying UniProt programmatically. No identifier
is transcribed, inferred from a pattern, or recalled.

Usage:
    python3 scripts/orthodb_fetch_alt_xrefs.py \
        --indir results/ranking/diagnostics/orthodb
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
BATCH = 40
FIELDS = "accession,id,gene_names,organism_id,sequence,xref_refseq,xref_ensembl,xref_embl,xref_geneid"


def fetch_batch(accs: list[str]) -> list[dict]:
    query = " OR ".join(f"accession:{a}" for a in accs)
    url = f"{UNIPROT}?" + urllib.parse.urlencode(
        {"query": query, "fields": FIELDS, "format": "json", "size": len(accs)}
    )
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                return json.load(fh).get("results", [])
        except Exception as exc:
            if attempt == 3:
                raise SystemExit(f"ERROR: UniProt query failed: {exc}")
            time.sleep(3 * (attempt + 1))
    return []


def extract(rec: dict) -> dict:
    acc = rec.get("primaryAccession", "")
    seq = (rec.get("sequence") or {}).get("value", "")
    genes = []
    for g in rec.get("genes", []) or []:
        nm = (g.get("geneName") or {}).get("value")
        if nm:
            genes.append(nm)
        for syn in g.get("synonyms", []) or []:
            if syn.get("value"):
                genes.append(syn["value"])
    xrefs: dict[str, list[str]] = {}
    for x in rec.get("uniProtKBCrossReferences", []) or []:
        db = x.get("database", "")
        xid = x.get("id", "")
        if db and xid:
            xrefs.setdefault(db, []).append(xid)
    org = str(((rec.get("organism") or {}).get("taxonId") or ""))
    return {
        "accession": acc,
        "uniprot_id": rec.get("uniProtkbId", ""),
        "taxid": org,
        "gene_names": ";".join(dict.fromkeys(genes)),
        "refseq": ";".join(dict.fromkeys(xrefs.get("RefSeq", []))),
        "ensembl": ";".join(dict.fromkeys(xrefs.get("Ensembl", []))),
        "embl": ";".join(dict.fromkeys(xrefs.get("EMBL", []))),
        "geneid": ";".join(dict.fromkeys(xrefs.get("GeneID", []))),
        "seq_len": str(len(seq)),
        "seq_md5": hashlib.md5(seq.encode()).hexdigest() if seq else "",
        "sequence": seq,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="results/ranking/diagnostics/orthodb")
    ap.add_argument(
        "--gap-type", default="identifier_gap",
        help="which unmapped anchors to fetch alternates for",
    )
    args = ap.parse_args(argv)
    indir = Path(args.indir)

    with open(indir / "mapping_audit.tsv", newline="") as fh:
        audit = list(csv.DictReader(fh, delimiter="\t"))
    targets = [
        r for r in audit
        if r["mapped"] != "True" and r["gap_type"] == args.gap_type
        and r["use_primary"] == "True"
    ]
    if not targets:
        raise SystemExit(f"ERROR: no unmapped anchors with gap_type={args.gap_type}")
    print(f"fetching UniProt records for {len(targets)} {args.gap_type} anchors",
          file=sys.stderr)

    accs = [r["accession"] for r in targets]
    out_rows = []
    for i in range(0, len(accs), BATCH):
        chunk = accs[i : i + BATCH]
        for rec in fetch_batch(chunk):
            out_rows.append(extract(rec))
        time.sleep(0.3)

    got = {r["accession"] for r in out_rows}
    missing = [a for a in accs if a not in got]
    if missing:
        print(f"WARNING: UniProt returned nothing for {len(missing)} accessions: "
              f"{missing[:8]}", file=sys.stderr)

    # Verify the organism UniProt reports matches the anchor table's taxid.
    # A mismatch means the anchor's recorded species is wrong, which would make
    # any recovered join meaningless.
    by_acc = {r["accession"]: r for r in targets}
    conflicts = [
        (r["accession"], r["taxid"], by_acc[r["accession"]]["taxid"])
        for r in out_rows
        if r["taxid"] and by_acc[r["accession"]]["taxid"]
        and r["taxid"] != by_acc[r["accession"]]["taxid"]
    ]
    if conflicts:
        print(f"WARNING: {len(conflicts)} accessions have a UniProt taxid differing "
              f"from the anchor table: {conflicts[:5]}", file=sys.stderr)

    cols = ["accession", "uniprot_id", "taxid", "gene_names", "refseq", "ensembl",
            "embl", "geneid", "seq_len", "seq_md5"]
    with open(indir / "identifier_gap_alt_xrefs.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(out_rows)

    with open(indir / "identifier_gap_sequences.fasta", "w") as fh:
        for r in out_rows:
            if r["sequence"]:
                fh.write(f">{r['accession']}\n{r['sequence']}\n")

    # Flat list of every candidate external id, for the OrthoDB-side scan.
    with open(indir / "identifier_gap_alt_ids.txt", "w") as fh:
        for r in out_rows:
            for field in ("refseq", "ensembl", "embl", "geneid", "gene_names"):
                for v in r[field].split(";"):
                    v = v.strip()
                    if v:
                        fh.write(f"{r['accession']}\t{field}\t{v}\n")

    counts = Counter()
    for r in out_rows:
        for field in ("refseq", "ensembl", "embl", "geneid", "gene_names"):
            if r[field]:
                counts[field] += 1
    print(json.dumps({
        "targets": len(targets),
        "uniprot_records": len(out_rows),
        "uniprot_missing": len(missing),
        "taxid_conflicts": len(conflicts),
        "anchors_with_alt_key": dict(counts),
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
