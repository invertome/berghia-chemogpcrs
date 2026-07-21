#!/usr/bin/env python3
"""Resolve the class-A anchors that carry no Pfam 7tm_1 (PF00001) hit.

A 7tm_1 miss is evidence of one of two very different things, and the whole
point of this script is that it refuses to guess which:

  GENUINELY DIVERGENT CLASS A -- ACKR1/DARC is the textbook case. UniProt
      curates it into the "G-protein coupled receptor 1 family", it is a real
      class-A receptor, and it simply does not score against the PF00001 HMM.
      Evicting it would remove a legitimate, informative corner of the class-A
      envelope.

  CLASS-B/C CONTAMINATION -- a receptor whose curated family statement says 2
      or 3. It has no business seeding a class-A family prototype.

The discriminator is ``curate_gpcr_references.gpcr_class_from_evidence``, in
which UniProt's own curated family statement is authoritative IN BOTH
DIRECTIONS and the Pfam signature is consulted only where curation names no
class. That matters here: an earlier pass used a domain-only test, which by
construction resolves every 7tm_1 miss as "not class A" and therefore erred in
both directions at once.

Every accession is re-queried LIVE. Cached values are not consulted, because
the cache was written by the same pass whose criterion is under review.

Outcome handling is EVICTION, never relabelling. The composite identifier
``ANCHOR_<class>_<tier>_<accession>`` embeds the class and is the join key into
the sequence store and the embedding npz, so relabelling an entry re-mints an
identifier that other artifacts already reference. Identifiers are write-once.

Usage:
    python3 scripts/orthodb_resolve_no7tm1.py \
        --ids <file of composite ids or accessions> \
        --anchors references/anchors/anchor_set_PROD.tsv \
        --out results/ranking/diagnostics/orthodb/no7tm1_resolution.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
from curate_gpcr_references import gpcr_class_from_evidence  # noqa: E402

STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = [
    "accession", "reviewed", "protein_name", "organism_name", "organism_id",
    "length", "cc_similarity", "xref_pfam", "protein_existence", "sequence",
]
CHUNK = 40

# UniProt's TSV header labels are NOT the field names requested. cc_similarity
# comes back as "Sequence similarities". Reading it under the wrong label
# yields an empty string for every record, which silently demotes the resolver
# to a Pfam-only test -- precisely the domain-only criterion this script exists
# to replace, and a failure that produces confident wrong evictions rather than
# an error. The label is therefore asserted against the live header.
COLUMN_LABELS = {
    "accession": "Entry",
    "curated_similarity": "Sequence similarities",
    "pfam": "Pfam",
    "reviewed": "Reviewed",
    "protein_name": "Protein names",
    "organism": "Organism",
    "taxid": "Organism (ID)",
    "length": "Length",
}


def fetch(accessions: list[str], max_retries: int = 4) -> dict[str, dict]:
    """Live UniProt fetch of the evidence fields, keyed by accession.

    Records are keyed by the accession UniProt RETURNS, and the caller checks
    for absences. A secondary accession that UniProt redirects to a primary one
    therefore shows up as a miss rather than silently binding to the wrong row.
    """
    out: dict[str, dict] = {}
    for i in range(0, len(accessions), CHUNK):
        chunk = accessions[i:i + CHUNK]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        url = f"{STREAM_URL}?{urlencode({'query': query, 'format': 'tsv', 'fields': ','.join(FIELDS)})}"
        for attempt in range(max_retries):
            try:
                req = Request(url, headers={"User-Agent": "berghia-chemogpcrs/1.0"})
                with urlopen(req, timeout=120) as fh:
                    text = fh.read().decode("utf-8")
                break
            except Exception as exc:
                if attempt == max_retries - 1:
                    raise SystemExit(f"ERROR: UniProt fetch failed: {exc}")
                time.sleep(2 ** attempt)
        reader = csv.DictReader(text.splitlines(), delimiter="\t")
        header = reader.fieldnames or []
        absent = [lbl for lbl in COLUMN_LABELS.values() if lbl not in header]
        if absent:
            raise SystemExit(
                f"ERROR: UniProt returned no column(s) {absent}. Header was "
                f"{header}. Refusing to resolve GPCR class against columns that "
                "did not arrive -- an absent curated-family column silently "
                "degrades this to a Pfam-only test."
            )
        for r in reader:
            out[r[COLUMN_LABELS["accession"]]] = r
        print(f"  fetched {len(out)}/{len(accessions)}", file=sys.stderr)
        time.sleep(0.3)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", required=True,
                    help="file of composite ids (ANCHOR_<class>_<tier>_<acc>) "
                         "or bare accessions, one per line")
    ap.add_argument("--anchors", default="references/anchors/anchor_set_PROD.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    raw = [l.strip() for l in open(args.ids) if l.strip()]
    acc_of = {}
    for r in raw:
        acc = r.rsplit("_", 1)[1] if r.startswith("ANCHOR_") else r
        acc_of[acc] = r

    with open(args.anchors, newline="") as fh:
        anchors = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

    accessions = sorted(acc_of)
    print(f"resolving {len(accessions)} accessions live against UniProt ...",
          file=sys.stderr)
    recs = fetch(accessions)

    missing = [a for a in accessions if a not in recs]
    if missing:
        print(f"WARNING: {len(missing)} accessions returned no record "
              f"(e.g. {missing[:5]})", file=sys.stderr)

    rows = []
    for acc in accessions:
        rec = recs.get(acc, {})
        sim = rec.get(COLUMN_LABELS["curated_similarity"], "") or ""
        pfam = rec.get(COLUMN_LABELS["pfam"], "") or ""
        resolved = gpcr_class_from_evidence(sim, pfam) if rec else "NO_RECORD"
        in_prod = acc in anchors
        prod_class = anchors[acc]["class"] if in_prod else ""
        if not rec:
            verdict = "UNRESOLVED_no_record"
        elif resolved == "A":
            verdict = "KEEP_divergent_class_A"
        elif resolved in ("B", "C"):
            verdict = f"EVICT_curated_class_{resolved}"
        else:
            verdict = "EVICT_class_UNKNOWN"
        rows.append({
            "accession": acc,
            "composite_id": acc_of[acc],
            "currently_in_prod": str(in_prod),
            "prod_class": prod_class,
            "prod_family": anchors[acc]["family"] if in_prod else "",
            "prod_tier": anchors[acc]["tier"] if in_prod else "",
            "resolved_class": resolved,
            "verdict": verdict,
            "reviewed": rec.get(COLUMN_LABELS["reviewed"], ""),
            "protein_name": rec.get(COLUMN_LABELS["protein_name"], ""),
            "organism": rec.get(COLUMN_LABELS["organism"], ""),
            "taxid": rec.get(COLUMN_LABELS["taxid"], ""),
            "length": rec.get(COLUMN_LABELS["length"], ""),
            "curated_similarity": sim,
            "pfam": pfam,
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(out)

    live = [r for r in rows if r["currently_in_prod"] == "True"]
    print(f"\n=== {len(rows)} accessions resolved; {len(live)} still in the anchor set ===")
    print("verdicts (ALL):", dict(Counter(r["verdict"] for r in rows)))
    print("verdicts (still in anchor set):",
          dict(Counter(r["verdict"] for r in live)))
    print(f"\n=== the {len(live)} live entries ===")
    for r in live:
        print(f"  {r['accession']:<10} {r['verdict']:<26} "
              f"family={r['prod_family']:<12} {r['protein_name'][:44]}")
        print(f"{'':<13}curated: {r['curated_similarity'][:100] or '(none)'}")
    print(f"\nwrote {out}")

    summary = {
        "n_resolved": len(rows),
        "n_still_in_anchor_set": len(live),
        "verdicts_all": dict(Counter(r["verdict"] for r in rows)),
        "verdicts_live": dict(Counter(r["verdict"] for r in live)),
        "evict_now": [r["accession"] for r in live if r["verdict"].startswith("EVICT")],
        "keep": [r["accession"] for r in live if r["verdict"].startswith("KEEP")],
    }
    out.with_name(out.stem + "_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
