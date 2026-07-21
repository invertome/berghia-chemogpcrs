#!/usr/bin/env python3
"""Reconcile ``anchor_set_PROD.tsv`` against the UniProt records it claims.

WHY
---
The anchor set's ``class`` column decides which anchors enter the class-A
envelope (``extract_prod_anchors`` selects ``class == "A"``), and that envelope
is the reference frame for both the A1 tree and the PLM novelty prototypes. Two
defects were measured in it:

* 54 rows carried ``class=A`` while UniProt curates them as GPCR family 2
  (class B, secretin/adhesion) or family 3 (class C, glutamate/GABA-B), plus one
  row (GPR143/OA1) whose only Pfam is PF02101 and whose curated family is the
  bespoke "OA family". None of the 55 carries PF00001. Sitting under a class-A
  family label, they silently contaminate that family's prototype.
* 26 of those rows also contradicted themselves, carrying ``class=A`` beside a
  ``family`` of ``class-B-secretin`` (17) or ``class-C`` (9).

Neither column was authoritative. The underlying UniProt record is, so both are
re-derived from it via ``curate_gpcr_references.gpcr_class_from_evidence``.

POLICY (uniform, applied to every row, no case-by-case judgement)
-----------------------------------------------------------------
1. Resolve each accession's class from its own UniProt record: the curated
   "Belongs to the G-protein coupled receptor N family" statement first, the
   specific Pfam family (PF00001/2/3) as fallback, else ``UNKNOWN``.
2. A row whose declared ``class`` is A but whose resolved class is not A is
   EVICTED from the anchor set. It is not relabelled, because the FASTA/npz key
   is ``ANCHOR_<class>_<tier>_<accession>`` -- rewriting ``class`` would re-mint
   that key, and identifiers here are write-once. Eviction leaves every retained
   row's key byte-identical.
3. ``UNKNOWN`` is evicted too. There is no default-to-A.
4. Rows outside the envelope under repair are left untouched.
5. Evicted rows are not destroyed: their sequences remain in the superset store
   ``anchor_set_FINAL_clean.fasta`` under their original ids, and every eviction
   is recorded with its evidence in the sidecar written by ``--evictions-out``.

The repair also backfills the ``taxid`` column (723 of 1094 rows were blank, so
every taxonomic count derived from it undercounted) and applies verified family
corrections. Neither changes any key.

USAGE
-----
    python3 scripts/anchorfix_class_reconcile.py \\
        --prod-tsv references/anchors/anchor_set_PROD.tsv \\
        --uniprot-cache references/anchors/anchor_set_PROD_uniprot.tsv \\
        --source-fasta references/anchors/anchor_set_FINAL_clean.fasta \\
        --evictions-out references/anchors/anchor_set_PROD_evictions.tsv \\
        --apply
"""
from __future__ import annotations

import csv
import os
import sys
from typing import Dict, Iterable, List, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from curate_gpcr_references import gpcr_class_from_evidence  # noqa: E402

PROD_COLUMNS = ["accession", "tier", "taxid", "species", "family", "class",
                "evidence"]

# Family corrections verified against UniProt's curated subfamily, each as
# accession -> (expected current family, corrected family). The expected current
# value is asserted so a correction cannot fire against a drifted file.
#
# All three are curated "Vasopressin/oxytocin receptor subfamily", which is a
# peptide-receptor subfamily. The two Platynereis entries were filed as
# glycoprotein-hormone through a substring collision on "gonadotropin" in their
# protein name; glycoprotein hormones are FSH/LH/TSH/CG, not GnRH.
#
# P46023 (GRL101, Lymnaea) is deliberately ABSENT: UniProt gives it the bare
# class-A family with no subfamily, so any promotion out of `orphan` would be
# inference rather than curation. It stays where it is.
# A0A0K0PUL2 is the same failure mode again, this time on "DH31": it was filed
# as class-B-secretin from its "DH31-like receptor 2" name, but InterPro matches
# it to IPR000276 (rhodopsin-like), IPR000405 (galanin receptor family),
# PF00001, PS50262 (family 1 profile) and PANTHER PTHR45695 (GPCR 1) with ZERO
# class-B signatures -- it is a class-A galanin-family peptide receptor. Its
# similarly named sibling A0A0K0PVL8 ("DH31-like receptor 1") really is class B
# (PF00002, PF02793, IPR003287 calcitonin receptor) and is evicted instead. All
# nine other galanin-family members in the set are labelled `peptide`.
VERIFIED_FAMILY_CORRECTIONS: Dict[str, Tuple[str, str]] = {
    "A0A6B9MRA0": ("glycoprotein-hormone", "peptide"),
    "A0A6B9MSD4": ("glycoprotein-hormone", "peptide"),
    "Q75W84": ("orphan", "peptide"),
    "A0A0K0PUL2": ("class-B-secretin", "peptide"),
}


# --- pure helpers -------------------------------------------------------------

def composite_id(row: dict) -> str:
    """The anchor's FASTA/npz key: ``ANCHOR_<class>_<tier>_<accession>``.

    Mirrors ``build_anchor_set.anchor_header``. Any edit that changes this
    string re-mints an identifier and is forbidden.
    """
    return f"ANCHOR_{row['class']}_{row['tier']}_{row['accession']}"


def partition_by_resolved_class(
    rows: Sequence[dict], resolved: Dict[str, str], envelope: str = "A",
) -> Tuple[List[dict], List[dict]]:
    """Split ``rows`` into (retained, evicted) for one class envelope.

    A row is evicted iff it declares ``class == envelope`` but its resolved
    class differs. Rows declaring another class pass through untouched. Order is
    preserved in both outputs. A row with no resolved class raises: a missing
    lookup is a broken join, never an implicit keep.
    """
    keep: List[dict] = []
    evict: List[dict] = []
    for row in rows:
        acc = row["accession"]
        if acc not in resolved:
            raise ValueError(
                f"no resolved class for accession {acc!r}: the UniProt cache "
                f"does not cover the anchor set")
        if row["class"] != envelope:
            keep.append(row)
        elif resolved[acc] == envelope:
            keep.append(row)
        else:
            evict.append(row)
    return keep, evict


def fill_taxids(rows: Iterable[dict], taxids: Dict[str, str]) -> int:
    """Backfill blank ``taxid`` cells from the source. Returns the fill count.

    A pre-existing taxid that contradicts the source raises rather than being
    silently overwritten or silently kept -- that disagreement means one of the
    two is wrong about which organism the row describes.
    """
    filled = 0
    for row in rows:
        acc = row["accession"]
        src = str(taxids.get(acc, "")).strip()
        if not src:
            raise ValueError(f"no source taxid for accession {acc!r}")
        cur = str(row.get("taxid", "")).strip()
        if not cur:
            row["taxid"] = src
            filled += 1
        elif cur != src:
            raise ValueError(
                f"taxid disagreement for {acc!r}: file has {cur!r}, "
                f"source has {src!r}")
    return filled


def apply_family_corrections(
    rows: Sequence[dict], corrections: Dict[str, Tuple[str, str]],
) -> int:
    """Apply verified ``accession -> (expected_family, new_family)`` moves.

    Returns the number of rows actually changed.

    IDEMPOTENT: a row already holding ``new_family`` is a no-op, so re-running
    the repair against an already-repaired file converges instead of raising.
    Any OTHER value still raises -- that means the file drifted to something
    neither expected nor already-corrected, and re-verification is required
    before overwriting it. Never touches ``class``, so no composite id moves.
    """
    by_acc = {r["accession"]: r for r in rows}
    applied = 0
    for acc, (expected, new_family) in corrections.items():
        row = by_acc.get(acc)
        if row is None:
            raise ValueError(f"correction target {acc!r} is not in the table")
        current = row["family"]
        if current == new_family:
            continue
        if current != expected:
            raise ValueError(
                f"correction target {acc!r} has family {current!r}, expected "
                f"{expected!r} (or the corrected {new_family!r}): the file "
                f"drifted, re-verify first")
        row["family"] = new_family
        applied += 1
    return applied


def write_table_atomic(path: str, rows: Sequence[dict],
                       columns: Sequence[str]) -> None:
    """Serialise to ``path`` via a temp file + ``os.replace``.

    A serialisation failure leaves the original file untouched.
    """
    tmp = f"{path}.tmp"
    try:
        with open(tmp, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(columns),
                                    delimiter="\t", extrasaction="raise")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def verify_written_table(path: str, expected_accessions: Sequence[str],
                         key_column: str = "accession") -> None:
    """Re-read ``path`` and assert count, uniqueness, and exact order."""
    back = list(csv.DictReader(open(path, newline=""), delimiter="\t"))
    if back and key_column not in back[0]:
        raise ValueError(
            f"key column {key_column!r} absent from {path}: "
            f"got {list(back[0])[:5]}")
    got = [r[key_column] for r in back]
    if len(got) != len(expected_accessions):
        raise ValueError(
            f"row count mismatch after write: {len(got)} on disk, "
            f"{len(expected_accessions)} expected")
    if len(set(got)) != len(got):
        dupes = sorted({a for a in got if got.count(a) > 1})
        raise ValueError(f"duplicate accessions after write: {dupes[:5]}")
    if got != list(expected_accessions):
        first = next(i for i, (a, b) in enumerate(zip(got, expected_accessions))
                     if a != b)
        raise ValueError(
            f"row order changed after write at index {first}: "
            f"{got[first]!r} != {expected_accessions[first]!r}")


def verify_sequence_integrity(rows: Sequence[dict], fasta: Dict[str, str],
                              sequences: Dict[str, str]) -> None:
    """Assert every retained row's stored sequence still matches its source.

    Guards the accumulate-by-accession failure mode, in which a merge doubles a
    sequence while every row count still reconciles. Compares the full string,
    not the length, because a doubled sequence changes both but a transposition
    changes only the string.
    """
    for row in rows:
        cid = composite_id(row)
        acc = row["accession"]
        if cid not in fasta:
            raise ValueError(f"retained anchor {cid!r} is absent from the "
                             f"sequence store: its key was re-minted")
        src = sequences.get(acc)
        if src and fasta[cid] != src:
            raise ValueError(
                f"sequence integrity failure for {acc!r}: store has "
                f"{len(fasta[cid])} aa, source has {len(src)} aa")


# --- io -----------------------------------------------------------------------

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = ("accession,protein_name,organism_name,organism_id,length,"
                  "cc_similarity,reviewed,xref_pfam,lineage_ids,sequence")
CACHE_COLUMNS = ["queried_accession", "Entry", "Protein names", "Organism",
                 "Organism (ID)", "Length", "Sequence similarities",
                 "Reviewed", "Pfam", "Taxonomic lineage (Ids)", "Sequence"]


def fetch_uniprot_cache(accessions: Sequence[str], out_path: str,
                        batch: int = 90) -> None:
    """Refresh the UniProt cache for ``accessions`` from the REST API.

    Every field in the cache comes verbatim from the API; nothing is recalled or
    hand-entered. An accession the API does not return is a hard error, because
    a silently absent record would drop a row from the class resolution and the
    partition would then treat it as unresolvable.
    """
    import time
    import urllib.error
    import urllib.parse
    import urllib.request

    found: Dict[str, dict] = {}
    for i in range(0, len(accessions), batch):
        chunk = accessions[i:i + batch]
        params = urllib.parse.urlencode({
            "query": " OR ".join(f"accession:{a}" for a in chunk),
            "fields": UNIPROT_FIELDS, "format": "tsv", "size": "500"})
        url = f"{UNIPROT_SEARCH}?{params}"
        for attempt in range(5):
            try:
                text = urllib.request.urlopen(url, timeout=120).read().decode()
                break
            except (urllib.error.URLError, OSError):
                if attempt == 4:
                    raise
                time.sleep(2 * (attempt + 1))
        for rec in csv.DictReader(text.splitlines(), delimiter="\t"):
            found[rec["Entry"]] = rec
        time.sleep(0.3)

    missing = [a for a in accessions if a not in found]
    if missing:
        raise ValueError(f"UniProt did not return {len(missing)} accession(s): "
                         f"{missing[:5]}")

    rows = [dict({"queried_accession": a},
                 **{c: found[a].get(c, "") for c in CACHE_COLUMNS[1:]})
            for a in accessions]
    write_table_atomic(out_path, rows, CACHE_COLUMNS)
    verify_written_table(out_path, list(accessions),
                         key_column="queried_accession")


def read_fasta(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    cur = None
    buf: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        out[cur] = "".join(buf)
    return out


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--prod-tsv", required=True)
    ap.add_argument("--uniprot-cache", required=True,
                    help="TSV from the UniProt REST API: queried_accession, "
                         "'Sequence similarities', Pfam, 'Organism (ID)', Sequence")
    ap.add_argument("--source-fasta", required=True,
                    help="anchor_set_FINAL_clean.fasta (composite-id store)")
    ap.add_argument("--evictions-out", required=True)
    ap.add_argument("--envelope", default="A")
    ap.add_argument("--apply", action="store_true",
                    help="write changes; otherwise report only")
    ap.add_argument("--fetch-cache", action="store_true",
                    help="refresh --uniprot-cache from the REST API first")
    args = ap.parse_args(argv)

    rows = list(csv.DictReader(open(args.prod_tsv, newline=""), delimiter="\t"))
    if args.fetch_cache:
        fetch_uniprot_cache([r["accession"] for r in rows], args.uniprot_cache)
        print(f"[anchorfix] refreshed {args.uniprot_cache} "
              f"({len(rows)} accessions)")
    uni = {r["queried_accession"]: r
           for r in csv.DictReader(open(args.uniprot_cache, newline=""),
                                   delimiter="\t")}
    fasta = read_fasta(args.source_fasta)

    resolved = {a: gpcr_class_from_evidence(r["Sequence similarities"], r["Pfam"])
                for a, r in uni.items()}
    taxids = {a: r["Organism (ID)"] for a, r in uni.items()}
    sequences = {a: r.get("Sequence", "") for a, r in uni.items()}

    keep, evict = partition_by_resolved_class(rows, resolved, args.envelope)
    filled = fill_taxids(keep, taxids)
    corrected = apply_family_corrections(keep, VERIFIED_FAMILY_CORRECTIONS)
    verify_sequence_integrity(keep, fasta, sequences)

    print(f"[anchorfix] {len(rows)} rows in -> {len(keep)} retained, "
          f"{len(evict)} evicted from envelope {args.envelope}")
    print(f"[anchorfix] taxids backfilled: {filled}")
    print(f"[anchorfix] family corrections applied: {corrected}")

    if not args.apply:
        print("[anchorfix] dry run; pass --apply to write")
        return 0

    write_table_atomic(args.prod_tsv, keep, PROD_COLUMNS)
    verify_written_table(args.prod_tsv, [r["accession"] for r in keep])

    ev_cols = ["accession", "prior_class", "prior_family", "tier",
               "composite_id", "resolved_class", "curated_similarity", "pfam",
               "organism", "taxid", "reason"]
    ev_rows = []
    for r in evict:
        u = uni[r["accession"]]
        ev_rows.append({
            "accession": r["accession"], "prior_class": r["class"],
            "prior_family": r["family"], "tier": r["tier"],
            "composite_id": composite_id(r),
            "resolved_class": resolved[r["accession"]],
            "curated_similarity": u["Sequence similarities"],
            "pfam": u["Pfam"], "organism": u["Organism"],
            "taxid": u["Organism (ID)"],
            "reason": f"resolved class {resolved[r['accession']]} != envelope "
                      f"{args.envelope}",
        })
    write_table_atomic(args.evictions_out, ev_rows, ev_cols)
    verify_written_table(args.evictions_out, [r["accession"] for r in evict])
    print(f"[anchorfix] wrote {args.prod_tsv} and {args.evictions_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
