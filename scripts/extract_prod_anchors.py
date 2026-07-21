#!/usr/bin/env python3
"""extract_prod_anchors.py — subset the anchor sequences to one GPCR class of
anchor_set_PROD (A1 dedicated tree build).

The PLM/embedding novelty is scored against ``anchor_set_PROD`` (the
``reference_<tag>_PROD.npz`` set). There is no ``anchor_set_PROD.fasta``, and no
single FASTA holds the whole reference set: ``anchor_set_FINAL_clean.fasta``
predates the 2026-07-21 widening and holds only the anchors that survived it,
while the widening's additions were written to a sidecar. The production sources
are therefore the PAIR

  * ``anchor_set_PROD_uniprot.tsv``    keyed by ``queried_accession``
  * ``anchor_set_PROD_additions.fasta`` keyed by ``ANCHOR_<class>_<tier>_<accession>``

resolved dump-first then sidecar. This selects the class-A PROD subset so the A1
tree's phylogenetic reference frame == the PLM novelty frame.

IDs are write-once: an anchor that resolves from neither source is an ERROR
(fail loud), never a silent drop. A present-but-blank sequence counts as
unresolved — an empty record parses cleanly and would corrupt the alignment
silently.
"""
from __future__ import annotations

import csv
import os
import tempfile
from typing import Dict, Set


def composite_anchor_ids(prod_tsv: str, gpcr_class: str = "A") -> Set[str]:
    """Composite ``ANCHOR_<class>_<tier>_<accession>`` ids for one class.

    Reads ``anchor_set_PROD.tsv`` (cols: accession, tier, ..., class, ...) and
    reconstructs the composite header verbatim (mirroring build_anchor_set.py's
    ``anchor_header``) for rows whose ``class`` equals ``gpcr_class``.
    """
    out: Set[str] = set()
    with open(prod_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["class"] == gpcr_class:
                out.add(f"ANCHOR_{row['class']}_{row['tier']}_{row['accession']}")
    return out


def load_uniprot_sequences(uniprot_tsv: str) -> Dict[str, str]:
    """``queried_accession`` -> sequence, skipping blank/whitespace-only cells."""
    out: Dict[str, str] = {}
    with open(uniprot_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            seq = (row.get("Sequence") or "").strip()
            if seq:
                out[row["queried_accession"]] = seq
    return out


def load_composite_fasta(fasta_path: str) -> Dict[str, str]:
    """First-token id -> sequence. A missing file yields an empty mapping.

    The sidecar is optional: the dump alone is a valid complete source. A file
    that exists but cannot be read is still an error.
    """
    out: Dict[str, str] = {}
    if not os.path.exists(fasta_path):
        return out
    cur = None
    buf: list = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            elif cur is not None:
                buf.append(line.strip())
    if cur is not None:
        out[cur] = "".join(buf)
    return {k: v for k, v in out.items() if v.strip()}


def resolve_anchor_sequences(prod_tsv: str, uniprot_tsv: str,
                             additions_fasta: str,
                             gpcr_class: str = "A") -> Dict[str, str]:
    """Composite id -> sequence for one class, from the (dump, sidecar) pair.

    Resolution order is dump-then-sidecar, mirroring what the probe-separation
    test treats as production truth. Raises ValueError naming the anchors that
    resolve from neither source.
    """
    by_acc = load_uniprot_sequences(uniprot_tsv)
    additions = load_composite_fasta(additions_fasta)

    resolved: Dict[str, str] = {}
    unresolved: list = []
    with open(prod_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["class"] != gpcr_class:
                continue
            cid = f"ANCHOR_{row['class']}_{row['tier']}_{row['accession']}"
            seq = by_acc.get(row["accession"]) or additions.get(cid) or ""
            if seq:
                resolved[cid] = seq
            else:
                unresolved.append(cid)

    if unresolved:
        raise ValueError(
            f"{len(unresolved)} class-{gpcr_class} anchor(s) resolve from neither "
            f"{uniprot_tsv} nor {additions_fasta}: {sorted(unresolved)[:5]}"
        )
    return resolved


def _publish_atomically(out_path: str, render) -> None:
    """Write via a sibling temp file then os.replace.

    Opening the destination with "w" truncates it at open and refills
    incrementally, so an interruption leaves a SHORT BUT VALID FASTA -- and a
    short FASTA parses, so every downstream consumer accepts it. That is not
    hypothetical: the anchor set this script reads was itself left half-written
    by an interrupted run, 55 of 63 rows applied, and only a failing test
    revealed it.

    The temp lives in the destination directory because os.replace is only
    atomic within one filesystem, and is removed on any failure so a partial
    file cannot be left behind for a glob to find.
    """
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    fd, tmp = tempfile.mkstemp(dir=out_dir,
                               prefix=os.path.basename(out_path) + ".",
                               suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as out:
            render(out)
            out.flush()
            os.fsync(out.fileno())
        os.replace(tmp, out_path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def write_anchor_fasta(records: Dict[str, str], out_path: str) -> int:
    """Write ``records`` as a composite-id FASTA, published atomically."""
    def _render(out):
        for cid in sorted(records):
            out.write(f">{cid}\n{records[cid]}\n")
    _publish_atomically(out_path, _render)
    return len(records)


def subset_fasta_by_ids(source_fasta: str, ids: Set[str], out_path: str) -> int:
    """Write the source-FASTA records whose first-token id is in ``ids``.

    Raises ValueError if any requested id is absent from the source (write-once
    id integrity). Returns the number of records written.
    """
    wanted = set(ids)
    kept: Dict[str, str] = {}
    cur_id = None
    buf: list = []

    def _flush():
        if cur_id is not None and cur_id in wanted:
            kept[cur_id] = "".join(buf)

    with open(source_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                cur_id = line[1:].split()[0]
                buf = [line]
            else:
                buf.append(line)
    _flush()

    missing = wanted - set(kept)
    if missing:
        raise ValueError(
            f"{len(missing)} requested anchor id(s) missing from {source_fasta}: "
            f"{sorted(missing)[:5]}"
        )
    def _render(out):
        for cid in sorted(kept):
            out.write(kept[cid])
    _publish_atomically(out_path, _render)
    return len(kept)


def main(argv=None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--prod-tsv", required=True, help="anchor_set_PROD.tsv")
    ap.add_argument("--uniprot-tsv", required=True,
                    help="anchor_set_PROD_uniprot.tsv (keyed by queried_accession)")
    ap.add_argument("--additions-fasta", required=True,
                    help="anchor_set_PROD_additions.fasta (keyed by composite id); "
                         "may be absent if the dump alone covers every anchor")
    ap.add_argument("--class", dest="gpcr_class", default="A", help="GPCR class (default A)")
    ap.add_argument("--out", required=True, help="output subset FASTA")
    args = ap.parse_args(argv)

    n_requested = len(composite_anchor_ids(args.prod_tsv, args.gpcr_class))
    records = resolve_anchor_sequences(args.prod_tsv, args.uniprot_tsv,
                                       args.additions_fasta, args.gpcr_class)
    n = write_anchor_fasta(records, args.out)
    # Coverage, not just completion: a count that silently differs from the
    # requested set is the failure this reports.
    print(f"[extract_prod_anchors] class {args.gpcr_class}: "
          f"requested={n_requested} resolved={len(records)} written={n} -> {args.out}")
    if n != n_requested:
        raise SystemExit(
            f"[extract_prod_anchors] ERROR: wrote {n} of {n_requested} class-"
            f"{args.gpcr_class} anchors")


if __name__ == "__main__":
    main()
