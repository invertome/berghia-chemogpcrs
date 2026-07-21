#!/usr/bin/env python3
"""extract_prod_anchors.py — subset a source anchor FASTA to one GPCR class of
anchor_set_PROD (A1 dedicated tree build).

The PLM/embedding novelty is scored against ``anchor_set_PROD`` (the
``reference_<tag>_PROD.npz`` set). There is no ``anchor_set_PROD.fasta``; the
sequences live in ``anchor_set_FINAL_clean.fasta`` keyed by the composite id
``ANCHOR_<class>_<tier>_<accession>``. This selects the class-A PROD subset (953
anchors) so the A1 tree's phylogenetic reference frame == the PLM novelty frame.

IDs are write-once: a requested anchor that is absent from the source is an
ERROR (fail loud), never a silent drop.
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
    # Publish by rename, never by writing the destination in place.
    #
    # This output is a production input: the reference npz is built from it and
    # the family prototypes are built from that. Opening the destination with
    # "w" truncates it at open and refills incrementally, so an interruption
    # leaves a SHORT BUT VALID FASTA -- and a short FASTA parses, so every
    # downstream consumer accepts it. That is not hypothetical: the anchor set
    # this script reads was itself left half-written by an interrupted run,
    # 55 of 63 rows applied, and only a failing test revealed it.
    #
    # The temp lives in the destination directory because os.replace is only
    # atomic within one filesystem, and is removed on any failure so a partial
    # file cannot be left behind for a glob to find.
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    fd, tmp = tempfile.mkstemp(dir=out_dir,
                               prefix=os.path.basename(out_path) + ".",
                               suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as out:
            for cid in sorted(kept):
                out.write(kept[cid])
            out.flush()
            os.fsync(out.fileno())
        os.replace(tmp, out_path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return len(kept)


def main(argv=None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--prod-tsv", required=True, help="anchor_set_PROD.tsv")
    ap.add_argument("--source-fasta", required=True,
                    help="FASTA holding the composite-id anchor sequences "
                         "(anchor_set_FINAL_clean.fasta)")
    ap.add_argument("--class", dest="gpcr_class", default="A", help="GPCR class (default A)")
    ap.add_argument("--out", required=True, help="output subset FASTA")
    args = ap.parse_args(argv)

    ids = composite_anchor_ids(args.prod_tsv, args.gpcr_class)
    n = subset_fasta_by_ids(args.source_fasta, ids, args.out)
    print(f"[extract_prod_anchors] class {args.gpcr_class}: {n} anchors -> {args.out}")


if __name__ == "__main__":
    main()
