#!/usr/bin/env python3
"""Build the combined 2,395-sequence FASTA for the TAAR embedding-placement job
(bead h0y0), from the four canonical source sets, with collision-free id
namespaces. Runs on the Unity checkout (paths are checkout-relative).

Sets and id conventions (see docs/plans/2026-07-22-taar-embedding-placement.md):
  anchors    (1410) references/anchors/derived/anchor_set_PROD_classA.fasta  ANCHOR_* (verbatim)
  candidates ( 790) results/chemogpcrs/chemogpcrs_berghia_classA.fa          BersteEVm* (verbatim)
  probes     ( 159) references/anchors/chemoreceptor_probe_set.fasta         PROBE_* (verbatim)
  TAAR2-9    (  36) references/non_chemo_gpcr/all_references.fasta            TAAR_<acc> (minted here)

Emits combined.fa AND membership.tsv (id<TAB>set) so the analysis never has to
guess a set from a prefix. Fails loud on any count mismatch, duplicate id, or
cross-set id collision.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

# TAAR2-9 selection = family aminergic / subfamily tyramine / gene ~ TAAR,
# MINUS the 3 TAAR1 orthologs and the 2 invertebrate tyramine receptors.
TAAR_EXCLUDE = {"Q96RJ0", "Q923Y8", "Q923Y9", "O02213", "Q19084"}
EXPECT = {"anchor": 1410, "candidate": 790, "probe": 159, "taar": 36}


def read_fasta(path: str) -> "list[tuple[str, str]]":
    if not os.path.exists(path):
        sys.exit(f"FATAL: missing FASTA {path}")
    recs, sid, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None:
                    recs.append((sid, "".join(chunks)))
                sid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if sid is not None:
            recs.append((sid, "".join(chunks)))
    return recs


def select_taars(all_ref_tsv: str, all_ref_fasta: str) -> "list[tuple[str, str]]":
    """The 36 TAAR2-9 accessions from the reference TSV, sequences from the FASTA."""
    accs = []
    with open(all_ref_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("family") == "aminergic" and row.get("subfamily") == "tyramine"
                    and "taar" in (row.get("gene") or "").lower()
                    and row["accession"] not in TAAR_EXCLUDE):
                accs.append(row["accession"])
    accs_set = set(accs)
    if len(accs) != len(accs_set):
        sys.exit(f"FATAL: duplicate TAAR accessions in {all_ref_tsv}")
    # all_references.fasta headers are `accession|family|subfamily|species` --
    # the accession is the first pipe-field, not the whole header token.
    seqs = {}
    for sid, s in read_fasta(all_ref_fasta):
        acc = sid.split("|")[0]
        if acc in accs_set:
            seqs[acc] = s
    missing = accs_set - set(seqs)
    if missing:
        sys.exit(f"FATAL: {len(missing)} TAAR accessions absent from {all_ref_fasta}: {sorted(missing)[:5]}")
    return [(f"TAAR_{acc}", seqs[acc]) for acc in accs]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="Unity checkout root (paths are relative to it)")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-membership", required=True)
    a = ap.parse_args()
    r = a.root

    sources = {
        "anchor": read_fasta(os.path.join(r, "references/anchors/derived/anchor_set_PROD_classA.fasta")),
        "candidate": read_fasta(os.path.join(r, "results/chemogpcrs/chemogpcrs_berghia_classA.fa")),
        "probe": read_fasta(os.path.join(r, "references/anchors/chemoreceptor_probe_set.fasta")),
        "taar": select_taars(
            os.path.join(r, "references/non_chemo_gpcr/all_references.tsv"),
            os.path.join(r, "references/non_chemo_gpcr/all_references.fasta")),
    }

    seen: dict[str, str] = {}
    combined: list[tuple[str, str, str]] = []  # (id, seq, set)
    for setname, recs in sources.items():
        if len(recs) != EXPECT[setname]:
            sys.exit(f"FATAL: {setname} has {len(recs)} seqs, expected {EXPECT[setname]}")
        for sid, seq in recs:
            if sid in seen:
                sys.exit(f"FATAL: id collision {sid!r} in {setname} (also in {seen[sid]})")
            if not seq or any(c not in "ACDEFGHIKLMNPQRSTVWYXBZUO*acdefghiklmnpqrstvwyxbzuo*" for c in seq):
                sys.exit(f"FATAL: empty/non-residue sequence for {sid} in {setname}")
            seen[sid] = setname
            combined.append((sid, seq.rstrip("*"), setname))

    total = len(combined)
    if total != sum(EXPECT.values()):
        sys.exit(f"FATAL: total {total} != {sum(EXPECT.values())}")

    tmp_fa, tmp_mem = a.out_fasta + ".tmp", a.out_membership + ".tmp"
    with open(tmp_fa, "w") as fa, open(tmp_mem, "w") as mem:
        mem.write("id\tset\n")
        for sid, seq, setname in combined:
            fa.write(f">{sid}\n{seq}\n")
            mem.write(f"{sid}\t{setname}\n")
    os.replace(tmp_fa, a.out_fasta)
    os.replace(tmp_mem, a.out_membership)

    counts = {s: sum(1 for _, _, ss in combined if ss == s) for s in EXPECT}
    print(f"[build_fasta] wrote {total} seqs -> {a.out_fasta}")
    print(f"[build_fasta] per-set: {counts}")
    print(f"[build_fasta] membership -> {a.out_membership}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
