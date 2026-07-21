#!/usr/bin/env python3
"""Merge the OrthoDB harvest and the literature-gated sweep into the anchor set.

Two independently-derived reference tiers are appended, each labelled so it can
be reverted ALONE:

  tier 9  ``orthodb-harvest``   -- taxonomic-breadth harvest from anchor-seeded
                                  OrthoDB orthogroups, gated at K>=1 seeds with
                                  strict single-family purity, verified class A
                                  per entry against Pfam 7tm_1 and filtered
                                  against the anchor population's own fragment
                                  thresholds.
  tier 10 ``literature-gated``  -- UniProt entries surviving the refexp2
                                  evidence gate: class A per entry, PE1/PE2,
                                  experimental ECO code, primary literature, and
                                  a name that commits to a function.

IDENTIFIERS ARE WRITE-ONCE. Every pre-existing row is copied through byte for
byte, in its original order. Nothing is renumbered, reordered or re-minted, and
gaps are correct. New rows take ``ANCHOR_<class>_<tier>_<accession>``, the same
composite the sequence store and the embedding npz are keyed on, so the class
must be right AT INSERT -- a later reclassification would re-mint a live join
key, which is why the earlier repair evicted rather than relabelled.

The OrthoDB tier is keyed on OrthoDB gene ids, whose colon is replaced by an
underscore. That is not cosmetic: anchors are consumed by Newick tree builders,
in which a colon terminates a taxon label and would silently corrupt every
downstream tree. The untransformed id is retained in the provenance sidecar.

The anchor table's seven-column schema is NOT extended. Readers that parse it
positionally would break, and the tier value alone already separates the tiers;
full per-entry provenance goes to a sidecar instead.

Usage:
    python3 scripts/orthodb_merge_reference_tiers.py \
        --anchors references/anchors/anchor_set_PROD.tsv \
        --harvest results/ranking/diagnostics/orthodb/harvest_final.tsv \
        --literature references/non_chemo_gpcr/refexp2_candidate_verdicts.tsv \
        --literature-fasta references/non_chemo_gpcr/refexp2_candidate_verdicts_passing.fasta \
        --outdir references/anchors
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from anchorfix_class_reconcile import PROD_COLUMNS, composite_id  # noqa: E402

TIER_ORTHODB = "9"
TIER_LITERATURE = "10"
TIER_LABELS = {TIER_ORTHODB: "orthodb-harvest", TIER_LITERATURE: "literature-gated"}

PROVENANCE_COLUMNS = [
    "composite_id", "accession", "tier", "tier_label", "class", "family",
    "taxid", "species", "source_db", "source_id", "evidence_basis",
    "sequence_length", "added_utc",
]


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def snapshot(paths: list[Path], outdir: Path) -> Path:
    """Copy the files about to change into a timestamped directory.

    references/ is gitignored, so there is no version history to recover from
    and an interrupted write has already left this file half-repaired once.
    """
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    snap = outdir / f".merge_snapshot_{stamp}"
    snap.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.exists():
            shutil.copy2(p, snap / p.name)
    return snap


def write_atomic(path: Path, rows: list[dict], columns: list[str]) -> None:
    """Temp file + os.replace, so a killed process cannot truncate the target."""
    if not rows:
        raise SystemExit(f"ERROR: refusing to write an empty table to {path}")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def write_fasta_atomic(path: Path, records: list[tuple[str, str]]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
    os.replace(tmp, path)


def build_orthodb_rows(harvest: list[dict], now: str) -> tuple[list[dict], list[dict], list[tuple[str, str]]]:
    """Rows, provenance and sequences for the OrthoDB tier.

    Only quality-PASS entries are admitted. Class is 'A' because every admitted
    entry carries its OWN verified 7tm_1 hit -- it is not inherited from the
    seeding orthogroup, which is the contamination route this tier is most
    exposed to.
    """
    rows, prov, seqs = [], [], []
    for r in harvest:
        if r.get("quality_verdict") != "PASS":
            continue
        if r.get("verdict") != "verified_class_A":
            raise SystemExit(
                f"ERROR: {r['gene_id']} is quality-PASS but its class-A verdict is "
                f"{r.get('verdict')!r}. Refusing to admit an entry whose class was "
                "not verified per entry."
            )
        acc = r["gene_id"].replace(":", "_")
        row = {
            "accession": acc,
            "tier": TIER_ORTHODB,
            "taxid": r["taxid"],
            "species": r["species"],
            "family": r["family"],
            "class": "A",
            "evidence": f"orthodb-harvest:og={r['og_id']}:seeds={r['og_seeds']}",
        }
        rows.append(row)
        seq = r.get("sequence", "")
        prov.append({
            "composite_id": composite_id(row),
            "accession": acc,
            "tier": TIER_ORTHODB,
            "tier_label": TIER_LABELS[TIER_ORTHODB],
            "class": "A",
            "family": r["family"],
            "taxid": r["taxid"],
            "species": r["species"],
            "source_db": "OrthoDB",
            "source_id": r["gene_id"],
            "evidence_basis": (
                f"og={r['og_id']}({r['og_name']}); seeds={r['og_seeds']}; "
                f"members={r['og_members']}; 7tm1_E={r['pf00001_evalue']}; "
                f"7tm1_cov={r['tm7_coverage']}; stratum={r.get('stratum', '')}"
            ),
            "sequence_length": str(len(seq)),
            "added_utc": now,
        })
        seqs.append((composite_id(row), seq))
    return rows, prov, seqs


def build_literature_rows(verdicts: list[dict], fasta: dict[str, str], now: str,
                          hold_families: frozenset[str] = frozenset()
                          ) -> tuple[list[dict], list[dict], list[tuple[str, str]], list[dict]]:
    """Rows, provenance and sequences for the literature-gated tier.

    Class comes from the entry's own class-A verification, recorded with its
    basis so an entry resting on a signature alone is distinguishable from one
    resting on curation.

    Entries whose family label is in ``hold_families`` are NOT merged. The
    family column is not decoration: it defines a novelty prototype, S_novel is
    a MIN over prototypes, so an extra prototype can only pull novelty DOWN.
    A label meaning "the reporting regex matched nothing" is not a family, and a
    centroid over its members would be a grab bag that suppresses novelty
    everywhere. Assigning them real families is a scientific call about what the
    production reference set contains, which refexp2 itself declines to make, so
    they are held for ratification rather than guessed at here.
    """
    rows, prov, seqs, held = [], [], [], []
    by_acc = {k.split("|")[0]: v for k, v in fasta.items()}
    for r in verdicts:
        if r.get("failed_at"):
            continue
        if r["family"] in hold_families:
            held.append(r)
            continue
        acc = r["accession"]
        seq = by_acc.get(acc, "")
        if not seq:
            raise SystemExit(
                f"ERROR: {acc} passed the gate but has no sequence in the passing "
                "FASTA. Refusing to admit an anchor with no sequence."
            )
        row = {
            "accession": acc,
            "tier": TIER_LITERATURE,
            "taxid": r["taxid"],
            "species": r["organism"],
            "family": r["family"],
            "class": "A",
            "evidence": f"literature-gated:{r['evidence_tier']}",
        }
        rows.append(row)
        prov.append({
            "composite_id": composite_id(row),
            "accession": acc,
            "tier": TIER_LITERATURE,
            "tier_label": TIER_LABELS[TIER_LITERATURE],
            "class": "A",
            "family": r["family"],
            "taxid": r["taxid"],
            "species": r["organism"],
            "source_db": "UniProtKB",
            "source_id": acc,
            "evidence_basis": (
                f"{r['evidence_tier']}; class_a={r.get('class_a_evidence', '')}; "
                f"pmids={r.get('primary_pmids', '')}; "
                f"reviewed={r.get('reviewed', '')}; name={r.get('protein_name', '')}"
            ),
            "sequence_length": r.get("sequence_length", str(len(seq))),
            "added_utc": now,
        })
        seqs.append((composite_id(row), seq))
    return rows, prov, seqs, held


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--anchors", default="references/anchors/anchor_set_PROD.tsv")
    ap.add_argument("--harvest", required=True)
    ap.add_argument("--literature", required=True)
    ap.add_argument("--literature-fasta", required=True)
    ap.add_argument("--outdir", default="references/anchors")
    ap.add_argument("--existing-sequences",
                    default="references/anchors/anchor_set_PROD_uniprot.tsv",
                    help="TSV carrying the CURRENT anchors' sequences, used to "
                         "detect a new entry that is the same protein as one "
                         "already present under a different identifier")
    ap.add_argument("--hold-family", action="append", default=None,
                    help="family label to HOLD out of the merge (repeatable); "
                         "defaults to 'unclassified'")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the merge and every check, write nothing")
    args = ap.parse_args(argv)

    anchors_path = Path(args.anchors)
    outdir = Path(args.outdir)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")

    existing = read_tsv(anchors_path)
    harvest = read_tsv(Path(args.harvest))
    verdicts = read_tsv(Path(args.literature))
    lit_fasta = read_fasta(Path(args.literature_fasta))

    # --- pre-merge state -------------------------------------------------
    existing_ids = [composite_id(r) for r in existing]
    if len(set(existing_ids)) != len(existing_ids):
        dupes = [k for k, v in Counter(existing_ids).items() if v > 1]
        raise SystemExit(
            f"ERROR: the EXISTING anchor table already contains duplicate "
            f"composite ids ({dupes[:5]}). Refusing to merge onto a table whose "
            "keys are not unique."
        )
    existing_acc = {r["accession"] for r in existing}

    hold = frozenset(args.hold_family if args.hold_family is not None
                     else ["unclassified"])
    new_odb, prov_odb, seq_odb = build_orthodb_rows(harvest, now)
    new_lit, prov_lit, seq_lit, held = build_literature_rows(
        verdicts, lit_fasta, now, hold)

    # --- collision checks, before anything is written --------------------
    problems: list[str] = []
    for label, rows in (("orthodb", new_odb), ("literature", new_lit)):
        clash = sorted({r["accession"] for r in rows} & existing_acc)
        if clash:
            problems.append(
                f"{label} tier reuses {len(clash)} accession(s) already in the "
                f"anchor set, e.g. {clash[:5]}"
            )
    cross = ({r["accession"] for r in new_odb} & {r["accession"] for r in new_lit})
    if cross:
        problems.append(f"the two new tiers share {len(cross)} accession(s): "
                        f"{sorted(cross)[:5]}")

    all_new_ids = [composite_id(r) for r in new_odb + new_lit]
    dupe_new = [k for k, v in Counter(all_new_ids).items() if v > 1]
    if dupe_new:
        problems.append(f"new rows carry duplicate composite ids: {dupe_new[:5]}")
    overlap_ids = sorted(set(all_new_ids) & set(existing_ids))
    if overlap_ids:
        problems.append(f"new composite ids collide with existing ones: "
                        f"{overlap_ids[:5]}")

    # Sequence-level duplication: a harvested OrthoDB gene can be the SAME
    # protein as an anchor already present under its UniProt accession. Counts
    # and key checks both pass in that case; only the sequence reveals it.
    new_seqs = dict(seq_odb + seq_lit)
    seq_to_ids: dict[str, list[str]] = defaultdict(list)
    for cid, s in new_seqs.items():
        if s:
            seq_to_ids[s].append(cid)
    internal_dupes = {s: ids for s, ids in seq_to_ids.items() if len(ids) > 1}

    # The same protein can arrive under an OrthoDB gene id having long been
    # present under its UniProt accession. Accession checks and row counts both
    # pass in that case; only comparing the SEQUENCES reveals it, and a
    # duplicated reference double-weights one point of the family envelope.
    existing_seq_ids: dict[str, str] = {}
    ex_path = Path(args.existing_sequences)
    if ex_path.exists():
        acc_to_row = {r["accession"]: r for r in existing}
        for r in read_tsv(ex_path):
            acc = r.get("queried_accession") or r.get("Entry", "")
            seq = (r.get("Sequence") or "").strip()
            if seq and acc in acc_to_row:
                existing_seq_ids[seq] = composite_id(acc_to_row[acc])
        print(f"cross-checking {len(new_seqs)} new sequences against "
              f"{len(existing_seq_ids)} existing anchor sequences")
    else:
        print(f"WARNING: {ex_path} absent -- cannot cross-check new sequences "
              "against the existing anchors", file=sys.stderr)
    redundant = {cid: existing_seq_ids[s] for cid, s in new_seqs.items()
                 if s and s in existing_seq_ids}

    if problems:
        for p in problems:
            print(f"ERROR: {p}", file=sys.stderr)
        raise SystemExit("refusing to merge")

    merged = existing + new_odb + new_lit

    # --- report ----------------------------------------------------------
    def spread(rows: list[dict]) -> dict:
        return {
            "rows": len(rows),
            "species": len({r["taxid"] for r in rows if r["taxid"]}),
        }

    print("=== MERGE ===")
    print(f"existing rows              {len(existing)}")
    print(f"+ tier {TIER_ORTHODB} orthodb-harvest    {len(new_odb)}")
    print(f"+ tier {TIER_LITERATURE} literature-gated  {len(new_lit)}")
    print(f"= merged rows              {len(merged)}")

    print("\n=== per-family, before -> after ===")
    fb, fa = Counter(r["family"] for r in existing), Counter(r["family"] for r in merged)
    print(f"{'family':<24}{'before':>8}{'after':>8}{'delta':>8}")
    for fam in sorted(set(fb) | set(fa)):
        print(f"{fam:<24}{fb.get(fam,0):>8}{fa.get(fam,0):>8}{fa.get(fam,0)-fb.get(fam,0):>+8}")

    print("\n=== per-class, before -> after ===")
    cb, ca = Counter(r["class"] for r in existing), Counter(r["class"] for r in merged)
    for k in sorted(set(cb) | set(ca)):
        print(f"  class {k:<4}{cb.get(k,0):>6} -> {ca.get(k,0):>6}")

    print("\n=== spread ===")
    print(f"  before: {spread(existing)}")
    print(f"  after : {spread(merged)}")

    if internal_dupes:
        print(f"\nWARNING: {len(internal_dupes)} sequence(s) appear under more than "
              f"one NEW composite id:")
        for s, ids in list(internal_dupes.items())[:5]:
            print(f"    len={len(s)}  {ids}")

    print(f"\n=== sequence-identical to an anchor already present: {len(redundant)} ===")
    for cid, old in list(redundant.items())[:10]:
        print(f"    {cid}  ==  {old}")

    if held:
        print(f"\n=== HELD, not merged: {len(held)} entr(ies) in families {sorted(hold)} ===")
        print("A family label is a novelty prototype. These need a ratified family "
              "before they can enter the reference set.")
        for r in held:
            print(f"  {r['accession']:<12}{r['evidence_tier']:<28}"
                  f"{r['organism'][:22]:<24}{r['protein_name'][:44]}")

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return 0

    # --- snapshot, then write --------------------------------------------
    prov_path = outdir / "anchor_set_PROD_tier_provenance.tsv"
    add_fa_path = outdir / "anchor_set_PROD_additions.fasta"
    snap = snapshot([anchors_path, prov_path, add_fa_path], outdir)
    print(f"\nsnapshot -> {snap}")

    write_atomic(anchors_path, merged, PROD_COLUMNS)
    write_atomic(prov_path, prov_odb + prov_lit, PROVENANCE_COLUMNS)
    write_fasta_atomic(add_fa_path, seq_odb + seq_lit)
    if held:
        write_atomic(outdir / "anchor_set_PROD_held_pending_family.tsv",
                     held, list(held[0].keys()))

    # --- re-read and verify ----------------------------------------------
    reread = read_tsv(anchors_path)
    reprov = read_tsv(prov_path)
    refa = read_fasta(add_fa_path)
    failures: list[str] = []

    if len(reread) != len(merged):
        failures.append(f"row count {len(reread)} != expected {len(merged)}")

    ids = [composite_id(r) for r in reread]
    if len(set(ids)) != len(ids):
        failures.append("composite ids are not unique after the write")

    # Every pre-existing id must still be present and byte-identical.
    lost = [i for i in existing_ids if i not in set(ids)]
    if lost:
        failures.append(f"{len(lost)} pre-existing composite id(s) vanished: {lost[:5]}")
    if ids[:len(existing_ids)] != existing_ids:
        failures.append("pre-existing rows were reordered or rewritten")

    prov_ids = {r["composite_id"] for r in reprov}
    dangling = sorted(prov_ids - set(ids))
    if dangling:
        failures.append(f"{len(dangling)} provenance row(s) name no anchor: {dangling[:5]}")
    unprovenanced = sorted(set(ids[len(existing_ids):]) - prov_ids)
    if unprovenanced:
        failures.append(f"{len(unprovenanced)} new anchor(s) have no provenance row: "
                        f"{unprovenanced[:5]}")

    # Sequence integrity: length recorded at source vs length actually written.
    for r in reprov:
        cid = r["composite_id"]
        if cid not in refa:
            failures.append(f"{cid} has no sequence in the additions FASTA")
            continue
        want, got = int(r["sequence_length"]), len(refa[cid])
        if want != got:
            failures.append(
                f"{cid} sequence length {got} != {want} recorded at source "
                "(an accumulate-by-accession bug has silently doubled sequences "
                "here before)"
            )

    print(f"\n=== POST-WRITE VERIFICATION ===")
    print(f"  anchor rows            {len(reread)}")
    print(f"  composite ids unique   {len(set(ids)) == len(ids)}")
    print(f"  pre-existing preserved {ids[:len(existing_ids)] == existing_ids}")
    print(f"  provenance rows        {len(reprov)}")
    print(f"  additions FASTA        {len(refa)}")
    print(f"  length mismatches      {sum(1 for f in failures if 'sequence length' in f)}")

    if failures:
        print("\nVERIFICATION FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  * {f}", file=sys.stderr)
        print(f"\nThe snapshot at {snap} holds the pre-merge state.", file=sys.stderr)
        return 1

    summary = {
        "merged_utc": now,
        "existing": len(existing),
        "added_orthodb": len(new_odb),
        "added_literature": len(new_lit),
        "total": len(merged),
        "snapshot": str(snap),
        "per_family_after": dict(fa),
        "per_class_after": dict(ca),
    }
    (outdir / "anchor_set_PROD_merge_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print("\nALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
