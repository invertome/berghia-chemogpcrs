#!/usr/bin/env python3
"""Ratify the families of the entries held out of the reference merge.

``orthodb_merge_reference_tiers.py`` refuses to merge any entry labelled
``unclassified`` and parks it in ``anchor_set_PROD_held_pending_family.tsv``.
That refusal is not fastidiousness: the family column defines a novelty
PROTOTYPE, S_novel is a MIN over prototypes, so an ``unclassified`` grab bag
would become one heterogeneous cloud whose nearest point is close to
everything -- it could only pull novelty DOWN, everywhere. The merge therefore
declines to guess and defers the call.

This module records the ratified call and applies it, in three buckets.

ADMITTED -- assigned a real family, merged into the reference set.
    The family is NOT taken from the protein name. Every admission names the
    curated InterPro FAMILY-type signature that has to be on the live UniProt
    record for the assignment to stand, and the run aborts if it is absent.
    A name is an author's label; an InterPro family entry is a curated
    classification, and only the second is allowed to place a prototype.

EXCLUDED -- evidence-gate survivors carrying no family-level signature. They
    pass the evidence gate but nothing curated places them in a family, so there
    is no family to assign. They stay out rather than being parked in a catch-all.

    A CORRECTED STANDARD. These were originally screened on their protein names
    ("Transmembrane receptor", "Multitransmembrane protein", "4656HH") while the
    admissions were required to show a curated InterPro family entry -- two
    different standards, and the looser one was wrong. Q9NJC9 was excluded as
    "RHO G-protein coupled receptor", a name that commits to nothing, while its
    record carried better opsin evidence than any of the nine opsins admitted.
    It is now admitted. All five exclusions were re-audited 2026-07-21 under the
    curated standard; the remaining four hold on that evidence, not on their
    names. Exclusion here means "no curated family signature", never "the name
    was uninformative".

PROBE -- the target class. Handled by refexp5_build_chemoreceptor_probe.py and
    named here only so this module is the single place that accounts for all
    held entries. Admitting a known chemosensory receptor to a NON-chemoreceptor
    reference set would teach the scorer that chemoreceptor-like space is
    already "known" and suppress the exact signal the pipeline hunts.

IDENTIFIERS ARE WRITE-ONCE. Admitted rows are appended and take the tier-10
composite ``ANCHOR_A_10_<accession>`` -- byte for byte the id the merge would
have minted had the family been assigned at merge time. Pre-existing rows are
copied through unchanged and the run aborts if any of them moves.

Usage:
    python3 scripts/refexp4_ratify_held_families.py --dry-run
    python3 scripts/refexp4_ratify_held_families.py
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import shutil
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from anchorfix_class_reconcile import PROD_COLUMNS, composite_id  # noqa: E402

TIER_LITERATURE = "10"
TIER_LABEL = "literature-gated"

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = (
    "accession,protein_name,organism_name,organism_id,length,reviewed,"
    "protein_existence,cc_similarity,xref_pfam,xref_interpro,go_id"
)

# --- the ratified call ---------------------------------------------------
#
# accession -> (family, organism, protein_name, required InterPro family ids)
#
# `required` is a set of InterPro accessions of which AT LEAST ONE must be on
# the live record. They are all InterPro type=family entries, so a hit is a
# curated family placement and not a bare 7TM signature. IPR000276
# (GPCR_Rhodpsn) and IPR017452 are deliberately NOT accepted anywhere here:
# every one of the 22 held entries carries them, so they discriminate nothing
# and would turn this check into a rubber stamp.
#
#   IPR050125 GPCR opsins            -- the opsin family entry
#   IPR002962 Peropsin               -- an opsin subfamily entry
#   IPR019427 7TM GPCR Srw           -- contains Drosophila SPR and its MIP ligand
#   IPR053219 GPCR Dmsr-1            -- RFamide/myosuppressin receptor family
#   IPR053071 GPCR1-related receptor -- SPR/gustatory-learning receptor family
#   IPR052954 GPCR_Ligand_Interaction-- peptide/pheromone receptor family
OPSIN_SIGNATURES = frozenset({"IPR050125", "IPR002962"})
PEPTIDE_SIGNATURES = frozenset({"IPR019427", "IPR053219", "IPR053071", "IPR052954"})

ADMIT: dict[str, tuple[str, str, str, frozenset[str]]] = {
    # --- opsin (9) -------------------------------------------------------
    "A0A0H5ANL2": ("opsin", "Idiosepius paradoxus", "Retinochrome1", OPSIN_SIGNATURES),
    "A0A1J0CN99": ("opsin", "Argopecten irradians", "Retinochrome", OPSIN_SIGNATURES),
    "A0A0H5ANK7": ("opsin", "Nautilus pompilius", "Retinochrome", OPSIN_SIGNATURES),
    "A0A0H5B8L2": ("opsin", "Idiosepius paradoxus", "Retinochrome2", OPSIN_SIGNATURES),
    "A0A2D3VF54": ("opsin", "Stylochus ellipticus", "XenopsinA", OPSIN_SIGNATURES),
    "A0A2D3VR06": ("opsin", "Stylochus ellipticus", "XenopsinB", OPSIN_SIGNATURES),
    "A0A2D3VJW5": ("opsin", "Prostheceraeus vittatus", "XenopsinB1", OPSIN_SIGNATURES),
    "A0A2D3VCG0": ("opsin", "Prostheceraeus vittatus", "XenopsinB2", OPSIN_SIGNATURES),
    "A0A2D3VP13": ("opsin", "Prostheceraeus vittatus", "XenopsinA", OPSIN_SIGNATURES),
    # --- peptide (5) -----------------------------------------------------
    "A0A0K0PUW0": ("peptide", "Platynereis dumerilii", "RGWamide receptor 1",
                   PEPTIDE_SIGNATURES),
    "A0A0K0PUF7": ("peptide", "Platynereis dumerilii", "NPY-4 receptor 1",
                   PEPTIDE_SIGNATURES),
    "A0A0K0PUE1": ("peptide", "Platynereis dumerilii", "Myomodulin receptor 1",
                   PEPTIDE_SIGNATURES),
    "A0A0K0PUL9": ("peptide", "Platynereis dumerilii", "MIP receptor 1",
                   PEPTIDE_SIGNATURES),
    "A0AAU6WVE3": ("peptide", "Haliotis discus hannai", "Sex peptide receptor",
                   PEPTIDE_SIGNATURES),
    # --- opsin, ratified later (1) ---------------------------------------
    # Admitted 2026-07-21, reversing an earlier exclusion. The exclusion was
    # made on the protein name ("RHO G-protein coupled receptor" commits to no
    # family) while the admissions were required to show curated evidence --
    # two different standards, and the looser one produced the error. Under the
    # SAME curated check the admissions face, this record is the strongest
    # opsin in the held set: IPR050125 (the opsin family entry) PLUS IPR027430
    # (visual-pigment retinal binding site), photoreceptor-activity and
    # phototransduction GO, and the only curated SIMILARITY statement of the 22
    # held entries. None of the nine opsins admitted above carries the retinal
    # binding site. Re-audited against live UniProt 2026-07-21.
    "Q9NJC9": ("opsin", "Schistosoma mansoni", "RHO G-protein coupled receptor",
               OPSIN_SIGNATURES),
}

# Admissions whose ruling CHANGED, so the disposition table carries the history
# rather than presenting a reversal as if it had always been the call.
ADMISSION_NOTES: dict[str, str] = {
    "Q9NJC9": ("; ratified 2026-07-21, REVERSING an earlier exclusion. The "
               "exclusion was made on the protein name; re-audited under the "
               "curated standard this record carries IPR050125 plus the "
               "IPR027430 retinal binding site, photoreceptor-activity and "
               "phototransduction GO, and the only SIMILARITY statement among "
               "the 22 held entries -- stronger opsin evidence than any of the "
               "nine opsins admitted first, none of which has the binding site"),
}

# accession -> (reason, protein_name). Kept out of the reference set.
EXCLUDE: dict[str, tuple[str, str]] = {
    "Q962I3": ("name commits to no family; no InterPro family-level signature "
               "beyond the generic rhodopsin-like 7TM entries", "Sj-Ts5"),
    "C9K4W2": ("name is an internal clone label; no Pfam hit at all and no "
               "InterPro family-level signature", "4656HH"),
    "F2VWU2": ("name commits to no family; no InterPro family-level signature "
               "beyond the generic rhodopsin-like 7TM entries",
               "Transmembrane receptor"),
    "C0M0N9": ("name commits to no family; pathway GO terms only, no InterPro "
               "family-level signature", "Multitransmembrane protein"),
}

# The four exclusions above were re-audited against live UniProt + the InterPro
# API on 2026-07-21 under the SAME curated standard the admissions face, because
# the fifth exclusion (Q9NJC9) turned out to be wrong and had been screened on
# its name. Result: all four hold. Not one carries an InterPro FAMILY-type entry
# beyond IPR000276/IPR017452, the generic pair every held entry has; none has a
# SIMILARITY statement; none has family-informative GO. C9K4W2 has no Pfam hit
# at all. The name-based screen happened to agree with the curated evidence on
# these four -- it was still the wrong standard to apply, and Q9NJC9 is what it
# cost. Measured entry types, not assumed: see tests/unit/test_refexp4_exclusion_audit.py.

# accession -> protein_name. Target class; see refexp5_build_chemoreceptor_probe.py.
PROBE: dict[str, str] = {
    "C5H877": "Chemosensory receptor A",
    "C5H675": "Chemosensory receptor B",
    "C5H674": "Chemosensory receptor C",
}

PROVENANCE_COLUMNS = [
    "composite_id", "accession", "tier", "tier_label", "class", "family",
    "taxid", "species", "source_db", "source_id", "evidence_basis",
    "sequence_length", "added_utc",
]

DISPOSITION_COLUMNS = [
    "accession", "protein_name", "organism", "taxid", "phylum", "disposition",
    "assigned_family", "family_basis", "reason",
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


def fetch_uniprot(accessions: list[str]) -> dict[str, dict]:
    """Live UniProt records, keyed by accession.

    Every field the ratification rests on is read from here, never from the
    held table: the held table is a derived artifact and a family assignment
    validated against a derived copy validates nothing.
    """
    query = " OR ".join(f"accession:{a}" for a in accessions)
    url = (f"{UNIPROT_SEARCH}?query={urllib.parse.quote(query)}"
           f"&format=json&size=500&fields={UNIPROT_FIELDS}")
    with urllib.request.urlopen(url, timeout=180) as resp:
        payload = json.load(resp)
    out: dict[str, dict] = {}
    for entry in payload.get("results", []):
        desc = entry.get("proteinDescription", {})
        name_block = desc.get("recommendedName") or (
            desc.get("submissionNames") or [{}])[0]
        xrefs = entry.get("uniProtKBCrossReferences", [])
        out[entry["primaryAccession"]] = {
            "accession": entry["primaryAccession"],
            "protein_name": name_block.get("fullName", {}).get("value", ""),
            "organism": entry.get("organism", {}).get("scientificName", ""),
            "taxid": str(entry.get("organism", {}).get("taxonId", "")),
            "length": entry["sequence"]["length"],
            "interpro": frozenset(x["id"] for x in xrefs
                                  if x["database"] == "InterPro"),
            "pfam": frozenset(x["id"] for x in xrefs if x["database"] == "Pfam"),
        }
    return out


def verify_admissions(records: dict[str, dict], held: dict[str, dict],
                      seqs: dict[str, str]) -> tuple[list[dict], list[str]]:
    """Check every admission against the live record. Returns (verified, errors).

    Four independent checks per entry, all of which must pass:
      identity   -- the returned record is the requested protein, by organism
                    AND by protein name (a resolving accession is not enough)
      family     -- a curated InterPro FAMILY entry corroborates the assignment
      length     -- live length == the length the held table recorded
      sequence   -- the passing FASTA carries a sequence of exactly that length
    """
    verified, errors = [], []
    for acc, (family, organism, protein_name, required) in sorted(ADMIT.items()):
        rec = records.get(acc)
        if rec is None:
            errors.append(f"{acc}: UniProt returned no record")
            continue
        if rec["organism"] != organism:
            errors.append(f"{acc}: organism is {rec['organism']!r}, expected "
                          f"{organism!r} -- the accession resolves to a "
                          "different record than the one ratified")
            continue
        if rec["protein_name"] != protein_name:
            errors.append(f"{acc}: protein name is {rec['protein_name']!r}, "
                          f"expected {protein_name!r}")
            continue
        hits = sorted(rec["interpro"] & required)
        if not hits:
            errors.append(
                f"{acc} ({protein_name}): assigned family {family!r} is NOT "
                f"corroborated -- none of {sorted(required)} is on the record, "
                f"which carries {sorted(rec['interpro'])}. Refusing to assign a "
                "family from the protein name alone")
            continue
        row = held.get(acc)
        if row is None:
            errors.append(f"{acc}: not present in the held table")
            continue
        if int(row["sequence_length"]) != rec["length"]:
            errors.append(f"{acc}: held length {row['sequence_length']} != live "
                          f"UniProt length {rec['length']}")
            continue
        seq = seqs.get(acc, "")
        if len(seq) != rec["length"]:
            errors.append(f"{acc}: FASTA sequence length {len(seq)} != live "
                          f"UniProt length {rec['length']}")
            continue
        verified.append({
            "accession": acc, "family": family, "organism": organism,
            "protein_name": protein_name, "taxid": rec["taxid"],
            "interpro_hits": hits, "length": rec["length"], "sequence": seq,
            "held": row,
        })
    return verified, errors


def snapshot(paths: list[Path], outdir: Path) -> Path:
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    snap = outdir / f".ratify_snapshot_{stamp}"
    snap.mkdir(parents=True, exist_ok=True)
    for p in paths:
        if p.exists():
            shutil.copy2(p, snap / p.name)
    return snap


def write_atomic(path: Path, rows: list[dict], columns: list[str]) -> None:
    if not rows:
        raise SystemExit(f"ERROR: refusing to write an empty table to {path}")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def append_fasta_atomic(path: Path, records: list[tuple[str, str]]) -> int:
    """Rewrite the additions FASTA with `records` appended. Returns the prior count.

    Read-modify-write through a temp file rather than an O_APPEND write: a
    killed append leaves a truncated final record that still parses.
    """
    existing = read_fasta(path) if path.exists() else {}
    for name, seq in records:
        if name in existing:
            raise SystemExit(f"ERROR: {name} is already in {path}")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for name, seq in existing.items():
            fh.write(f">{name}\n{seq}\n")
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
    os.replace(tmp, path)
    return len(existing)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--anchors", default="references/anchors/anchor_set_PROD.tsv")
    ap.add_argument("--held",
                    default="references/anchors/anchor_set_PROD_held_pending_family.tsv")
    ap.add_argument("--literature-fasta",
                    default="references/non_chemo_gpcr/refexp2_candidate_verdicts_passing.fasta")
    ap.add_argument("--outdir", default="references/anchors")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    anchors_path, outdir = Path(args.anchors), Path(args.outdir)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")

    existing = read_tsv(anchors_path)
    held_rows = read_tsv(Path(args.held))
    held = {r["accession"]: r for r in held_rows}
    fasta = read_fasta(Path(args.literature_fasta))
    seqs = {k.split("|")[0]: v for k, v in fasta.items()}

    # Every held entry must be accounted for exactly once. A held entry that
    # this module forgets would silently stay out of the reference set with no
    # record of why, which is the failure this file exists to prevent.
    ruled = set(ADMIT) | set(EXCLUDE) | set(PROBE)
    unaccounted = sorted(set(held) - ruled)
    phantom = sorted(ruled - set(held))
    if unaccounted:
        raise SystemExit(f"ERROR: held entries with no ruling: {unaccounted}")
    if phantom:
        raise SystemExit(f"ERROR: ruling names entries not in the held table: "
                         f"{phantom}")

    records = fetch_uniprot(sorted(ruled))
    verified, errors = verify_admissions(records, held, seqs)

    print("=== VERIFICATION (live UniProt) ===")
    for v in verified:
        print(f"  OK   {v['accession']:<12}{v['family']:<10}"
              f"{v['protein_name'][:30]:<32}{v['organism'][:24]:<26}"
              f"len={v['length']:<5}via {','.join(v['interpro_hits'])}")
    for e in errors:
        print(f"  FAIL {e}", file=sys.stderr)
    if errors:
        raise SystemExit(f"refusing to ratify: {len(errors)} check(s) failed")
    if len(verified) != len(ADMIT):
        raise SystemExit(f"verified {len(verified)} != {len(ADMIT)} admissions")

    existing_ids = [composite_id(r) for r in existing]
    existing_by_acc = {r["accession"]: r for r in existing}

    # A ratified entry already in the reference set is NOT a clash -- it is a
    # prior run of this same module. Ratifying a later entry (Q9NJC9) must not
    # require re-admitting, renumbering or rewriting the fourteen that landed
    # first: identifiers are write-once and pre-existing rows are copied through
    # untouched. So the admissions are partitioned rather than rejected.
    #
    # The already-present ones are still VERIFIED above against live UniProt and
    # are re-checked here for family agreement. Silently skipping them would let
    # the ratified family drift away from the recorded one with no error.
    already, to_add = [], []
    for v in verified:
        prior = existing_by_acc.get(v["accession"])
        if prior is None:
            to_add.append(v)
            continue
        if prior["family"] != v["family"]:
            raise SystemExit(
                f"ERROR: {v['accession']} is already in the reference set with "
                f"family {prior['family']!r}, but this module ratifies it as "
                f"{v['family']!r}. Refusing to rewrite a write-once row; "
                "reconcile the ruling first.")
        already.append(v)

    if not to_add:
        print(f"\n[idempotent] all {len(already)} ratified admission(s) are "
              "already in the reference set; nothing to add")

    new_rows, prov_rows, new_seqs = [], [], []
    for v in to_add:
        h = v["held"]
        row = {
            "accession": v["accession"], "tier": TIER_LITERATURE,
            "taxid": v["taxid"], "species": v["organism"],
            "family": v["family"], "class": "A",
            "evidence": f"literature-gated:{h['evidence_tier']}",
        }
        cid = composite_id(row)
        new_rows.append(row)
        new_seqs.append((cid, v["sequence"]))
        prov_rows.append({
            "composite_id": cid, "accession": v["accession"],
            "tier": TIER_LITERATURE, "tier_label": TIER_LABEL, "class": "A",
            "family": v["family"], "taxid": v["taxid"], "species": v["organism"],
            "source_db": "UniProtKB", "source_id": v["accession"],
            "evidence_basis": (
                f"{h['evidence_tier']}; family ratified from InterPro "
                f"{','.join(v['interpro_hits'])} (curated family entry, not the "
                f"protein name); class_a={h.get('class_a_evidence','')}; "
                f"pmids={h.get('primary_pmids','')}; name={v['protein_name']}"),
            "sequence_length": str(v["length"]), "added_utc": now,
        })

    cids = [composite_id(r) for r in new_rows]
    if len(set(cids)) != len(cids):
        raise SystemExit("ERROR: admissions carry duplicate composite ids")
    if set(cids) & set(existing_ids):
        raise SystemExit(f"ERROR: new composite ids collide with existing: "
                         f"{sorted(set(cids) & set(existing_ids))}")

    merged = existing + new_rows
    before = Counter(r["family"] for r in existing)
    after = Counter(r["family"] for r in merged)

    print(f"\n=== RATIFICATION ===")
    print(f"  held entries        {len(held)}")
    print(f"  ratified admissions {len(verified)}")
    print(f"    already present   {len(already)} (write-once, untouched)")
    print(f"    newly added       {len(new_rows)}")
    print(f"  excluded            {len(EXCLUDE)}")
    print(f"  probe (target class){len(PROBE):>5}")
    print(f"  anchors {len(existing)} -> {len(merged)}")
    for fam in sorted(set(before) | set(after)):
        d = after.get(fam, 0) - before.get(fam, 0)
        if d:
            print(f"    {fam:<22}{before.get(fam,0):>6} -> {after.get(fam,0):>6} "
                  f"({d:+d})")

    disposition = []
    for v in verified:
        disposition.append({
            "accession": v["accession"], "protein_name": v["protein_name"],
            "organism": v["organism"], "taxid": v["taxid"],
            "phylum": v["held"]["phylum"], "disposition": "admitted",
            "assigned_family": v["family"],
            "family_basis": "InterPro:" + ",".join(v["interpro_hits"]),
            "reason": ("curated InterPro family entry corroborates the assignment"
                       + ADMISSION_NOTES.get(v["accession"], "")),
        })
    for acc, (reason, name) in sorted(EXCLUDE.items()):
        h = held[acc]
        disposition.append({
            "accession": acc, "protein_name": name, "organism": h["organism"],
            "taxid": h["taxid"], "phylum": h["phylum"], "disposition": "excluded",
            "assigned_family": "", "family_basis": "", "reason": reason,
        })
    for acc, name in sorted(PROBE.items()):
        h = held[acc]
        disposition.append({
            "accession": acc, "protein_name": name, "organism": h["organism"],
            "taxid": h["taxid"], "phylum": h["phylum"], "disposition": "probe",
            "assigned_family": "", "family_basis": "",
            "reason": "target class -- a molluscan chemosensory receptor. Scored "
                      "alongside candidates by refexp5_build_chemoreceptor_probe.py; "
                      "never a reference, prototype or envelope member",
        })

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return 0

    prov_path = outdir / "anchor_set_PROD_tier_provenance.tsv"
    add_fa_path = outdir / "anchor_set_PROD_additions.fasta"
    disp_path = outdir / "anchor_set_PROD_held_disposition.tsv"
    snap = snapshot([anchors_path, prov_path, add_fa_path], outdir)
    print(f"\nsnapshot -> {snap}")

    write_atomic(anchors_path, merged, PROD_COLUMNS)
    prov_before = read_tsv(prov_path) if prov_path.exists() else []
    write_atomic(prov_path, prov_before + prov_rows, PROVENANCE_COLUMNS)
    fa_before = append_fasta_atomic(add_fa_path, new_seqs)
    write_atomic(disp_path, disposition, DISPOSITION_COLUMNS)

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
    if ids[:len(existing_ids)] != existing_ids:
        failures.append("pre-existing rows were reordered or rewritten")
    if len(refa) != fa_before + len(new_seqs):
        failures.append(f"additions FASTA holds {len(refa)}, expected "
                        f"{fa_before + len(new_seqs)}")
    if len(reprov) != len(prov_before) + len(prov_rows):
        failures.append(f"provenance holds {len(reprov)}, expected "
                        f"{len(prov_before) + len(prov_rows)}")
    # `cids` covers the NEWLY ADDED rows only, so it is parallel to `to_add`,
    # NOT to `verified` (which also holds the already-present admissions).
    # Zipping it against `verified` silently pairs the wrong entries.
    if len(cids) != len(to_add):
        failures.append(f"composite ids {len(cids)} != admissions to add "
                        f"{len(to_add)} -- the id list is not parallel")
    for v, cid in zip(to_add, cids):
        if cid not in refa:
            failures.append(f"{cid} has no sequence in the additions FASTA")
        elif len(refa[cid]) != v["length"]:
            failures.append(f"{cid} sequence length {len(refa[cid])} != "
                            f"{v['length']} at source")
    # The already-present admissions are re-checked against their OWN ids,
    # recomputed from the on-disk row rather than assumed by position.
    for v in already:
        cid = composite_id(existing_by_acc[v["accession"]])
        if cid not in refa:
            failures.append(f"{cid} (already present) has no sequence in the "
                            "additions FASTA")
        elif len(refa[cid]) != v["length"]:
            failures.append(f"{cid} (already present) sequence length "
                            f"{len(refa[cid])} != {v['length']} at source")

    print("\n=== POST-WRITE VERIFICATION ===")
    print(f"  anchor rows            {len(reread)}")
    print(f"  composite ids unique   {len(set(ids)) == len(ids)}")
    print(f"  pre-existing preserved {ids[:len(existing_ids)] == existing_ids}")
    print(f"  additions FASTA        {len(refa)}")
    print(f"  provenance rows        {len(reprov)}")
    print(f"  length mismatches      {sum(1 for f in failures if 'length' in f)}")

    if failures:
        print("\nVERIFICATION FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  * {f}", file=sys.stderr)
        print(f"\nThe snapshot at {snap} holds the pre-ratification state.",
              file=sys.stderr)
        return 1
    print("\nALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
