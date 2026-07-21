#!/usr/bin/env python3
"""Build the chemoreceptor PROBE set -- scored like a candidate, never a reference.

The three Aplysia californica "Chemosensory receptor A/B/C" entries (C5H877,
C5H675, C5H674; Cummins et al., PMID 19493360) survived the refexp2 evidence
gate, but they are TARGET class. Admitting them to a NON-chemoreceptor
reference set would place a prototype in chemoreceptor-like space and teach the
scorer that the very region the pipeline hunts is already "known" -- novelty is
a MIN over prototypes, so every genuinely chemoreceptor-like candidate would be
pulled DOWN by the one reference that most resembles it. The reference set must
stay a set of things a candidate is NOT.

So they are scored the same way a candidate is, and read as a behavioural probe
on the exclusion machinery.

WHAT THE PROBE CAN AND CANNOT SHOW
----------------------------------
Exclusion is RARE. With an exclusion rate p among the 790 class-A candidates,
the chance that a probe set of n behaving like average candidates produces even
one exclusion is 1-(1-p)^n -- a few percent at the rates this pipeline calls.
"Not excluded" is therefore the overwhelmingly likely outcome under BOTH the
hypothesis that the machinery works and the hypothesis that it is blind, and
carries almost no evidential weight. This is not a pass/fail test and must not
be reported as one.

Two readouts ARE informative:

  (a) EXCLUSION IS A RED FLAG. If the machinery labels a published molluscan
      chemosensory receptor a non-chemoreceptor, that is a genuine falsification
      of the exclusion rule that produced it -- one-directional, in the same
      sense as the project's Apgr check. It cannot confirm the machinery; it can
      break it.

  (b) NOVELTY RANK. Where the probes fall on the novelty axis RELATIVE TO the
      790-candidate distribution is a graded, always-available reading, unlike
      the near-degenerate binary. Report the probes' percentile in the candidate
      novelty distribution. Note the direction is not predicted a priori: these
      are molluscan class-A receptors whose family is not in the reference set,
      so high novelty is expected, and a LOW novelty percentile would say the
      reference set already explains chemoreceptor-like space -- which is the
      contamination this file exists to prevent.

A MEASURED INTERACTION WITH THE RATIFIED PEPTIDE ADMISSIONS
-----------------------------------------------------------
Readout (b) is not measured against a neutral reference set. All FIVE peptide
entries admitted by refexp4_ratify_held_families.py share a family-level
InterPro entry with two of the three probes:

    IPR019427 (7TM GPCR Srw)  -- C5H877, C5H675 and the admitted Myomodulin
                                 receptor 1, MIP receptor 1, Sex peptide receptor
    IPR052954 (GPCR_Ligand_Interaction) -- C5H877, C5H675 and the admitted
                                 RGWamide receptor 1, NPY-4 receptor 1

InterPro describes Srw as "a solo family amongst the superfamilies of
CHEMORECEPTORS", and it contains both the nematode chemoreceptor-class Srw
genes and the deorphanized peptide receptors SPR/MIP-R. Those admissions are
still correct -- their curated family entries and GO:0008528 peptide-receptor
activity place them in `peptide`, not in chemoreceptor space -- but they put
prototype mass in the signature neighbourhood of probes A and B, which can only
push those two probes' novelty DOWN. Measured on the pre-merge anchors carrying
Pfam data (1039 of 1550), PF10324/Srw anchors went from 4 (3 peptide, 1 orphan)
to 7. Probe C (C5H674) carries no signature beyond the generic rhodopsin-like
7TM entries and is unaffected.

So a low novelty percentile on probes A/B is confounded with this admission and
must not be read as "the reference set explains chemoreceptor space". Probe C is
the clean reading. If the difference matters, measure the probes against the
reference set with and without the five peptide admissions.

STRUCTURAL SEPARATION (four independent layers)
-----------------------------------------------
  1. Probes live in their OWN files. Nothing in the reference path reads them:
     `build_embedding_channel.load_ref_labels` reads anchor_set_PROD.tsv, which
     never contains a probe -- `verify_separation` asserts it.
  2. Probe ids carry the `PROBE_` prefix. `load_ref_labels` MINTS its keys as
     f"ANCHOR_{cls}_{tier}_{acc}", so it is not merely unlikely but impossible
     for it to emit a probe id.
  3. The probe TSV deliberately has NO `family` and no `tier` column. Pointing
     a reference loader at this file raises KeyError instead of silently
     building a prototype out of chemoreceptors -- the failure is loud.
  4. `assert_no_probe_in_reference` is exported for the scoring path to call on
     its actual reference/prototype inputs, and is asserted in the unit tests.
  5. `assert_no_probe_sequence_in_reference` catches the leak the first four
     layers cannot see: the SAME protein entering the reference set under a
     DIFFERENT accession. This is not hypothetical. UniProt currently holds
     A0ABM0GH42 and A0ABM0GH43, which are byte-identical to probes C5H675 and
     C5H674 but carry neither their accessions nor their literature -- a
     taxonomic-breadth harvest would admit them as ordinary molluscan class-A
     GPCRs and every id-based and accession-based check would pass. They are
     listed in `KNOWN_PROBE_TWINS` so they can also be rejected by accession.

Usage:
    python3 scripts/refexp5_build_chemoreceptor_probe.py --dry-run
    python3 scripts/refexp5_build_chemoreceptor_probe.py
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable

PROBE_ID_PREFIX = "PROBE_"
REFERENCE_ID_PREFIX = "ANCHOR_"

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = (
    "accession,protein_name,organism_name,organism_id,length,sequence,reviewed,"
    "protein_existence,xref_pfam,xref_interpro,lit_pubmed_id"
)

# The probe roster: accession -> (organism, protein_name, evidence_note).
#
# Deliberately short. A probe whose chemoreceptor claim is itself weak cannot
# falsify anything -- if the machinery excludes it, the honest reading is that
# the probe was never a chemoreceptor, not that the machinery is wrong. Only
# entries with a primary publication proposing them as chemosensory receptors
# belong here.
PROBE_ROSTER: dict[str, tuple[str, str, str]] = {
    "C5H877": ("Aplysia californica", "Chemosensory receptor A",
               "Cummins et al. PMID:19493360; rhinophore/tentacle-expressed "
               "class-A GPCR proposed as a chemosensory receptor; not "
               "deorphanized"),
    "C5H675": ("Aplysia californica", "Chemosensory receptor B",
               "Cummins et al. PMID:19493360; rhinophore/tentacle-expressed "
               "class-A GPCR proposed as a chemosensory receptor; not "
               "deorphanized"),
    "C5H674": ("Aplysia californica", "Chemosensory receptor C",
               "Cummins et al. PMID:19493360; rhinophore/tentacle-expressed "
               "class-A GPCR proposed as a chemosensory receptor; not "
               "deorphanized"),
}

# NOTE: the probe TSV has no `family` and no `tier` column ON PURPOSE. See
# layer 3 of the structural separation above. Do not add them.
PROBE_COLUMNS = [
    "probe_id", "accession", "class", "organism", "taxid", "protein_name",
    "sequence_length", "pfam", "interpro", "pmids", "evidence_note", "added_utc",
]


def probe_id(accession: str, gpcr_class: str = "A") -> str:
    """The probe's FASTA/npz key: ``PROBE_<class>_<accession>``.

    Distinct by construction from the anchor composite
    ``ANCHOR_<class>_<tier>_<accession>``, so a probe and an anchor can never
    collide even for the same accession.
    """
    return f"{PROBE_ID_PREFIX}{gpcr_class}_{accession}"


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


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- the separation guards ------------------------------------------------

def assert_no_probe_in_reference(reference_ids: Iterable[str],
                                 probe_ids: Iterable[str],
                                 context: str = "reference set") -> None:
    """Raise if any probe reaches a reference/prototype input.

    Call this from the scoring path on its ACTUAL inputs -- the reference npz
    keys, the reference label map, the per-family prototype members. Checking a
    fixture proves only that the fixture is clean.

    Two independent checks, because either alone can be defeated: exact id
    overlap catches a probe copied in verbatim, and the prefix check catches a
    probe re-keyed to look like an anchor.
    """
    reference_ids = list(reference_ids)
    probes = set(probe_ids)
    overlap = sorted(probes & set(reference_ids))
    if overlap:
        raise AssertionError(
            f"{len(overlap)} chemoreceptor probe(s) reached the {context}: "
            f"{overlap[:5]}. A probe in a novelty prototype teaches the scorer "
            "that chemoreceptor-like space is already known and suppresses the "
            "signal this pipeline exists to find.")
    smuggled = sorted(r for r in reference_ids if r.startswith(PROBE_ID_PREFIX))
    if smuggled:
        raise AssertionError(
            f"{len(smuggled)} id(s) in the {context} carry the probe prefix "
            f"{PROBE_ID_PREFIX!r}: {smuggled[:5]}")


def assert_no_probe_sequence_in_reference(
        reference_sequences: dict[str, str],
        probe_sequences: dict[str, str],
        context: str = "reference set") -> None:
    """Raise if any reference sequence is byte-identical to a probe sequence.

    The leak the id and accession checks cannot see. A harvest that admits the
    same protein under a different accession passes every key-based check --
    row counts, id uniqueness, accession disjointness all hold -- while placing
    a chemoreceptor prototype at distance zero from the probe it duplicates.
    Only comparing the sequences reveals it.

    Verified against UniProt 2026-07-21: A0ABM0GH42 and A0ABM0GH43 are
    byte-identical to C5H675 and C5H674.
    """
    by_seq = {seq: pid for pid, seq in probe_sequences.items() if seq}
    hits = sorted((ref_id, by_seq[seq])
                  for ref_id, seq in reference_sequences.items()
                  if seq and seq in by_seq)
    if hits:
        raise AssertionError(
            f"{len(hits)} reference sequence(s) in the {context} are "
            f"byte-identical to a chemoreceptor probe: {hits[:5]}. The "
            "accession differs, so every key-based check passed -- this is the "
            "duplicate-under-another-accession leak.")


# Accessions known to be the SAME protein as a probe. Reject these in any
# reference harvest; see assert_no_probe_sequence_in_reference.
KNOWN_PROBE_TWINS: dict[str, str] = {
    "A0ABM0GH42": "C5H675",
    "A0ABM0GH43": "C5H674",
}


def verify_separation(probe_tsv: Path, anchors_tsv: Path,
                      anchor_sequences: dict[str, str] | None = None) -> list[str]:
    """Check the on-disk probe set against the on-disk reference set.

    Returns a list of failures (empty when separation holds).
    """
    failures: list[str] = []
    probes = read_tsv(probe_tsv)
    anchors = read_tsv(anchors_tsv)
    probe_acc = {r["accession"] for r in probes}
    anchor_acc = {r["accession"] for r in anchors}

    shared = sorted(probe_acc & anchor_acc)
    if shared:
        failures.append(f"{len(shared)} probe accession(s) are ALSO anchors: "
                        f"{shared[:5]}")
    twins = sorted(set(KNOWN_PROBE_TWINS) & anchor_acc)
    if twins:
        failures.append(
            f"{len(twins)} anchor(s) are known duplicates of a probe under a "
            f"different accession: "
            f"{[(t, KNOWN_PROBE_TWINS[t]) for t in twins[:5]]}")
    if anchor_sequences is not None:
        probe_seqs = read_fasta(probe_tsv.with_suffix(".fasta"))
        try:
            assert_no_probe_sequence_in_reference(
                anchor_sequences, probe_seqs, context="anchor set")
        except AssertionError as exc:
            failures.append(str(exc))
    try:
        assert_no_probe_in_reference(
            # the ids load_ref_labels would mint from the anchor table
            [f"ANCHOR_{r['class']}_{r['tier']}_{r['accession']}" for r in anchors],
            [r["probe_id"] for r in probes],
            context="anchor table",
        )
    except AssertionError as exc:
        failures.append(str(exc))

    if any("family" in r for r in probes):
        failures.append("the probe TSV carries a `family` column -- that is the "
                        "column a reference loader turns into a prototype, and "
                        "its absence is a deliberate structural guard")
    for r in probes:
        if not r["probe_id"].startswith(PROBE_ID_PREFIX):
            failures.append(f"probe id {r['probe_id']!r} lacks the "
                            f"{PROBE_ID_PREFIX!r} prefix")
        if r["probe_id"].startswith(REFERENCE_ID_PREFIX):
            failures.append(f"probe id {r['probe_id']!r} carries the reference "
                            f"prefix {REFERENCE_ID_PREFIX!r}")
    return failures


# --- the readout ----------------------------------------------------------

def exclusion_power(n_probes: int, exclusion_rate: float) -> float:
    """P(at least one probe excluded) if probes behave like average candidates.

    This is the probability the probe produces its own red flag under the NULL.
    It is small, which is exactly why "not excluded" is weak evidence: the
    outcome is near-certain whether or not the machinery works.
    """
    return 1.0 - (1.0 - exclusion_rate) ** n_probes


def fetch_uniprot(accessions: list[str]) -> dict[str, dict]:
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
            "sequence": entry["sequence"]["value"],
            "pfam": ",".join(sorted(x["id"] for x in xrefs
                                    if x["database"] == "Pfam")),
            "interpro": ",".join(sorted(x["id"] for x in xrefs
                                        if x["database"] == "InterPro")),
            "pmids": ",".join(sorted({x["id"] for x in xrefs
                                      if x["database"] == "PubMed"})),
        }
    return out


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


def write_fasta_atomic(path: Path, records: list[tuple[str, str]]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
    os.replace(tmp, path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--anchors", default="references/anchors/anchor_set_PROD.tsv")
    ap.add_argument("--outdir", default="references/anchors")
    ap.add_argument("--candidate-fasta",
                    default="reports/plm_report/data/chemogpcrs_berghia_classA.fa",
                    help="candidate FASTA, read ONLY to report the population "
                         "size the probe is scored against")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    tsv_path = outdir / "chemoreceptor_probe_set.tsv"
    fa_path = outdir / "chemoreceptor_probe_set.fasta"

    records = fetch_uniprot(sorted(PROBE_ROSTER))
    rows, seqs, errors = [], [], []
    for acc, (organism, protein_name, note) in sorted(PROBE_ROSTER.items()):
        rec = records.get(acc)
        if rec is None:
            errors.append(f"{acc}: UniProt returned no record")
            continue
        if rec["organism"] != organism:
            errors.append(f"{acc}: organism {rec['organism']!r} != {organism!r}")
            continue
        if rec["protein_name"] != protein_name:
            errors.append(f"{acc}: protein name {rec['protein_name']!r} != "
                          f"{protein_name!r}")
            continue
        if len(rec["sequence"]) != rec["length"]:
            errors.append(f"{acc}: sequence length {len(rec['sequence'])} != "
                          f"declared {rec['length']}")
            continue
        pid = probe_id(acc)
        rows.append({
            "probe_id": pid, "accession": acc, "class": "A",
            "organism": rec["organism"], "taxid": rec["taxid"],
            "protein_name": rec["protein_name"],
            "sequence_length": str(rec["length"]), "pfam": rec["pfam"],
            "interpro": rec["interpro"], "pmids": rec["pmids"],
            "evidence_note": note, "added_utc": now,
        })
        seqs.append((pid, rec["sequence"]))

    print("=== PROBE ROSTER (verified against live UniProt) ===")
    for r in rows:
        print(f"  {r['probe_id']:<16}{r['accession']:<10}"
              f"{r['protein_name'][:30]:<32}{r['organism'][:24]:<26}"
              f"len={r['sequence_length']}")
    for e in errors:
        print(f"  FAIL {e}", file=sys.stderr)
    if errors:
        raise SystemExit(f"refusing to build the probe: {len(errors)} failed")

    # Separation against the reference set, checked BEFORE anything is written.
    anchors = read_tsv(Path(args.anchors))
    anchor_acc = {r["accession"] for r in anchors}
    leaked = sorted({r["accession"] for r in rows} & anchor_acc)
    if leaked:
        raise SystemExit(
            f"ERROR: probe accession(s) {leaked} are ALREADY in the reference "
            "set. A chemoreceptor in the reference set is the contamination "
            "this probe exists to keep out -- evict it before building.")

    n_cand = 0
    cand_path = Path(args.candidate_fasta)
    if cand_path.exists():
        n_cand = len(read_fasta(cand_path))

    print(f"\n=== READOUT DESIGN ===")
    print(f"  probes                        {len(rows)}")
    print(f"  candidate population          {n_cand or 'unknown (FASTA absent)'}")
    print("  primary readout (a) exclusion = RED FLAG (one-directional; can "
          "falsify the\n      exclusion rule, can never confirm it)")
    print("  primary readout (b) novelty percentile of each probe within the "
          "candidate\n      distribution (graded, always available)")
    print("\n  power of readout (a), if probes behave like average candidates:")
    for p in (0.010, 0.0177, 0.030):
        n_exc = round(p * n_cand) if n_cand else "?"
        print(f"      exclusion rate {p:6.2%} (~{n_exc} of {n_cand or '?'} "
              f"candidates) -> P(>=1 probe flagged) = "
              f"{exclusion_power(len(rows), p):6.2%}")
    print("  => 'not excluded' is the expected outcome under BOTH hypotheses "
          "and is WEAK\n     evidence. Do not report the probe as a pass.")

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return 0

    write_atomic(tsv_path, rows, PROBE_COLUMNS)
    write_fasta_atomic(fa_path, seqs)

    # --- re-read and verify ----------------------------------------------
    reread = read_tsv(tsv_path)
    refa = read_fasta(fa_path)
    failures: list[str] = []
    if len(reread) != len(rows):
        failures.append(f"probe rows {len(reread)} != {len(rows)}")
    ids = [r["probe_id"] for r in reread]
    if len(set(ids)) != len(ids):
        failures.append("probe ids are not unique")
    for r in reread:
        pid = r["probe_id"]
        if pid not in refa:
            failures.append(f"{pid} has no sequence in the probe FASTA")
        elif len(refa[pid]) != int(r["sequence_length"]):
            failures.append(f"{pid} sequence length {len(refa[pid])} != "
                            f"{r['sequence_length']} at source")
    failures.extend(verify_separation(tsv_path, Path(args.anchors)))

    print("\n=== POST-WRITE VERIFICATION ===")
    print(f"  probe rows             {len(reread)}")
    print(f"  probe ids unique       {len(set(ids)) == len(ids)}")
    print(f"  probe FASTA records    {len(refa)}")
    print(f"  separation failures    {len(failures)}")
    if failures:
        print("\nVERIFICATION FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  * {f}", file=sys.stderr)
        return 1
    print("\nALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
