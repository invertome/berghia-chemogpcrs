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

THE SET WAS EXPANDED 3 -> 159 ON 2026-07-21
-------------------------------------------
Three probes could not falsify anything: at the exclusion rates this pipeline
calls (1.0-3.0%), three probes had a 3.0-8.7% chance of producing a single
exclusion, so surviving was the expected outcome either way. Cummins et al.
reported 90 genes (AcCRa 28, AcCRb 38, AcCRc 24 -- read from the full text,
PMC2700072) and deposited three. `refexp6_recover_apgr_family.py` recovers the
rest from the AplCal3.0 RefSeq proteome by seed-calibrated HMM search; 155
survive. With C6FGJ8 (a different species and publication) the set is 159 in
three strata of UNEQUAL evidential weight -- see `stratum`.

CORRELATED EVIDENCE, NOT INDEPENDENT EVIDENCE
---------------------------------------------
The recovered members were found BY HOMOLOGY to the seeds. They resemble each
other by construction and do NOT multiply information the way 155 independent
probes would. Measured, not assumed:

    median pairwise identity            31.6%   (the set is not a blob)
    median nearest-neighbour identity   75.9%   (but it IS clustered)
    single-linkage clusters       146 at >=90%, 86 at >=70%, 39 at >=50%

Quoting only the 31.6% would badly understate the redundancy: 136 of the 159
have a >=50%-identical partner and 22 have a >=90%-identical one. What follows
differs sharply between the two readouts.

  (a) FALSIFICATION -- genuinely helped. Every member is another chance to trip
      the failure, and correlated draws still each get that chance. Deflating
      n by the cluster count, P(>=1 exclusion) rises from 3.0-8.7% at n=3 to
      32-77% (p=1.0%) or 70-99% (p=3.0%) over n_eff = 39..146. That is the
      difference between a test that cannot speak and one that can. Do NOT
      quote the naive n=159 figure; it is an upper bound.

  (b) NOVELTY PERCENTILE -- barely helped. 155 homologs of three receptors say
      where ONE FAMILY sits on the novelty axis. They say nothing about where
      chemoreceptors in general sit, and the expansion does not make the
      percentile more representative -- it makes one family's position more
      precisely measured. Worse, the recovery is BIASED across the family: it
      returns 2.89x Cummins's count for the tight AcCRa subfamily, 1.63x for
      AcCRb, but only 0.50x for the divergent AcCRc (whose members share as
      little as 19% identity, per Cummins). The set is enriched for exactly the
      most redundant, least novel part of the family.

WHAT THE PROBE CAN AND CANNOT SHOW
----------------------------------
Exclusion is RARE. With an exclusion rate p among the 790 class-A candidates,
the chance that a probe set of n_eff behaving like average candidates produces
even one exclusion is 1-(1-p)^n_eff. At n=3 this was a few percent, which is
why "not excluded" carried almost no evidential weight. The expansion changes
the magnitude of readout (a) but NOT its direction: exclusion still falsifies,
non-exclusion still cannot confirm. This is not a pass/fail test.

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

DOES THE EXPANDED SET INHERIT THE ENTANGLEMENT? YES, COMPLETELY.
----------------------------------------------------------------
Measured 2026-07-21, and the answer is worse than "partly":

  * All 155 recovered members pass PF10324/Srw above Pfam's gathering threshold
    BY CONSTRUCTION -- Srw membership is one of the selection criteria. Their
    Srw bitscores run 25.0 to 82.8 (median 51.7), every one of them ABOVE probe
    C's 24.6. Not a single recovered member is as Srw-marginal as probe C.
  * So the expansion adds NO clean members. Probe C remains the sole reading
    unconfounded by the peptide admissions -- it is now 1 of 159 rather than
    1 of 3, which makes the clean stratum a SMALLER fraction of the set, not a
    larger one.
  * C6FGJ8 does not help either: it carries IPR019427 (Srw), the same
    family-level entry that entangles probes A and B.

What the expansion does NOT do is put the family assignment in doubt. Every one
of the 155 scores decisively closer to a seed than to any of the 460
characterized peptide anchors -- median margin 212.4 bits, minimum 53.2, and
zero members where a peptide reference wins. The 5 ratified admissions score at
most 64.8 bits against any member versus a median 244.1 for the seeds. The
entanglement is a SIGNATURE-SHARING problem that biases novelty, not a
misassignment problem.

Practical consequence: report probe C separately, always. Reporting a single
pooled novelty percentile over 159 members would bury the one clean reading
under 158 confounded ones.

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
# accession -> (organism, protein_name, stratum, is_fragment, evidence_note)
PROBE_ROSTER: dict[str, tuple[str, str, str, bool, str]] = {
    "C5H877": ("Aplysia californica", "Chemosensory receptor A",
               "published-deposited", False,
               "Cummins et al. PMID:19493360 (BMC Biol 7:28, "
               "doi:10.1186/1741-7007-7-28, verified via NCBI 2026-07-21); "
               "rhinophore/tentacle-expressed class-A GPCR proposed as a "
               "chemosensory receptor; not deorphanized. AcCRa subfamily "
               "representative. NOT in the AplCal3.0 RefSeq proteome"),
    "C5H675": ("Aplysia californica", "Chemosensory receptor B",
               "published-deposited", False,
               "Cummins et al. PMID:19493360 (BMC Biol 7:28, "
               "doi:10.1186/1741-7007-7-28, verified via NCBI 2026-07-21); "
               "rhinophore/tentacle-expressed class-A GPCR proposed as a "
               "chemosensory receptor; not deorphanized. AcCRb subfamily "
               "representative; == RefSeq NP_001191494.1"),
    "C5H674": ("Aplysia californica", "Chemosensory receptor C",
               "published-deposited", False,
               "Cummins et al. PMID:19493360 (BMC Biol 7:28, "
               "doi:10.1186/1741-7007-7-28, verified via NCBI 2026-07-21); "
               "rhinophore/tentacle-expressed class-A GPCR proposed as a "
               "chemosensory receptor; not deorphanized. AcCRc subfamily "
               "representative; == RefSeq NP_001191495.1. The ONLY member "
               "carrying no curated Srw family signature -- the clean reading"),
    # --- the one member not derived from the seeds by homology ------------
    "C6FGJ8": ("Aplysia dactylomela", "Putative chemosensory receptor",
               "published-independent", True,
               "Cummins et al. PMID:19525430 (J Exp Biol 212:2037-44, "
               "doi:10.1242/jeb.026427, verified via NCBI 2026-07-21) -- a "
               "DIFFERENT species and a DIFFERENT publication from the three "
               "seeds, so the only cross-species, cross-publication evidence "
               "in the set. FRAGMENT: 246 aa vs 354 for the seeds (UniProt "
               "carries the Fragment flag). Truncation affects alignment- and "
               "embedding-based scoring alike, so its novelty percentile may "
               "report partly on length rather than on chemoreceptor identity "
               "and must never be pooled with the full-length members without "
               "that being visible. Note the authorship overlaps the seeds' "
               "paper (Cummins, Nagle, Degnan): independent publication and "
               "species, NOT an independent research group. AND ITS SEQUENCE "
               "INDEPENDENCE IS LIMITED: measured 77.6% identical to seed A "
               "(C5H877), vs 42.0% to B and 25.3% to C, and the three seeds "
               "are only 26.4-38.7% identical to one another. So it is CLOSER "
               "to seed A than any seed is to any other seed -- effectively an "
               "A. dactylomela ortholog of ApGPCR-A. It joins seed A's cluster "
               "at both the 50% and 70% identity thresholds and adds no new "
               "cluster. It is independent in PUBLICATION and SPECIES, not in "
               "sequence, and must not be counted as an independent draw"),
}

# NOTE: the probe TSV has no `family` and no `tier` column ON PURPOSE. See
# layer 3 of the structural separation above. Do not add them.
#
# `stratum` and `is_fragment` ARE present, and they are load-bearing: the
# strata do not carry equal evidential weight (see STRATA below) and a
# truncated member's scores are not comparable to a full-length one's. Pooling
# them into a single number is the misreading these columns exist to block.
PROBE_COLUMNS = [
    "probe_id", "accession", "class", "organism", "taxid", "protein_name",
    "sequence_length", "stratum", "is_fragment", "source_db", "subfamily",
    "best_seed", "best_seed_bitscore", "pfam", "interpro", "pmids",
    "evidence_note", "added_utc",
]

# Which seed each recovered member is nearest, mapped to the subfamily that
# seed represents in Cummins et al. Verified from the paper's full text
# (PMC2700072): AcCRa 28 genes, AcCRb 38, AcCRc 24, "a total of 90".
SEED_SUBFAMILY = {
    "C5H877": "AcCRa", "C5H675": "AcCRb", "C5H674": "AcCRc",
    "PROBE_A_C5H877": "AcCRa", "PROBE_A_C5H675": "AcCRb",
    "PROBE_A_C5H674": "AcCRc",
}

# What Cummins et al. reported, read from the full text and not from memory.
CUMMINS_SUBFAMILY_COUNTS = {"AcCRa": 28, "AcCRb": 38, "AcCRc": 24}
CUMMINS_TOTAL = 90

# --- the measured correlation structure of the expanded set ---------------
#
# MAFFT L-INS-i over all 159 members, pairwise identity over columns where both
# sequences have a residue. These are pinned so the power statement cannot
# drift away from the data it was computed from.
#
# The two numbers say different things and BOTH are needed:
#   global pairwise identity  median 31.6%  -- the set is not a uniform blob
#   nearest-neighbour identity median 75.9% -- but it IS clustered into groups
#                                              of close relatives
# Quoting only the first would understate the redundancy badly.
IDENTITY_MEDIAN_PAIRWISE = 31.6
IDENTITY_MEDIAN_NEAREST_NEIGHBOUR = 75.9

# Single-linkage clusters at a given % identity. This is the effective sample
# size used for the power statement -- interpretable, and computed from the
# data rather than assumed. n_eff is bracketed by these, not asserted as one
# number, because the right identity threshold for "would share an exclusion
# outcome" is not known.
IDENTITY_CLUSTERS = {90: 146, 80: 105, 70: 86, 60: 64, 50: 39, 40: 21}
N_EFF_BRACKET = (IDENTITY_CLUSTERS[50], IDENTITY_CLUSTERS[90])  # (39, 146)


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
#
# The UniProt pair was found by comparing sequences, not ids. The RefSeq pair
# was found the same way, while recovering the family: the AplCal3.0 proteome
# carries two of the three seeds under RefSeq accessions, so a recovery keyed
# on accession would have re-admitted C5H675 and C5H674 as if they were newly
# found members and inflated the count. Sequence comparison caught both.
# ApGPCR-A has no RefSeq twin -- it is absent from the annotation entirely.
KNOWN_PROBE_TWINS: dict[str, str] = {
    "A0ABM0GH42": "C5H675",
    "A0ABM0GH43": "C5H674",
    "NP_001191494.1": "C5H675",
    "NP_001191495.1": "C5H674",
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
    """P(at least one probe excluded) if probes behave like INDEPENDENT candidates.

    This is the probability the probe produces its own red flag under the NULL.
    It is small at n=3, which is exactly why "not excluded" was weak evidence:
    the outcome is near-certain whether or not the machinery works.

    THIS IS AN UPPER BOUND for the expanded set and must not be quoted for it.
    The recovered members were found BY HOMOLOGY to the seeds, so they are not
    independent draws; use `effective_exclusion_power`.
    """
    return 1.0 - (1.0 - exclusion_rate) ** n_probes


def effective_sample_size(n: int, mean_correlation: float) -> float:
    """Effective number of INDEPENDENT probes given mean pairwise correlation.

    The standard design-effect deflation, n_eff = n / (1 + (n-1)*rho). It is
    the honest denominator for a set built by homology: 155 members that
    resemble each other by construction do not carry 155 members' worth of
    information about whether the machinery excludes chemoreceptors.

    rho is a CORRELATION, not a sequence identity -- they are not the same
    quantity and identity must not be substituted for it directly. What
    identity gives is a defensible bracket: members sharing high identity will
    almost always share an exclusion outcome, so their outcome correlation is
    at least as high as one would guess from a weak-correlation model. Report a
    RANGE across plausible rho rather than a single number.
    """
    if n <= 1:
        return float(n)
    if not 0.0 <= mean_correlation <= 1.0:
        raise ValueError(f"correlation must be in [0,1], got {mean_correlation}")
    return n / (1.0 + (n - 1) * mean_correlation)


def effective_exclusion_power(n_probes: int, exclusion_rate: float,
                              mean_correlation: float) -> float:
    """P(>=1 exclusion) deflated for the correlation among a homologous set.

    At rho=0 this is `exclusion_power`. As rho -> 1 the whole set behaves like
    a SINGLE probe and the power collapses to the one-probe value, which is the
    correct limit: a set of near-identical sequences either all trip the
    exclusion or none of them do.
    """
    n_eff = effective_sample_size(n_probes, mean_correlation)
    return 1.0 - (1.0 - exclusion_rate) ** n_eff


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
            # UniProt's OWN fragment flag, so the roster's declaration is
            # checked against the source rather than trusted.
            "is_fragment": "fragment" in str(desc.get("flag", "")).lower(),
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
    ap.add_argument("--recovered",
                    default="references/anchors/apgr_recovered_members.tsv",
                    help="manifest written by refexp6_recover_apgr_family.py")
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
    for acc, (organism, protein_name, stratum, is_fragment, note) in sorted(
            PROBE_ROSTER.items()):
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
        # The fragment flag is scoring-relevant, so it is verified, not asserted.
        if rec["is_fragment"] != is_fragment:
            errors.append(
                f"{acc}: roster declares is_fragment={is_fragment} but UniProt "
                f"says {rec['is_fragment']}. A truncated probe scored as though "
                "full-length reports partly on length, not on identity")
            continue
        pid = probe_id(acc)
        rows.append({
            "probe_id": pid, "accession": acc, "class": "A",
            "organism": rec["organism"], "taxid": rec["taxid"],
            "protein_name": rec["protein_name"],
            "sequence_length": str(rec["length"]),
            "stratum": stratum, "is_fragment": "1" if is_fragment else "0",
            "source_db": "UniProtKB",
            "subfamily": SEED_SUBFAMILY.get(acc, ""),
            "best_seed": "", "best_seed_bitscore": "",
            "pfam": rec["pfam"],
            "interpro": rec["interpro"], "pmids": rec["pmids"],
            "evidence_note": note, "added_utc": now,
        })
        seqs.append((pid, rec["sequence"]))

    print("=== PROBE ROSTER (verified against live UniProt) ===")
    for r in rows:
        print(f"  {r['probe_id']:<16}{r['accession']:<10}"
              f"{r['protein_name'][:30]:<32}{r['organism'][:24]:<26}"
              f"len={r['sequence_length']:<5}{r['stratum']}"
              f"{' FRAGMENT' if r['is_fragment'] == '1' else ''}")
    for e in errors:
        print(f"  FAIL {e}", file=sys.stderr)
    if errors:
        raise SystemExit(f"refusing to build the probe: {len(errors)} failed")

    # --- the HMM-recovered stratum ---------------------------------------
    rec_tsv = Path(args.recovered)
    rec_fa = rec_tsv.with_suffix(".fasta")
    if rec_tsv.exists() and rec_fa.exists():
        rec_rows = read_tsv(rec_tsv)
        rec_seqs = read_fasta(rec_fa)
        uniprot_seq_ids = {s: r["probe_id"] for r, s in
                           zip(rows, (x[1] for x in seqs))}
        rec_errors = []
        for r in rec_rows:
            acc = r["refseq_accession"]
            seq = rec_seqs.get(acc, "")
            if len(seq) != int(r["sequence_length"]):
                rec_errors.append(f"{acc}: FASTA length {len(seq)} != "
                                  f"{r['sequence_length']}")
                continue
            # A recovered member that duplicates a published probe would
            # inflate the count and fake independence. Checked on SEQUENCE.
            if seq in uniprot_seq_ids:
                rec_errors.append(
                    f"{acc}: byte-identical to {uniprot_seq_ids[seq]} -- a "
                    "recovered member duplicating a published probe")
                continue
            pid = probe_id(acc)
            rows.append({
                "probe_id": pid, "accession": acc, "class": "A",
                "organism": r["organism"], "taxid": r["taxid"],
                "protein_name": r["protein_name"],
                "sequence_length": r["sequence_length"],
                "stratum": "hmm-recovered", "is_fragment": "0",
                "source_db": "RefSeq",
                "subfamily": SEED_SUBFAMILY.get(r["best_seed"], ""),
                "best_seed": r["best_seed"],
                "best_seed_bitscore": r["best_seed_bitscore"],
                "pfam": "PF00001,PF10324", "interpro": "", "pmids": "",
                "evidence_note": (
                    "recovered from the AplCal3.0 RefSeq proteome "
                    f"({r['assembly']}) by seed-calibrated HMM search; see "
                    "refexp6_recover_apgr_family.py. Evidence is homology to "
                    "the three published seeds plus PF10324/PF00001 family "
                    "membership -- NOT an individual publication. Exclusion of "
                    "this member is a WEAKER red flag than exclusion of a "
                    "published-deposited one"),
                "added_utc": now,
            })
            seqs.append((pid, seq))
        for e in rec_errors:
            print(f"  FAIL {e}", file=sys.stderr)
        if rec_errors:
            raise SystemExit(f"refusing to build: {len(rec_errors)} recovered "
                             "member(s) failed")
        print(f"\n  + {len(rec_rows)} hmm-recovered members merged")
    else:
        print(f"\n  [no recovered-member manifest at {rec_tsv}] -- probe set "
              "will hold only the published strata")

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

    strata = {}
    for r in rows:
        strata[r["stratum"]] = strata.get(r["stratum"], 0) + 1

    print(f"\n=== READOUT DESIGN ===")
    print(f"  probes                        {len(rows)}")
    for s in sorted(strata):
        print(f"      {s:<24}{strata[s]:>5}")
    print(f"  candidate population          {n_cand or 'unknown (FASTA absent)'}")
    print("  primary readout (a) exclusion = RED FLAG (one-directional; can "
          "falsify the\n      exclusion rule, can never confirm it)")
    print("  primary readout (b) novelty percentile of each probe within the "
          "candidate\n      distribution (graded, always available)")

    lo, hi = N_EFF_BRACKET
    print(f"\n  power of readout (a). The members were recovered BY HOMOLOGY to "
          "the seeds,\n  so they are NOT independent draws. n_eff is bracketed "
          f"by single-linkage\n  identity clustering: {hi} clusters at >=90% "
          f"identity, {lo} at >=50%.")
    print(f"\n      {'rate':>7}  {'n=3 (was)':>11}  {'n_eff=' + str(lo):>11}  "
          f"{'n_eff=' + str(hi):>11}  {'n=' + str(len(rows)) + ' (naive)':>15}")
    for p in (0.010, 0.0177, 0.030):
        print(f"      {p:6.2%}  {exclusion_power(3, p):10.1%}  "
              f"{exclusion_power(lo, p):10.1%}  {exclusion_power(hi, p):10.1%}  "
              f"{exclusion_power(len(rows), p):14.1%}")
    print("\n  The naive column is an UPPER BOUND and must not be quoted. The "
          "honest\n  reading is the bracketed pair: the expansion moves "
          "readout (a) from 'expected\n  under both hypotheses' to genuinely "
          "informative, but by roughly the factor a\n  few dozen independent "
          "probes would buy, not a hundred and fifty.")
    print(f"\n  Within-set identity: median pairwise "
          f"{IDENTITY_MEDIAN_PAIRWISE}%, median nearest-neighbour "
          f"{IDENTITY_MEDIAN_NEAREST_NEIGHBOUR}%.")
    print("  For readout (b) the expansion buys much less: a tight homologous "
          "cluster\n  locates ONE family on the novelty axis, not "
          "chemoreceptors in general.")

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
