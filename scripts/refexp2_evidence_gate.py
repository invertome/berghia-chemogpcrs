#!/usr/bin/env python3
"""refexp2_evidence_gate.py — expand the class-A GPCR reference set beyond
Swiss-Prot reviewed entries, admitting ONLY entries with positive evidence of
genuine characterization.

Why this exists
---------------
The embedding-novelty envelope scores a Berghia candidate against the spread of
a family's OWN members. The anchor set is overwhelmingly vertebrate (38
molluscan / 57 lophotrochozoan of 1094), so a molluscan query is measured across
~550 My against a bar calibrated on vertebrate-to-vertebrate variation. Swiss-Prot
reviewed holds only 19 class-A GPCRs for all of Lophotrochozoa and all 19 are
already anchors: the reviewed ceiling is exhausted. Expansion therefore has to
reach into TrEMBL, and the entire value of that relaxation rests on the gate
below being strict. A loose gate reproduces exactly the failure the reviewed-only
rule prevented -- propagating misidentifications onto invertebrate proteins,
which is the worst case for chemoreceptor expansions.

The gate
--------
S0  universe      class-A membership VERIFIED per entry, not inferred from a
                  family name (see verify_class_a).
S1  existence     Protein Existence 1 (protein level) or 2 (transcript level).
                  PE3 (inferred from homology) and PE4 (predicted) are exactly
                  the automated inference this gate exists to exclude.
S2  evidence code Reviewed entries must carry ECO:0000269 (experimental evidence
                  from a publication). Entries resting only on ECO:0000256
                  (automatic assertion by rule), ECO:0000259 (signature match) or
                  ECO:0000313 (imported) do not qualify on annotation alone.

                  MEASURED CAVEAT, and the reason S3/S4 exist: ECO:0000269 is
                  applied by Swiss-Prot curators and is STRUCTURALLY ABSENT from
                  TrEMBL. Of 363 unreviewed lophotrochozoan PE1/PE2 class-A
                  entries, exactly 0 carry it. Requiring it of TrEMBL is
                  therefore identical to requiring reviewed status, which would
                  make the relaxation a no-op. For unreviewed entries the
                  curator code is replaced by -- not waived in favour of -- the
                  literature and characterization gates, which are stricter in
                  practice than an annotation code.
S3  literature    At least one PubMed-indexed JOURNAL ARTICLE. Bare EMBL
                  submissions (no peer review) and bulk genome / transcriptome /
                  proteome announcement papers are rejected: a sequence reported
                  in a genome paper is not a characterized receptor.
S4  characterized The entry's own name must assign a specific ligand or
                  function. Submitter placeholders ("Orphan G-protein coupled
                  receptor 34", "GCR002"), hedges ("Putative", "Probable") and
                  uncharacterized labels are rejected.

S4 is not a formality. Of the 47 new Platynereis entries carried by the
large-scale deorphanization study, all 47 are named "Orphan G-protein coupled
receptor N" by the authors themselves -- they are the receptors that were tested
and NOT deorphanized; the ones that were are already anchors. Of the 82 new
planarian entries, all are "GCR###" RNAi-screen placeholders. The submitters'
own naming does the discrimination the gate needs.

Modes
-----
discover     find and gate lophotrochozoan candidates; write candidates + funnel
audit        apply the same gate to the EXISTING anchors; measurement only
corrections  verify the known family misassignments against UniProt

This script never edits the anchor set, never mints and never renumbers an
identifier. Identifiers are write-once; every artifact it writes is additive and
keyed by UniProt accession.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from typing import Dict, Iterable, List, Optional, Sequence

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# The class call is delegated, not reimplemented. `gpcr_class_from_evidence` is
# the resolver that decided anchor-set membership, so sharing it is what keeps
# this audit and the anchor set from disagreeing about which accessions are
# class A. The two module-private names are imported alongside it for the same
# reason: re-typing the family regex here would create a second copy of the rule,
# free to drift from the one that is authoritative.
from curate_gpcr_references import (  # noqa: E402
    _CURATED_FAMILY_TO_CLASS as CURATED_FAMILY_TO_CLASS,
    _CURATED_GPCR_FAMILY_RE as CURATED_GPCR_FAMILY_RE,
    gpcr_class_from_evidence,
)

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "berghia-chemogpcrs/1.0"

# Lophotrochozoa. Molluscs are the target, but annelids / flatworms / brachiopods
# are vastly closer to Berghia than any vertebrate and widen the envelope in the
# right direction.
LOPHOTROCHOZOA_TAXID = "1206795"

# Class-A membership (see verify_class_a). UniProt's curated family statement is
# the authority; these signatures corroborate it, and stand in for it only where
# curation is silent.
PFAM_7TM_1 = "PF00001"          # 7tm_1, Pfam clan CL0192 (rhodopsin-like)
INTERPRO_RHODOPSIN = "IPR000276"  # GPCR, rhodopsin-like superfamily

# SEARCH-SIDE ONLY: the UniProt `family:` query term used to cast the candidate
# net wide (refexp3 builds its universe query from it). It is deliberately NOT
# the verdict rule -- verify_class_a resolves class from the record's own curated
# statement, and matching a query string is not evidence of anything. UniProt's
# `family:` field is hyphenation-insensitive, measured 2026-07-20: this spelling
# and "g protein-coupled receptor 1 family" both return 396250 entries globally
# and both return 5700 within Mollusca, so the net is not narrowed by the choice.
CURATED_CLASS_A_FAMILY = "g-protein coupled receptor 1 family"

# The family string ALONE is not the universe. Unreviewed records frequently
# carry the domain signature but no curated family statement, so a family-only
# net never examines them at all. Measured live 2026-07-21 for Lophotrochozoa
# PE1/PE2: 381 entries by family string, 572 by the union below -- a third of
# the eligible records were invisible to the narrower query.
#
# This widens the SEARCH ONLY. verify_class_a still re-derives class per entry
# from the record's own curated statement, so admitting a record to the
# candidate pool is not admitting it to the reference set. refexp3 builds the
# same union from these same two constants.
CLASS_A_UNIVERSE = (
    f'(family:"{CURATED_CLASS_A_FAMILY}" '
    f"OR xref:pfam-{PFAM_7TM_1} "
    f"OR xref:interpro-{INTERPRO_RHODOPSIN})"
)

# Pfam clan CL0176 (Chemosens_recp) -- NOT the rhodopsin clan. A hit on these
# does not establish class A, and silently treating it as if it did is the
# specific trap this constant exists to name.
PFAM_NON_CLASS_A = {
    "PF08395": "7tm_7, clan CL0176 Chemosens_recp",
    "PF02949": "7tm_6, clan CL0176 Chemosens_recp",
}

EXPERIMENTAL_ECO = {
    "ECO:0000269",  # experimental evidence from a publication
    "ECO:0007744",  # combinatorial evidence with experimental support
    "ECO:0007829",  # combinatorial evidence, X-ray/structure derived
}
AUTOMATED_ECO = {
    "ECO:0000256",  # automatic assertion according to rule
    "ECO:0000259",  # automatic assertion, signature match
    "ECO:0000313",  # imported from another database
    "ECO:0000255",  # sequence-model evidence
    "ECO:0000250",  # sequence similarity
}

LOPHOTROCHOZOAN_PHYLA = (
    "Mollusca", "Annelida", "Platyhelminthes", "Brachiopoda",
    "Phoronida", "Nemertea", "Rotifera", "Bryozoa",
)

# ---------------------------------------------------------------------------
# S3: literature
# ---------------------------------------------------------------------------

# Bulk sequence-survey papers. A receptor that appears only in one of these was
# sequenced, not characterized.
#
# 'cDNA library' / 'expressed genes' / 'preliminary' were added after reading all
# 26 supporting abstracts: an oyster hemocyte cDNA-library survey and an Octopus
# gastric-ganglion transcriptome survey (whose own title says "Preliminary
# Characterization ... Putative") were assigning receptor names purely by
# homology. Those are sequence surveys wearing a characterization title, and the
# entries they carry are exactly what this gate must exclude.
_SURVEY_PAPER = re.compile(
    r"\bgenome\b|\bgenomic\b|\btranscriptome\b|\btranscriptomic\b|"
    r"\bdraft\b|\bassembly\b|\bassemblies\b|\bsequencing\b|\bsequenced\b|"
    r"\bproteome\b|\bEST\b|\bexpressed sequence tag|\bcDNA librar|"
    r"\bexpressed genes\b|\bpreliminary\b|"
    r"\bgene (?:catalog|catalogue|repertoire|content|inventory)\b|"
    r"\bwhole[- ]genome\b|\bde novo assembl",
    re.I,
)

# Language marking an experimental receptor-function study: a ligand was applied
# and a response measured. This is the HIGHER bar (tier
# 'functionally-characterized'), distinct from merely being named in a paper.
_FUNCTIONAL_ASSAY = re.compile(
    r"\bdeorphani[sz]|\bagonist\b|\bantagonist\b|\bpharmacolog|"
    r"\bEC50\b|\bIC50\b|\bdose[- ]response\b|\bcognate (?:ligand|peptide)\b|"
    r"\breceptor[- ](?:ligand|peptide)[- ]pair\b|\bligand[- ]receptor pair\b|"
    r"\bactivated by\b|\bactivates? the receptor\b|\bheterologous(?:ly)? express|"
    r"\bcalcium (?:imaging|mobili[sz]ation)\b|\bcAMP (?:assay|accumulation)\b|"
    r"\bfunctional(?:ly)? charact",
    re.I,
)

# ---------------------------------------------------------------------------
# S4: characterization, judged from the entry's own name
# ---------------------------------------------------------------------------

# Submitter placeholders and hedges. These are positive statements that the
# protein was NOT characterized, made by the people who deposited it.
_UNCHARACTERIZED_NAME = re.compile(
    r"\borphan\b|\buncharacteri[sz]ed\b|\bputative\b|\bprobable\b|"
    r"\bpredicted\b|\bhypothetical\b|\bunknown\b|\bunnamed\b|"
    r"\blike protein\b|\bfragment\b|"
    # '...receptor-like' / '...-like receptor' is a homology assignment, not an
    # identification. "Serotonin receptor-like planarian receptor 1" says the
    # submitter matched it to serotonin receptors, not that serotonin was shown
    # to activate it.
    r"receptor[- ]like\b|\blike[- ]receptor\b|[- ]like\b",
    re.I,
)

# Names that carry no ligand or function at all: bare GPCR labels and locus-style
# placeholders (GCR002, A1, LOC123456, Smed_...).
_PLACEHOLDER_NAME = re.compile(
    r"^(?:G[- ]?protein[- ]coupled receptor|GPCR|G[- ]?protein receptor|"
    r"seven transmembrane[- ]?\w*|7TM \w*)\s*[-_ ]?\d*$"
    r"|^GCR\d+$|^GPR\d+$|^[A-Z]{1,3}\d{1,4}$|^LOC\d+$|"
    r"^(?:GPCR|GCR|Receptor)[-_ ]?\d+$",
    re.I,
)

# A bare gene-family abbreviation carrying only a serial number -- NPYR-10,
# GCR002, TKR2. The family label in such a name is a HOMOLOGY call that was then
# numbered, which is precisely the automated inference this gate exists to
# exclude; it is not an experimentally established ligand.
#
# The discriminator is a spelled-out word. A genuine characterization reads
# "Adipokinetic hormone receptor 1A" or "5-HT1 receptor" -- both carry a real
# term alongside any abbreviation -- whereas "NPYR-10" is an abbreviation and a
# serial number and nothing else.
#
# Added after 14 Schmidtea mediterranea "NPYR-N" entries cleared every other
# stage. They come from a 100-plus-receptor RNAi screen that does not
# individually deorphanize them, and their family label is a homology call.
# Admitting them would have been exactly the "lower the bar to produce a number"
# failure this gate is meant to prevent.
_ABBREVIATION_INDEX_NAME = re.compile(r"^[A-Za-z]{2,8}[-_ ]?\d{1,4}[A-Za-z]?$")

# The rule above matches on SHAPE, and that shape is also what a genuine
# spelled-out name takes once it carries an isoform number: "Opsin-3",
# "Peropsin 1", "Xenopsin1" and "Acropsin 1" are all letters-then-digits and
# were all being rejected as placeholders. That contradicts the rule's own
# stated discriminator -- a spelled-out word -- because "opsin" IS one.
#
# So the shape test is kept and a spelled-out functional word exempts a name
# from it. The vocabulary is spelled-out terms only: an acronym plus a serial
# number ("NPYR-10", "GCR002", "TKR2", "NPFR1") carries no such word and is
# still rejected, which is the case the rule was added for.
_SPELLED_OUT_FUNCTION_WORD = re.compile(
    r"opsin|rhodopsin|receptor|hormone|serotonin|dopamine|adrenergic|"
    r"histamine|octopamine|tyramine|muscarinic|melatonin|corazonin|"
    r"tachykinin|kinin|orexin|gastrin|somatostatin|opioid|purinoceptor|"
    r"prostaglandin|prostanoid|thromboxane|cannabinoid|leukotriene|"
    r"chemokine|allatotropin|allatostatin|adipokinetic|sulfakinin",
    re.I,
)


def protein_existence_level(entry: dict) -> int:
    """Numeric PE level from a UniProt JSON entry ('2: Evidence at transcript
    level' -> 2). Raises on an unparseable value rather than guessing."""
    raw = entry.get("proteinExistence", "")
    match = re.match(r"\s*([1-5])\s*:", raw)
    if not match:
        raise ValueError(f"unparseable proteinExistence: {raw!r}")
    return int(match.group(1))


def is_reviewed(entry: dict) -> bool:
    return "reviewed (Swiss-Prot)" in entry.get("entryType", "")


def evidence_codes(entry: dict) -> set:
    """Every ECO code attached anywhere in the entry's annotation."""
    codes: set = set()

    def walk(node) -> None:
        if isinstance(node, dict):
            code = node.get("evidenceCode")
            if code:
                codes.add(code)
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    for key in ("comments", "features", "proteinDescription", "genes",
                "keywords", "references"):
        walk(entry.get(key))
    return codes


def has_experimental_evidence(entry: dict) -> bool:
    """True when the entry carries a curator-assigned experimental ECO code.

    Structurally unavailable to TrEMBL entries -- see the module docstring.
    """
    return bool(evidence_codes(entry) & EXPERIMENTAL_ECO)


def rests_only_on_automated_evidence(entry: dict) -> bool:
    """True when every ECO code present is an automated assertion."""
    codes = evidence_codes(entry)
    if not codes:
        return True
    return not (codes - AUTOMATED_ECO)


def protein_name(entry: dict) -> str:
    """The entry's own name: recommended (curated) first, then the submitter's."""
    description = entry.get("proteinDescription", {})
    recommended = description.get("recommendedName") or {}
    if recommended:
        return recommended.get("fullName", {}).get("value", "")
    submitted = description.get("submissionNames") or []
    if submitted:
        return submitted[0].get("fullName", {}).get("value", "")
    alternative = description.get("alternativeNames") or []
    if alternative:
        return alternative[0].get("fullName", {}).get("value", "")
    return ""


def lineage(entry: dict) -> List[str]:
    return list(entry.get("organism", {}).get("lineage", []))


def phylum_of(entry: dict) -> str:
    names = lineage(entry)
    for phylum in LOPHOTROCHOZOAN_PHYLA:
        if phylum in names:
            return phylum
    if "Vertebrata" in names:
        return "Vertebrata"
    for phylum in ("Arthropoda", "Nematoda", "Echinodermata", "Tunicata",
                   "Cephalochordata", "Cnidaria", "Placozoa", "Porifera"):
        if phylum in names:
            return phylum
    return "other"


def cross_reference_ids(entry: dict, database: str) -> List[str]:
    return [x["id"] for x in entry.get("uniProtKBCrossReferences", [])
            if x.get("database") == database]


def curated_family_text(entry: dict) -> str:
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SIMILARITY":
            for text in comment.get("texts", []):
                return text.get("value", "")
    return ""


def verify_class_a(entry: dict) -> tuple:
    """Verify class-A membership from the record's own evidence.

    Returns ``(is_class_a, basis, evidence_string)``.

    The determination is a CLASS test, not a domain-detection test. Detecting
    PF00001/IPR000276 and reporting the result as though it were a class call
    errs in both directions, and both directions were measured against live
    UniProt:

    * False positives on the negative side -- eight anchors (ACKR1/DARC
      Q5Y7A3, Q8IWP5, Q53ZP8 and Q09554, Q09964, Q09965, Q09966, Q9UJ42) are
      curated "Belongs to the G protein-coupled receptor 1 family", which IS
      class A, yet are divergent enough that no rhodopsin profile fires. Domain
      detection fails all eight; the class they are curated into passes them.
    * False negatives on the positive side -- A0A0A8JZN4 and P46091 are curated
      into the CMKLR family, which names no GPCR class at all, but both carry
      PF00001 and both are genuinely class A. A curation-only test would drop
      them.

    So curation decides where it speaks, and the signature decides where it does
    not. Precedence, delegated to ``gpcr_class_from_evidence`` -- the SAME
    resolver that decided anchor-set membership -- so that this audit and the
    anchor set cannot drift apart:

    1. The curated "Belongs to the G[- ]protein[- ]coupled receptor N family"
       statement. N=1 is class A, N=2 class B, N=3 class C. It is decisive in
       BOTH directions: a curated family-2/3 entry is rejected even when a
       rhodopsin signature fires on it.
    2. Where curation names no class family (absent, or a family such as CMKLR
       or Fz/Smo that carries no number), the specific Pfam family PF00001.
    3. Failing that, InterPro IPR000276 -- the weakest basis, recorded as such.

    A chemosensory-clan hit (CL0176, NOT the rhodopsin clan CL0192) is rejected
    with its clan named, but only in the signature-fallback branch: it cannot
    override an explicit curated class.
    """
    pfam = set(cross_reference_ids(entry, "Pfam"))
    interpro = set(cross_reference_ids(entry, "InterPro"))
    curated = curated_family_text(entry)

    signatures = []
    if PFAM_7TM_1 in pfam:
        signatures.append(f"Pfam:{PFAM_7TM_1}")
    if INTERPRO_RHODOPSIN in interpro:
        signatures.append(f"InterPro:{INTERPRO_RHODOPSIN}")
    corroboration = ("corroborated by " + "+".join(signatures) if signatures
                     else "no rhodopsin-clan signature")

    # 1. Curated class statement, decisive in both directions.
    resolved = gpcr_class_from_evidence(curated, ";".join(sorted(pfam)))
    match = CURATED_GPCR_FAMILY_RE.search(curated or "")
    curated_class = CURATED_FAMILY_TO_CLASS.get(match.group(1)) if match else None
    if curated_class == "A":
        return True, "curated-family-1", f"curated:GPCR-1-family; {corroboration}"
    if curated_class:
        return (False, f"curated-family-{match.group(1)}",
                f"curated:GPCR-{match.group(1)}-family = class {curated_class}, "
                f"not A; {corroboration}")

    # 2/3. Curation names no class family: fall back to the signature, and say so.
    context = f"curation non-decisive ({curated or 'no curated family'})"
    wrong_clan = pfam & set(PFAM_NON_CLASS_A)
    if wrong_clan and PFAM_7TM_1 not in pfam:
        names = "; ".join(f"{p} ({PFAM_NON_CLASS_A[p]})" for p in sorted(wrong_clan))
        return False, "chemosensory-clan", f"chemosensory clan, not rhodopsin: {names}; {context}"
    if PFAM_7TM_1 in pfam:
        assert resolved == "A", f"resolver disagreement: PF00001 present but resolved {resolved}"
        return True, "signature-pfam", f"Pfam:{PFAM_7TM_1}; {context}"
    if INTERPRO_RHODOPSIN in interpro:
        return True, "signature-interpro", f"InterPro:{INTERPRO_RHODOPSIN} only; {context}"
    return False, "no-evidence", f"no curated class family and no rhodopsin-clan signature ({curated or 'none'})"


def citations(entry: dict) -> List[dict]:
    return [reference.get("citation", {}) for reference in entry.get("references", [])]


def pubmed_id(citation: dict) -> str:
    for xref in citation.get("citationCrossReferences", []):
        if xref.get("database") == "PubMed":
            return xref.get("id", "")
    return ""


def is_journal_article(citation: dict) -> bool:
    """A peer-reviewed journal article with a PubMed ID.

    Bare EMBL 'submission' records are rejected: a database deposit has had no
    peer review at all.
    """
    return citation.get("citationType") == "journal article" and bool(pubmed_id(citation))


def is_survey_paper(text: str) -> bool:
    """A sequence survey rather than a characterization of this receptor."""
    return bool(_SURVEY_PAPER.search(text or ""))


def shows_functional_assay(text: str) -> bool:
    """True when the paper reports a measured receptor response, not just a name."""
    return bool(_FUNCTIONAL_ASSAY.search(text or ""))


def primary_literature(entry: dict,
                       abstracts: Optional[Dict[str, str]] = None) -> List[dict]:
    """Peer-reviewed journal articles that are not sequence surveys.

    Survey detection reads the TITLE only, deliberately. Scanning abstracts as
    well was tried and over-rejected badly: a functional study routinely mentions
    "genome" or "sequenced" in passing, which is not the same as being a genome
    announcement. Titles name what a paper IS. Both real survey papers in this
    dataset are caught by their titles alone -- an oyster hemocyte study titled
    "Identification of expressed genes in cDNA library ..." and an Octopus
    gastric-ganglion study titled "Preliminary Characterization of Gene- and
    Putative Neurochemical-Complexity ...".

    ``abstracts`` is accepted for signature symmetry with the tiering step and is
    intentionally unused here.
    """
    return [c for c in citations(entry)
            if is_journal_article(c) and not is_survey_paper(c.get("title", ""))]


def name_asserts_function(name: str) -> bool:
    """True when the entry's own name commits to a ligand or function.

    Rejects submitter placeholders and hedges. 'Orphan G-protein coupled
    receptor 34' and 'GCR002' are their authors stating the protein was not
    characterized.
    """
    if not name or not name.strip():
        return False
    stripped = name.strip()
    if _UNCHARACTERIZED_NAME.search(stripped):
        return False
    if _PLACEHOLDER_NAME.match(stripped):
        return False
    if (_ABBREVIATION_INDEX_NAME.match(stripped)
            and not _SPELLED_OUT_FUNCTION_WORD.search(stripped)):
        return False
    # A name that is only a bare GPCR label carries no ligand.
    bare = re.sub(r"[\s\-_,\.]+", " ", stripped).strip().lower()
    if bare in {"g protein coupled receptor", "g-protein coupled receptor",
                "gpcr", "receptor", "g protein coupled receptor 1",
                "seven transmembrane receptor", "rhodopsin like receptor"}:
        return False
    return True


# ---------------------------------------------------------------------------
# Family assignment, for reporting the per-family gain only
# ---------------------------------------------------------------------------

# Deliberately LOCAL to this script rather than reusing
# curate_gpcr_references.classify_family. That production classifier misses
# every one of these real records, measured on this run's survivors:
#   * "Xenopsin", "Peropsin"   -- its needle is \bopsin\b, and there is no word
#                                 boundary inside a compound like Xenopsin
#   * "Rhabdomeric opsin1"     -- same needle, defeated by the digit suffix
#   * "5-HT1 receptor"         -- its needle is \b5-HT\b, defeated by "5-HT1"
#   * "Serotonergic GPCR"      -- it matches "serotonin receptor", not the
#                                 adjectival form
# Widening those patterns would change what the production non-chemoreceptor
# reference set contains, which is a pipeline behaviour change and the user's
# call, so it is REPORTED rather than made here.
_FAMILY_FOR_REPORT = [
    ("opsin", re.compile(r"opsin\b|opsin\d|\brhodopsin", re.I)),
    ("aminergic", re.compile(
        r"5-?HT\d?|5-?hydroxytryptamine|serotonin|serotonergic|serotoninergic|"
        r"\bdopamine\b|adrenergic|histamine receptor|octopamine|tyramine|"
        r"trace amine|muscarinic", re.I)),
    ("peptide", re.compile(
        r"\bFLRFamide\b|\bFMRFamide\b|allatotropin|allatostatin|melatonin|"
        r"adipokinetic|\bAKH\b|corazonin|tachykinin|neuropeptide|vasopressin|"
        r"oxytocin|annetocin|cholecystokinin|gastrin|orexin|sulfakinin|"
        r"thyrotropin-releasing|myoinhibitory|kinin|opioid|somatostatin", re.I)),
    ("glycoprotein-hormone", re.compile(
        r"glycoprotein hormone|thyrotropin receptor|follicle-stimulating|"
        r"lutropin|luteinizing hormone receptor", re.I)),
    ("nucleotide", re.compile(r"adenosine receptor|purinoceptor|\bP2[RY]", re.I)),
    ("lipid", re.compile(
        r"prostaglandin|prostanoid|thromboxane|cannabinoid|"
        r"sphingosine 1-phosphate|lysophosphatidic|leukotriene", re.I)),
    ("chemokine", re.compile(r"chemokine receptor|\bCCR\d|\bCXCR\d|\bACKR\d", re.I)),
]


def family_for_report(name: str) -> str:
    """Coarse family label for the per-family summary. Reporting only -- it
    never feeds the gate and never decides whether an entry is admitted."""
    for family, pattern in _FAMILY_FOR_REPORT:
        if pattern.search(name or ""):
            return family
    return "unclassified"


# ---------------------------------------------------------------------------
# Gate
# ---------------------------------------------------------------------------

GATE_STAGES = ("class_a", "existence", "evidence_code", "literature",
               "characterized")


def gate_entry(entry: dict, pubmed_abstracts: Optional[Dict[str, str]] = None) -> dict:
    """Apply the full gate to one entry. Returns a verdict dict; never raises
    on a failing entry, so the funnel can be reported honestly."""
    pubmed_abstracts = pubmed_abstracts or {}
    accession = entry.get("primaryAccession", "")
    name = protein_name(entry)
    reviewed = is_reviewed(entry)

    verdict = {
        "accession": accession,
        "protein_name": name,
        "organism": entry.get("organism", {}).get("scientificName", ""),
        "taxid": entry.get("organism", {}).get("taxonId", ""),
        "phylum": phylum_of(entry),
        "reviewed": "reviewed" if reviewed else "unreviewed",
        "sequence_length": entry.get("sequence", {}).get("length", ""),
        "pass_class_a": False, "pass_existence": False,
        "pass_evidence_code": False, "pass_literature": False,
        "pass_characterized": False, "deorphanized": False,
        "evidence_tier": "", "family": family_for_report(name),
        "class_a_evidence": "", "class_a_basis": "", "protein_existence": "",
        "evidence_codes": "", "primary_pmids": "", "failed_at": "",
    }

    is_class_a, class_a_basis, class_a_evidence = verify_class_a(entry)
    verdict["class_a_evidence"] = class_a_evidence
    verdict["class_a_basis"] = class_a_basis
    verdict["pass_class_a"] = is_class_a
    if not is_class_a:
        verdict["failed_at"] = "class_a"
        return verdict

    try:
        level = protein_existence_level(entry)
    except ValueError:
        verdict["failed_at"] = "existence"
        return verdict
    verdict["protein_existence"] = level
    verdict["pass_existence"] = level in (1, 2)
    if not verdict["pass_existence"]:
        verdict["failed_at"] = "existence"
        return verdict

    codes = evidence_codes(entry)
    verdict["evidence_codes"] = ";".join(sorted(codes))
    if reviewed:
        # A curator has looked at this entry, so the curator code is available
        # and is required.
        verdict["pass_evidence_code"] = has_experimental_evidence(entry)
    else:
        # ECO:0000269 cannot exist here. Record the fact and defer to S3/S4.
        verdict["pass_evidence_code"] = rests_only_on_automated_evidence(entry)
    if not verdict["pass_evidence_code"]:
        verdict["failed_at"] = "evidence_code"
        return verdict

    papers = primary_literature(entry, pubmed_abstracts)
    verdict["primary_pmids"] = ";".join(pubmed_id(c) for c in papers)
    verdict["pass_literature"] = bool(papers)
    if not papers:
        verdict["failed_at"] = "literature"
        return verdict

    verdict["pass_characterized"] = name_asserts_function(name)
    if not verdict["pass_characterized"]:
        verdict["failed_at"] = "characterized"
        return verdict

    # Two honest tiers among the survivors. Clearing the gate means the receptor
    # is named, published and not homology-inferred; it does NOT by itself mean a
    # ligand was applied and a response measured. Only the higher tier claims
    # that, and the two are reported separately rather than merged into one
    # flattering number.
    for citation in papers:
        blob = (f"{citation.get('title', '')} "
                f"{pubmed_abstracts.get(pubmed_id(citation), '')}")
        if shows_functional_assay(blob):
            verdict["deorphanized"] = True
            break
    verdict["evidence_tier"] = ("functionally-characterized"
                               if verdict["deorphanized"]
                               else "published-not-deorphanized")
    return verdict


def funnel(verdicts: Sequence[dict]) -> dict:
    """Survivor count at each gate stage, plus the deorphanized subset."""
    counts = {"input": len(verdicts)}
    surviving = list(verdicts)
    for stage in GATE_STAGES:
        surviving = [v for v in surviving if v[f"pass_{stage}"]]
        counts[stage] = len(surviving)
    counts["deorphanized"] = sum(1 for v in surviving if v["deorphanized"])
    return counts


# ---------------------------------------------------------------------------
# UniProt / PubMed access
# ---------------------------------------------------------------------------

def _get(url: str, retries: int = 4) -> tuple:
    last = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=180) as response:
                total = int(response.headers.get("x-total-results", "-1"))
                return (response.read().decode("utf-8"),
                        response.headers.get("Link", ""), total)
        except Exception as exc:  # network flake; retry with backoff
            last = exc
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    raise last


def fetch_json_paginated(query: str, size: int = 200) -> List[dict]:
    """Run a UniProt query in JSON, following Link pagination to exhaustion.

    Refuses an empty result and asserts the row count matches the server's own
    x-total-results, so neither a renamed field nor a truncated page can be
    mistaken for the complete answer.
    """
    params = {"query": query, "format": "json", "size": str(size)}
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"
    results: List[dict] = []
    total = None
    while url:
        body, link, page_total = _get(url)
        if total is None:
            total = page_total
        results.extend(json.loads(body).get("results", []))
        match = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = match.group(1) if match else None

    if not results:
        raise RuntimeError(
            "UniProt returned zero rows -- refusing to treat an empty result as "
            f"success. A field or keyword was probably renamed. Query: {query}")
    if total is not None and total >= 0 and len(results) != total:
        raise RuntimeError(
            f"pagination lost rows: parsed {len(results)} but server reported "
            f"{total} for query: {query}")
    return results


def fetch_entries_by_accession(accessions: Iterable[str], batch: int = 100) -> Dict[str, dict]:
    """Resolve accessions to full JSON entries. Unresolvable accessions raise --
    a silent drop would let a stale accession vanish without a trace."""
    wanted = list(dict.fromkeys(a for a in accessions if a))
    resolved: Dict[str, dict] = {}
    for start in range(0, len(wanted), batch):
        chunk = wanted[start:start + batch]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        params = {"query": query, "format": "json", "size": "500"}
        body, _, _ = _get(f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}")
        for entry in json.loads(body).get("results", []):
            resolved[entry["primaryAccession"]] = entry
            for secondary in entry.get("secondaryAccessions", []):
                resolved.setdefault(secondary, entry)
    missing = [a for a in wanted if a not in resolved]
    if missing:
        raise RuntimeError(
            f"{len(missing)} accession(s) did not resolve against UniProt: {missing[:10]}")
    return resolved


def fetch_pubmed_abstracts(pmids: Iterable[str], batch: int = 100) -> Dict[str, str]:
    """Title + abstract per PMID, from NCBI E-utilities."""
    wanted = [p for p in dict.fromkeys(pmids) if p]
    out: Dict[str, str] = {}
    for start in range(0, len(wanted), batch):
        chunk = wanted[start:start + batch]
        params = {"db": "pubmed", "id": ",".join(chunk), "retmode": "xml",
                  "rettype": "abstract"}
        body, _, _ = _get(f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(params)}")
        for block in re.split(r"</PubmedArticle>", body):
            pmid_match = re.search(r"<PMID[^>]*>(\d+)</PMID>", block)
            if not pmid_match:
                continue
            text = " ".join(re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>",
                                       block, re.S))
            title = " ".join(re.findall(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>",
                                        block, re.S))
            out[pmid_match.group(1)] = re.sub(r"<[^>]+>", " ", f"{title} {text}")
        time.sleep(0.4)  # NCBI rate limit
    return out


def verify_sequence_integrity(entries: Dict[str, dict], records: Sequence[dict]) -> None:
    """Assert every emitted sequence matches its source record exactly.

    Accumulating by accession has previously DOUBLED sequences in this project
    while passing every count check, so length equality is checked per record
    against the freshly-fetched source, and duplicate accessions are refused.
    """
    seen = set()
    for record in records:
        accession = record["accession"]
        if accession in seen:
            raise RuntimeError(f"duplicate accession in output: {accession}")
        seen.add(accession)
        source = entries.get(accession)
        if source is None:
            raise RuntimeError(f"record {accession} has no source entry")
        expected = source.get("sequence", {}).get("value", "")
        if record["sequence"] != expected:
            raise RuntimeError(
                f"{accession}: sequence differs from source "
                f"(len {len(record['sequence'])} vs {len(expected)})")
        declared = source.get("sequence", {}).get("length")
        if declared is not None and len(record["sequence"]) != declared:
            raise RuntimeError(
                f"{accession}: length {len(record['sequence'])} != declared {declared}")


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_tsv(rows: Sequence[dict], columns: Sequence[str], path: str) -> None:
    """Atomic write. Refuses an empty table so a failed query can never replace
    a good artifact with headers only."""
    if not rows:
        raise RuntimeError("refusing to write an empty table")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp"
    with open(tmp, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp, path)


VERDICT_COLUMNS = [
    "accession", "protein_name", "family", "organism", "taxid", "phylum", "reviewed",
    "sequence_length", "protein_existence", "class_a_basis", "class_a_evidence",
    "pass_class_a", "pass_existence", "pass_evidence_code", "pass_literature",
    "pass_characterized", "deorphanized", "evidence_tier", "failed_at",
    "primary_pmids", "evidence_codes",
]


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------

def discovery_query() -> str:
    return (f"(taxonomy_id:{LOPHOTROCHOZOA_TAXID}) AND "
            f"{CLASS_A_UNIVERSE} AND (existence:1 OR existence:2)")


def run_discover(args) -> int:
    print("[refexp2] fetching lophotrochozoan class-A PE1/PE2 entries...",
          file=sys.stderr)
    entries = fetch_json_paginated(discovery_query())
    print(f"[refexp2] universe: {len(entries)} entries", file=sys.stderr)

    existing = set()
    if args.anchor_tsv:
        with open(args.anchor_tsv, newline="") as handle:
            existing = {r["accession"] for r in csv.DictReader(handle, delimiter="\t")}
    candidates = [e for e in entries if e["primaryAccession"] not in existing]
    print(f"[refexp2] {len(entries) - len(candidates)} already anchors; "
          f"{len(candidates)} new candidates", file=sys.stderr)

    pmids = {pubmed_id(c) for e in candidates for c in citations(e) if pubmed_id(c)}
    print(f"[refexp2] fetching {len(pmids)} PubMed abstracts...", file=sys.stderr)
    abstracts = fetch_pubmed_abstracts(pmids)

    verdicts = [gate_entry(e, abstracts) for e in candidates]
    write_tsv(verdicts, VERDICT_COLUMNS, args.output)

    counts = funnel(verdicts)
    print("\n[refexp2] GATE FUNNEL (new candidates)", file=sys.stderr)
    for stage in ("input",) + GATE_STAGES + ("deorphanized",):
        print(f"  {stage:<16} {counts[stage]:>5}", file=sys.stderr)

    survivors = [v for v in verdicts if not v["failed_at"]]
    if survivors:
        by_accession = {e["primaryAccession"]: e for e in candidates}
        records = [{"accession": v["accession"],
                    "sequence": by_accession[v["accession"]]["sequence"]["value"]}
                   for v in survivors]
        verify_sequence_integrity(by_accession, records)
        fasta = os.path.splitext(args.output)[0] + "_passing.fasta"
        with open(fasta + ".tmp", "w") as handle:
            for verdict, record in zip(survivors, records):
                handle.write(f">{verdict['accession']}|{verdict['phylum']}|"
                             f"{verdict['organism'].replace(' ', '_')}\n"
                             f"{record['sequence']}\n")
        os.replace(fasta + ".tmp", fasta)
        print(f"[refexp2] wrote {len(records)} verified sequences -> {fasta}",
              file=sys.stderr)

    write_family_summary(verdicts, args.anchor_tsv,
                         os.path.splitext(args.output)[0] + "_family_summary.tsv")
    return 0


SUMMARY_COLUMNS = [
    "family", "anchor_total", "anchor_mollusca", "anchor_lophotrochozoa",
    "gated_additions", "additions_mollusca", "additions_lophotrochozoa",
    "additions_functionally_characterized", "lophotrochozoa_after",
    "fold_change_lophotrochozoa",
]


def write_family_summary(verdicts: Sequence[dict], anchor_tsv: str, path: str) -> None:
    """Per-family gain: what the envelope actually gets, family by family.

    The point of the expansion is molluscan/lophotrochozoan representation, so
    that is what is counted -- not a raw total, which would flatter the result.
    """
    with open(anchor_tsv, newline="") as handle:
        anchors = [r for r in csv.DictReader(handle, delimiter="\t")
                   if r.get("class") == "A"]
    anchor_entries = fetch_entries_by_accession(a["accession"] for a in anchors)

    families = sorted({a["family"] for a in anchors} |
                      {v["family"] for v in verdicts if not v["failed_at"]})
    rows = []
    for family in families:
        in_family = [a for a in anchors if a["family"] == family]
        lineages = [lineage(anchor_entries[a["accession"]]) for a in in_family]
        anchor_moll = sum(1 for l in lineages if "Mollusca" in l)
        anchor_lopho = sum(1 for l in lineages
                           if any(p in l for p in LOPHOTROCHOZOAN_PHYLA))
        added = [v for v in verdicts if not v["failed_at"] and v["family"] == family]
        add_moll = sum(1 for v in added if v["phylum"] == "Mollusca")
        add_lopho = len(added)  # every survivor is lophotrochozoan by construction
        after = anchor_lopho + add_lopho
        rows.append({
            "family": family,
            "anchor_total": len(in_family),
            "anchor_mollusca": anchor_moll,
            "anchor_lophotrochozoa": anchor_lopho,
            "gated_additions": len(added),
            "additions_mollusca": add_moll,
            "additions_lophotrochozoa": add_lopho,
            "additions_functionally_characterized": sum(
                1 for v in added if v["deorphanized"]),
            "lophotrochozoa_after": after,
            "fold_change_lophotrochozoa": (
                f"{after / anchor_lopho:.2f}x" if anchor_lopho else
                ("new" if after else "0")),
        })
    write_tsv(rows, SUMMARY_COLUMNS, path)
    print(f"\n[refexp2] per-family summary -> {path}", file=sys.stderr)
    header = (f"  {'family':<22}{'anchLopho':>10}{'added':>7}{'addMoll':>9}"
              f"{'funcChar':>10}{'after':>7}{'fold':>8}")
    print(header, file=sys.stderr)
    for row in rows:
        print(f"  {row['family']:<22}{row['anchor_lophotrochozoa']:>10}"
              f"{row['gated_additions']:>7}{row['additions_mollusca']:>9}"
              f"{row['additions_functionally_characterized']:>10}"
              f"{row['lophotrochozoa_after']:>7}"
              f"{row['fold_change_lophotrochozoa']:>8}", file=sys.stderr)


def run_audit(args) -> int:
    """Apply the same gate to the EXISTING anchors. Measurement only: nothing is
    removed, no identifier is touched."""
    with open(args.anchor_tsv, newline="") as handle:
        anchors = list(csv.DictReader(handle, delimiter="\t"))
    class_a = [a for a in anchors if a.get("class") == "A"]
    print(f"[refexp2] auditing {len(class_a)} class-A anchors "
          f"(of {len(anchors)} total)...", file=sys.stderr)

    entries = fetch_entries_by_accession(a["accession"] for a in class_a)
    pmids = {pubmed_id(c) for a in class_a
             for c in citations(entries[a["accession"]]) if pubmed_id(c)}
    print(f"[refexp2] fetching {len(pmids)} PubMed abstracts...", file=sys.stderr)
    abstracts = fetch_pubmed_abstracts(pmids)

    verdicts = []
    for anchor in class_a:
        verdict = gate_entry(entries[anchor["accession"]], abstracts)
        verdict["anchor_family"] = anchor.get("family", "")
        verdict["anchor_tier"] = anchor.get("tier", "")
        verdict["anchor_evidence"] = anchor.get("evidence", "")
        verdicts.append(verdict)
    write_tsv(verdicts, VERDICT_COLUMNS + ["anchor_family", "anchor_tier",
                                           "anchor_evidence"], args.output)

    counts = funnel(verdicts)
    print("\n[refexp2] GATE FUNNEL (existing class-A anchors)", file=sys.stderr)
    for stage in ("input",) + GATE_STAGES + ("deorphanized",):
        print(f"  {stage:<16} {counts[stage]:>5}", file=sys.stderr)
    return 0


# Known family misassignments in the anchor set, each stated as a HYPOTHESIS to
# be checked against UniProt rather than a conclusion. run_corrections resolves
# every accession and reports whether the record supports the change; nothing is
# written on the strength of this table alone.
#
# The first two are a substring collision on "gonadotropin": GnRH (gonadotropin-
# RELEASING hormone) is a peptide receptor, whereas the glycoprotein hormones are
# FSH / LH / TSH / CG. P46023 is deliberately listed as unresolved -- UniProt
# gives GRL101 only the bare class-A family with no FSH/LSH/TSH subfamily, so
# promoting it would be inference, not curation, and that is the user's call.
PROPOSED_CORRECTIONS = [
    ("A0A6B9MRA0", "glycoprotein-hormone", "peptide"),
    ("A0A6B9MSD4", "glycoprotein-hormone", "peptide"),
    ("Q75W84", "orphan", "peptide"),
    ("P46023", "orphan", "UNRESOLVED-needs-user-call"),
]

CORRECTION_COLUMNS = [
    "accession", "current_family", "proposed_family", "verdict", "reviewed",
    "protein_name", "organism", "taxid", "sequence_length", "class_a_evidence",
    "uniprot_curated_family", "supporting_pmids", "rationale",
]


def run_corrections(args) -> int:
    """Resolve each proposed correction against UniProt and report support.

    Read-only: writes a proposal artifact, never the anchor set.
    """
    with open(args.anchor_tsv, newline="") as handle:
        anchors = {r["accession"]: r for r in csv.DictReader(handle, delimiter="\t")}

    accessions = [a for a, _, _ in PROPOSED_CORRECTIONS]
    entries = fetch_entries_by_accession(accessions)

    rows = []
    for accession, current, proposed in PROPOSED_CORRECTIONS:
        entry = entries[accession]
        name = protein_name(entry)
        curated = curated_family_text(entry)
        _, _, class_a_evidence = verify_class_a(entry)

        in_anchor = anchors.get(accession, {}).get("family", "(not an anchor)")
        if in_anchor != current:
            verdict = f"STALE-PROPOSAL: anchor family is {in_anchor!r}, not {current!r}"
            rationale = "re-check before acting; the anchor set moved"
        elif proposed.startswith("UNRESOLVED"):
            verdict = "UNRESOLVED"
            rationale = ("UniProt assigns only the bare class-A family with no "
                         "subfamily, so any promotion is inference not curation")
        else:
            # The curated subfamily is the authority. It must actually name the
            # subfamily the proposal claims.
            supports = "vasopressin/oxytocin" in curated.lower()
            verdict = "SUPPORTED" if supports else "NOT-SUPPORTED"
            rationale = (f"UniProt curated subfamily {curated!r} is a peptide-receptor "
                         "subfamily; glycoprotein hormones are FSH/LH/TSH/CG, not GnRH"
                         if supports else
                         f"curated family {curated!r} does not name the proposed subfamily")

        rows.append({
            "accession": accession,
            "current_family": current,
            "proposed_family": proposed,
            "verdict": verdict,
            "reviewed": "reviewed" if is_reviewed(entry) else "unreviewed",
            "protein_name": name,
            "organism": entry.get("organism", {}).get("scientificName", ""),
            "taxid": entry.get("organism", {}).get("taxonId", ""),
            "sequence_length": entry.get("sequence", {}).get("length", ""),
            "class_a_evidence": class_a_evidence,
            "uniprot_curated_family": curated,
            "supporting_pmids": ";".join(
                pubmed_id(c) for c in citations(entry) if pubmed_id(c)),
            "rationale": rationale,
        })

    write_tsv(rows, CORRECTION_COLUMNS, args.output)
    print(f"\n[refexp2] correction proposals -> {args.output}\n", file=sys.stderr)
    for row in rows:
        print(f"  {row['accession']:12} {row['current_family']:>21} -> "
              f"{row['proposed_family']:<26} {row['verdict']}", file=sys.stderr)
        print(f"    {row['organism']} | {row['protein_name']}", file=sys.stderr)
        print(f"    curated: {row['uniprot_curated_family']}", file=sys.stderr)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--mode", required=True,
                        choices=("discover", "audit", "corrections"))
    parser.add_argument("--anchor-tsv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    return {"discover": run_discover, "audit": run_audit,
            "corrections": run_corrections}[args.mode](args)


if __name__ == "__main__":
    sys.exit(main())
