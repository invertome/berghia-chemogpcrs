#!/usr/bin/env python3
"""curate_gpcr_references.py — Build the curated non-chemoreceptor GPCR
reference set from UniProt Swiss-Prot reviewed entries.

Phase 1 / Task 1.1 of the non-chemoreceptor classification pipeline.

Design (revised after API probing):
    Query UniProt's REST API ONCE for all reviewed GPCRs (keyword
    "G-protein coupled receptor") in the curated taxa (human, mouse,
    rat, Drosophila, C. elegans), paginating via the Link header for
    full coverage. Then partition each entry into our 10 coarse
    families (with medium drill-down for aminergic and peptide) by
    protein-name + gene-name pattern matching.

    This avoids brittle per-family UniProt queries (the `family:` field
    requires exact in-database family names that don't match our
    taxonomy) and gives us full control over how each GPCR maps to our
    classification scheme.

    Excludes ALL invertebrate auto-annotated entries — the curated taxa
    list (5 species, all manually-curated by Swiss-Prot) ensures we only
    see vetted entries. Tasks 1.2-1.4 add invertebrate-specific
    receptors from FlyBase / WormBase / literature where coverage is
    needed.

Usage:
    python3 curate_gpcr_references.py --output-dir references/non_chemo_gpcr/

Output:
    01_swissprot_backbone.fasta  — sequences with `>accession|family|subfamily|species` headers
    01_swissprot_backbone.tsv    — provenance metadata (one row per entry)
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

UNIPROT_REST_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_FIELDS = [
    "accession", "reviewed", "protein_name",
    "organism_name", "gene_names", "length", "sequence",
]
PAGE_SIZE = 500  # UniProt's hard maximum

# Curated taxa: human, mouse, rat, Drosophila melanogaster, C. elegans.
# Restricted to these species because (a) Swiss-Prot has manually reviewed
# them at depth, (b) auto-propagated invertebrate annotations from RefSeq
# are unreliable for chemoreceptor LSE expansions (see project memory
# feedback_no_track_plans.md / non-chemoreceptor classification design).
CURATED_TAXA = {
    "9606": "Homo sapiens",
    "10090": "Mus musculus",
    "10116": "Rattus norvegicus",
    "7227": "Drosophila melanogaster",
    "6239": "Caenorhabditis elegans",
}


# Stable UniProt keyword ACCESSION for "G protein-coupled receptor".
#
# 2026-07-20: the previous query used the keyword's DISPLAY NAME,
# keyword:"g-protein coupled receptor". UniProt renamed KW-0297 to
# "G protein-coupled receptor" (no hyphen after G), so the old string
# returned HTTP 200 with x-total-results: 0 — curate_all() produced an empty
# record list and the writers happily overwrote the reference FASTA/TSV with
# headers only, exit 0. Verified the same day: the display name gives 0 rows
# for human, the accession gives 839. Accessions cannot be renamed; display
# names can, so always query the accession.
GPCR_KEYWORD = "keyword:KW-0297"


def base_query() -> str:
    """The unified UniProt query: reviewed GPCRs in curated taxa."""
    taxa = " OR ".join(f"organism_id:{tid}" for tid in CURATED_TAXA)
    return f"reviewed:true AND {GPCR_KEYWORD} AND ({taxa})"


# ---- Family / subfamily classification by protein-name pattern ---------

# Family-level patterns. First match wins; order matters (more-specific
# patterns first, generic Class-A catchall last).
#
# The patterns use word-keyword matching (not adjacent phrases) because
# protein names have variable adjective insertions: "Tachykinin-like
# PEPTIDES receptor" is still a tachykinin family member, "Cholinergic
# receptor MUSCARINIC 1" is still aminergic, etc.
_FAMILY_PATTERNS: list[tuple[re.Pattern, str]] = [
    # Opsin (specific; matches "rhodopsin" too)
    (re.compile(r"\bopsin\b|\brhodopsin\b", re.I), "opsin"),
    # Class F — Frizzled / Smoothened
    (re.compile(r"\bfrizzled\b|\bsmoothened\b", re.I),
     "class-F-frizzled"),
    # Class C — metabotropic glutamate, GABA-B, calcium-sensing, taste-2
    (re.compile(r"metabotropic glutamate|\bmGlu|\bGABA-B|gamma-aminobutyric.*B|"
                r"calcium-sensing receptor|\bCASR\b|taste receptor type 1",
                re.I), "class-C"),
    # Class B secretin/adhesion — secretin-related peptide hormones AND
    # adhesion-class GPCRs (latrophilin, methuselah, cadherin/stan, boss).
    (re.compile(r"\bsecretin receptor\b|\bglucagon receptor\b|"
                r"\bcalcitonin receptor\b|parathyroid hormone receptor|"
                r"corticotropin-releasing factor|gastric inhibitory polypep|"
                r"glucagon-like peptide|pituitary adenylate cyclase|"
                r"vasoactive intestinal pep|"
                r"\blatrophilin\b|\bmethuselah\b|\badhesion G[- ]?protein|"
                r"\bcadherin\b.*EGF.*LAG|"
                r"\bbride of sevenless\b", re.I),
     "class-B-secretin"),
    # Glycoprotein hormone receptors
    (re.compile(r"thyroid-stimulating hormone|follicle-stimulating hormone|"
                r"luteinizing hormone receptor|chorionic gonadotropin|"
                r"\bLHCGR\b|\bTSHR\b|\bFSHR\b", re.I),
     "glycoprotein-hormone"),
    # Aminergic — keep BEFORE generic peptide patterns. Includes muscarinic
    # acetylcholine (Class A aminergic group per GPCRdb convention) and
    # invertebrate dopamine receptors (Drosophila Dop1R/Dop2R, C. elegans
    # dop-1 .. dop-6, ser-1 .. ser-7, gar-1/2/3 muscarinic, octr-*, tyra-*).
    (re.compile(r"5-?hydroxytryptamine|serotonin receptor|\b5-HT\b|"
                r"\bdopamine\b|\bDRD\d|\bDop\d?R\b|\bdop-\d+\b|"
                r"adrenergic|noradrenergic|norepinephrine receptor|"
                r"histamine receptor|\bHRH\d|"
                r"octopamine|\boctr-\d+\b|"
                r"tyramine|\btyra-\d+\b|"
                r"trace amine|"
                r"muscarinic.*acetylcholine|cholinergic.*muscarinic|"
                r"\bmAChR\b|\bgar-\d\b|"
                r"\bser-\d+\b|\bmod-\d+\b", re.I), "aminergic"),
    # Lipid (prostaglandin, lysophosphatidic, sphingosine 1-phosphate,
    # cannabinoid, leukotriene, FFA, oxysterol).
    (re.compile(r"prostaglandin.*receptor|lysophosphatidic acid|"
                r"sphingosine 1-phosphate|cannabinoid receptor|"
                r"leukotriene.*receptor|free fatty acid receptor|"
                r"oxysterol receptor|thromboxane.*receptor|"
                r"platelet-activating factor receptor", re.I), "lipid"),
    # Nucleotide (adenosine, P2Y purinergic)
    (re.compile(r"adenosine receptor|\bADORA\d|P2Y receptor|"
                r"purinergic receptor P2Y|\bP2RY\d|\bP2Y\d{1,2}\b", re.I),
     "nucleotide"),
    # Peptide — broad keyword matching, not phrase matching. Order: most
    # specific subfamilies first, generic "neuropeptide" / "peptide" last.
    (re.compile(
        r"\btachykinin\b|\bTkR\d|"
        r"\bneuropeptide [YFS]\b|\bNPY\d?R\b|\bNPF[FR]?\b|\bsNPF\b|"
        r"\bvasopressin\b|\boxytocin\b|\bAVPR\d|"
        r"\ballatostatin\b|\bAst[ABC]-?R\b|"
        r"\bmelatonin\b|\bMTNR\d|"
        r"\bopioid\b|\bOPR[KMD]\b|"
        r"\bsomatostatin\b|\bSSTR\d|"
        r"\bcholecystokinin\b|\bCCKR\b|\bCCK[12]R?\b|"
        r"\bkinin\b|\bbradykinin\b|\bBDKRB\d\b|\bLkr\b|"
        r"\bpyrokinin\b|\bPK1-?R\b|"
        r"\bAKHR\b|\badipokinetic\b|"
        r"\bcorazonin\b|\bCrzR\b|"
        r"\bdiuretic hormone\b|\bDh\d+-R\b|"
        r"\bsex peptide\b|"
        r"\bpigment-dispersing\b|\bPdfr\b|"
        r"\bFMRFamide\b|"
        r"\borexin\b|\bhypocretin\b|"
        r"\bghrelin\b|"
        r"\bgalanin\b|"
        r"\bmotilin\b|"
        r"\brelaxin\b|"
        r"\bkisspeptin\b|"
        r"\bangiotensin\b|"
        r"\bendothelin\b|"
        r"\bneurotensin\b|"
        r"\bneuromedin\b|"
        r"\bcrustacean cardioactive\b|\bCCAP\b|"
        r"\bCCHamide\b|\bCCHa\d?-R\b|"
        r"\bcapability receptor\b|\bCapaR\b|"
        r"\bproctolin\b|"
        r"\bSIFamide\b|\bSIFR\b|"
        r"\bmyoinhibitory peptide\b|\bMip-?R\b|"
        r"\bmyosuppressin\b|"
        r"\bnatriuretic peptide\b|"
        r"\bglycoprotein hormone-like\b|"
        r"\bDmsR\b|\bdmsr-\d+\b|"
        r"\bnpr-\d+\b|\bntr-\d+\b|"
        r"\bdaf-3[78]\b|\baex-2\b|"
        r"\btrissin\b|"
        r"\bECTH\b|\bETHR\b|"
        r"\bproctolin\b|"
        # --- invertebrate neuropeptide receptors (2026-07-14 granularity pass) ---
        r"\bpheromone biosynthesis[- ]activating\b|\bPBAN\b|"
        r"\bACP receptor\b|adipokinetic hormone/corazonin|"
        r"\ballatotropin\b|"
        r"\bCAPA\b|\bperiviscerokinin\b|"
        r"\bdiapause hormone\b|"
        r"\belevenin\b|"
        r"\bcardioacceleratory\b|"
        r"\bthyrotropin-releasing hormone\b|\bTRH-?R\b|"
        r"\bachatin\b|"
        r"\bleucokinin\b|"
        r"\bNPYLR\b|\bneuropeptide Y-like\b|"
        r"\bsulfakinin\b|"
        r"\bRYamide\b|\bluqin\b|\bnatalisin\b|\bFLRFamide\b|"
        r"\bTK receptor\b",
        re.I), "peptide"),
]

# Aminergic medium-granularity sub-patterns (only used if family is aminergic).
_AMINERGIC_SUBFAMILY_PATTERNS = [
    (re.compile(r"\b(5-?HT|5-hydroxytryptamine|serotonin)\b", re.I), "5HT"),
    (re.compile(r"\b(dopamine|DRD\d)\b", re.I), "dopamine"),
    (re.compile(r"\b(adrenergic|noradrenergic|norepinephrine|ADR[AB]\d?)\b", re.I),
     "norepinephrine"),
    (re.compile(r"\bhistamine\b|\bHRH\d\b", re.I), "histamine"),
    (re.compile(r"\boctopamine\b|\bOct[A-Za-z]?\d?R?\b", re.I), "octopamine"),
    # Trace amine-associated receptors — MUST be tested before tyramine.
    #
    # 2026-07-20: the tyramine needle used to include `\bTAR\d?\b`, which
    # matches the alternative names UniProt packs into the `protein_name` TSV
    # field: "Trace amine-associated receptor 5 (TaR-5) ...". The shipped
    # aminergic_tyramine.aln had 41 members, 39 of them TAARs and only 2
    # genuine tyramine receptors. The needle contributed zero unique true
    # positives on the curated taxa (all real tyramine receptors match on the
    # word "tyramine" or on \bTyrR\b), so it was pure false-positive surface.
    #
    # TAARs keep an honest label rather than being dropped: TAAR1 is a
    # bona-fide non-olfactory aminergic receptor, while TAAR2-9 are
    # vertebrate olfactory chemoreceptors. Whether the olfactory ones belong
    # in a NON-chemoreceptor reference set at all is a scientific call for
    # the user, not something to decide silently here.
    (re.compile(r"\btrace amine\b", re.I), "trace-amine"),
    (re.compile(r"\btyramine\b|\bTyrR\b", re.I), "tyramine"),
]

# Peptide medium-granularity sub-patterns.
_PEPTIDE_SUBFAMILY_PATTERNS = [
    (re.compile(r"\b(neuropeptide [YF]|NPY\d?R?|NPFR?|sNPF)\b", re.I), "NPY-NPF"),
    (re.compile(r"\btachykinin\b|\bTACR\d\b|\bTkR\d", re.I), "tachykinin"),
    (re.compile(r"\bvasopressin\b|\boxytocin\b|\bAVPR\d", re.I),
     "vasopressin-oxytocin"),
    (re.compile(r"\ballatostatin\b|\bAst[ABC]-?R\b", re.I), "allatostatin"),
    (re.compile(r"\bkinin receptor\b|\bbradykinin\b|\bBDKRB\d\b|\bLkr\b", re.I),
     "kinin"),
    (re.compile(r"\bopioid\b|\bOPR[KMD]\b", re.I), "opioid"),
    (re.compile(r"\bmelatonin\b|\bMTNR\d", re.I), "melatonin"),
    (re.compile(r"\bsomatostatin\b|\bSSTR\d", re.I), "somatostatin"),
    # invertebrate neuropeptide subfamilies (2026-07-14 granularity pass)
    (re.compile(r"pheromone biosynthesis|\bPBAN\b|\bpyrokinin\b", re.I), "PBAN-pyrokinin"),
    (re.compile(r"\bACP receptor\b|adipokinetic hormone/corazonin", re.I), "ACP"),
    (re.compile(r"\bAKH\b|\badipokinetic\b", re.I), "AKH"),
    (re.compile(r"\ballatotropin\b", re.I), "allatotropin"),
    (re.compile(r"\bCAPA\b|\bperiviscerokinin\b|\bcapability\b", re.I), "CAPA"),
    (re.compile(r"\bdiapause hormone\b", re.I), "diapause-hormone"),
    (re.compile(r"\belevenin\b", re.I), "elevenin"),
    (re.compile(r"\bthyrotropin-releasing\b|\bTRH\b", re.I), "TRH"),
    (re.compile(r"\bleucokinin\b|\bLkr\b|\bkinin\b", re.I), "kinin"),
    (re.compile(r"\bproctolin\b", re.I), "proctolin"),
    (re.compile(r"\bSIFamide\b", re.I), "SIFamide"),
    (re.compile(r"\bcorazonin\b|\bCrzR\b", re.I), "corazonin"),
    (re.compile(r"\bmyosuppressin\b", re.I), "myosuppressin"),
    (re.compile(r"\bFMRFamide\b|\bFLRFamide\b", re.I), "FMRFamide"),
    (re.compile(r"\bsulfakinin\b", re.I), "sulfakinin"),
    (re.compile(r"\bachatin\b", re.I), "achatin"),
    (re.compile(r"\bcrustacean cardioactive\b|\bCCAP\b", re.I), "CCAP"),
    (re.compile(r"\bCCHamide\b", re.I), "CCHamide"),
    (re.compile(r"\ballatostatin\b|\bAst[ABC]", re.I), "allatostatin"),
    (re.compile(r"\bcholecystokinin\b|\bCCKR\b|\bCCK[12]R?\b", re.I),
     "cholecystokinin"),
]


def classify_family(protein_name: str, gene_names: str) -> str:
    """Assign a coarse family from the 10-family taxonomy. Returns
    'unclassified-gpcr' if no pattern matches (still a reviewed GPCR
    per the input query, just not in our scoped families)."""
    text = f"{protein_name} {gene_names}"
    for pat, family in _FAMILY_PATTERNS:
        if pat.search(text):
            return family
    return "unclassified-gpcr"


# ---- GPCR class from the underlying record (NOT from the protein name) ------

# UniProt curates GPCR class as "Belongs to the G[-]protein coupled receptor N
# family": 1 = rhodopsin (A), 2 = secretin/adhesion (B), 3 = glutamate (C).
# The hyphenation is inconsistent across entries, so both spellings are matched.
_CURATED_GPCR_FAMILY_RE = re.compile(
    r"Belongs to the G[- ]?protein[- ]coupled receptor (\d+) family", re.I)
_CURATED_FAMILY_TO_CLASS = {"1": "A", "2": "B", "3": "C"}

# Fallback: the SPECIFIC Pfam family. Clan membership is useless here --
# PF00001/PF00002/PF00003/PF02101 all sit in clan CL0192 (GPCR_A), verified
# against InterPro -- so only the family accession discriminates class.
_PFAM_TO_CLASS = {"PF00001": "A", "PF00002": "B", "PF00003": "C"}
_PFAM_TOKEN_RE = re.compile(r"\bPF\d{5}\b")


def gpcr_class_from_evidence(curated_similarity: str, pfam: str) -> str:
    """Resolve a GPCR class (A/B/C) from a UniProt record's own evidence.

    This is the single, uniform rule that decides anchor-set membership of the
    class-A envelope. Precedence:

    1. UniProt's curated "Belongs to the G-protein coupled receptor N family"
       statement (expert-assigned; outranks a profile hit).
    2. The specific Pfam family (PF00001/PF00002/PF00003).
    3. ``UNKNOWN`` -- there is deliberately no default to "A", because a silent
       default is exactly how a class-B receptor ends up polluting a class-A
       family prototype.

    Both arguments come straight from the API (``cc_similarity``, ``xref_pfam``);
    neither the protein name nor any existing label is consulted, so a mislabelled
    row cannot launder itself through this function.
    """
    m = _CURATED_GPCR_FAMILY_RE.search(curated_similarity or "")
    if m:
        klass = _CURATED_FAMILY_TO_CLASS.get(m.group(1))
        if klass:
            return klass
    for token in _PFAM_TOKEN_RE.findall(pfam or ""):
        if token in _PFAM_TO_CLASS:
            return _PFAM_TO_CLASS[token]
    return "UNKNOWN"


def classify_subfamily(family: str, protein_name: str,
                       gene_names: str) -> str:
    """Drill-down for aminergic / peptide; empty string for other families."""
    text = f"{protein_name} {gene_names}"
    if family == "aminergic":
        for pat, label in _AMINERGIC_SUBFAMILY_PATTERNS:
            if pat.search(text):
                return label
    elif family == "peptide":
        for pat, label in _PEPTIDE_SUBFAMILY_PATTERNS:
            if pat.search(text):
                return label
    return ""


# ---- API + parsing -----------------------------------------------------

def query_uniprot(query: str, fields: list[str] | None = None,
                  size: int = PAGE_SIZE, max_retries: int = 3) -> str:
    """Query UniProt REST API; return raw TSV text (single page)."""
    fields = fields or DEFAULT_FIELDS
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
        "size": str(size),
    }
    url = f"{UNIPROT_REST_URL}?{urlencode(params)}"
    return _http_get_with_retry(url, max_retries)


def _http_get_with_retry(url: str, max_retries: int = 3) -> str:
    for attempt in range(max_retries):
        try:
            req = Request(url, headers={"User-Agent": "berghia-chemogpcrs/1.0"})
            with urlopen(req, timeout=60) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            if attempt < max_retries - 1:
                wait = 2 ** attempt
                print(f"WARN: UniProt query failed (attempt {attempt + 1}/{max_retries}): {e}; "
                      f"retrying in {wait}s", file=sys.stderr)
                time.sleep(wait)
            else:
                raise


def query_uniprot_paginated(query: str,
                            fields: list[str] | None = None) -> str:
    """Query UniProt and follow Link header pagination until exhausted.
    Returns concatenated TSV (header from first page only)."""
    fields = fields or DEFAULT_FIELDS
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
        "size": str(PAGE_SIZE),
    }
    url = f"{UNIPROT_REST_URL}?{urlencode(params)}"

    pages: list[str] = []
    while url:
        req = Request(url, headers={"User-Agent": "berghia-chemogpcrs/1.0"})
        with urlopen(req, timeout=60) as resp:
            text = resp.read().decode("utf-8")
            pages.append(text)
            link = resp.headers.get("Link", "")
        # Parse Link header for next page: <URL>; rel="next"
        m = re.search(r"<([^>]+)>;\s*rel=\"next\"", link)
        url = m.group(1) if m else None
    if not pages:
        return ""
    # Keep header from first page; strip headers from subsequent pages
    header = pages[0].split("\n", 1)[0]
    bodies: list[str] = []
    for i, p in enumerate(pages):
        rest = p.split("\n", 1)[1] if "\n" in p else ""
        bodies.append(rest)
    return header + "\n" + "".join(bodies)


def parse_uniprot_tsv(text: str) -> list[dict[str, str]]:
    """Parse UniProt TSV response into list of row dicts.

    Skips rows where Reviewed != 'reviewed' (defensive — query should be
    reviewed-only but this protects against any leakage)."""
    if not text or not text.strip():
        return []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    out: list[dict[str, str]] = []
    for row in reader:
        if row.get("Reviewed", "").strip().lower() != "reviewed":
            continue
        out.append(row)
    return out


def make_curated_record(uniprot_row: dict[str, str], family: str,
                        subfamily: str) -> dict[str, str]:
    """Convert a parsed UniProt row to our canonical record format."""
    organism_full = uniprot_row.get("Organism", "")
    species = re.sub(r"\s*\(.*\)\s*$", "", organism_full).strip()
    gene_names = uniprot_row.get("Gene Names", "").strip()
    primary_gene = gene_names.split()[0] if gene_names else ""
    return {
        "accession": uniprot_row.get("Entry", "").strip(),
        "family": family,
        "subfamily": subfamily,
        "species": species,
        "gene": primary_gene,
        "source": "swissprot",
        "sequence": uniprot_row.get("Sequence", "").strip(),
    }


def curate_all() -> list[dict[str, str]]:
    """Run the full curation: query UniProt, classify each entry, return
    canonical records."""
    query = base_query()
    text = query_uniprot_paginated(query)
    rows = parse_uniprot_tsv(text)
    if not rows:
        # A malformed keyword/field returns HTTP 200 with zero results, which
        # is indistinguishable from success until the reference artifacts have
        # already been overwritten with headers only. Fail loudly instead.
        raise RuntimeError(
            "UniProt returned zero rows for the curation query — refusing to "
            "overwrite the reference set with an empty result. This usually "
            "means a query term stopped matching (e.g. a renamed keyword). "
            f"Query was: {query}"
        )
    records: list[dict[str, str]] = []
    for row in rows:
        family = classify_family(
            row.get("Protein names", ""), row.get("Gene Names", ""))
        if family == "unclassified-gpcr":
            # Keep these but mark — useful for debugging coverage gaps;
            # they're filtered out at the consolidation step (Task 1.6).
            sub = ""
        else:
            sub = classify_subfamily(
                family, row.get("Protein names", ""),
                row.get("Gene Names", ""))
        records.append(make_curated_record(row, family, sub))
    return records


# ---- Backward-compat shims for tests written against the old design ----

def family_queries() -> dict[str, str]:
    """Legacy interface: return per-family queries. With the unified-query
    design, each family's "query" is the base query plus a per-family
    name filter — unused at runtime but kept so existing tests still
    exercise the taxon-restriction invariant."""
    base = base_query()
    families = [
        "aminergic", "peptide", "opsin", "lipid", "nucleotide",
        "metabotropic-neurotransmitter", "glycoprotein-hormone",
        "class-B-secretin", "class-C", "class-F-frizzled",
    ]
    return {f: base for f in families}


def curate_family(family: str, subfamily_hint: str,
                  query: str) -> list[dict[str, str]]:
    """Legacy interface: run a single family query, return canonical
    records. With the unified-query design we just query and filter;
    kept for existing tests."""
    text = query_uniprot(query)
    rows = parse_uniprot_tsv(text)
    out: list[dict[str, str]] = []
    for row in rows:
        if family in ("aminergic", "peptide"):
            sub = classify_subfamily(family, row.get("Protein names", ""),
                                     row.get("Gene Names", ""))
        else:
            sub = subfamily_hint
        out.append(make_curated_record(row, family, sub))
    return out


def dedup_records(records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Drop duplicate accessions; first occurrence wins."""
    seen: set[str] = set()
    out: list[dict[str, str]] = []
    for r in records:
        acc = r["accession"]
        if acc and acc not in seen:
            seen.add(acc)
            out.append(r)
    return out


# ---- Output -----------------------------------------------------------

def write_fasta(records: list[dict[str, str]], path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for r in records:
            header = f">{r['accession']}|{r['family']}|{r['subfamily']}|{r['species']}"
            f.write(f"{header}\n{r['sequence']}\n")


def write_tsv(records: list[dict[str, str]], path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    cols = ["accession", "family", "subfamily", "species", "gene", "source", "length"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in records:
            w.writerow([
                r["accession"], r["family"], r["subfamily"], r["species"],
                r["gene"], r["source"], len(r.get("sequence", "")),
            ])


# ---- Main -------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--output-dir", required=True,
                    help="Directory to write 01_swissprot_backbone.{fasta,tsv}")
    ap.add_argument("--keep-unclassified", action="store_true",
                    help="Keep entries where no family pattern matched "
                         "(default: drop 'unclassified-gpcr' from output).")
    args = ap.parse_args()

    print("[curate] Querying UniProt for all reviewed GPCRs in curated taxa "
          "(human, mouse, rat, D.mel, C.elegans)...", file=sys.stderr)
    records = curate_all()
    print(f"[curate] Retrieved {len(records)} entries before family filter",
          file=sys.stderr)

    # Family breakdown
    by_family: dict[str, int] = {}
    for r in records:
        by_family[r["family"]] = by_family.get(r["family"], 0) + 1
    print("\n[curate] Per-family counts:", file=sys.stderr)
    for fam in sorted(by_family.keys()):
        print(f"  {fam:<35} {by_family[fam]:>4}", file=sys.stderr)

    if not args.keep_unclassified:
        records = [r for r in records if r["family"] != "unclassified-gpcr"]
        print(f"\n[curate] After dropping unclassified-gpcr: {len(records)} records",
              file=sys.stderr)

    records = dedup_records(records)
    out_dir = Path(args.output_dir)
    # Write directly to all_references.{fasta,tsv} — the canonical filename
    # downstream consumers (build_classification_hmms.py, select_backbone_reps.py,
    # classify_via_og_vote.py, 06c orchestrator) expect. Phase 1 was originally
    # designed to layer multiple sources (1.1 SwissProt backbone, 1.2 FlyBase,
    # 1.3 WormBase, 1.4 mollusc literature) into a separate consolidator step
    # producing all_references.*; with 1.2/1.3/1.4 deferred to bead -qgs,
    # SwissProt IS the consolidated set, so we write to the canonical name
    # directly. This eliminates the manual `cp` step that was needed before.
    write_fasta(records, str(out_dir / "all_references.fasta"))
    write_tsv(records, str(out_dir / "all_references.tsv"))
    print(f"\n[curate] DONE: {len(records)} unique records -> {out_dir}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
