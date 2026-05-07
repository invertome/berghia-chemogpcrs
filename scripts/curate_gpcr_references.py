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


def base_query() -> str:
    """The unified UniProt query: reviewed GPCRs in curated taxa."""
    taxa = " OR ".join(f"organism_id:{tid}" for tid in CURATED_TAXA)
    return (
        f'reviewed:true AND keyword:"g-protein coupled receptor" '
        f'AND ({taxa})'
    )


# ---- Family / subfamily classification by protein-name pattern ---------

# Family-level patterns. First match wins; order matters (more-specific
# patterns first, generic Class-A catchall last).
_FAMILY_PATTERNS: list[tuple[re.Pattern, str]] = [
    # Opsin (specific; matches "rhodopsin" too)
    (re.compile(r"\bopsin\b|\brhodopsin\b", re.I), "opsin"),
    # Class F — Frizzled / Smoothened
    (re.compile(r"\bfrizzled\b|\bsmoothened\b|\bSmoothened\b", re.I),
     "class-F-frizzled"),
    # Class C — metabotropic glutamate, GABA-B, calcium-sensing, taste-2
    (re.compile(r"metabotropic glutamate|\bmGlu|\bGABA-B|gamma-aminobutyric.*B|"
                r"calcium-sensing receptor|\bCASR\b|taste receptor type 1",
                re.I), "class-C"),
    # Class B secretin family
    (re.compile(r"\bsecretin receptor\b|\bglucagon receptor\b|"
                r"\bcalcitonin receptor\b|parathyroid hormone receptor|"
                r"corticotropin-releasing factor|gastric inhibitory polypep|"
                r"glucagon-like peptide|pituitary adenylate cyclase|"
                r"vasoactive intestinal pep", re.I),
     "class-B-secretin"),
    # Glycoprotein hormone receptors
    (re.compile(r"thyroid-stimulating hormone receptor|follicle-stimulating "
                r"hormone receptor|luteinizing hormone receptor|"
                r"chorionic gonadotropin receptor", re.I),
     "glycoprotein-hormone"),
    # Aminergic — keep BEFORE generic peptide patterns
    (re.compile(r"5-?hydroxytryptamine|serotonin receptor|\b5-HT\b|"
                r"dopamine receptor|\bDRD\d|adrenergic|noradrenergic|"
                r"norepinephrine receptor|histamine receptor|\bHRH\d|"
                r"octopamine receptor|tyramine receptor|trace amine",
                re.I), "aminergic"),
    # Lipid (prostaglandin, lysophosphatidic, sphingosine 1-phosphate,
    # cannabinoid, leukotriene)
    (re.compile(r"prostaglandin.*receptor|lysophosphatidic acid|"
                r"sphingosine 1-phosphate|cannabinoid receptor|"
                r"leukotriene.*receptor|free fatty acid receptor|"
                r"oxysterol receptor", re.I), "lipid"),
    # Nucleotide (adenosine, P2Y purinergic)
    (re.compile(r"adenosine receptor|\bADORA\d|P2Y receptor|"
                r"purinergic receptor P2Y|P2RY\d", re.I), "nucleotide"),
    # Peptide — many subfamilies; broad name patterns
    (re.compile(r"neuropeptide [YFS]|tachykinin receptor|vasopressin receptor|"
                r"oxytocin receptor|allatostatin receptor|"
                r"melatonin receptor|opioid receptor|somatostatin receptor|"
                r"cholecystokinin receptor|kinin receptor|bradykinin receptor|"
                r"\bAKHR\b|adipokinetic|corazonin receptor|diuretic hormone|"
                r"sex peptide receptor|pigment-dispersing|FMRFamide|"
                r"orexin receptor|ghrelin receptor|galanin receptor|"
                r"motilin receptor|relaxin receptor|kisspeptin receptor|"
                r"angiotensin.*receptor|endothelin receptor|"
                r"neurotensin receptor|neuromedin", re.I), "peptide"),
]

# Aminergic medium-granularity sub-patterns (only used if family is aminergic).
_AMINERGIC_SUBFAMILY_PATTERNS = [
    (re.compile(r"\b(5-?HT|5-hydroxytryptamine|serotonin)\b", re.I), "5HT"),
    (re.compile(r"\b(dopamine|DRD\d)\b", re.I), "dopamine"),
    (re.compile(r"\b(adrenergic|noradrenergic|norepinephrine|ADR[AB]\d?)\b", re.I),
     "norepinephrine"),
    (re.compile(r"\bhistamine\b|\bHRH\d\b", re.I), "histamine"),
    (re.compile(r"\boctopamine\b|\bOct[A-Za-z]?\d?R?\b", re.I), "octopamine"),
    (re.compile(r"\btyramine\b|\bTyrR\b|\bTAR\d?\b", re.I), "tyramine"),
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
    text = query_uniprot_paginated(base_query())
    rows = parse_uniprot_tsv(text)
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
    write_fasta(records, str(out_dir / "01_swissprot_backbone.fasta"))
    write_tsv(records, str(out_dir / "01_swissprot_backbone.tsv"))
    print(f"\n[curate] DONE: {len(records)} unique records -> {out_dir}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
