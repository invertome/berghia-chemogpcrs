#!/usr/bin/env python3
"""census_reviewed_classA_gpcr_families.py — measure the CEILING of curated
class-A GPCR references available per family, and the taxonomic composition
of what is available.

Why this exists
---------------
The embedding-novelty channel scores each Berghia candidate against a
per-family reference envelope, and the envelope percentile is calibrated on
each family's OWN members' pairwise distances. When a family's members are
all vertebrate, that calibration spread is vertebrate-to-vertebrate, which is
far narrower than the mollusc-to-vertebrate distance every real query has to
cross. The result is a systematic, directional inflation of the percentile for
every candidate regardless of true family membership.

Before spending effort expanding a family, we therefore need a MEASURED answer
to a prior question: how many genuinely reviewed class-A references EXIST for
that family at all, and how many of them are non-vertebrate? If the reviewed
ceiling is zero outside Vertebrata, the family cannot be fixed by expansion and
its envelope is uninformative for a molluscan query — which is a finding, not a
shortfall.

What it measures
----------------
For each family, against UniProtKB:
  * ``reviewed_total``      reviewed class-A entries in all of Swiss-Prot
  * ``vertebrate`` / ``non_vertebrate`` / ``mollusca`` / ``lophotrochozoa``
  * the same breakdown for what the anchor set currently holds, with review
    status verified per accession rather than assumed

Curation rules enforced here (violating these has already cost this project):
  * Swiss-Prot **reviewed only** for ceiling counts. TrEMBL entries are counted
    separately and never folded into a "reviewed" number.
  * **Strictly class A**, established from UniProt's curated
    ``G protein-coupled receptor 1 family`` annotation, never inferred from a
    family name.
  * Every accession is **resolved against UniProt**; an accession that fails to
    resolve is an error, never a silent drop.
  * A query returning **zero rows is refused**, never written out. UniProt
    answers a stale field name with HTTP 200 and an empty body, which is
    indistinguishable from success until the artifact has been overwritten.

Usage:
    python3 scripts/census_reviewed_classA_gpcr_families.py \
        --anchor-tsv references/anchors/anchor_set_PROD.tsv \
        --output references/non_chemo_gpcr/reviewed_classA_family_census.tsv

This script is READ-ONLY with respect to every existing reference artifact. It
never edits the anchor set and never mints or renumbers an identifier.
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
import urllib.parse
import urllib.request
from typing import Dict, Iterable, List

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
USER_AGENT = "berghia-chemogpcrs/1.0"

# Class-A membership comes from UniProt's curated protein-family annotation.
# This is a verified property of each record, not an inference from its name.
CLASS_A_QUERY = 'reviewed:true AND family:"g-protein coupled receptor 1 family"'

FIELDS = [
    "accession", "reviewed", "protein_name", "organism_name", "organism_id",
    "gene_names", "length", "protein_families", "lineage",
]

# Family assignment patterns.
#
# PRECEDENCE MATTERS, and it is protein-name-first by design. UniProt packs
# legacy gene aliases into the gene-name field, and several of those aliases
# name a DIFFERENT family than the receptor actually belongs to:
#   * LTB4R  carries the alias ``CMKRL1`` ("chemokine receptor-like 1") but is
#     a leukotriene (lipid) receptor.
#   * LPAR4 / LPAR6 carry the aliases ``P2RY9`` / ``P2RY5`` but are
#     lysophosphatidic acid (lipid) receptors, not nucleotide receptors.
# Matching gene aliases with equal weight to the curated protein name puts all
# of these in the wrong family, silently. So the recommended protein name is
# tested first and wins outright; gene names are consulted only if the name is
# uninformative.
#
# The 'gonadotropin' trap: "Gonadotropin-releasing hormone receptor" (GnRH) is
# a PEPTIDE receptor. The glycoprotein hormones are FSH / LH / TSH / CG. A
# substring match on 'gonadotropin' pulls GnRH receptors into
# glycoprotein-hormone, so the glycoprotein-hormone needle requires the
# specific hormone names and explicitly excludes the releasing-hormone form.
_GNRH = re.compile(r"gonadotropin[- ]releasing|gonadoliberin|\bGNRHR?\b", re.I)

FAMILY_PATTERNS: List[tuple[str, re.Pattern]] = [
    ("chemokine", re.compile(
        r"chemokine receptor|chemokine.*receptor|\bCCR\d|\bCXCR\d|"
        r"\bCX3CR1\b|\bXCR1\b|\bACKR\d", re.I)),
    ("glycoprotein-hormone", re.compile(
        r"thyrotropin receptor|thyroid-stimulating hormone receptor|"
        r"follicle-stimulating hormone receptor|follitropin|"
        r"lutropin|luteinizing hormone receptor|choriogonadotropin|"
        r"glycoprotein hormone", re.I)),
    ("nucleotide", re.compile(
        r"adenosine receptor|purinoceptor|purinergic receptor P2Y|"
        r"\bADORA\d|\bP2RY\d|\bP2Y\d", re.I)),
    ("lipid", re.compile(
        r"prostaglandin|prostacyclin|prostanoid|thromboxane|"
        r"cannabinoid receptor|sphingosine 1-phosphate|"
        r"lysophosphatidic acid|free fatty acid receptor|"
        r"leukotriene|oxoeicosanoid|platelet-activating factor|"
        r"oxysterol|hydroxycarboxylic acid receptor", re.I)),
]

LOPHOTROCHOZOAN_PHYLA = (
    "Mollusca", "Annelida", "Platyhelminthes", "Brachiopoda",
    "Phoronida", "Nemertea", "Rotifera", "Bryozoa",
)


def classify_family(protein_name: str, gene_names: str) -> str:
    """Assign one of the four target families, or '' if none match.

    Protein name is authoritative and is tested alone first; gene names are
    only consulted when the protein name matches nothing. See the module note
    on the CMKRL1 / P2RY5 alias traps for why equal weighting is unsafe.
    """
    if _GNRH.search(protein_name):
        # GnRH is a peptide receptor; never let it reach glycoprotein-hormone.
        return ""
    for family, pattern in FAMILY_PATTERNS:
        if pattern.search(protein_name):
            return family
    for family, pattern in FAMILY_PATTERNS:
        if pattern.search(gene_names):
            return family
    return ""


def phylum_of(lineage: str) -> str:
    """Coarse clade label from a UniProt taxonomic lineage string."""
    for phylum in LOPHOTROCHOZOAN_PHYLA:
        if phylum in lineage:
            return phylum
    if "Vertebrata" in lineage:
        return "Vertebrata"
    for phylum in ("Arthropoda", "Nematoda", "Echinodermata", "Tunicata",
                   "Cephalochordata", "Cnidaria", "Placozoa", "Porifera"):
        if phylum in lineage:
            return phylum
    return "other"


def is_lophotrochozoan(lineage: str) -> bool:
    return any(p in lineage for p in LOPHOTROCHOZOAN_PHYLA)


def _get(url: str) -> tuple[str, str, int]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        total = int(response.headers.get("x-total-results", "-1"))
        return response.read().decode("utf-8"), response.headers.get("Link", ""), total


def fetch_paginated(query: str, fields: Iterable[str] = FIELDS) -> List[dict]:
    """Run a UniProt query, following Link pagination to exhaustion.

    Refuses an empty result and asserts the row count matches the server's own
    ``x-total-results``, so a truncated page can never be mistaken for the
    complete answer.
    """
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
        "size": "500",
    }
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"
    pages: List[str] = []
    total = None
    while url:
        text, link, page_total = _get(url)
        if total is None:
            total = page_total
        pages.append(text)
        match = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = match.group(1) if match else None

    header = pages[0].split("\n", 1)[0]
    body = "".join(p.split("\n", 1)[1] if "\n" in p else "" for p in pages)
    rows = list(csv.DictReader(io.StringIO(header + "\n" + body), delimiter="\t"))

    if not rows:
        raise RuntimeError(
            "UniProt returned zero rows — refusing to treat an empty result as "
            "success. This usually means a field or keyword was renamed. "
            f"Query was: {query}"
        )
    if total is not None and total >= 0 and len(rows) != total:
        raise RuntimeError(
            f"pagination lost rows: parsed {len(rows)} but server reported "
            f"{total} for query: {query}"
        )
    return rows


def fetch_by_accession(accessions: Iterable[str], batch: int = 100) -> Dict[str, dict]:
    """Resolve each accession against UniProt. Unresolvable accessions raise."""
    wanted = list(dict.fromkeys(accessions))
    resolved: Dict[str, dict] = {}
    for start in range(0, len(wanted), batch):
        chunk = wanted[start:start + batch]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        params = {
            "query": query, "format": "tsv",
            "fields": ",".join(FIELDS), "size": "500",
        }
        text, _, _ = _get(f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}")
        for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
            resolved[row["Entry"]] = row
    missing = [a for a in wanted if a not in resolved]
    if missing:
        raise RuntimeError(
            f"{len(missing)} accession(s) did not resolve against UniProt: "
            f"{missing[:10]}"
        )
    return resolved


def census(class_a_rows: List[dict], anchor_rows: List[dict],
           anchor_info: Dict[str, dict]) -> List[dict]:
    """Build the per-family census rows."""
    families = [name for name, _ in FAMILY_PATTERNS]
    out: List[dict] = []

    for family in families:
        available = [
            r for r in class_a_rows
            if classify_family(r["Protein names"], r["Gene Names"]) == family
        ]
        anchors = [
            r for r in anchor_rows
            if r.get("family") == family and r.get("class") == "A"
        ]
        anchor_lineages = [
            anchor_info[r["accession"]]["Taxonomic lineage"] for r in anchors
        ]
        anchor_reviewed = sum(
            1 for r in anchors
            if anchor_info[r["accession"]]["Reviewed"] == "reviewed"
        )
        non_vert = [r for r in available
                    if "Vertebrata" not in r["Taxonomic lineage"]]
        out.append({
            "family": family,
            "swissprot_reviewed_classA_total": len(available),
            "swissprot_vertebrate": len(available) - len(non_vert),
            "swissprot_non_vertebrate": len(non_vert),
            "swissprot_mollusca": sum(
                1 for r in available if "Mollusca" in r["Taxonomic lineage"]),
            "swissprot_lophotrochozoa": sum(
                1 for r in available
                if is_lophotrochozoan(r["Taxonomic lineage"])),
            "anchor_current_total": len(anchors),
            "anchor_reviewed": anchor_reviewed,
            "anchor_trembl": len(anchors) - anchor_reviewed,
            "anchor_mollusca": sum(
                1 for lin in anchor_lineages if "Mollusca" in lin),
            "anchor_lophotrochozoa": sum(
                1 for lin in anchor_lineages if is_lophotrochozoan(lin)),
        })
    return out


COLUMNS = [
    "family",
    "swissprot_reviewed_classA_total", "swissprot_vertebrate",
    "swissprot_non_vertebrate", "swissprot_mollusca",
    "swissprot_lophotrochozoa",
    "anchor_current_total", "anchor_reviewed", "anchor_trembl",
    "anchor_mollusca", "anchor_lophotrochozoa",
]


def write_tsv(rows: List[dict], path: str) -> None:
    """Write atomically so a partial write can never replace a good file."""
    if not rows:
        raise RuntimeError("refusing to write an empty census")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp"
    with open(tmp, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp, path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--anchor-tsv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    print("[census] fetching all reviewed class-A GPCRs from UniProt...",
          file=sys.stderr)
    class_a_rows = fetch_paginated(CLASS_A_QUERY)
    print(f"[census] reviewed class-A universe: {len(class_a_rows)}",
          file=sys.stderr)

    with open(args.anchor_tsv, newline="") as handle:
        anchor_rows = list(csv.DictReader(handle, delimiter="\t"))
    print(f"[census] resolving {len(anchor_rows)} anchor accessions...",
          file=sys.stderr)
    anchor_info = fetch_by_accession(r["accession"] for r in anchor_rows)

    rows = census(class_a_rows, anchor_rows, anchor_info)
    write_tsv(rows, args.output)

    print(f"\n[census] wrote {args.output}\n", file=sys.stderr)
    header = (f"{'family':<24}{'revd':>6}{'vert':>6}{'nonvert':>9}"
              f"{'moll':>6}{'lopho':>7}{'|':>3}{'anch':>6}{'aRevd':>7}"
              f"{'aTrEMBL':>9}{'aMoll':>7}")
    print(header, file=sys.stderr)
    for row in rows:
        print(f"{row['family']:<24}"
              f"{row['swissprot_reviewed_classA_total']:>6}"
              f"{row['swissprot_vertebrate']:>6}"
              f"{row['swissprot_non_vertebrate']:>9}"
              f"{row['swissprot_mollusca']:>6}"
              f"{row['swissprot_lophotrochozoa']:>7}"
              f"{'|':>3}"
              f"{row['anchor_current_total']:>6}"
              f"{row['anchor_reviewed']:>7}"
              f"{row['anchor_trembl']:>9}"
              f"{row['anchor_mollusca']:>7}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
