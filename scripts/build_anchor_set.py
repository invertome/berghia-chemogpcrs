#!/usr/bin/env python3
"""build_anchor_set.py — Source the three-tier, reliability-graded GPCR
anchor set used to seed per-class reference trees (C2 of the per-class
pool overhaul, bead berghia-chemogpcrs-521.3).

Anchors are well-characterized GPCRs injected into the per-class pools as
must-keep landmarks (functional reference frame) so a Berghia clade can be
read as chemo vs non-chemo by phylogenetic placement. Three tiers, graded
by reliability (bar = experimental characterization); deuterostomes are
excluded by construction.

  Tier 1 — Mollusca (in-group functional landmarks, all reviewed/experimental)
    1a  reviewed Swiss-Prot molluscan receptors (genuine 7TM receptors only —
        the keyword:Transmembrane / NOT keyword:Toxin filter removes cone-snail
        venom peptides that carry the GPCR keyword by target annotation).
    1c  unreviewed-but-experimentally-characterized molluscan receptors:
        transcript-level (existence:2) named-family receptors that cite a
        VERIFIED functional/pharmacology paper (the PMID allow-list in
        references/anchors/tier1c_functional_pmids.tsv). Includes Berghia's own
        Heterobranchia (Aplysia/Lymnaea/Planorbella serotonin & octopamine
        receptors) absent from the reviewed set.
  Tier 2 — Annelida: Platynereis dumerilii deorphanized neuropeptide GPCRs
        (Bauknecht & Jekely 2015, PMID 26190115). The only systematic
        lophotrochozoan deorphanization dataset; nothing phylogenetically
        closer to Mollusca is functionally characterized.
  Tier 3 — Ecdysozoa: reviewed Drosophila + C. elegans GPCRs, read from the
        existing curated non-chemo reference TSV. The only tier with
        Class-B/C/F members; out-group, so eval-gated downstream (C3).

family -> class rule: class-B*->B, class-C*->C, class-F*->F, else A.
Output: references/anchors/anchor_set.{fasta,tsv}; FASTA headers are
ANCHOR_<class>_<tier>_<accession>.

UniProt data: CC BY 4.0 — https://www.uniprot.org/help/license

Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

# Reuse the curated protein-name -> family classifier (no network at import).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from curate_gpcr_references import _FAMILY_PATTERNS  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_REST_URL = "https://rest.uniprot.org/uniprotkb/search"
PAGE_SIZE = 500

MOLLUSCA_TAXID = 6447
PLATYNEREIS_TAXID = 6359
BERGHIA_TAXID = 1287507
PLATYNEREIS_DEORPH_PMID = "26190115"  # Bauknecht & Jekely 2015

# Stable UniProt keyword ACCESSION for "G protein-coupled receptor".
# 2026-07-20: the display-name form keyword:"g-protein coupled receptor"
# returns HTTP 200 with zero results — UniProt renamed KW-0297 to
# "G protein-coupled receptor" (no hyphen after G). Verified the same day:
# display name -> 0 rows for human, KW-0297 -> 839. Accessions cannot be
# renamed, so always query the accession.
GPCR_KEYWORD = "keyword:KW-0297"

# Tier-3 source taxa (reviewed, from the curated non-chemo TSV).
TIER3_TAXA = {
    "Drosophila melanogaster": 7227,
    "Caenorhabditis elegans": 6239,
}

# Named functional families used to gate tier-1c (a specific receptor name, not
# an uncharacterized/predicted ORF). Same families surveyed during sourcing.
TIER1C_NAMED_FAMILIES = [
    "opsin", "rhodopsin", "serotonin", "5-HT", "cholecystokinin", "gastrin",
    "prostaglandin", "dopamine", "octopamine", "GnRH", "neuropeptide F",
]

# Verified tier-1c functional-paper allow-list. Each PMID was obtained
# programmatically (UniProt lit_pubmed cross-refs of named molluscan receptors)
# and verified by reading its title (NCBI esummary) to confirm it is a genuine
# functional/pharmacology/cloning characterization — NOT a transcriptome/genome
# announcement or a toxicology mRNA-expression survey (those were dropped).
# Embedded here (not only in references/anchors/tier1c_functional_pmids.tsv) so
# the script is self-contained on Unity, where references/ is not git-tracked.
TIER1C_FUNCTIONAL_PMIDS = frozenset({
    "8891606",   # 1996 Functional characterisation of a 5-HT2 receptor (Lymnaea)
    "10677541",  # 2000 Heterologously expressed octopamine receptor (Aplysia)
    "18310116",  # 2008 Two cloned serotonin receptors (Helisoma/Planorbella)
    "19509343",  # 2009 Light perception in a bioluminescent organ (opsin)
    "19706550",  # 2009 Serotonin receptor coupled to adenylyl cyclase (Aplysia)
    "20100484",  # 2010 Serotonin receptor role in spawning (scallop)
    "20392722",  # 2010 Distributed light sensing in cuttlefish skin (opsin)
    "20964689",  # 2010 PKC regulation by serotonin receptors (Aplysia)
    "23274282",  # 2013 Cloning/characterization of dopamine-like receptor (oyster)
    "25463295",  # 2015 PTGER4 in host immune protection (oyster)
    "25994633",  # 2015 Light-activated chromatophore expansion / opsins (Octopus)
    "26002349",  # 2015 Cephalopod visual-system evolution (opsins)
    "26702352",  # 2015 Chiton larval eyespot opsins
    "27992549",  # 2016 Functional analysis of octopamine/tyramine receptor (oyster)
    "30401878",  # 2018 CCK/sulfakinin signalling system in Lophotrochozoa
    "30849410",  # 2019 Functional properties of 5-HT receptors (abalone)
    "31696430",  # 2020 Serotonin receptor characterization (abalone)
    "32545589",  # 2020 GnRH receptor characterization (abalone)
})


# ---------------------------------------------------------------------------
# Pure logic
# ---------------------------------------------------------------------------

def family_to_class(family: str) -> str:
    """Map a family label to a GPCR class (class-B*->B, class-C*->C,
    class-F*->F, everything else -> A)."""
    f = family or ""
    if f.startswith("class-B"):
        return "B"
    if f.startswith("class-C"):
        return "C"
    if f.startswith("class-F"):
        return "F"
    return "A"


def anchor_header(klass: str, tier: str, accession: str) -> str:
    """FASTA header for an anchor: ANCHOR_<class>_<tier>_<accession>."""
    return f"ANCHOR_{klass}_{tier}_{accession}"


def classify_family(protein_name: str) -> str:
    """Classify a protein name into a family via the curated patterns;
    names matching no specific family resolve to a generic Class-A bucket."""
    for pat, fam in _FAMILY_PATTERNS:
        if pat.search(protein_name or ""):
            return fam
    return "class-A-other"


def load_functional_pmids(manifest_path: str = None) -> set[str]:
    """Return the verified tier-1c functional-paper PMID allow-list.

    With no ``manifest_path``, returns the embedded TIER1C_FUNCTIONAL_PMIDS
    (the authoritative, version-controlled list). A manifest TSV (with a 'pmid'
    column) overrides it — useful for auditing or extending the set.
    """
    if not manifest_path:
        return set(TIER1C_FUNCTIONAL_PMIDS)
    pmids: set[str] = set()
    with open(manifest_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            p = str(row.get("pmid", "")).strip()
            if p.isdigit():
                pmids.add(p)
    return pmids


def _evidence_priority(rec: dict) -> tuple[int, int]:
    """Lower is better: reviewed beats experimental, then lower tier number."""
    rev = 0 if rec.get("evidence") == "reviewed" else 1
    try:
        tier = int(rec.get("tier", "9"))
    except ValueError:
        tier = 9
    return (rev, tier)


def dedup_by_accession(records: list[dict]) -> list[dict]:
    """Collapse duplicate accessions, preferring reviewed over experimental
    then the lower tier number. Order follows first appearance."""
    best: dict[str, dict] = {}
    order: list[str] = []
    for r in records:
        acc = r["accession"]
        if acc not in best:
            best[acc] = r
            order.append(acc)
        elif _evidence_priority(r) < _evidence_priority(best[acc]):
            best[acc] = r
    return [best[a] for a in order]


# ---------------------------------------------------------------------------
# UniProt search (position-based TSV parse; robust to display-header names)
# ---------------------------------------------------------------------------

def uniprot_search(query: str, fields: list[str], *,
                   page_size: int = PAGE_SIZE, max_retries: int = 3) -> list[dict]:
    """Query UniProtKB and return rows as dicts keyed by the requested
    programmatic field names. Columns come back in field order, so we zip by
    position rather than relying on UniProt's display headers."""
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
        "size": str(page_size),
    }
    url = f"{UNIPROT_REST_URL}?{urlencode(params)}"
    out: list[dict] = []
    first = True
    while url:
        text = ""
        for attempt in range(max_retries):
            try:
                req = Request(url, headers={"User-Agent": "berghia-chemogpcrs/1.0"})
                with urlopen(req, timeout=60) as resp:
                    text = resp.read().decode("utf-8")
                    link = resp.headers.get("Link", "")
                break
            except Exception:
                if attempt == max_retries - 1:
                    raise
        lines = text.split("\n")
        body = lines[1:] if first else lines[1:]  # always skip the header row
        first = False
        for line in body:
            if not line.strip():
                continue
            vals = line.split("\t")
            out.append({f: (vals[i] if i < len(vals) else "")
                        for i, f in enumerate(fields)})
        m = re.search(r"<([^>]+)>;\s*rel=\"next\"", link)
        url = m.group(1) if m else None
    return out


# ---------------------------------------------------------------------------
# Tier sourcing
# ---------------------------------------------------------------------------

_SOURCE_FIELDS = ["accession", "protein_name", "organism_id", "organism_name",
                  "sequence"]
_SOURCE_FIELDS_LIT = _SOURCE_FIELDS + ["lit_pubmed_id"]


def _make_record(row: dict, *, tier: str, evidence: str) -> dict:
    family = classify_family(row.get("protein_name", ""))
    try:
        taxid = int(row.get("organism_id") or 0)
    except ValueError:
        taxid = 0
    return {
        "accession": row.get("accession", ""),
        "tier": tier,
        "taxid": taxid,
        "species": row.get("organism_name", ""),
        "family": family,
        "klass": family_to_class(family),
        "evidence": evidence,
        "sequence": (row.get("sequence", "") or "").strip(),
    }


def source_tier1a(search_fn=uniprot_search) -> list[dict]:
    """Tier 1a — reviewed genuine molluscan 7TM receptors (venom peptides
    excluded via keyword:Transmembrane / NOT keyword:Toxin)."""
    query = (
        f"taxonomy_id:{MOLLUSCA_TAXID} AND {GPCR_KEYWORD} "
        f"AND keyword:Transmembrane AND NOT keyword:Toxin "
        f"AND reviewed:true AND NOT taxonomy_id:{BERGHIA_TAXID}"
    )
    return [_make_record(r, tier="1", evidence="reviewed")
            for r in search_fn(query, _SOURCE_FIELDS)]


def source_tier1c(functional_pmids: set[str],
                  search_fn=uniprot_search) -> list[dict]:
    """Tier 1c — unreviewed-but-experimentally-characterized molluscan
    receptors: transcript-level named-family receptors citing a verified
    functional-paper PMID."""
    if not functional_pmids:
        return []
    named = " OR ".join(f'protein_name:"{n}"' for n in TIER1C_NAMED_FAMILIES)
    pmid_filter = " OR ".join(f"lit_pubmed:{p}" for p in sorted(functional_pmids))
    query = (
        f"taxonomy_id:{MOLLUSCA_TAXID} AND {GPCR_KEYWORD} "
        f"AND reviewed:false AND existence:2 "
        f"AND ({named}) AND ({pmid_filter}) "
        f"AND NOT taxonomy_id:{BERGHIA_TAXID}"
    )
    out: list[dict] = []
    for r in search_fn(query, _SOURCE_FIELDS_LIT):
        cited = {p for p in re.split(r"[;,\s]+", r.get("lit_pubmed_id", ""))
                 if p.isdigit()}
        hit = cited & functional_pmids
        if not hit:
            continue
        evidence = f"experimental:{min(hit, key=int)}"
        out.append(_make_record(r, tier="1", evidence=evidence))
    return out


def source_tier2(search_fn=uniprot_search) -> list[dict]:
    """Tier 2 — Platynereis dumerilii deorphanized neuropeptide GPCRs."""
    query = (
        f"taxonomy_id:{PLATYNEREIS_TAXID} AND {GPCR_KEYWORD} "
        f"AND lit_pubmed:{PLATYNEREIS_DEORPH_PMID}"
    )
    evidence = f"experimental:{PLATYNEREIS_DEORPH_PMID}"
    return [_make_record(r, tier="2", evidence=evidence)
            for r in search_fn(query, _SOURCE_FIELDS)]


def parse_tier3_from_curated_tsv(tsv_path: str, fasta_path: str,
                                 taxa: dict[str, int] = None) -> list[dict]:
    """Tier 3 — reviewed Drosophila + C. elegans GPCRs from the curated
    non-chemo reference set (already family-labelled)."""
    taxa = taxa or TIER3_TAXA
    seqs: dict[str, str] = {}
    with open(fasta_path) as fh:
        acc = None
        buf: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf)
                acc = line[1:].split("|")[0].split()[0]
                buf = []
            else:
                buf.append(line)
        if acc is not None:
            seqs[acc] = "".join(buf)

    out: list[dict] = []
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            species = row.get("species", "")
            if species not in taxa:
                continue
            family = row.get("family", "")
            acc = row.get("accession", "")
            out.append({
                "accession": acc,
                "tier": "3",
                "taxid": taxa[species],
                "species": species,
                "family": family,
                "klass": family_to_class(family),
                "evidence": "reviewed",
                "sequence": seqs.get(acc, ""),
            })
    return out


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_anchor_set(records: list[dict], out_dir: str) -> None:
    """Write anchor_set.fasta (ANCHOR_<class>_<tier>_<accession> headers) and
    anchor_set.tsv (provenance)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "anchor_set.fasta", "w") as fh:
        for r in records:
            header = anchor_header(r["klass"], r["tier"], r["accession"])
            fh.write(f">{header}\n{r['sequence']}\n")
    with open(out / "anchor_set.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "tier", "taxid", "species", "family", "class",
                    "evidence"])
        for r in records:
            w.writerow([r["accession"], r["tier"], r["taxid"], r["species"],
                        r["family"], r["klass"], r["evidence"]])


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def build_anchor_set(out_dir: str, curated_tsv: str, curated_fasta: str,
                     pmid_manifest: str = None, *,
                     search_fn=uniprot_search) -> list[dict]:
    """Source all tiers, de-dup by accession, write outputs, return records."""
    pmids = load_functional_pmids(pmid_manifest)
    records: list[dict] = []
    records += source_tier1a(search_fn)
    records += source_tier1c(pmids, search_fn)
    records += source_tier2(search_fn)
    records += parse_tier3_from_curated_tsv(curated_tsv, curated_fasta)
    records = [r for r in records if r["accession"] and r["sequence"]]
    records = dedup_by_accession(records)
    write_anchor_set(records, out_dir)
    return records


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_args_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Source the three-tier GPCR anchor set "
                    "(references/anchors/anchor_set.{fasta,tsv}).")
    p.add_argument("--out-dir", default="references/anchors",
                   help="Output directory (default: references/anchors)")
    p.add_argument("--curated-tsv",
                   default="references/non_chemo_gpcr/all_references.tsv",
                   help="Curated non-chemo reference TSV (tier-3 source)")
    p.add_argument("--curated-fasta",
                   default="references/non_chemo_gpcr/all_references.fasta",
                   help="Curated non-chemo reference FASTA (tier-3 source)")
    p.add_argument("--pmid-manifest", default=None,
                   help="Optional tier-1c functional-paper PMID allow-list TSV "
                        "(overrides the embedded TIER1C_FUNCTIONAL_PMIDS default)")
    return p


def main(argv=None) -> None:
    args = build_args_parser().parse_args(argv)
    records = build_anchor_set(
        out_dir=args.out_dir,
        curated_tsv=args.curated_tsv,
        curated_fasta=args.curated_fasta,
        pmid_manifest=args.pmid_manifest,
    )
    from collections import Counter
    by_tier = Counter(r["tier"] for r in records)
    by_class = Counter(r["klass"] for r in records)
    print(f"[build_anchor_set] {len(records)} anchors -> {args.out_dir}/anchor_set.{{fasta,tsv}}",
          file=sys.stderr)
    print(f"  by tier: {dict(sorted(by_tier.items()))}", file=sys.stderr)
    print(f"  by class: {dict(sorted(by_class.items()))}", file=sys.stderr)


if __name__ == "__main__":
    main()
