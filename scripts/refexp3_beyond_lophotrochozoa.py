#!/usr/bin/env python3
"""refexp3_beyond_lophotrochozoa.py -- widen the class-A reference envelope
BEYOND Lophotrochozoa, into the rest of the non-vertebrate Metazoa.

Why this exists, measured
------------------------
The novelty envelope scores a Berghia candidate against the spread of a
family's OWN members. That spread is 86% vertebrate and arthropod, so it
encodes one lineage's idiosyncrasy rather than the family's real variation. A
calibration run put a number on the consequence: a random molluscan class-A
receptor sits at or beyond the 100th percentile for six of seven families, and
the current rule calls 22% of ordinary molluscs "inside" peptide.

refexp2 attacked this by adding lophotrochozoan references and is exhausted:
296 candidates yielded 36, and the three families with the most unreviewed
molluscan material (LGR, adenosine, prostanoid) yielded ZERO, because 97.4% of
molluscan class-A entries are PE3 -- inferred from homology, not observed.

The remaining move is not more references, it is broader NON-VERTEBRATE
phylogenetic coverage, so each family's spread reflects the family rather than
the vertebrates in it. Adding more mammals would deepen the exact bias being
corrected. Vertebrata is therefore OUT OF SCOPE here by construction, and an
assertion enforces it rather than trusting the taxonomy query.

Clades swept, and why each is worth the query
---------------------------------------------
Ecdysozoa        sister clade to Lophotrochozoa; shares protostome-specific
                 receptor features with molluscs. ~220 ecdysozoan entries are
                 ALREADY anchors, so this script reports new-vs-present rather
                 than a gross count that would double-count them.
Cnidaria +       pre-bilaterian, so informative precisely for the ancient
non-bilaterians  families (LGR, adenosine, prostanoid, opsin) where the
                 vertebrate radiation is derived.
non-vertebrate   echinoderms, hemichordates, amphioxus, tunicates.
deuterostomes    Topologically no closer to molluscs than vertebrates are, but
                 far less derived, so they may retain ancestral receptor
                 character the vertebrate radiation lost.
Xenacoelomorpha  early-branching bilaterian; swept for completeness.

The gate is NOT reinvented here
-------------------------------
Every pass/fail decision is made by refexp2_evidence_gate, imported as a
module. This script contributes the taxonomic sweep, the clade attribution and
the breadth accounting -- nothing else. If the gate's strictness is to change,
it changes in ONE place, in the file the user already approved.

Two measured calibrations carried forward from that run, restated because they
are what keep this from being a bar-lowering exercise:
  * ECO:0000269 is applied by Swiss-Prot curators and is structurally absent
    from TrEMBL (0 of 363 unreviewed candidates carried it). Requiring it of
    unreviewed entries is arithmetically identical to requiring reviewed
    status. For those entries it is replaced by -- not waived in favour of --
    the literature and naming gates.
  * Titles are not evidence; PubMed indexing is. Q6HA06 reads "Identification
    and characterization of a glycoprotein hormone receptor from the oyster
    Crassostrea gigas" and was never published.

One deliberate difference from refexp2: the SEARCH universe
-----------------------------------------------------------
refexp2 drew its universe from the curated family string alone. Measured
against UniProt on 2026-07-20, that string misses a third of the eligible
records: Lophotrochozoa PE1/PE2 returns 381 entries by family string but 572
by (family string OR PF00001 OR IPR000276). The curated string is a SIMILARITY
comment that TrEMBL entries often lack entirely.

This script therefore searches the union. That widens what is LOOKED AT, not
what is ADMITTED: verify_class_a still re-derives class-A membership per entry
from the rhodopsin clan CL0192 domains, and still rejects a CL0176
Chemosens_recp hit outright. A wider net feeding an unchanged gate is the
correct direction; a narrower net feeding a looser gate would be the failure.

This script never edits the anchor set, never merges anything, never mints and
never renumbers an identifier. It writes a verdicts table and a FASTA of what
passes; merging is scoring-affecting and is the user's call.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import refexp2_evidence_gate as gate  # noqa: E402

UNIPROT_TAXONOMY = "https://rest.uniprot.org/taxonomy"

# Every taxid below was resolved against the UniProt taxonomy API rather than
# recalled, and verify_clade_taxids re-queries each one at runtime and ABORTS on
# any mismatch. A wrong taxid here would silently sweep the wrong clade and
# report a confident number for it -- the LSE-clade incident in this project was
# exactly that failure (two "clade" taxids were bacteria).
CLADES: List[Tuple[str, str, int]] = [
    # (clade label, expected scientific name, taxid)
    ("Ecdysozoa",        "Ecdysozoa",        1206794),
    ("Cnidaria",         "Cnidaria",         6073),
    ("Placozoa",         "Placozoa",         10226),
    ("Porifera",         "Porifera",         6040),
    ("Ctenophora",       "Ctenophora",       10197),
    ("Echinodermata",    "Echinodermata",    7586),
    ("Hemichordata",     "Hemichordata",     10219),
    ("Cephalochordata",  "Cephalochordata",  7735),
    ("Tunicata",         "Tunicata",         7712),
    ("Xenacoelomorpha",  "Xenacoelomorpha",  1312402),
]

# Out of scope, and asserted against rather than assumed absent.
VERTEBRATA_TAXID = 7742

# The clade whose non-empty result proves the query TEMPLATE still works, so a
# zero from Porifera can be reported as a real biological zero rather than a
# silently broken keyword. Measured 2026-07-20: 3226 entries.
CONTROL_CLADE = "Ecdysozoa"

PE_FILTER = "(existence:1 OR existence:2)"
CLASS_A_UNIVERSE = (
    f'(family:"{gate.CURATED_CLASS_A_FAMILY}" '
    f"OR xref:pfam-{gate.PFAM_7TM_1} "
    f"OR xref:interpro-{gate.INTERPRO_RHODOPSIN})"
)


def universe_query(taxid: int) -> str:
    return f"(taxonomy_id:{taxid}) AND {CLASS_A_UNIVERSE} AND {PE_FILTER}"


def diagnostic_query(taxid: int) -> str:
    """Same universe WITHOUT the PE filter.

    A clade that yields nothing is a legitimate result, but 'no entries at all'
    and 'entries exist and every one is PE3 homology-inference' are different
    findings and the report must not conflate them.
    """
    return f"(taxonomy_id:{taxid}) AND {CLASS_A_UNIVERSE}"


def verify_clade_taxids() -> None:
    """Re-query every taxid and abort on any name mismatch.

    Obtain-and-verify, not recall: this is the guard that would have caught the
    bacterial-taxid incident.
    """
    for label, expected, taxid in CLADES + [("Vertebrata", "Vertebrata",
                                             VERTEBRATA_TAXID)]:
        body, _, _ = gate._get(f"{UNIPROT_TAXONOMY}/{taxid}?format=json")
        import json as _json
        actual = _json.loads(body).get("scientificName", "")
        if actual != expected:
            raise RuntimeError(
                f"taxid {taxid} resolves to {actual!r}, expected {expected!r} "
                f"for clade {label} -- refusing to sweep the wrong clade")
    print(f"[refexp3] verified {len(CLADES) + 1} clade taxids against UniProt",
          file=sys.stderr)


def fetch_json_allow_empty(query: str, size: int = 200) -> Tuple[List[dict], int]:
    """Paginated JSON fetch that permits a genuine zero.

    gate.fetch_json_paginated refuses an empty result outright, which is right
    for a query that MUST match but wrong here: Porifera legitimately has no
    PE1/PE2 class-A entry, and turning that into a crash would prevent the
    report from stating a real biological zero.

    The safety it gives up is restored two ways: the server's own
    x-total-results must equal the parsed row count (so a truncated page cannot
    pass as complete), and the caller runs a CONTROL clade whose non-empty
    result proves the query template still resolves.
    """
    params = {"query": query, "format": "json", "size": str(size)}
    url = f"{gate.UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"
    results: List[dict] = []
    total: Optional[int] = None
    while url:
        body, link, page_total = gate._get(url)
        if total is None:
            total = page_total
        import json as _json
        results.extend(_json.loads(body).get("results", []))
        match = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = match.group(1) if match else None
    if total is not None and total >= 0 and len(results) != total:
        raise RuntimeError(
            f"pagination lost rows: parsed {len(results)} but server reported "
            f"{total} for query: {query}")
    return results, (total if total is not None else len(results))


def count_only(query: str) -> int:
    """Row count without downloading the rows."""
    params = {"query": query, "format": "json", "size": "1", "fields": "accession"}
    _, _, total = gate._get(f"{gate.UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}")
    return total


def is_vertebrate(entry: dict) -> bool:
    return "Vertebrata" in gate.lineage(entry)


def clade_of(entry: dict, queried_clade: str) -> str:
    """Which swept clade an entry belongs to.

    Falls back to the clade that was queried, so an entry can never be silently
    attributed to a clade it was not found under.
    """
    names = gate.lineage(entry)
    for label, expected, _ in CLADES:
        if expected in names:
            return label
    return queried_clade


def phylum_within_clade(entry: dict) -> str:
    """A finer label for reporting inside the big clades.

    Ecdysozoa spans Arthropoda / Nematoda / Tardigrada / Priapulida, and
    reporting 'Ecdysozoa: N' alone would hide whether the yield is all
    Drosophila again.
    """
    names = gate.lineage(entry)
    for phylum in ("Arthropoda", "Nematoda", "Tardigrada", "Priapulida",
                   "Onychophora", "Nematomorpha", "Kinorhyncha", "Loricifera",
                   "Cnidaria", "Placozoa", "Porifera", "Ctenophora",
                   "Echinodermata", "Hemichordata", "Chordata",
                   "Xenacoelomorpha"):
        if phylum in names:
            return phylum
    return "unresolved"


# ---------------------------------------------------------------------------
# Reporting layer
#
# Everything below LABELS survivors. None of it admits or rejects an entry --
# the gate has already spoken by this point. It exists because refexp2's
# reporting vocabulary was built on lophotrochozoan and vertebrate material and
# is measurably blind to this universe, and because a raw survivor count would
# hide two things the user must see before deciding whether to merge.
# ---------------------------------------------------------------------------

# refexp2's family_for_report needles miss a large, systematic slice of this
# universe. Measured on this run's survivors:
#   * 58 opsins named "Ultraviolet-sensitive visual pigment" / "Blue-sensitive
#     visual pigment" -- its needle is \bopsin\b|\brhodopsin, and the entire
#     arthropod visual literature says "visual pigment" instead.
#   * ~60 arthropod peptide receptors (SIFamide, myosuppressin, AstA/AstC,
#     ecdysis-triggering hormone, CCAP, inotocin, sex peptide, ACP) -- its
#     peptide needle is a vertebrate/lophotrochozoan vocabulary that simply
#     does not contain the ecdysozoan peptide names.
#   * LGR entries named "Leucine-rich repeat-containing G-protein-coupled
#     receptor" -- the glycoprotein-hormone needle looks for FSH/LH/TSH.
#
# Those misses land on exactly the families the sweep is meant to inform (LGR,
# peptide, opsin), so leaving them in 'unclassified' would understate the
# result. This extension is LOCAL and REPORTING-ONLY, mirroring the precedent
# refexp2 set: widening the production classifier changes what the
# non-chemoreceptor reference set contains, which is a pipeline behaviour
# change and the user's call, so it is reported here rather than made there.
_FAMILY_EXTENSION = [
    ("opsin", re.compile(
        r"visual pigment|photopigment|\bopsin|rhodopsin|"
        r"(?:ultraviolet|blue|violet|green|red|long[- ]wavelength|"
        r"short[- ]wavelength|middle[- ]wavelength)[- ]sensitive", re.I)),
    ("glycoprotein-hormone", re.compile(
        r"leucine[- ]rich repeat[- ]?(?:containing)? ?G[- ]?protein|\bLGR\d?\b|"
        r"bursicon|glycoprotein hormone", re.I)),
    ("peptide", re.compile(
        r"SIFamide|SIFa\b|myosuppressin|allatostatin|\bAst[ABC]\b|"
        r"ecdysis[- ]triggering|\bETH\b|crustacean cardioactive|\bCCAP\b|"
        r"cardioacceleratory|inotocin|sex peptide|\bACP receptor|"
        r"adipokinetic|\bAKH\b|corazonin|proctolin|pyrokinin|"
        r"diuretic hormone|CAPA|leucokinin|\bNPF\b|neuropeptide F|"
        r"bombyxin|eclosion hormone|trissin|natalisin|elevenin|"
        r"\bTK receptor|tachykinin|sulfakinin|FMRFamide|FLRFamide|"
        r"gonadotropin[- ]releasing|\bGnRH\b|\bNPS/CCAP|relaxin|"
        r"neuroparsin|orcokinin|RYamide|CCHamide|pigment[- ]dispersing", re.I)),
    ("aminergic", re.compile(r"octopamine|tyramine|\bOAMB\b|dopamine|serotonin", re.I)),
]


def family_for_report_extended(name: str) -> str:
    """refexp2's label, falling back to the ecdysozoan/visual vocabulary above.

    refexp2's classifier is consulted FIRST so that any entry it already labels
    keeps the same label it would have had in the previous sweep, and the two
    runs stay comparable.
    """
    base = gate.family_for_report(name)
    if base != "unclassified":
        return base
    for family, pattern in _FAMILY_EXTENSION:
        if pattern.search(name or ""):
            return family
    return "unclassified"


# Names that cleared S4 while still committing to no ligand and no function.
# These are a measured LEAK, not a design: refexp2's placeholder patterns were
# tuned on 'Orphan GPCR N' and 'GCR###' and never saw a UniProt automatic name
# like 'G-protein coupled receptors family 1 profile domain-containing
# protein'. 10 of 1133 survivors match. They are FLAGGED rather than silently
# dropped, because silently dropping them would be this script editing the
# approved gate's verdict from the outside.
_RESIDUAL_PLACEHOLDER = re.compile(
    r"^(?:G[- ]?protein[- ]coupled receptors? family \d profile "
    r"domain[- ]containing protein"
    r"|7[- ]?transmembrane protein.*"
    r"|Peptide receptor GPCR"
    r"|Rhodopsin[- ]?like .*domain[- ]containing.*"
    r"|G[- ]?protein[- ]coupled receptor.*domain[- ]containing.*"
    r"|[A-Z]{2}\d{4,}p)$",
    re.I,
)


def is_residual_placeholder(name: str) -> bool:
    return bool(_RESIDUAL_PLACEHOLDER.match((name or "").strip()))


VERDICT_COLUMNS = list(gate.VERDICT_COLUMNS) + [
    "clade", "phylum_within_clade", "already_anchor",
    "family_extended", "residual_placeholder", "min_paper_entry_count",
]


def load_anchor_accessions(path: str) -> Tuple[set, List[dict]]:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {r["accession"] for r in rows}, rows


def sweep(anchor_tsv: str, output: str) -> int:
    verify_clade_taxids()
    anchor_accessions, anchor_rows = load_anchor_accessions(anchor_tsv)
    print(f"[refexp3] anchor set: {len(anchor_accessions)} accessions "
          f"({sum(1 for r in anchor_rows if r.get('class') == 'A')} class A)",
          file=sys.stderr)

    per_clade: List[dict] = []
    all_verdicts: List[dict] = []
    all_entries: Dict[str, dict] = {}
    control_seen = 0

    for label, _expected, taxid in CLADES:
        entries, total = fetch_json_allow_empty(universe_query(taxid))
        no_pe_filter = count_only(diagnostic_query(taxid))
        if label == CONTROL_CLADE:
            control_seen = len(entries)

        vertebrate_contamination = [e for e in entries if is_vertebrate(e)]
        if vertebrate_contamination:
            raise RuntimeError(
                f"{label}: {len(vertebrate_contamination)} entries carry "
                f"Vertebrata in their lineage -- Vertebrata is out of scope and "
                f"the taxonomy filter did not hold")

        already = [e for e in entries
                   if e["primaryAccession"] in anchor_accessions]
        new = [e for e in entries
               if e["primaryAccession"] not in anchor_accessions]

        verdicts = [gate.gate_entry(e) for e in new]
        for verdict, entry in zip(verdicts, new):
            verdict["clade"] = clade_of(entry, label)
            verdict["phylum_within_clade"] = phylum_within_clade(entry)
            verdict["already_anchor"] = "no"
            all_entries[verdict["accession"]] = entry
        all_verdicts.extend(verdicts)

        counts = gate.funnel(verdicts)
        per_clade.append({
            "clade": label,
            "taxid": taxid,
            "universe_any_pe": no_pe_filter,
            "universe_pe1_pe2": total,
            "already_anchor": len(already),
            "new_candidates": len(new),
            **{f"pass_{stage}": counts[stage] for stage in gate.GATE_STAGES},
            "passing": counts[gate.GATE_STAGES[-1]],
        })
        print(f"[refexp3] {label:<18} anyPE={no_pe_filter:>5} PE1/2={total:>5} "
              f"anchored={len(already):>4} new={len(new):>5} "
              f"pass={counts[gate.GATE_STAGES[-1]]:>4}", file=sys.stderr)
        time.sleep(0.3)

    if control_seen == 0:
        raise RuntimeError(
            f"control clade {CONTROL_CLADE} returned zero entries -- the query "
            f"template is broken and every zero in this run is meaningless")
    print(f"[refexp3] control clade {CONTROL_CLADE} returned {control_seen} "
          f"entries; zeros elsewhere are real", file=sys.stderr)

    # Abstracts drive ONLY the deorphanized tier, never a pass/fail decision
    # (gate.primary_literature reads titles by design). Fetching them for the
    # survivors alone is therefore exactly equivalent to fetching them for the
    # whole universe, and avoids thousands of needless E-utilities calls.
    survivors = [v for v in all_verdicts if not v["failed_at"]]
    pmids = {p for v in survivors for p in v["primary_pmids"].split(";") if p}
    print(f"[refexp3] fetching {len(pmids)} abstracts for {len(survivors)} "
          f"survivors (tiering only)...", file=sys.stderr)
    abstracts = gate.fetch_pubmed_abstracts(pmids) if pmids else {}
    for verdict in survivors:
        retiered = gate.gate_entry(all_entries[verdict["accession"]], abstracts)
        if retiered["failed_at"]:
            raise RuntimeError(
                f"{verdict['accession']}: re-gating with abstracts changed the "
                f"verdict to {retiered['failed_at']!r}; abstracts must not "
                f"affect pass/fail")
        verdict["deorphanized"] = retiered["deorphanized"]
        verdict["evidence_tier"] = retiered["evidence_tier"]

    annotate_paper_burden(survivors)
    for verdict in all_verdicts:
        verdict["family_extended"] = family_for_report_extended(
            verdict["protein_name"])
        verdict["residual_placeholder"] = (
            "yes" if is_residual_placeholder(verdict["protein_name"]) else "no")

    gate.write_tsv(all_verdicts, VERDICT_COLUMNS, output)
    print(f"[refexp3] verdicts -> {output}", file=sys.stderr)

    write_clade_funnel(per_clade, _sibling(output, "_clade_funnel.tsv"))
    write_family_breadth(survivors, anchor_rows,
                         _sibling(output, "_family_breadth.tsv"))
    write_clade_family(survivors, _sibling(output, "_clade_family.tsv"))
    write_passing_fasta(survivors, all_entries, _sibling(output, "_passing.fasta"))

    leaks = [v for v in survivors if v["residual_placeholder"] == "yes"]
    if leaks:
        print(f"\n[refexp3] WARNING: {len(leaks)} survivors carry a name that "
              f"commits to no ligand and no function -- the approved gate's "
              f"placeholder patterns did not anticipate them. Flagged in "
              f"column residual_placeholder, NOT removed:", file=sys.stderr)
        for verdict in leaks:
            print(f"    {verdict['accession']:12} {verdict['clade']:<16} "
                  f"{verdict['protein_name']}", file=sys.stderr)
    return 0


def annotate_paper_burden(survivors: Sequence[dict]) -> None:
    """How many survivors share each survivor's smallest supporting paper.

    This is the sweep's single most important qualifier and the reason the
    headline count must not be read on its own. A paper contributing ONE
    entry characterized that receptor. A paper contributing 198 entries -- and
    one beetle opsin-duplication study here contributes exactly 198 -- mined a
    comparative dataset across dozens of species and named its sequences by
    phylogenetic position.

    Both kinds of entry clear the approved gate, and for opsins the second kind
    is arguably still legitimate family material, because opsin membership is
    diagnosable from sequence in a way that peptide-receptor ligand assignment
    is not. But they are NOT the same evidence, the difference is invisible in
    a survivor count, and the decision to merge belongs to the user. So it is
    measured and carried as a column rather than argued about here.

    The MINIMUM across an entry's supporting papers is used, so an entry that
    also appears in a focused study is credited to that study.
    """
    burden: Dict[str, int] = {}
    for verdict in survivors:
        for pmid in {p for p in verdict["primary_pmids"].split(";") if p}:
            burden[pmid] = burden.get(pmid, 0) + 1
    for verdict in survivors:
        pmids = [p for p in {p for p in verdict["primary_pmids"].split(";") if p}]
        verdict["min_paper_entry_count"] = (
            min(burden[p] for p in pmids) if pmids else "")


def _sibling(output: str, suffix: str) -> str:
    return os.path.splitext(output)[0] + suffix


CLADE_FUNNEL_COLUMNS = [
    "clade", "taxid", "universe_any_pe", "universe_pe1_pe2", "already_anchor",
    "new_candidates", "pass_class_a", "pass_existence", "pass_evidence_code",
    "pass_literature", "pass_characterized", "passing",
]


def write_clade_funnel(rows: Sequence[dict], path: str) -> None:
    gate.write_tsv(rows, CLADE_FUNNEL_COLUMNS, path)
    print(f"[refexp3] clade funnel -> {path}", file=sys.stderr)


FAMILY_BREADTH_COLUMNS = [
    "family", "anchor_total", "anchor_vertebrate", "anchor_nonvertebrate",
    "anchor_vertebrate_fraction", "refexp3_additions",
    "additions_focused_evidence", "additions_functionally_characterized",
    "additions_clades", "nonvertebrate_after", "vertebrate_fraction_after",
]

CLADE_FAMILY_COLUMNS = ["clade", "family", "passing", "focused_evidence",
                        "functionally_characterized", "distinct_species"]


def write_clade_family(survivors: Sequence[dict], path: str) -> None:
    """Passing entries broken out by clade AND family.

    Reported because a clade total can be one species sequenced many times.
    distinct_species is the honest denominator for 'phylogenetic breadth'.
    """
    keys = sorted({(v["clade"], v["family_extended"]) for v in survivors})
    rows = []
    for clade, family in keys:
        subset = [v for v in survivors
                  if v["clade"] == clade and v["family_extended"] == family]
        rows.append({
            "clade": clade,
            "family": family,
            "passing": len(subset),
            "focused_evidence": sum(
                1 for v in subset
                if isinstance(v.get("min_paper_entry_count"), int)
                and v["min_paper_entry_count"] <= 3),
            "functionally_characterized": sum(1 for v in subset if v["deorphanized"]),
            "distinct_species": len({v["organism"] for v in subset}),
        })
    gate.write_tsv(rows, CLADE_FAMILY_COLUMNS, path)
    print(f"[refexp3] clade x family -> {path}", file=sys.stderr)


def write_family_breadth(survivors: Sequence[dict], anchor_rows: Sequence[dict],
                         path: str) -> None:
    """Per family: does it gain real non-vertebrate breadth, or is it
    vertebrate-only no matter how wide the net goes?

    That question is the whole point of the sweep -- it decides whether a
    family's envelope can ever be trusted for a molluscan query -- so the
    vertebrate FRACTION is reported before and after, not a raw addition count
    that would flatter a family with 300 mammals in it.
    """
    class_a = [r for r in anchor_rows if r.get("class") == "A"]
    entries = gate.fetch_entries_by_accession(r["accession"] for r in class_a)

    families = sorted({r["family"] for r in class_a} |
                      {v["family_extended"] for v in survivors})
    rows = []
    for family in families:
        in_family = [r for r in class_a if r["family"] == family]
        vert = sum(1 for r in in_family
                   if "Vertebrata" in gate.lineage(entries[r["accession"]]))
        nonvert = len(in_family) - vert
        added = [v for v in survivors if v["family_extended"] == family]
        after_nonvert = nonvert + len(added)
        after_total = len(in_family) + len(added)
        rows.append({
            "family": family,
            "anchor_total": len(in_family),
            "anchor_vertebrate": vert,
            "anchor_nonvertebrate": nonvert,
            "anchor_vertebrate_fraction": (
                f"{vert / len(in_family):.3f}" if in_family else "n/a"),
            "refexp3_additions": len(added),
            "additions_focused_evidence": sum(
                1 for v in added
                if isinstance(v.get("min_paper_entry_count"), int)
                and v["min_paper_entry_count"] <= 3),
            "additions_functionally_characterized": sum(
                1 for v in added if v["deorphanized"]),
            "additions_clades": ",".join(sorted({v["clade"] for v in added})) or "-",
            "nonvertebrate_after": after_nonvert,
            "vertebrate_fraction_after": (
                f"{vert / after_total:.3f}" if after_total else "n/a"),
        })
    gate.write_tsv(rows, FAMILY_BREADTH_COLUMNS, path)
    print(f"[refexp3] family breadth -> {path}", file=sys.stderr)


def write_passing_fasta(survivors: Sequence[dict], entries: Dict[str, dict],
                        path: str) -> None:
    """Sequences that cleared the gate, verified byte-for-byte against source.

    Accumulating by accession has previously DOUBLED sequences in this project
    while passing every count check, so gate.verify_sequence_integrity re-checks
    every emitted sequence against its freshly-fetched source.
    """
    if not survivors:
        print("[refexp3] no survivors; no FASTA written", file=sys.stderr)
        return
    records = [{"accession": v["accession"],
                "sequence": entries[v["accession"]]["sequence"]["value"]}
               for v in survivors]
    gate.verify_sequence_integrity(entries, records)
    tmp = f"{path}.tmp"
    with open(tmp, "w") as handle:
        for verdict, record in zip(survivors, records):
            handle.write(f">{verdict['accession']}|{verdict['clade']}|"
                         f"{verdict['phylum_within_clade']}|"
                         f"{verdict['organism'].replace(' ', '_')}\n"
                         f"{record['sequence']}\n")
    os.replace(tmp, path)
    print(f"[refexp3] wrote {len(records)} verified sequences -> {path}",
          file=sys.stderr)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--anchor-tsv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    return sweep(args.anchor_tsv, args.output)


if __name__ == "__main__":
    sys.exit(main())
