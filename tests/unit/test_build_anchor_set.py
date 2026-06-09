"""Unit tests for build_anchor_set.py (C2 — three-tier anchor sourcing).

Covers the pure logic (family->class, header format, family classification,
tier-3 from the curated TSV, accession dedup, output writing) and the
network-backed tier sourcing via an injected fake UniProt search function.

No network / no compute: the UniProt calls are dependency-injected and the
curated-TSV input is a tiny in-test fixture.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import build_anchor_set as bas  # noqa: E402


# ---------------------------------------------------------------------------
# family -> class rule (class-B*->B, class-C*->C, class-F*->F, else A)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("family,expected", [
    ("class-B-secretin", "B"),
    ("class-C", "C"),
    ("class-F-frizzled", "F"),
    ("peptide", "A"),
    ("aminergic", "A"),
    ("opsin", "A"),
    ("lipid", "A"),
    ("nucleotide", "A"),
    ("glycoprotein-hormone", "A"),
    ("", "A"),
    ("something-unknown", "A"),
])
def test_family_to_class(family, expected):
    assert bas.family_to_class(family) == expected


# ---------------------------------------------------------------------------
# anchor header: ANCHOR_<class>_<tier>_<accession>
# ---------------------------------------------------------------------------

def test_anchor_header_format():
    assert bas.anchor_header("A", "1", "P12345") == "ANCHOR_A_1_P12345"
    assert bas.anchor_header("F", "3", "Q9XYZ1") == "ANCHOR_F_3_Q9XYZ1"


# ---------------------------------------------------------------------------
# classify_family: reuse the curate_gpcr_references protein-name patterns
# ---------------------------------------------------------------------------

def test_classify_family_known_names():
    # class-B / C / F discriminators
    assert bas.family_to_class(bas.classify_family("Frizzled-7")) == "F"
    assert bas.family_to_class(
        bas.classify_family("Metabotropic glutamate receptor 5")) == "C"
    assert bas.family_to_class(
        bas.classify_family("Secretin receptor")) == "B"
    # Class A families
    assert bas.family_to_class(bas.classify_family("Rhodopsin")) == "A"
    assert bas.family_to_class(
        bas.classify_family("5-hydroxytryptamine receptor 2")) == "A"


def test_classify_family_unknown_defaults_to_class_a():
    # A name matching no specific family still resolves to Class A.
    assert bas.family_to_class(bas.classify_family("Hypothetical receptor X")) == "A"


# ---------------------------------------------------------------------------
# tier-3 from the curated non-chemo TSV: Drosophila + C. elegans only
# ---------------------------------------------------------------------------

@pytest.fixture
def curated_tsv(tmp_path):
    p = tmp_path / "all_references.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "family", "subfamily", "species", "gene",
                    "source", "length"])
        w.writerow(["P1", "peptide", "", "Drosophila melanogaster", "g1", "sp", "400"])
        w.writerow(["P2", "class-B-secretin", "", "Drosophila melanogaster", "g2", "sp", "450"])
        w.writerow(["P3", "class-C", "", "Caenorhabditis elegans", "g3", "sp", "900"])
        w.writerow(["P4", "class-F-frizzled", "", "Caenorhabditis elegans", "g4", "sp", "500"])
        # Non-target taxa must be excluded:
        w.writerow(["P5", "aminergic", "", "Homo sapiens", "g5", "sp", "420"])
        w.writerow(["P6", "opsin", "", "Mus musculus", "g6", "sp", "350"])
    # a matching FASTA so sequences are attached
    fa = tmp_path / "all_references.fasta"
    fa.write_text(
        ">P1|x\nAAAA\n>P2|x\nCCCC\n>P3|x\nDDDD\n>P4|x\nEEEE\n>P5|x\nFFFF\n>P6|x\nGGGG\n")
    return p, fa


def test_parse_tier3_keeps_only_fly_and_worm(curated_tsv):
    tsv, fa = curated_tsv
    recs = bas.parse_tier3_from_curated_tsv(str(tsv), str(fa))
    accs = {r["accession"] for r in recs}
    assert accs == {"P1", "P2", "P3", "P4"}  # human + mouse excluded


def test_parse_tier3_assigns_class_tier_and_evidence(curated_tsv):
    tsv, fa = curated_tsv
    recs = {r["accession"]: r for r in bas.parse_tier3_from_curated_tsv(str(tsv), str(fa))}
    assert recs["P1"]["klass"] == "A"
    assert recs["P2"]["klass"] == "B"
    assert recs["P3"]["klass"] == "C"
    assert recs["P4"]["klass"] == "F"
    for r in recs.values():
        assert r["tier"] == "3"
        assert r["evidence"] == "reviewed"
    assert recs["P1"]["sequence"] == "AAAA"


# ---------------------------------------------------------------------------
# dedup by accession: collision prefers reviewed over experimental
# ---------------------------------------------------------------------------

def test_dedup_by_accession_prefers_reviewed():
    recs = [
        {"accession": "P1", "tier": "1", "evidence": "experimental:111", "klass": "A"},
        {"accession": "P1", "tier": "1", "evidence": "reviewed", "klass": "A"},
        {"accession": "P2", "tier": "2", "evidence": "experimental:222", "klass": "A"},
    ]
    out = {r["accession"]: r for r in bas.dedup_by_accession(recs)}
    assert len(out) == 2
    assert out["P1"]["evidence"] == "reviewed"
    assert out["P2"]["evidence"] == "experimental:222"


# ---------------------------------------------------------------------------
# write_anchor_set: FASTA headers + TSV columns
# ---------------------------------------------------------------------------

def test_write_anchor_set_outputs(tmp_path):
    recs = [
        {"accession": "P1", "tier": "1", "taxid": 6500, "species": "Aplysia californica",
         "family": "aminergic", "klass": "A", "evidence": "reviewed", "sequence": "MKT"},
        {"accession": "Q9", "tier": "3", "taxid": 7227, "species": "Drosophila melanogaster",
         "family": "class-F-frizzled", "klass": "F", "evidence": "reviewed", "sequence": "MWW"},
    ]
    bas.write_anchor_set(recs, str(tmp_path))
    fasta = (tmp_path / "anchor_set.fasta").read_text()
    assert ">ANCHOR_A_1_P1\nMKT" in fasta
    assert ">ANCHOR_F_3_Q9\nMWW" in fasta
    with open(tmp_path / "anchor_set.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {row["accession"] for row in rows} == {"P1", "Q9"}
    cols = set(rows[0].keys())
    assert {"accession", "tier", "taxid", "species", "family", "class",
            "evidence"} <= cols


# ---------------------------------------------------------------------------
# Sourcing tiers via an injected fake UniProt search
# ---------------------------------------------------------------------------

def _fake_search(responses):
    """Return a search_fn that yields canned rows keyed by a query substring."""
    def search_fn(query, fields):
        for key, rows in responses.items():
            if key in query:
                return rows
        return []
    return search_fn


def test_source_tier1a_reviewed_tags_and_classifies():
    rows = [
        {"accession": "P31356", "protein_name": "Rhodopsin",
         "organism_id": "6500", "organism_name": "Aplysia californica",
         "sequence": "MKTII"},
    ]
    search = _fake_search({"reviewed:true": rows})
    recs = bas.source_tier1a(search)
    assert len(recs) == 1
    r = recs[0]
    assert r["tier"] == "1"
    assert r["evidence"] == "reviewed"
    assert r["klass"] == "A"
    assert r["accession"] == "P31356"
    assert r["taxid"] == 6500


def test_source_tier1c_restricted_to_pmid_allowlist():
    # Two entries: one cites an allow-listed PMID, one does not.
    rows_good = [
        {"accession": "O61232", "protein_name": "Serotonin receptor 5-HT2",
         "organism_id": "6523", "organism_name": "Lymnaea stagnalis",
         "sequence": "MABC", "lit_pubmed_id": "8891606"},
    ]
    # The query must include the OR-joined PMID filter; the fake matches on that.
    search = _fake_search({"8891606": rows_good})
    recs = bas.source_tier1c({"8891606", "10677541"}, search)
    assert len(recs) == 1
    r = recs[0]
    assert r["tier"] == "1"
    assert r["evidence"] == "experimental:8891606"
    assert r["klass"] == "A"


def test_functional_pmids_embedded_default():
    """With no manifest, the verified 18-PMID functional allow-list is embedded
    in the module (so the script is self-contained on Unity, where references/
    is not git-tracked)."""
    pmids = bas.load_functional_pmids()
    assert len(pmids) == 18
    # a couple of anchor pharmacology papers must be present
    assert "8891606" in pmids    # Lymnaea 5-HT2 functional characterisation
    assert "26190115" not in pmids  # that's the tier-2 Platynereis paper, not tier-1c


def test_functional_pmids_manifest_overrides(tmp_path):
    """An explicit manifest TSV overrides the embedded default."""
    p = tmp_path / "pmids.tsv"
    p.write_text("pmid\tyear\ttitle\n111\t2020\tfoo\n222\t2021\tbar\n")
    pmids = bas.load_functional_pmids(str(p))
    assert pmids == {"111", "222"}


def test_source_tier2_platynereis_evidence():
    rows = [
        {"accession": "PLAT1", "protein_name": "Neuropeptide receptor",
         "organism_id": "6359", "organism_name": "Platynereis dumerilii",
         "sequence": "MPLT"},
    ]
    search = _fake_search({"26190115": rows})
    recs = bas.source_tier2(search)
    assert len(recs) == 1
    assert recs[0]["tier"] == "2"
    assert recs[0]["evidence"] == "experimental:26190115"
    assert recs[0]["klass"] == "A"
