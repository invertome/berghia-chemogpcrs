"""Tests for the anchor GPCR-class policy (``gpcr_class_from_evidence``).

The anchor set's ``class`` column decides which anchors enter the class-A
envelope (``extract_prod_anchors``), which in turn is the reference frame for
both the A1 tree and the PLM novelty prototypes. A class-B or class-C receptor
sitting under ``class=A`` silently corrupts a family prototype, so class must be
derived from the underlying UniProt record by ONE uniform rule rather than from
the protein name or a hand edit.

The rule, in precedence order:

1. UniProt's curated ``SIMILARITY: Belongs to the G[-]protein coupled receptor
   N family`` statement, where N in {1, 2, 3} maps to class {A, B, C}.
2. Failing that, the specific Pfam family: PF00001 (7tm_1) -> A,
   PF00002 (7tm_2) -> B, PF00003 (7tm_3) -> C.
3. Failing both, ``UNKNOWN`` -- never a silent default to A.

Pfam CLAN membership is deliberately NOT used: PF00001/PF00002/PF00003/PF02101
all sit in clan CL0192 (GPCR_A), so the clan cannot discriminate class. Only
the specific family can. Verified against InterPro 2026-07-20.

The real-data cases below are pinned to records fetched from the UniProt REST
API, not recalled.
"""
import pytest

from curate_gpcr_references import gpcr_class_from_evidence


# --- rule 1: the curated family statement ------------------------------------

@pytest.mark.parametrize("n,expected", [("1", "A"), ("2", "B"), ("3", "C")])
def test_curated_family_number_maps_to_class(n, expected):
    sim = f"SIMILARITY: Belongs to the G-protein coupled receptor {n} family."
    assert gpcr_class_from_evidence(sim, "") == expected


@pytest.mark.parametrize("prefix", ["G-protein coupled", "G protein-coupled"])
def test_both_uniprot_hyphenation_variants_are_recognised(prefix):
    """UniProt is inconsistent about the hyphen; both spellings occur in the
    anchor set and both must parse."""
    sim = f"SIMILARITY: Belongs to the {prefix} receptor 1 family."
    assert gpcr_class_from_evidence(sim, "") == "A"


def test_subfamily_suffix_does_not_break_the_match():
    sim = ("SIMILARITY: Belongs to the G-protein coupled receptor 3 family. "
           "GABA-B receptor subfamily. {ECO:0000256|ARBA:ARBA00008991}.")
    assert gpcr_class_from_evidence(sim, "") == "C"


def test_curated_family_outranks_pfam():
    """Curation is expert-assigned; Pfam is a profile hit. When they disagree
    the curated statement wins, and the disagreement is the caller's to audit."""
    sim = "SIMILARITY: Belongs to the G-protein coupled receptor 2 family."
    assert gpcr_class_from_evidence(sim, "PF00001;") == "B"


# --- rule 2: the Pfam fallback -----------------------------------------------

@pytest.mark.parametrize("pfam,expected", [
    ("PF00001;", "A"),
    ("PF00002;", "B"),
    ("PF00003;PF01094;", "C"),
])
def test_pfam_family_used_when_no_curated_statement(pfam, expected):
    assert gpcr_class_from_evidence("", pfam) == expected


def test_pfam_match_is_exact_not_substring():
    """PF00001 must not be matched by PF000012 or similar."""
    assert gpcr_class_from_evidence("", "PF000012;") == "UNKNOWN"


# --- rule 3: no silent default ------------------------------------------------

def test_no_evidence_is_unknown_never_class_a():
    assert gpcr_class_from_evidence("", "") == "UNKNOWN"


def test_unrelated_similarity_and_pfam_is_unknown():
    sim = "SIMILARITY: Belongs to the G protein-coupled receptor OA family."
    assert gpcr_class_from_evidence(sim, "PF02101;") == "UNKNOWN"


def test_clan_sibling_pf02101_does_not_confer_class_a():
    """PF02101 shares clan CL0192 with PF00001 but is a different family;
    clan membership must not leak into the class call."""
    assert gpcr_class_from_evidence("", "PF02101;") == "UNKNOWN"


# --- real records, pinned to the live UniProt API (verified 2026-07-20) -------

# (accession, curated similarity, pfam, expected class, why it matters)
_REAL = [
    # ACKR1/DARC: curated GPCR family 1, but NO Pfam at all. A domain-only
    # audit calls this a class-A failure; it is a genuine, divergent class-A
    # atypical chemokine receptor and must stay in the class-A envelope.
    ("Q5Y7A3",
     "SIMILARITY: Belongs to the G-protein coupled receptor 1 family. "
     "Atypical chemokine receptor subfamily.", "", "A"),
    # C. elegans putative GPCR: curated family 1, no Pfam. Same shape.
    ("Q09965",
     "SIMILARITY: Belongs to the G protein-coupled receptor 1 family. "
     "B0244 subfamily.", "", "A"),
    # Drosophila GABA-B R2: curated family 3 -> class C, NOT class B.
    ("Q9BML5",
     "SIMILARITY: Belongs to the G-protein coupled receptor 3 family. "
     "GABA-B receptor subfamily.", "PF00003;PF01094;", "C"),
    # Manduca diuretic hormone receptor: curated family 2 -> class B.
    ("P35464",
     "SIMILARITY: Belongs to the G protein-coupled receptor 2 family.",
     "PF00002;PF02793;", "B"),
    # CMKLR1: no "receptor N family" statement, but carries PF00001 -> class A.
    ("P46091",
     "SIMILARITY: Belongs to the chemokine-like receptor (CMKLR) family.",
     "PF00001;", "A"),
    # GPR143/OA1: its own curated family, Pfam PF02101 only -> not class A.
    ("P51810",
     "SIMILARITY: Belongs to the G protein-coupled receptor OA family.",
     "PF02101;", "UNKNOWN"),
]


@pytest.mark.parametrize("acc,sim,pfam,expected", _REAL,
                         ids=[r[0] for r in _REAL])
def test_real_uniprot_records(acc, sim, pfam, expected):
    assert gpcr_class_from_evidence(sim, pfam) == expected, acc
