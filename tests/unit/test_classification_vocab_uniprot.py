"""Pin the UniProt query vocabulary and the aminergic subfamily needles.

Audit finding #8: both `curate_gpcr_references.base_query()` and
`build_anchor_set.GPCR_KEYWORD` queried `keyword:"g-protein coupled
receptor"`. UniProt renamed KW-0297's display name to "G protein-coupled
receptor" (no hyphen after G), so the old string returns HTTP 200 with
x-total-results: 0. `curate_all()` then returned [] and OVERWROTE the
reference FASTA/TSV with headers only, exiting 0. Verified 2026-07-20:

    keyword:"g-protein coupled receptor" AND organism_id:9606  -> 0
    keyword:KW-0297                      AND organism_id:9606  -> 839

The fix pins the stable keyword ACCESSION (KW-0297), which cannot be
renamed, plus a fail-loud guard so a zero-row result can never silently
overwrite a shipped artifact again.

Audit finding #13: the tyramine subfamily needle `\\bTAR\\d?\\b` matches
the alternative names UniProt packs into the `protein_name` TSV field --
"Trace amine-associated receptor 5 (TaR-5) ...". Verified 2026-07-20
against the shipped aminergic_tyramine.aln: 41 members, 39 of them TAARs,
only 2 genuine tyramine receptors.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import build_anchor_set as bas
import curate_gpcr_references as cgr


# Stable UniProt keyword accession for "G protein-coupled receptor".
GPCR_KEYWORD_ACCESSION = "KW-0297"

# Real UniProt `protein_name` TSV values, copied verbatim 2026-07-20.
TAAR_PROTEIN_NAMES = [
    "Trace amine-associated receptor 1 (TaR-1) (Trace amine receptor 1)",
    "Trace amine-associated receptor 5 (TaR-5) (Trace amine receptor 5) (hTaar5)",
    "Trace amine-associated receptor 8c (TaR-8c) (Trace amine receptor 8c) (mTaar8c)",
    "Trace amine-associated receptor 7e (TaR-7e) (Trace amine receptor 7e) (mTaar7e)",
    "Putative trace amine-associated receptor 3 (TaR-3) (Trace amine receptor 3)",
]

# P22270 "Tyramine/octopamine receptor" is deliberately excluded: it is a
# genuine dual-specificity receptor and the octopamine pattern (which is
# ordered first) claims it. That is pre-existing, correct behaviour.
TRUE_TYRAMINE_PROTEIN_NAMES = [
    "Tyramine receptor Ser-2",
    "Tyramine receptor tyra-2",
]


# ---------------------------------------------------------------------------
# 1. Keyword accession, not display name
# ---------------------------------------------------------------------------

def test_curate_base_query_uses_keyword_accession():
    q = cgr.base_query()
    assert GPCR_KEYWORD_ACCESSION in q


def test_anchor_set_uses_keyword_accession():
    assert GPCR_KEYWORD_ACCESSION in bas.GPCR_KEYWORD


@pytest.mark.parametrize("module_query", [
    "curate_base_query",
    "anchor_gpcr_keyword",
])
def test_renamed_display_name_is_not_used(module_query):
    """The renamed display name silently returns zero rows — it must not
    appear in any query string."""
    q = (cgr.base_query() if module_query == "curate_base_query"
         else bas.GPCR_KEYWORD)
    assert "g-protein coupled receptor" not in q.lower()


def test_curate_base_query_still_restricts_to_curated_taxa():
    """Guard against the fix widening the query."""
    q = cgr.base_query()
    assert "reviewed:true" in q
    for taxid in cgr.CURATED_TAXA:
        assert f"organism_id:{taxid}" in q


# ---------------------------------------------------------------------------
# 2. Zero rows must fail loudly, never overwrite artifacts
# ---------------------------------------------------------------------------

def test_curate_all_raises_when_uniprot_returns_no_rows(monkeypatch):
    """A 200-with-zero-results must not produce an empty reference set."""
    monkeypatch.setattr(cgr, "query_uniprot_paginated",
                        lambda *a, **k: "Entry\tProtein names\tGene Names\n")
    with pytest.raises(RuntimeError, match="(?i)zero|no rows|empty"):
        cgr.curate_all()


def test_curate_all_succeeds_on_a_normal_result(monkeypatch):
    """The guard must not fire on a healthy response."""
    header = ("Entry\tReviewed\tProtein names\tGene Names\tOrganism\t"
              "Organism (ID)\tSequence\n")
    row = ("P08908\treviewed\t5-hydroxytryptamine receptor 1A\tHTR1A\t"
           "Homo sapiens\t9606\tMDVLSPGQGNNTT\n")
    monkeypatch.setattr(cgr, "query_uniprot_paginated",
                        lambda *a, **k: header + row)
    records = cgr.curate_all()
    assert len(records) == 1
    assert records[0]["accession"] == "P08908"


# ---------------------------------------------------------------------------
# 3. TAARs are not tyramine receptors
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", TAAR_PROTEIN_NAMES)
def test_taar_is_not_labelled_tyramine(name):
    """39 of the 41 shipped 'tyramine' reference sequences were TAARs."""
    family = cgr.classify_family(name, "")
    sub = cgr.classify_subfamily(family, name, "")
    assert sub != "tyramine", f"{name!r} was mislabelled as tyramine"


@pytest.mark.parametrize("name", TAAR_PROTEIN_NAMES)
def test_taar_gets_its_own_truthful_subfamily(name):
    """TAARs stay in the aminergic family but carry an honest label,
    so the mislabelled sequences are traceable rather than merely gone."""
    family = cgr.classify_family(name, "")
    assert family == "aminergic"
    assert cgr.classify_subfamily(family, name, "") == "trace-amine"


@pytest.mark.parametrize("name", TRUE_TYRAMINE_PROTEIN_NAMES)
def test_real_tyramine_receptors_still_classify_as_tyramine(name):
    """Guard against over-correction: the 2 genuine members must survive."""
    family = cgr.classify_family(name, "")
    assert cgr.classify_subfamily(family, name, "") == "tyramine"


def test_ambiguous_tar_needle_is_gone():
    """`\\bTAR\\d?\\b` matches 'TaR-1'..'TaR-9' in UniProt alternative
    names. It contributed zero unique true positives on the curated taxa
    (all 3 real tyramine receptors match on the word 'tyramine'), so it is
    pure false-positive surface."""
    tyramine_patterns = [p.pattern for p, label in
                         cgr._AMINERGIC_SUBFAMILY_PATTERNS if label == "tyramine"]
    assert tyramine_patterns, "tyramine subfamily pattern disappeared"
    for pattern in tyramine_patterns:
        assert "TAR" not in pattern.upper().replace("TYRAMINE", ""), (
            f"ambiguous TAR needle still present in {pattern!r}"
        )


def test_octopamine_is_unaffected():
    """Adjacent aminergic subfamilies must not shift."""
    name = "Octopamine receptor Oamb"
    assert cgr.classify_subfamily(cgr.classify_family(name, ""), name, "") == "octopamine"
