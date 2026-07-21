"""The refexp2 naming rule and the candidate universe it draws from.

Two defects are pinned here, both of which produced a plausible number rather
than an error:

  * the abbreviation-plus-index rule rejected genuine spelled-out names that
    merely carried an isoform number ("Opsin-3", "Xenopsin1"), contradicting
    its own documented discriminator;
  * the discovery query drew on the curated family string alone, which
    unreviewed records frequently lack, so roughly a third of eligible records
    were never examined.
"""

import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import refexp2_evidence_gate as gate  # noqa: E402


# --- the naming rule -------------------------------------------------------

# Spelled-out names carrying an isoform number. Every one of these has the
# letters-then-digits SHAPE the rule matches on, and every one names a real
# receptor, so the rule must not be what decides them.
SPELLED_OUT_WITH_INDEX = [
    "Opsin-3",
    "Opsin 1",
    "Peropsin 1",
    "Acropsin 1",
    "Xenopsin1",
    "Xenopsin2",
]

# True placeholders: an acronym and a serial number and nothing else. These are
# the case the rule was added for and must keep rejecting.
ABBREVIATION_PLACEHOLDERS = [
    "NPYR-10",
    "NPYR-5",
    "GCR002",
    "GCR484",
    "TKR2",
    "NPFR1",
]


@pytest.mark.parametrize("name", SPELLED_OUT_WITH_INDEX)
def test_spelled_out_names_with_an_index_are_not_rejected_as_placeholders(name):
    """A spelled-out word plus a number is a named receptor, not a placeholder."""
    assert gate.name_asserts_function(name), (
        f"{name!r} carries a spelled-out functional word and must not be "
        "rejected by the abbreviation-plus-index rule"
    )


@pytest.mark.parametrize("name", ABBREVIATION_PLACEHOLDERS)
def test_abbreviation_placeholders_are_still_rejected(name):
    """The rule's original purpose survives the exemption."""
    assert not gate.name_asserts_function(name)


@pytest.mark.parametrize("name", SPELLED_OUT_WITH_INDEX + ABBREVIATION_PLACEHOLDERS)
def test_every_case_here_actually_has_the_shape_the_rule_matches(name):
    """Guard the test itself.

    If these names stopped matching the shape regex, the parametrisations above
    would pass for the wrong reason -- they would be asserting about names the
    rule never touches.
    """
    assert gate._ABBREVIATION_INDEX_NAME.match(name), (
        f"{name!r} no longer has the abbreviation-plus-index shape, so it no "
        "longer exercises the exemption this module exists to pin"
    )


def test_exemption_requires_a_spelled_out_word_not_merely_lowercase():
    """The discriminator is a known spelled-out term, not letter case.

    A lowercase acronym must not launder itself through the exemption.
    """
    assert not gate._SPELLED_OUT_FUNCTION_WORD.search("npyr-10")
    assert gate._SPELLED_OUT_FUNCTION_WORD.search("Xenopsin1")


def test_rejection_of_true_placeholders_is_not_delegated_to_other_rules():
    """The abbreviation rule itself must do the rejecting for these names.

    Without this, the placeholder test above could pass because some other
    filter happened to catch them, leaving the abbreviation rule silently
    defanged.
    """
    for name in ABBREVIATION_PLACEHOLDERS:
        assert not gate._SPELLED_OUT_FUNCTION_WORD.search(name)


# --- the candidate universe ------------------------------------------------

def test_discovery_query_is_the_union_not_the_family_string_alone():
    """The family string alone misses records that carry only the signature."""
    query = gate.discovery_query()
    assert f'family:"{gate.CURATED_CLASS_A_FAMILY}"' in query
    assert f"xref:pfam-{gate.PFAM_7TM_1}" in query
    assert f"xref:interpro-{gate.INTERPRO_RHODOPSIN}" in query


def test_discovery_query_keeps_its_clade_and_existence_restrictions():
    """Widening the family term must not widen the clade or evidence filters."""
    query = gate.discovery_query()
    assert f"taxonomy_id:{gate.LOPHOTROCHOZOA_TAXID}" in query
    assert "existence:1" in query and "existence:2" in query


def test_universe_terms_are_or_ed_not_and_ed():
    """AND would NARROW the net to records carrying every term at once."""
    assert " OR " in gate.CLASS_A_UNIVERSE
    assert " AND " not in gate.CLASS_A_UNIVERSE


def test_widening_the_search_does_not_widen_the_class_verdict():
    """The class call stays per-entry and curation-led.

    The universe is a search convenience; a Pfam hit alone must never resolve a
    record to class A, because that is how a class-B/C receptor would enter the
    class-A envelope.
    """
    assert "xref:pfam" in gate.CLASS_A_UNIVERSE
    # PF00002 is class B. Curation naming family 2 outranks any signature.
    from curate_gpcr_references import gpcr_class_from_evidence
    assert gpcr_class_from_evidence(
        "SIMILARITY: Belongs to the G-protein coupled receptor 2 family.",
        "PF00001;",
    ) == "B"
