"""bead 7uby — an ambiguous species code must never be coloured as one phylum.

``species_code_lookup.species_for`` already refuses to guess for the ``phau``
code, which is claimed by BOTH *Physella acuta* (Gastropoda, 1078 headers) and
*Phoronis australis* (Phoronida, 428 headers). But ``get_taxon_group`` resolved
the same code through ``PREFIX_TO_GROUP``, which had no ambiguity register at
all, so those 428 phoronid sequences were coloured and legend-counted as
Gastropoda in published tree figures.

Both functions run on the same leaf three lines apart in
``visualize_gpcr_tree_rect.draw_lse``; they must agree that the code is
unresolvable.
"""
from __future__ import annotations

import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("Bio")

import species_code_lookup as scl
import visualize_gpcr_tree as vgt


def test_phau_is_registered_ambiguous():
    assert "phau" in vgt.AMBIGUOUS_PREFIXES


def test_ambiguous_code_is_not_silently_resolved_to_a_phylum():
    """The whole defect in one line: 'phau' must not come back 'Gastropoda'."""
    assert vgt.get_taxon_group("ref_phau_NMRA01000193.1_28409") == "Ambiguous"
    assert vgt.get_taxon_group("ref_phau_JAPYMB010000312.1_29984") == "Ambiguous"


def test_ambiguous_codes_are_absent_from_the_group_map():
    """A code that cannot be resolved must not sit in the map as a landmine."""
    for code in vgt.AMBIGUOUS_PREFIXES:
        assert code not in vgt.PREFIX_TO_GROUP


def test_ambiguous_group_is_renderable():
    """'Ambiguous' needs its own colour, or it silently falls back to grey."""
    assert "Ambiguous" in vgt.TAXON_COLORS


def test_both_seams_agree_the_code_is_unresolvable():
    """get_taxon_group and species_for must not disagree on the same leaf."""
    leaf = "ref_phau_NMRA01000193.1_28409"
    code_map = scl.CodeSpeciesMap(
        {"phau": ("109671", "Physella acuta")},
        ambiguous={"phau": ["Phoronis australis", "Physella acuta"]},
    )
    with pytest.raises(scl.AmbiguousSpeciesCodeError):
        scl.species_for(leaf, code_map)
    assert vgt.get_taxon_group(leaf) == "Ambiguous"


def test_code_map_ambiguity_is_honoured_when_supplied():
    """A collision discovered in data (not yet hardcoded) must also be refused."""
    code_map = scl.CodeSpeciesMap(
        {"alvmar": ("1198043", "Alviniconcha marisindica")},
        ambiguous={"alvmar": ["Alviniconcha marisindica", "Other genus species"]},
    )
    assert vgt.get_taxon_group("ref_alvmar_g1") == "Gastropoda"       # no map: as mapped
    assert vgt.get_taxon_group("ref_alvmar_g1", code_map) == "Ambiguous"


def test_unambiguous_codes_are_unaffected():
    assert vgt.get_taxon_group("ref_alvmar_g1") == "Gastropoda"
    assert vgt.get_taxon_group("BersteEVm009528t6") == "Berghia"
    assert vgt.get_taxon_group("") == "Unknown"


def test_prefix_twin_codes_still_resolve_independently():
    """'baar' is a strict prefix of 'baare'; longest-match must keep them apart."""
    assert vgt.get_taxon_group("ref_baare_g1") == vgt.PREFIX_TO_GROUP["baare"]
    assert vgt.get_taxon_group("ref_baar_g1") == vgt.PREFIX_TO_GROUP["baar"]


def test_berghia_display_name_refuses_to_fabricate_a_species():
    """It used to return 'Berghia stephanieae <anything>' for ANY input."""
    assert scl.berghia_display_name("BersteEVm009528t6") == "Berghia stephanieae EVm009528t6"
    assert scl.berghia_display_name("TRINITY_DN1").startswith("Berghia stephanieae ")
    for alien in ("ref_alvmar_x", "ref_phau_g1", "sp|P12345|XYZ", "", None):
        with pytest.raises(ValueError):
            scl.berghia_display_name(alien)
