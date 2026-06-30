"""Tests for species_code_lookup — code → species resolution.

Verifies the merge of the two authoritative sources (SPECIES_MAP + reference
proteome filenames) and the leaf-name parsing used to label tips.
"""
from __future__ import annotations

import species_code_lookup as scl


def test_parse_leaf():
    assert scl.parse_leaf("BersteEVm011558t16") == ("berghia", None, "BersteEVm011558t16")
    assert scl.parse_leaf("TRINITY_DN1_c0") == ("berghia", None, "TRINITY_DN1_c0")
    assert scl.parse_leaf("ref_seph_CAHIK_3") == ("ref", "seph", "CAHIK_3")
    assert scl.parse_leaf("ref_laant_scaffold1") == ("ref", "laant", "scaffold1")
    assert scl.parse_leaf("sp|P12345|XYZ") == ("other", None, "sp|P12345|XYZ")


def test_build_map_scans_filenames_and_merges_species_map(tmp_path):
    refs = tmp_path / "references"
    refs.mkdir()
    (refs / "9999_Genus_species.faa").write_text(">xyz_locus_1\nMACDEF\n")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)

    # picked up from the fake proteome filename + its header prefix
    assert code_map["xyz"] == ("9999", "Genus species")
    # SPECIES_MAP entries are always merged in (curated source)
    assert code_map["alvmar"][1] == "Alviniconcha marisindica"


def test_species_for(tmp_path):
    refs = tmp_path / "references"
    refs.mkdir()
    (refs / "9999_Genus_species.faa").write_text(">xyz_locus_1\nMACDEF\n")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)

    assert scl.species_for("ref_xyz_abc_1", code_map) == ("Genus species", "abc_1")
    assert scl.species_for("ref_unknowncode_5", code_map) == (None, "5")
    assert scl.species_for("BersteEVm1", code_map) == (None, "BersteEVm1")


def test_berghia_display_name():
    assert scl.berghia_display_name("BersteEVm009528t6") == "Berghia stephanieae EVm009528t6"
    assert scl.berghia_display_name("BersteEVm011558t16") == "Berghia stephanieae EVm011558t16"
    # non-Berste-prefixed Berghia ids (e.g. TRINITY_) still get the binomial
    assert scl.berghia_display_name("TRINITY_DN1").startswith("Berghia stephanieae ")


def test_species_map_wins_on_conflict(tmp_path):
    # A filename claiming a different species for an existing SPECIES_MAP code
    # must not override the curated entry.
    refs = tmp_path / "references"
    refs.mkdir()
    (refs / "1_Wrong_species.faa").write_text(">alvmar_x_1\nMACDEF\n")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)
    assert code_map["alvmar"][1] == "Alviniconcha marisindica"
