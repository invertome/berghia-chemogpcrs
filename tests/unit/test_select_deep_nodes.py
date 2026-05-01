"""Unit tests for select_deep_nodes.py focal-leaf detection (bead -mqt)."""
import pytest

from select_deep_nodes import leaf_belongs_to_focal


class TestLeafBelongsToFocal:
    def test_exact_match(self):
        assert leaf_belongs_to_focal("1287507", "1287507")
        assert leaf_belongs_to_focal("1287507_berghia_stephanieae",
                                    "1287507_berghia_stephanieae")

    def test_prefix_with_underscore(self):
        assert leaf_belongs_to_focal(
            "1287507_berghia_stephanieae_TRINITY_DN1", "1287507_berghia_stephanieae")
        assert leaf_belongs_to_focal("1287507_GENE1", "1287507")

    def test_ref_taxid_in_middle(self):
        # Reference header convention: ref_TAXID_N
        assert leaf_belongs_to_focal("ref_1287507_42", "1287507")
        assert leaf_belongs_to_focal("ref_outgroup_1287507_5", "1287507")

    def test_no_match(self):
        # Substring without underscore boundaries should NOT match
        # (otherwise '12' would match '128' etc.)
        assert not leaf_belongs_to_focal("ref_12_42", "1287507")
        assert not leaf_belongs_to_focal("9876_genus_species", "1287507")
        assert not leaf_belongs_to_focal("Orexin_human", "or")

    def test_empty_safety(self):
        assert not leaf_belongs_to_focal("", "1287507")
