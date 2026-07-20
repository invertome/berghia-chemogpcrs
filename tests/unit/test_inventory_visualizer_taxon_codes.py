"""Taxon-group prefixes in the tree visualizer must match real header codes.

2026-07 audit finding: ``visualize_gpcr_tree.PREFIX_TO_GROUP`` listed ``phau2``
under 'Other Lophotrochozoa' for *Phoronis australis*, but no sequence anywhere
carries a ``phau2_`` prefix -- the Phoronis headers are ``phau_`` and the
documented rename target (``recover_cds_from_assemblies.SPECIES_MAP``) is
``phaust``. So the entry matched nothing and every phoronid leaf fell through
to the 'phau' entry, colouring it as a gastropod.

2026-07 follow-up (bead 7uby): that fix swapped one dead key for another. The
``s/^>phau_/>phaust_/`` rename was never applied, so ``phaust`` also matched
ZERO sequences and all 428 real Phoronis headers still resolved through the
gastropod ``phau`` entry. The tests below asserted the rename as fact -- using
a fictional ``ref_phaust_g1`` input that no sequence carries -- and so passed
the whole time the figures were wrong. They now assert what the DATA says, and
``test_7uby_prefix_map_validation.py`` checks the map against the real
proteomes so a fixture can never again stand in for the real key.
"""
from __future__ import annotations

import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("Bio")

import visualize_gpcr_tree as vgt


def test_dead_phau2_code_is_gone():
    assert "phau2" not in vgt.PREFIX_TO_GROUP


def test_dead_phaust_code_is_gone():
    """The rename target is not a real code until the rename is applied."""
    assert "phaust" not in vgt.PREFIX_TO_GROUP


def test_phau_is_ambiguous_not_gastropod():
    """'phau' is carried by both Physella acuta AND Phoronis australis."""
    assert "phau" not in vgt.PREFIX_TO_GROUP
    assert "phau" in vgt.AMBIGUOUS_PREFIXES


def test_phoronis_leaf_is_not_coloured_as_a_gastropod():
    """A real Phoronis header (NMRA01* = its WGS project) must not read Gastropoda."""
    assert vgt.get_taxon_group("ref_phau_NMRA01000193.1_28409") == "Ambiguous"


def test_physella_leaf_is_not_confidently_resolved_either():
    """The collision is symmetric: neither species can be claimed from the code."""
    assert vgt.get_taxon_group("ref_phau_JAPYMB010000312.1_29984") == "Ambiguous"


def test_registered_codes_are_unique():
    """A code may belong to exactly one taxonomic group."""
    seen: dict[str, str] = {}
    for group, prefixes in vgt._groups.items():
        for p in prefixes:
            assert p not in seen, f"code {p!r} in both {seen.get(p)!r} and {group!r}"
            seen[p] = group


def test_prefix_twin_codes_resolve_independently():
    """'baar' is a strict prefix of 'baare'; longest-match must keep them apart."""
    assert vgt.get_taxon_group("ref_baare_g1") == vgt.PREFIX_TO_GROUP["baare"]
    assert vgt.get_taxon_group("ref_baar_g1") == vgt.PREFIX_TO_GROUP["baar"]
