"""Taxon-group prefixes in the tree visualizer must match real header codes.

2026-07 audit finding: ``visualize_gpcr_tree.PREFIX_TO_GROUP`` listed ``phau2``
under 'Other Lophotrochozoa' for *Phoronis australis*, but no sequence anywhere
carries a ``phau2_`` prefix -- the Phoronis headers are ``phau_`` and the
documented rename target (``recover_cds_from_assemblies.SPECIES_MAP``) is
``phaust``. So the entry matched nothing and every phoronid leaf fell through
to the 'phau' entry, colouring it as a gastropod.
"""
from __future__ import annotations

import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("Bio")

import visualize_gpcr_tree as vgt


def test_phaust_is_the_registered_phoronis_code():
    assert vgt.PREFIX_TO_GROUP.get("phaust") == "Other Lophotrochozoa"


def test_dead_phau2_code_is_gone():
    assert "phau2" not in vgt.PREFIX_TO_GROUP


def test_phau_remains_the_gastropod_code():
    """After the data rename, 'phau' is unambiguously *Physella acuta*."""
    assert vgt.PREFIX_TO_GROUP.get("phau") == "Gastropoda"


def test_phoronis_leaf_groups_as_lophotrochozoan():
    assert vgt.get_taxon_group("ref_phaust_g1") == "Other Lophotrochozoa"


def test_physella_leaf_groups_as_gastropod():
    assert vgt.get_taxon_group("ref_phau_g1") == "Gastropoda"


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
