"""Unit tests for OrthoDB clade assignment used by the expansion ceiling.

The clade breakdown is the number the whole exercise is meant to produce, and
it keys entirely on parsing an organism id out of an OrthoDB gene id and then
reading that organism's level path. Both are cross-namespace steps, so both are
tested against the real OrthoDB formats rather than assumed.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "orthodb_expansion_ceiling",
    Path(__file__).resolve().parents[2] / "scripts" / "orthodb_expansion_ceiling.py",
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


# ------------------------------------------------------- gene id -> organism --

@pytest.mark.parametrize(
    "gene_id,expected",
    [
        ("9606_0:001a2b", "9606_0"),      # human, real odb12v2 shape
        ("6500_1:0000ff", "6500_1"),      # Aplysia californica, second assembly
        ("225164_0:00012", "225164_0"),   # Lottia gigantea
    ],
)
def test_parse_org_id_extracts_organism(gene_id, expected):
    assert mod.parse_org_id(gene_id) == expected


@pytest.mark.parametrize(
    "bad", ["9606:001", "abc_0:001", "9606_0", "", "9606_0_1:x"]
)
def test_parse_org_id_rejects_unexpected_formats(bad):
    """Silently mis-parsing an id would mis-assign a clade with no error."""
    with pytest.raises(ValueError):
        mod.parse_org_id(bad)


# --------------------------------------------------- level path -> clade --

def _write_l2s(tmp_path, rows):
    p = tmp_path / "odb_level2species.tsv"
    p.write_text("".join("\t".join(r) + "\n" for r in rows))
    return p


def test_clade_assignment_follows_the_level_path(tmp_path):
    p = _write_l2s(
        tmp_path,
        [
            ("2759", "9606_0", "1", "{2759,33208,7742,32523,40674,9606}"),
            ("2759", "6500_0", "1", "{2759,33208,1206795,6447,6500}"),
            ("2759", "7227_0", "1", "{2759,33208,6656,6960,50557,7227}"),
            ("2759", "6239_0", "1", "{2759,33208,6231,119089,6239}"),
            ("2759", "6085_0", "1", "{2759,33208,6073,6101,6085}"),
            ("2759", "4932_0", "1", "{2759,4890,4932}"),
        ],
    )
    c = mod.load_org_clades(p)
    assert c["9606_0"] == "Vertebrata"
    assert c["6500_0"] == "Mollusca"
    assert c["7227_0"] == "Arthropoda"
    assert c["6239_0"] == "Nematoda"
    assert c["6085_0"] == "Cnidaria"
    assert c["4932_0"] == "non_Metazoa"


def test_mollusca_wins_over_the_lophotrochozoa_it_sits_inside(tmp_path):
    """Mollusca is nested in Lophotrochozoa; order must not double-count it."""
    p = _write_l2s(
        tmp_path, [("2759", "6500_0", "1", "{2759,33208,1206795,6447,6500}")]
    )
    assert mod.load_org_clades(p)["6500_0"] == "Mollusca"


def test_non_molluscan_lophotrochozoan_is_labelled_separately(tmp_path):
    p = _write_l2s(
        tmp_path, [("2759", "6358_0", "1", "{2759,33208,1206795,6358}")]
    )
    assert mod.load_org_clades(p)["6358_0"] == "other_Lophotrochozoa"


def test_metazoan_outside_every_modelled_clade_falls_to_other_metazoa(tmp_path):
    """Tunicates/echinoderms have no OrthoDB level, so they must land here
    rather than being silently dropped or mislabelled as non-Metazoa."""
    p = _write_l2s(
        tmp_path, [("2759", "7719_0", "1", "{2759,33208,7719}")]  # Ciona
    )
    assert mod.load_org_clades(p)["7719_0"] == "other_Metazoa"


def test_empty_level2species_is_fatal(tmp_path):
    p = tmp_path / "odb_level2species.tsv"
    p.write_text("")
    with pytest.raises(SystemExit):
        mod.load_org_clades(p)


def test_clade_order_taxids_are_the_verified_orthodb_levels():
    """Guards against someone adding a taxid that is not an OrthoDB level.

    Bilateria (33213), Protostomia (33317), Deuterostomia (33511), Chordata
    (7711) and Ecdysozoa (1206794) are NOT levels in odb12v2 and must never
    appear here, because no orthogroup is ever built on them.
    """
    taxids = {t for _, t in mod.CLADE_ORDER}
    forbidden = {"33213", "33317", "33511", "7711", "1206794"}
    assert not (taxids & forbidden)
    assert taxids == {"7742", "6447", "1206795", "6656", "6231", "6073"}
