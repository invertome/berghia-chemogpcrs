"""bead 7uby — PREFIX_TO_GROUP must be provably in sync with the real headers.

``PREFIX_TO_GROUP`` is hand-maintained against auto-generated header prefixes
and had drifted badly: 230 map keys vs 238 real header codes, four keys matching
ZERO sequences (each a half-applied rename) and twelve real codes never mapped
at all, so 1,317 sequences fell through ``.get(first, 'Unknown')`` and rendered
as "Unknown". Because the fallback is a ``.get`` default, no run ever failed.

``validate_prefix_map`` closes that: it derives the truth from the reference
proteomes themselves (directory = taxonomic group, header prefix = code) and
refuses a map that has drifted. The last test runs it against the REAL
reference directory, which is what actually catches a half-applied rename the
moment it happens — a fixture cannot, because a fixture encodes the same
assumption the map does.
"""
from __future__ import annotations

import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("Bio")

import visualize_gpcr_tree as vgt

REFS = vgt.PROJECT_ROOT / "references" / "nath_et_al"


def _proteome(root, group_dir, taxid, species, code, n=3):
    d = root / group_dir
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{taxid}_{species}.faa"
    f.write_text("".join(f">{code}_locus{i}\nMACDEF\n" for i in range(n)))
    return f


def test_scan_derives_codes_and_groups_from_the_data(tmp_path):
    _proteome(tmp_path, "gastropoda", "1", "Genus_one", "gone")
    _proteome(tmp_path, "other_lophotrochozoan_phyla", "2", "Genus_two", "gtwo")
    obs = vgt.scan_reference_codes(tmp_path)
    assert obs["gone"].group == "Gastropoda"
    assert obs["gtwo"].group == "Other Lophotrochozoa"
    assert obs["gone"].count == 3


def test_validation_passes_on_a_map_that_matches(tmp_path):
    _proteome(tmp_path, "gastropoda", "1", "Genus_one", "gone")
    vgt.validate_prefix_map(tmp_path, prefix_to_group={"gone": "Gastropoda"},
                            ambiguous=set())


def test_validation_rejects_a_key_matching_zero_sequences(tmp_path):
    """The 'phaust'/'lian2'/'meme'/'scin' failure: a half-applied rename."""
    _proteome(tmp_path, "gastropoda", "1", "Genus_one", "gone")
    with pytest.raises(vgt.PrefixMapDriftError) as e:
        vgt.validate_prefix_map(tmp_path, ambiguous=set(),
                                prefix_to_group={"gone": "Gastropoda",
                                                 "goneRENAMED": "Gastropoda"})
    assert "goneRENAMED" in str(e.value)


def test_validation_rejects_an_unmapped_code(tmp_path):
    """The 'liant'/'memmem'/'taso'... failure: 1,317 sequences as 'Unknown'."""
    _proteome(tmp_path, "gastropoda", "1", "Genus_one", "gone")
    _proteome(tmp_path, "bivalvia", "2", "Genus_two", "gtwo", n=772)
    with pytest.raises(vgt.PrefixMapDriftError) as e:
        vgt.validate_prefix_map(tmp_path, prefix_to_group={"gone": "Gastropoda"},
                                ambiguous=set())
    assert "gtwo" in str(e.value)


def test_validation_rejects_a_wrong_group(tmp_path):
    _proteome(tmp_path, "bivalvia", "1", "Genus_one", "gone")
    with pytest.raises(vgt.PrefixMapDriftError) as e:
        vgt.validate_prefix_map(tmp_path, prefix_to_group={"gone": "Gastropoda"},
                                ambiguous=set())
    assert "Bivalvia" in str(e.value)


def test_validation_rejects_an_unregistered_collision(tmp_path):
    """Two proteomes in different phyla sharing a code must be declared."""
    _proteome(tmp_path, "gastropoda", "1", "Genus_one", "dup")
    _proteome(tmp_path, "other_lophotrochozoan_phyla", "2", "Genus_two", "dup")
    with pytest.raises(vgt.PrefixMapDriftError) as e:
        vgt.validate_prefix_map(tmp_path, prefix_to_group={"dup": "Gastropoda"},
                                ambiguous=set())
    assert "dup" in str(e.value)


def test_validation_accepts_a_declared_collision(tmp_path):
    _proteome(tmp_path, "gastropoda", "1", "Genus_one", "dup")
    _proteome(tmp_path, "other_lophotrochozoan_phyla", "2", "Genus_two", "dup")
    vgt.validate_prefix_map(tmp_path, prefix_to_group={}, ambiguous={"dup"})


def test_validation_rejects_a_declared_collision_left_in_the_group_map(tmp_path):
    _proteome(tmp_path, "gastropoda", "1", "Genus_one", "dup")
    _proteome(tmp_path, "other_lophotrochozoan_phyla", "2", "Genus_two", "dup")
    with pytest.raises(vgt.PrefixMapDriftError):
        vgt.validate_prefix_map(tmp_path, prefix_to_group={"dup": "Gastropoda"},
                                ambiguous={"dup"})


@pytest.mark.skipif(not REFS.is_dir(), reason="reference proteomes not present")
def test_shipped_map_matches_the_real_reference_proteomes():
    """The assertion that would have caught BOTH defects the day they landed.

    No fixture can do this job: a fixture encodes the same assumption the map
    does, so it agrees with the map even when the map is wrong about the data.
    """
    vgt.validate_prefix_map()


@pytest.mark.skipif(not REFS.is_dir(), reason="reference proteomes not present")
def test_no_real_sequence_falls_through_to_unknown():
    """Every real header code must render as a real group, not 'Unknown'."""
    obs = vgt.scan_reference_codes()
    unknown = {c: o.count for c, o in obs.items()
               if vgt.get_taxon_group(f"ref_{c}_locus1") == "Unknown"}
    assert not unknown, f"{sum(unknown.values())} sequences render as Unknown: {unknown}"


@pytest.mark.skipif(not REFS.is_dir(), reason="reference proteomes not present")
def test_phoronis_headers_are_not_coloured_as_gastropods():
    """428 real Phoronis australis headers, measured, must not read Gastropoda."""
    obs = vgt.scan_reference_codes()
    assert obs["phau"].ambiguous_over == {"Gastropoda", "Other Lophotrochozoa"}
    assert vgt.get_taxon_group("ref_phau_NMRA01000193.1_28409") != "Gastropoda"
