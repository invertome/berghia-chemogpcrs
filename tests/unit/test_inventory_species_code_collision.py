"""Collision-safe reference species-code resolution.

2026-07 audit finding: the header-prefix code ``phau`` is claimed by two
proteomes in *different phyla* --

    references/nath_et_al/lse/gastropoda/109671_Physella_acuta.faa
        848 headers starting ">phau_"          (Mollusca, Gastropoda)
    references/nath_et_al/lse/other_lophotrochozoan_phyla/115415_Phoronis_australis.faa
        251 headers starting ">phau_"          (Phoronida)

``recover_cds_from_assemblies.SPECIES_MAP`` documents a one-time
``sed 's/^>phau_/>phaust_/'`` rename for the Phoronis file, but it was never
applied (no ">phaust_" header exists anywhere under references/). The curated
SPECIES_MAP entry for ``phau`` then silently "resolved" the collision toward
*Physella acuta*, so every phoronid sequence was reported as a snail.

The fix here is in the LOOKUP, not the data (the data is not ours to edit):
a code claimed by two different proteome files is AMBIGUOUS and must raise
rather than return a confident wrong answer.
"""
from __future__ import annotations

import pytest

import species_code_lookup as scl


def _write(refs, relpath, filename, code, n=1):
    d = refs / relpath
    d.mkdir(parents=True, exist_ok=True)
    body = "".join(f">{code}_g{i}\nMACDEF\n" for i in range(n))
    (d / filename).write_text(body)


# ------------------------------------------------------- ambiguity detect

def test_code_claimed_by_two_files_is_flagged_ambiguous(tmp_path):
    refs = tmp_path / "references"
    _write(refs, "lse/gastropoda", "109671_Physella_acuta.faa", "phau")
    _write(refs, "lse/other_phyla", "115415_Phoronis_australis.faa", "phau")

    code_map = scl.build_code_species_map(references_dir=refs, warn=False)
    assert "phau" in code_map.ambiguous
    assert code_map.ambiguous["phau"] == ["Phoronis australis", "Physella acuta"]


def test_species_for_raises_on_ambiguous_code(tmp_path):
    refs = tmp_path / "references"
    _write(refs, "lse/gastropoda", "109671_Physella_acuta.faa", "phau")
    _write(refs, "lse/other_phyla", "115415_Phoronis_australis.faa", "phau")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)

    with pytest.raises(scl.AmbiguousSpeciesCodeError) as exc:
        scl.species_for("ref_phau_g0", code_map)
    msg = str(exc.value)
    assert "phau" in msg
    assert "Physella acuta" in msg and "Phoronis australis" in msg


def test_species_for_non_strict_returns_none_not_wrong_species(tmp_path):
    """Opt-out degrades to 'unknown', never to a confident wrong answer."""
    refs = tmp_path / "references"
    _write(refs, "lse/gastropoda", "109671_Physella_acuta.faa", "phau")
    _write(refs, "lse/other_phyla", "115415_Phoronis_australis.faa", "phau")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)

    assert scl.species_for("ref_phau_g0", code_map, strict=False) == (None, "g0")


def test_curated_species_map_cannot_mask_a_two_file_collision(tmp_path):
    """SPECIES_MAP wins ordinary disagreements but must NOT silence a real
    two-file collision -- that is exactly how the phoronid mislabel happened."""
    refs = tmp_path / "references"
    # 'alvmar' is a real curated SPECIES_MAP code.
    _write(refs, "a", "1_Alviniconcha_marisindica.faa", "alvmar")
    _write(refs, "b", "2_Some_other.faa", "alvmar")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)
    assert "alvmar" in code_map.ambiguous
    with pytest.raises(scl.AmbiguousSpeciesCodeError):
        scl.species_for("ref_alvmar_x", code_map)


# --------------------------------------------------------- no regressions

def test_single_file_code_still_resolves(tmp_path):
    refs = tmp_path / "references"
    _write(refs, "a", "999_Genus_species.faa", "xyz")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)
    assert code_map.ambiguous == {}
    assert scl.species_for("ref_xyz_abc_1", code_map) == ("Genus species", "abc_1")


def test_species_map_override_of_single_file_still_wins(tmp_path):
    """Pre-existing contract: a lone filename disagreeing with the curated
    map is overridden by the curated map, silently and without ambiguity."""
    refs = tmp_path / "references"
    _write(refs, "a", "1_Wrong_species.faa", "alvmar")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)
    assert code_map.ambiguous == {}
    assert code_map["alvmar"][1] == "Alviniconcha marisindica"


def test_code_map_is_still_a_plain_mapping(tmp_path):
    """Downstream callers (visualize_gpcr_tree_rect) treat the result as a
    dict and call .get() on it -- keep that contract."""
    refs = tmp_path / "references"
    _write(refs, "a", "999_Genus_species.faa", "xyz")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)
    assert isinstance(code_map, dict)
    assert code_map.get("xyz") == ("999", "Genus species")
    assert code_map.get("no_such_code") is None
    assert dict(code_map)["xyz"] == ("999", "Genus species")


def test_curated_prefix_twin_codes_do_not_collide():
    """Latent twin flagged by the audit: 'baar' is a strict prefix of 'baare'.

    Both are curated SPECIES_MAP codes; the pairing below is read from that
    source, not recalled. Safe only because leaf parsing anchors on the '_'
    separator, so pin the behaviour.
    """
    from recover_cds_from_assemblies import SPECIES_MAP
    assert SPECIES_MAP["baar"][1] == "Batillaria_attramentaria"
    assert SPECIES_MAP["baare"][1] == "Babylonia_areolata"

    code_map = scl.build_code_species_map(warn=False)
    assert scl.species_for("ref_baar_g1", code_map)[0] == "Batillaria attramentaria"
    assert scl.species_for("ref_baare_g1", code_map)[0] == "Babylonia areolata"


def test_prefix_twin_codes_from_filenames_do_not_collide(tmp_path):
    """Same guarantee for codes that exist only as proteome filenames."""
    refs = tmp_path / "references"
    _write(refs, "a", "1_Genus_shortcode.faa", "zzq")
    _write(refs, "b", "2_Genus_longcode.faa", "zzqe")
    code_map = scl.build_code_species_map(references_dir=refs, warn=False)
    assert code_map.ambiguous == {}
    assert scl.species_for("ref_zzq_g1", code_map)[0] == "Genus shortcode"
    assert scl.species_for("ref_zzqe_g1", code_map)[0] == "Genus longcode"
