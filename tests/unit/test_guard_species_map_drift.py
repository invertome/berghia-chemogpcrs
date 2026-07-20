"""Guard: ``SPECIES_MAP`` must stay honest about the reference proteomes.

``recover_cds_from_assemblies.SPECIES_MAP`` is a hand-curated
``code -> (taxid, species, assembly)`` table. Its sibling ``PREFIX_TO_GROUP``
gained a real-data validator (``visualize_gpcr_tree.validate_prefix_map``);
``SPECIES_MAP`` had none, and it had drifted:

  ``"phaust": ("115415", "Phoronis_australis", "GCA_055505105.1")``

No sequence anywhere carries a ``phaust_`` prefix -- the documented
``sed 's/^>phau_/>phaust_/'`` rename was never applied. The entry was therefore
a PHANTOM key: ``species_code_lookup.build_code_species_map()`` injected it into
the code map, so ``species_for('ref_phaust_x')`` confidently returned
*Phoronis australis* for an identifier that cannot exist, while the 428 real
phoronid sequences still carried ``phau_`` and collided with *Physella acuta*.

These tests derive truth from the proteomes under ``references/`` rather than
from a fixture. A fixture would encode the ASSUMED key and pass while the real
map was wrong -- the failure mode this repo has hit repeatedly.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import recover_cds_from_assemblies as rca  # noqa: E402

REFERENCES = PROJECT_ROOT / "references"
needs_references = pytest.mark.skipif(
    not REFERENCES.is_dir(),
    reason="references/ proteomes not present in this checkout",
)


# ---------------------------------------------------------------------------
# 1. The phantom key itself
# ---------------------------------------------------------------------------

def test_phaust_phantom_key_is_gone():
    """'phaust' is a rename TARGET, not a code any sequence carries.

    Keeping it in SPECIES_MAP made ``species_for('ref_phaust_...')`` resolve to
    a species for an identifier that does not exist.
    """
    assert "phaust" not in rca.SPECIES_MAP


def test_phau_still_maps_to_the_genome_it_really_names():
    """Removing the phantom must not disturb the real, live entry.

    The rename was never applied, so 'phau' is still the only phoronid/physellid
    code in the data. Its SPECIES_MAP entry names the *Physella* assembly.
    """
    assert rca.SPECIES_MAP["phau"][1] == "Physella_acuta"


def test_phau_is_declared_ambiguous():
    """The collision must be REGISTERED, not resolved by curation fiat."""
    assert "phau" in rca.AMBIGUOUS_SPECIES_CODES
    assert set(rca.AMBIGUOUS_SPECIES_CODES["phau"]) == {
        "Physella_acuta", "Phoronis_australis"}


# ---------------------------------------------------------------------------
# 2. The validator, against the REAL proteomes
# ---------------------------------------------------------------------------

@needs_references
def test_species_map_validates_against_real_proteomes():
    """The whole point: the shipped map must match the data on disk."""
    rca.validate_species_map()


@needs_references
def test_scan_finds_the_real_phau_collision():
    """Truth is derived, not declared: both species really do claim 'phau'."""
    observed = rca.scan_reference_species()
    assert "phau" in observed
    assert observed["phau"].ambiguous_over == {"Physella_acuta",
                                               "Phoronis_australis"}
    assert observed["phau"].species is None, "an ambiguous code has no species"


@needs_references
def test_scan_confirms_phaust_matches_zero_sequences():
    """The measurement that proves the phantom: zero sequences, not 'few'."""
    assert "phaust" not in rca.scan_reference_species()


@needs_references
def test_every_mapped_code_matches_at_least_one_real_sequence():
    """No dead keys. This is the check that would have caught 'phaust'."""
    observed = rca.scan_reference_species()
    dead = sorted(set(rca.SPECIES_MAP) - set(observed))
    assert dead == [], f"SPECIES_MAP code(s) match zero sequences: {dead}"


# ---------------------------------------------------------------------------
# 3. The validator must actually REFUSE drift (not just pass today)
# ---------------------------------------------------------------------------

def _proteome(tmp_path: Path, rel: str, taxid: str, species: str, code: str,
              n: int = 3) -> None:
    p = tmp_path / rel / f"{taxid}_{species}.faa"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("".join(f">{code}_g{i}\nMAAA\n" for i in range(n)))


def test_validator_refuses_a_dead_key(tmp_path: Path):
    """A reintroduced phantom must fail loudly."""
    _proteome(tmp_path, "lse/gastropoda", "109671", "Physella_acuta", "phau")
    with pytest.raises(rca.SpeciesMapDriftError, match="ZERO sequences"):
        rca.validate_species_map(
            references_dir=tmp_path,
            species_map={"phau": ("109671", "Physella_acuta", "GCF_x"),
                         "ghost": ("1", "Not_real", "GCA_x")},
            ambiguous={})


def test_validator_refuses_a_wrong_species(tmp_path: Path):
    """A code whose proteome names a different species than the map claims."""
    _proteome(tmp_path, "lse/gastropoda", "6500", "Aplysia_californica", "aplcal")
    with pytest.raises(rca.SpeciesMapDriftError, match="mis-named"):
        rca.validate_species_map(
            references_dir=tmp_path,
            species_map={"aplcal": ("6500", "Aplysia_kurodai", "GCA_x")},
            ambiguous={})


def test_validator_refuses_an_undeclared_collision(tmp_path: Path):
    """Two species sharing one code, with nothing declaring the ambiguity."""
    _proteome(tmp_path, "lse/gastropoda", "109671", "Physella_acuta", "phau")
    _proteome(tmp_path, "lse/other", "115415", "Phoronis_australis", "phau")
    with pytest.raises(rca.SpeciesMapDriftError, match="not declared ambiguous"):
        rca.validate_species_map(
            references_dir=tmp_path,
            species_map={"phau": ("109671", "Physella_acuta", "GCF_x")},
            ambiguous={})


def test_validator_refuses_a_stale_ambiguity_register(tmp_path: Path):
    """If the rename IS applied one day, the register must be updated, not left
    behind to suppress lookups for a code that no longer collides."""
    _proteome(tmp_path, "lse/gastropoda", "109671", "Physella_acuta", "phau")
    with pytest.raises(rca.SpeciesMapDriftError, match="no longer collide"):
        rca.validate_species_map(
            references_dir=tmp_path,
            species_map={"phau": ("109671", "Physella_acuta", "GCF_x")},
            ambiguous={"phau": ("Physella_acuta", "Phoronis_australis")})


def test_validator_refuses_a_changed_collision_membership(tmp_path: Path):
    """A THIRD species joining the collision must not be silently absorbed."""
    _proteome(tmp_path, "a", "1", "Physella_acuta", "phau")
    _proteome(tmp_path, "b", "2", "Phoronis_australis", "phau")
    _proteome(tmp_path, "c", "3", "Physa_fontinalis", "phau")
    with pytest.raises(rca.SpeciesMapDriftError, match="declared ambiguous over"):
        rca.validate_species_map(
            references_dir=tmp_path,
            species_map={"phau": ("1", "Physella_acuta", "GCF_x")},
            ambiguous={"phau": ("Physella_acuta", "Phoronis_australis")})


# ---------------------------------------------------------------------------
# 4. Guard against over-correction
# ---------------------------------------------------------------------------

def test_validator_does_not_require_total_coverage(tmp_path: Path):
    """SPECIES_MAP is a curated subset of species that have recoverable genome
    assemblies -- measured 2026-07-20: 57 entries vs 238 observed header codes,
    so 182 unmapped codes are CORRECT, not drift.

    Porting PREFIX_TO_GROUP's 'missing codes' check verbatim would fire 182
    false positives and make the validator useless.
    """
    _proteome(tmp_path, "lse/gastropoda", "6500", "Aplysia_californica", "aplcal")
    _proteome(tmp_path, "lse/annelida", "6358", "Capitella_teleta", "cate")
    rca.validate_species_map(
        references_dir=tmp_path,
        species_map={"aplcal": ("6500", "Aplysia_californica", "GCA_x")},
        ambiguous={})


def test_clean_map_passes(tmp_path: Path):
    _proteome(tmp_path, "lse/gastropoda", "6500", "Aplysia_californica", "aplcal")
    observed = rca.validate_species_map(
        references_dir=tmp_path,
        species_map={"aplcal": ("6500", "Aplysia_californica", "GCA_x")},
        ambiguous={})
    assert observed["aplcal"].count == 3
