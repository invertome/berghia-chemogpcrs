"""The CDS-recovery path must REFUSE any species code declared ambiguous.

``AMBIGUOUS_SPECIES_CODES`` registers header-prefix codes carried by more than
one species' proteome. ``SPECIES_MAP`` resolves each code to exactly ONE genome
assembly, so for a registered code that mapping is a guess:
``get_missing_proteins`` collects every sequence whose id starts with
``<code>_`` and ``process_species`` aligns the whole batch against that single
assembly.

Live instance (measured 2026-07-20): ``phau`` is the prefix of 1,078 *Physella
acuta* headers (Gastropoda) and 428 *Phoronis australis* headers (Phoronida).
``process_species('phau', ...)`` pulled all 1,506 and ran miniprot against
``GCF_028476545.2`` -- the *Physella acuta* genome -- so 428 phoronid proteins
had their CDS "recovered" from a snail. Nothing failed; the alignments simply
came back worse, and the wrong-species CDS flowed on into the dN/dS axis.

The register and ``species_code_lookup.species_for()`` already DECLARED this
collision. Nothing in the recovery path consulted the declaration, which is
this repo's recurring defect shape: a correct guard sitting at the wrong seam.

These tests pin the refusal at every seam that resolves a code, and they drive
the guard from the register rather than from the literal string 'phau', so a
future collision is refused the first time it is declared.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import recover_cds_from_assemblies as rca  # noqa: E402


# A synthetic register + map, so the tests prove the guard is DATA-driven and
# do not depend on 'phau' remaining the live collision after remediation.
SYNTHETIC_AMBIGUOUS = {"zzcol": ("Genus_alpha", "Genus_beta")}


@pytest.fixture()
def og_dir(tmp_path):
    """An orthogroup dir whose 'zzcol_' prefix is carried by two species."""
    d = tmp_path / "og"
    d.mkdir()
    (d / "1_Genus_alpha.faa").write_text(">zzcol_a1\nMKV\n>zzcol_a2\nMKW\n")
    (d / "2_Genus_beta.faa").write_text(">zzcol_b1\nMKY\n")
    (d / "3_Genus_clean.faa").write_text(">clean_c1\nMKF\n")
    return d


@pytest.fixture()
def explode(monkeypatch):
    """Make every downstream step fatal, so a refusal that is not the FIRST
    thing to happen shows up as the wrong exception."""
    def boom(name):
        def _f(*_a, **_k):
            raise AssertionError(
                f"{name}() reached for an ambiguous code -- the refusal did "
                f"not fire before the alignment path")
        return _f

    monkeypatch.setattr(rca, "download_assembly", boom("download_assembly"))
    monkeypatch.setattr(rca, "run_miniprot", boom("run_miniprot"))
    monkeypatch.setattr(rca, "extract_cds_from_gff", boom("extract_cds_from_gff"))


# ---------------------------------------------------------------------------
# 1. The register is consulted, not a hardcoded code
# ---------------------------------------------------------------------------

def test_guard_is_driven_by_the_register_not_a_hardcoded_code():
    """A collision declared for the first time must be refused immediately."""
    with pytest.raises(rca.SpeciesCodeCollisionError) as exc:
        rca.assert_species_code_unambiguous(
            "zzcol", ambiguous=SYNTHETIC_AMBIGUOUS)
    assert "Genus alpha" in str(exc.value) or "Genus_alpha" in str(exc.value)
    assert "Genus beta" in str(exc.value) or "Genus_beta" in str(exc.value)


def test_unambiguous_code_passes_the_guard():
    rca.assert_species_code_unambiguous("losc", ambiguous=SYNTHETIC_AMBIGUOUS)
    rca.assert_species_code_unambiguous("clean", ambiguous=SYNTHETIC_AMBIGUOUS)


def test_every_registered_collision_is_refused_by_process_species(monkeypatch, tmp_path):
    """Data-driven over the LIVE register: no declared code may proceed."""
    assert rca.AMBIGUOUS_SPECIES_CODES, "register is empty; nothing to assert"
    for code in rca.AMBIGUOUS_SPECIES_CODES:
        with pytest.raises(rca.SpeciesCodeCollisionError):
            rca.process_species(code, str(tmp_path), str(tmp_path / "cds.fna"),
                                str(tmp_path), str(tmp_path))


# ---------------------------------------------------------------------------
# 2. process_species refuses -- loudly, and BEFORE any alignment work
# ---------------------------------------------------------------------------

def test_process_species_refuses_ambiguous_code(monkeypatch, tmp_path, og_dir, explode):
    monkeypatch.setattr(rca, "AMBIGUOUS_SPECIES_CODES", SYNTHETIC_AMBIGUOUS)
    monkeypatch.setitem(rca.SPECIES_MAP, "zzcol",
                        ("111", "Genus_alpha", "GCF_TEST.1"))
    with pytest.raises(rca.SpeciesCodeCollisionError) as exc:
        rca.process_species("zzcol", str(og_dir), str(tmp_path / "cds.fna"),
                            str(tmp_path), str(tmp_path))
    msg = str(exc.value)
    # Names BOTH colliding species...
    assert "Genus_alpha" in msg and "Genus_beta" in msg
    # ...names the single assembly that would have been (mis)used...
    assert "GCF_TEST.1" in msg
    # ...and says why the output cannot be trusted.
    assert "wrong" in msg.lower()


def test_process_species_does_not_return_a_success_like_zero(monkeypatch, tmp_path, og_dir):
    """A skip that returns 0 is indistinguishable from 'no missing proteins'.

    Pre-fix, ``process_species`` had exactly two non-raising outcomes -- a real
    count, or ``return 0`` for the several skip paths. If the refusal were
    implemented as another ``return 0`` the caller could not tell a refused
    species from a complete one, and ``main()`` would print a clean SUMMARY and
    exit 0. Assert the refusal is an exception, never a value.
    """
    monkeypatch.setattr(rca, "AMBIGUOUS_SPECIES_CODES", SYNTHETIC_AMBIGUOUS)
    monkeypatch.setitem(rca.SPECIES_MAP, "zzcol",
                        ("111", "Genus_alpha", "GCF_TEST.1"))
    try:
        result = rca.process_species("zzcol", str(og_dir),
                                     str(tmp_path / "cds.fna"),
                                     str(tmp_path), str(tmp_path))
    except rca.SpeciesCodeCollisionError:
        return
    pytest.fail(f"ambiguous code returned {result!r} instead of raising")


def test_reextract_only_path_also_refuses(monkeypatch, tmp_path, og_dir, explode):
    """--reextract-only skips miniprot but still resolves code -> genome."""
    monkeypatch.setattr(rca, "AMBIGUOUS_SPECIES_CODES", SYNTHETIC_AMBIGUOUS)
    monkeypatch.setitem(rca.SPECIES_MAP, "zzcol",
                        ("111", "Genus_alpha", "GCF_TEST.1"))
    (tmp_path / "zzcol_miniprot.gff").write_text("##gff-version 3\n")
    (tmp_path / "GCF_TEST.1.fa").write_text(">chr1\nACGT\n")
    with pytest.raises(rca.SpeciesCodeCollisionError):
        rca.process_species("zzcol", str(og_dir), str(tmp_path / "cds.fna"),
                            str(tmp_path), str(tmp_path), reextract_only=True)


def test_unambiguous_code_still_processes(monkeypatch, tmp_path, og_dir):
    """Regression: the guard must not block ordinary species."""
    monkeypatch.setattr(rca, "AMBIGUOUS_SPECIES_CODES", SYNTHETIC_AMBIGUOUS)
    monkeypatch.setitem(rca.SPECIES_MAP, "clean",
                        ("222", "Genus_clean", "GCF_CLEAN.1"))
    monkeypatch.setattr(rca, "download_assembly",
                        lambda *a, **k: str(tmp_path / "g.fa"))
    monkeypatch.setattr(rca, "run_miniprot", lambda *a, **k: True)
    monkeypatch.setattr(rca, "extract_cds_from_gff",
                        lambda *a, **k: ({"clean_c1": "ATGAAATTT"},
                                         {"found": 1, "no_hit": 0,
                                          "bad_cds": 0, "multi_hit": 0}))
    assert rca.process_species("clean", str(og_dir), str(tmp_path / "cds.fna"),
                               str(tmp_path), str(tmp_path)) == 1


# ---------------------------------------------------------------------------
# 3. get_missing_proteins -- the seam where the contamination is assembled
# ---------------------------------------------------------------------------

def test_get_missing_proteins_refuses_ambiguous_prefix(monkeypatch, tmp_path, og_dir):
    """THE defect, at its source.

    Pre-fix this returned {'zzcol_a1','zzcol_a2','zzcol_b1'} -- both species'
    sequences fused into one query batch destined for one genome. The function
    is importable and callable independently of ``process_species``, so it
    carries its own refusal rather than relying on its caller's.
    """
    monkeypatch.setattr(rca, "AMBIGUOUS_SPECIES_CODES", SYNTHETIC_AMBIGUOUS)
    with pytest.raises(rca.SpeciesCodeCollisionError):
        rca.get_missing_proteins(str(og_dir), str(tmp_path / "cds.fna"), "zzcol")


def test_get_missing_proteins_still_works_for_clean_prefix(monkeypatch, tmp_path, og_dir):
    monkeypatch.setattr(rca, "AMBIGUOUS_SPECIES_CODES", SYNTHETIC_AMBIGUOUS)
    got = rca.get_missing_proteins(str(og_dir), str(tmp_path / "cds.fna"), "clean")
    assert set(got) == {"clean_c1"}


# ---------------------------------------------------------------------------
# 4. main() -- fails fast, and fails VISIBLY at the process boundary
# ---------------------------------------------------------------------------

def test_main_refuses_before_processing_any_species(tmp_path):
    """Pre-flight: a full run must abort at second 0, not after N species.

    Run as a subprocess so the assertion covers the real exit status -- a
    refusal that raises inside ``main()`` but exits 0 would be the same
    silent-success defect in a new costume.
    """
    script = PROJECT_ROOT / "scripts" / "recover_cds_from_assemblies.py"
    code = sorted(rca.AMBIGUOUS_SPECIES_CODES)[0]
    r = subprocess.run(
        [sys.executable, str(script), "--species", code,
         "--og-dir", str(tmp_path), "--cds-file", str(tmp_path / "cds.fna"),
         "--output-dir", str(tmp_path / "out")],
        capture_output=True, text=True, timeout=120)
    assert r.returncode != 0, "ambiguous code exited 0 (silent success)"
    combined = r.stdout + r.stderr
    for species in rca.AMBIGUOUS_SPECIES_CODES[code]:
        assert species in combined, f"error text does not name {species}"
    assert "Downloading" not in combined, "download started despite refusal"


def test_main_preflight_lists_every_offending_code(monkeypatch, tmp_path):
    """The pre-flight reports ALL offending codes, not just the first."""
    monkeypatch.setattr(rca, "AMBIGUOUS_SPECIES_CODES",
                        {"zzcol": ("Genus_alpha", "Genus_beta"),
                         "zzcolb": ("Genus_gamma", "Genus_delta")})
    with pytest.raises(rca.SpeciesCodeCollisionError) as exc:
        rca.assert_species_list_unambiguous(["losc", "zzcol", "zzcolb"])
    msg = str(exc.value)
    assert "zzcol" in msg and "zzcolb" in msg


def test_main_exit_status_is_propagated():
    """``main()`` returns 1 for a failed audit; ``__main__`` must not drop it."""
    source = (PROJECT_ROOT / "scripts" / "recover_cds_from_assemblies.py").read_text()
    tail = source[source.index("if __name__ =="):]
    assert "sys.exit(main())" in tail, (
        "main()'s return code is discarded, so --validate-species-map exits 0 "
        "even when the audit fails")
