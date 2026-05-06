"""Tests for ``scripts/recover_cds_from_assemblies.get_missing_proteins``.

Bug context (2026-05-06, follow-up to bead -lfy):
    The Unity B2 run (job 56724467) reported "No missing proteins for X,
    skipping" for all 55 species and recovered 0 CDS. Root cause: the
    glob in ``get_missing_proteins`` is::

        glob.glob(os.path.join(og_dir, '*.fa'))

    The Nath et al. reference tree contains:
      - 426 ``.faa`` files (NOT ``.fa``)
      - all of them in nested subdirectories (``lse/<phylum>/...``,
        ``one_to_one_ortholog/...``)
      - zero ``.fa`` files at the top level

    So the glob returns ``[]`` for every species, ``missing`` is always
    empty, and miniprot recovery is silently skipped.

These tests pin the contract:
    1. The protein scan must recurse into subdirectories.
    2. It must accept ``.faa`` (the real Nath extension), and ideally
       also ``.fa`` / ``.fasta`` for forward compatibility with other
       reference packages that use those extensions.
    3. Proteins that already have a CDS (in ``cds_file``) are excluded
       — regression guard for the existing behaviour.
"""
from __future__ import annotations

from pathlib import Path

# conftest.py adds scripts/ to sys.path
import recover_cds_from_assemblies as rca


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for sid, seq in records.items():
            fh.write(f">{sid}\n{seq}\n")


def test_finds_faa_in_nested_subdirectory(tmp_path: Path) -> None:
    """The Nath et al. tree puts protein FASTAs at
    ``<og_dir>/lse/<phylum>/<taxid>_<genus>_<species>.faa``. The scanner
    must descend into those subdirectories AND accept the ``.faa``
    extension."""
    og_dir = tmp_path / "nath_et_al"
    nested = og_dir / "lse" / "annelida"
    _write_fasta(
        nested / "283909_Capitella_teleta.faa",
        {
            "cate_KB300474.1_27": "MKLLRSV",
            "cate_KB292229.1_1519": "MAQTPLM",
        },
    )
    cds_file = tmp_path / "all_references_cds.fna"
    cds_file.write_text("")  # empty: nothing has a CDS yet

    missing = rca.get_missing_proteins(str(og_dir), str(cds_file), "cate")

    assert set(missing.keys()) == {"cate_KB300474.1_27", "cate_KB292229.1_1519"}, (
        "scanner must descend into subdirectories and accept .faa")


def test_finds_proteins_in_one_to_one_ortholog_subdir(tmp_path: Path) -> None:
    """The conserved (non-LSE) Nath proteins live in
    ``<og_dir>/one_to_one_ortholog/`` — also a subdirectory."""
    og_dir = tmp_path / "nath_et_al"
    sub = og_dir / "one_to_one_ortholog"
    _write_fasta(sub / "OG0000123.faa", {"aplcal_NM_001234_1": "MERTQ"})
    cds_file = tmp_path / "all_references_cds.fna"
    cds_file.write_text("")

    missing = rca.get_missing_proteins(str(og_dir), str(cds_file), "aplcal")

    assert "aplcal_NM_001234_1" in missing


def test_excludes_proteins_already_in_cds_file(tmp_path: Path) -> None:
    """Regression guard: proteins that already have a CDS in
    ``cds_file`` must NOT be returned as missing (otherwise miniprot
    re-recovers everything every run)."""
    og_dir = tmp_path / "nath_et_al"
    nested = og_dir / "lse" / "annelida"
    _write_fasta(
        nested / "283909_Capitella_teleta.faa",
        {
            "cate_present_1": "MKLLRSV",
            "cate_missing_2": "MAQTPLM",
        },
    )
    cds_file = tmp_path / "all_references_cds.fna"
    cds_file.write_text(">cate_present_1\nATGAAACTG\n")

    missing = rca.get_missing_proteins(str(og_dir), str(cds_file), "cate")

    assert "cate_present_1" not in missing, "already-present CDS must be excluded"
    assert "cate_missing_2" in missing


def test_filters_by_species_prefix(tmp_path: Path) -> None:
    """Multi-species ``.faa`` files exist in the conserved tree (one
    file per orthogroup, many species). The prefix filter must keep
    only the requested species' proteins."""
    og_dir = tmp_path / "nath_et_al"
    sub = og_dir / "one_to_one_ortholog"
    _write_fasta(
        sub / "OG0000456.faa",
        {
            "cate_X1": "MAAA",
            "aplcal_Y2": "MBBB",
            "phaust_Z3": "MCCC",
        },
    )
    cds_file = tmp_path / "all_references_cds.fna"
    cds_file.write_text("")

    missing = rca.get_missing_proteins(str(og_dir), str(cds_file), "cate")

    assert set(missing.keys()) == {"cate_X1"}


def test_accepts_legacy_dot_fa_extension(tmp_path: Path) -> None:
    """If a future / legacy reference package uses ``.fa`` (the original
    glob pattern), it must still work."""
    og_dir = tmp_path / "nath_et_al"
    nested = og_dir / "lse" / "mollusca"
    _write_fasta(nested / "1287507_Berghia.fa", {"bste_test_1": "MNNN"})
    cds_file = tmp_path / "all_references_cds.fna"
    cds_file.write_text("")

    missing = rca.get_missing_proteins(str(og_dir), str(cds_file), "bste")

    assert "bste_test_1" in missing


def test_returns_empty_dict_when_no_proteins_for_species(tmp_path: Path) -> None:
    """If the requested species has no proteins anywhere in the tree,
    return an empty dict (current contract — used by the caller's
    ``if not missing: skip`` shortcut)."""
    og_dir = tmp_path / "nath_et_al"
    nested = og_dir / "lse" / "annelida"
    _write_fasta(nested / "283909_Capitella_teleta.faa", {"cate_X1": "MAAA"})
    cds_file = tmp_path / "all_references_cds.fna"
    cds_file.write_text("")

    missing = rca.get_missing_proteins(str(og_dir), str(cds_file), "nonexistent")

    assert missing == {}
