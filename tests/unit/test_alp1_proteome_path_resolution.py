"""alp1: per-phase proteome DIRECTORY resolution + loud zero-match guard.

Every path and filename in these fixtures was verified against the live Unity
checkout ``/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05``
on 2026-07-20; the per-test comments record what was measured where.

Why this file exists instead of extra cases in
``test_consolidate_proteomes_for_genome_wide_og.py``: those tests build their
Phase-1a fixture at ``tmp_path/cache/proteomes`` — the directory the *consumer*
assumed — so they passed green while all 131 Phase-1a species on Unity resolved
to nothing (``consolidation_report.tsv``: 538 ``missing_proteome``, 0 ``ok``,
job exit 0). A fixture that encodes the assumed path cannot detect a wrong one.
These fixtures use the directory the *producer* actually writes, taken from
``download_species_tree_phase1a.py`` (``<cache_dir>/<taxid>_<binomial>.faa``)
with ``cache_dir`` = ``species_tree_data/ncbi_proteomes`` per
``scripts/unity/run_extension_proteome_harvest.sh`` (``NCBI_CACHE_DIR``).
"""
from __future__ import annotations

from pathlib import Path

import consolidate_proteomes_for_genome_wide_og as cons

# Column headers copied from the real manifests on Unity (verified 2026-07-20:
# `head -1 references/species_tree/{proteome_manifest,genome_inventory}.tsv`).
_PHASE1A_COLS = (
    "taxid", "binomial", "clade", "source", "accession", "assembly_level",
    "annotation_status", "est_protein_count", "submission_date", "drop_reason",
)
_INVENTORY_COLS = (
    "taxid", "binomial", "clade", "policy_class", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "contig_n50", "total_length_bp", "drop_reason",
    "source_batch",
)

# Real records, copied from the Unity manifests rather than invented.
#   6183 Schistosoma mansoni GCF_000237925.1 -> ncbi_proteomes/6183_Schistosoma_mansoni.faa
#     (species_tree_data/ncbi_proteomes/download_report.tsv, status=ok)
#   102321 Hyotissa hyotis -> proteomes_braker4/102321_Hyotissa_hyotis.aa.fna
#     (first entry of `ls references/species_tree/cache/proteomes_braker4`)
_P1A_TAXID, _P1A_BINOMIAL = 6183, "Schistosoma mansoni"
_P1A_LEAF = "6183_Schistosoma_mansoni"
_P1F_TAXID, _P1F_BINOMIAL = 102321, "Hyotissa hyotis"
_P1F_LEAF = "102321_Hyotissa_hyotis"


def _write_tsv(path: Path, cols, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in cols) + "\n")


def _phase1a_manifest(path: Path) -> Path:
    _write_tsv(path, _PHASE1A_COLS, [{
        "taxid": str(_P1A_TAXID), "binomial": _P1A_BINOMIAL,
        "clade": "platyhelminthes", "accession": "GCF_000237925.1",
    }])
    return path


def _inventory_manifest(path: Path) -> Path:
    _write_tsv(path, _INVENTORY_COLS, [{
        "taxid": str(_P1F_TAXID), "binomial": _P1F_BINOMIAL,
        "clade": "bivalvia", "accession": "GCA_000000001.1",
    }])
    return path


def _put(path: Path, text: str = ">p1\nMKV\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
# Directory resolution, per phase
# ---------------------------------------------------------------------------

class TestPhaseDirectoryResolution:

    def test_phase1a_resolves_from_ncbi_proteomes_dir(self, tmp_path: Path) -> None:
        """Phase 1a lives in species_tree_data/ncbi_proteomes/<leaf>.faa.

        Verified on Unity 2026-07-20: that directory holds 133 `.faa` files and
        130 of the 131 accession-bearing proteome_manifest.tsv leaves match one.
        `references/species_tree/cache/proteomes/` (what the code read) does not
        exist there at all — `ls references/species_tree/cache/` returns only
        the proteomes_braker4* entries.
        """
        want = _put(tmp_path / "species_tree_data" / "ncbi_proteomes" / f"{_P1A_LEAF}.faa")
        got = cons._locate_proteome(
            _P1A_TAXID, _P1A_BINOMIAL, _P1A_LEAF, "1a",
            tmp_path, tmp_path / "braker4_output",
        )
        assert got == want

    def test_phase1a_canonical_dir_wins_over_legacy_cache_dir(self, tmp_path: Path) -> None:
        """When both exist, the producer's directory must win.

        Guards a stale/partial legacy tree from shadowing live downloads.
        """
        canonical = _put(
            tmp_path / "species_tree_data" / "ncbi_proteomes" / f"{_P1A_LEAF}.faa",
            ">canonical\nMKV\n",
        )
        _put(
            tmp_path / "references" / "species_tree" / "cache" / "proteomes" / f"{_P1A_LEAF}.faa",
            ">legacy\nMKV\n",
        )
        got = cons._locate_proteome(
            _P1A_TAXID, _P1A_BINOMIAL, _P1A_LEAF, "1a",
            tmp_path, tmp_path / "braker4_output",
        )
        assert got == canonical

    def test_phase1f_resolves_from_braker4_cache(self, tmp_path: Path) -> None:
        """Regression guard: Phase 1f's directory is already correct.

        Verified on Unity 2026-07-20: `references/species_tree/cache/
        proteomes_braker4/` holds 132 `.aa.fna` + 132 `.cds.fna`, and all 132
        `.aa.fna` stems match a genome_inventory.tsv leaf (132/132 overlap).
        Repairing Phase 1a must not move this one.
        """
        want = _put(
            tmp_path / "references" / "species_tree" / "cache"
            / "proteomes_braker4" / f"{_P1F_LEAF}.aa.fna"
        )
        got = cons._locate_proteome(
            _P1F_TAXID, _P1F_BINOMIAL, _P1F_LEAF, "1f",
            tmp_path, tmp_path / "braker4_output",
        )
        assert got == want

    def test_missing_phase1a_reports_canonical_directory(self, tmp_path: Path) -> None:
        """A miss must name the directory an operator should go look in."""
        got = cons._locate_proteome(
            _P1A_TAXID, _P1A_BINOMIAL, _P1A_LEAF, "1a",
            tmp_path, tmp_path / "braker4_output",
        )
        assert got == tmp_path / "species_tree_data" / "ncbi_proteomes" / f"{_P1A_LEAF}.faa"


# ---------------------------------------------------------------------------
# Loud per-phase zero-match guard
# ---------------------------------------------------------------------------

class TestPerPhaseZeroMatchGuard:

    def test_main_fails_when_one_phase_resolves_nothing(self, tmp_path: Path, capsys) -> None:
        """THE regression test for bead alp1.

        Phase 1f resolves; Phase 1a resolves nothing because its directory is
        wrong. The old global `n_ok == 0` guard cannot see this — one healthy
        phase masks a totally broken one — which is exactly how 538
        `missing_proteome` rows were reported as a clean exit 0.
        """
        _put(
            tmp_path / "references" / "species_tree" / "cache"
            / "proteomes_braker4" / f"{_P1F_LEAF}.aa.fna"
        )
        # Phase 1a proteome deliberately absent.
        rc = cons.main([
            "--out-dir", str(tmp_path / "of_input"),
            "--base-dir", str(tmp_path),
            "--braker4-output-dir", str(tmp_path / "braker4_output"),
            "--phase1a-manifest", str(_phase1a_manifest(tmp_path / "proteome_manifest.tsv")),
            "--manifest", str(_inventory_manifest(tmp_path / "genome_inventory.tsv")),
        ])
        assert rc == 2
        err = capsys.readouterr().err
        assert "1a" in err
        assert "0/1" in err or "0 of 1" in err

    def test_main_succeeds_when_every_phase_resolves(self, tmp_path: Path) -> None:
        _put(tmp_path / "species_tree_data" / "ncbi_proteomes" / f"{_P1A_LEAF}.faa")
        _put(
            tmp_path / "references" / "species_tree" / "cache"
            / "proteomes_braker4" / f"{_P1F_LEAF}.aa.fna"
        )
        rc = cons.main([
            "--out-dir", str(tmp_path / "of_input"),
            "--base-dir", str(tmp_path),
            "--braker4-output-dir", str(tmp_path / "braker4_output"),
            "--phase1a-manifest", str(_phase1a_manifest(tmp_path / "proteome_manifest.tsv")),
            "--manifest", str(_inventory_manifest(tmp_path / "genome_inventory.tsv")),
        ])
        assert rc == 0
        assert (tmp_path / "of_input" / f"{_P1A_LEAF}.fa").exists()
        assert (tmp_path / "of_input" / f"{_P1F_LEAF}.fa").exists()

    def test_partial_gap_within_a_resolving_phase_stays_a_warning(self, tmp_path: Path) -> None:
        """Only a TOTAL phase failure is fatal.

        Phase 1f is mid-flight on Unity (132 of 453 species annotated as of
        2026-07-20), so a phase that resolves some-but-not-all species must
        keep exiting 0 — 03c holds the hard cleanliness assertion.
        """
        _write_tsv(tmp_path / "genome_inventory.tsv", _INVENTORY_COLS, [
            {"taxid": str(_P1F_TAXID), "binomial": _P1F_BINOMIAL,
             "clade": "bivalvia", "accession": "GCA_000000001.1"},
            {"taxid": "999999", "binomial": "Absent species",
             "clade": "bivalvia", "accession": "GCA_000000002.1"},
        ])
        _put(
            tmp_path / "references" / "species_tree" / "cache"
            / "proteomes_braker4" / f"{_P1F_LEAF}.aa.fna"
        )
        rc = cons.main([
            "--out-dir", str(tmp_path / "of_input"),
            "--base-dir", str(tmp_path),
            "--braker4-output-dir", str(tmp_path / "braker4_output"),
            "--manifest", str(tmp_path / "genome_inventory.tsv"),
        ])
        assert rc == 0
