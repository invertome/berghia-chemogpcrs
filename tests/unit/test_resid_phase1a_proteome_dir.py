"""RESID REPAIR 2 — the P5 Phase 1a manifest must look where proteomes are.

``scripts/build_p5_phase1a_manifest.py`` defaulted ``--proteomes-dir`` to
``references/species_tree/cache/proteomes``. No producer writes there.

MEASURED on the cluster (Unity checkout
``/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05``,
metadata listing only, 2026-07-20):

    references/species_tree/cache/proteomes        -> No such file or directory
    references/species_tree/cache/                  -> proteomes_braker4,
                                                      proteomes_braker4_anisus_*
    species_tree_data/ncbi_proteomes/              -> 133 *.faa
                                                      133 *.cds.fna
    e.g. 100268_Paragonimus_heterotremus.faa
         100452_Candidula_unifasciata.faa

So the real Phase 1a proteomes live in ``species_tree_data/ncbi_proteomes``
as ``<taxid>_<sanitized_binomial>.faa`` -- exactly the layout
``consolidate_proteomes_for_genome_wide_og.py`` already encodes in its
declarative ``PHASE_PROTEOME_DIRS`` table after the same defect was fixed
there (where it produced 538 ``missing_proteome`` rows and exit 0).

Two failure modes are tested:

  1. WRONG DIRECTORY: every species resolves to nothing, and the operator is
     handed N individually "missing" species instead of one bad path.
  2. SILENT SUCCESS: that run exited 0 and wrote a header-only manifest, so
     the downstream SLURM scan array simply had no work and also succeeded.

Note the ``.cds.fna`` NUCLEOTIDE files sitting beside the proteins in the
same directory: any suffix probing must not reach them.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import build_p5_phase1a_manifest as bpm  # noqa: E402
import consolidate_proteomes_for_genome_wide_og as cons  # noqa: E402


MANIFEST_COLUMNS = (
    "taxid", "binomial", "clade", "source", "accession",
    "assembly_level", "annotation_status", "est_protein_count",
    "submission_date", "drop_reason",
)

# Real rows, copied from the format the cluster listing showed.
REAL_ROWS = [
    {"taxid": "100268", "binomial": "Paragonimus heterotremus",
     "accession": "GCA_012934845.1"},
    {"taxid": "100452", "binomial": "Candidula unifasciata",
     "accession": "GCA_905116865.2"},
]


def _write_manifest(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r.get(c, "") for c in MANIFEST_COLUMNS) + "\n")


def _write_fasta(path: Path, n_seqs: int = 3) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for i in range(n_seqs):
            f.write(f">seq{i}\nMKTLL\n")


def _real_layout(base_dir: Path) -> Path:
    """Recreate the measured cluster layout under a tmp base dir."""
    d = base_dir / "species_tree_data" / "ncbi_proteomes"
    for r in REAL_ROWS:
        leaf = f"{r['taxid']}_{bpm.sanitize_sample_name(r['binomial'])}"
        _write_fasta(d / f"{leaf}.faa")
        # the nucleotide companion that really sits next to it
        _write_fasta(d / f"{leaf}.cds.fna")
    return d


# --------------------------------------------------------------------------
# the phase -> directory table
# --------------------------------------------------------------------------
def test_phase1a_canonical_dir_is_where_the_producer_actually_writes():
    dirs = bpm.PHASE_PROTEOME_DIRS["1a"]
    assert dirs[0] == "species_tree_data/ncbi_proteomes", (
        "the canonical (producer) Phase 1a directory must lead; measured on "
        "the cluster as the only directory holding the 133 *.faa proteomes"
    )


def test_phase1a_table_is_shared_with_the_consolidator_not_re_declared():
    """One source of truth. Two copies of a path table drift, and the drift
    is silent -- which is how this defect survived being fixed next door."""
    assert bpm.PHASE_PROTEOME_DIRS is cons.PHASE_PROTEOME_DIRS, (
        "build_p5_phase1a_manifest must reuse "
        "consolidate_proteomes_for_genome_wide_og.PHASE_PROTEOME_DIRS, "
        "not re-declare its own copy"
    )
    assert bpm.PROTEOME_SUFFIXES is cons.PROTEOME_SUFFIXES


def test_the_dead_default_is_not_the_canonical_directory():
    """references/species_tree/cache/proteomes does not exist on the cluster.
    It may survive as a trailing legacy fallback, never as the lead."""
    assert bpm.PHASE_PROTEOME_DIRS["1a"][0] != "references/species_tree/cache/proteomes"


def test_bare_fna_is_not_probed_so_cds_nucleotides_cannot_be_picked_up():
    """`<leaf>.cds.fna` sits beside `<leaf>.faa` in the real directory."""
    assert ".fna" not in bpm.PROTEOME_SUFFIXES, (
        "a bare .fna probe would feed CDS nucleotides in as if they were "
        "protein; the real ncbi_proteomes dir holds 133 of them"
    )


# --------------------------------------------------------------------------
# default resolution
# --------------------------------------------------------------------------
def test_default_run_finds_the_real_layout_with_no_proteomes_dir_flag(tmp_path):
    """The whole defect in one test: with NO --proteomes-dir (and nothing in
    the tree passes one -- verified by grep across *.sh), the default must
    resolve the proteomes that are actually on disk."""
    _real_layout(tmp_path)
    mf = tmp_path / "proteome_manifest.tsv"
    _write_manifest(mf, REAL_ROWS)
    out = tmp_path / "p5_manifest.tsv"

    rc = bpm.main([
        "--manifest", str(mf),
        "--base-dir", str(tmp_path),
        "--out", str(out),
    ])

    assert rc == 0
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert len(rows) == len(REAL_ROWS), (
        "the default proteomes dir did not resolve the real cluster layout"
    )
    for r in rows:
        p = Path(r["proteome_path"])
        assert p.is_absolute() and p.exists() and p.suffix == ".faa"


def test_explicit_proteomes_dir_still_overrides_the_table(tmp_path):
    """An operator pointing somewhere explicit must win."""
    custom = tmp_path / "elsewhere"
    _write_fasta(custom / "100268_Paragonimus_heterotremus.faa")
    mf = tmp_path / "m.tsv"
    _write_manifest(mf, REAL_ROWS[:1])
    out = tmp_path / "out.tsv"

    rc = bpm.main([
        "--manifest", str(mf),
        "--proteomes-dir", str(custom),
        "--out", str(out),
    ])
    assert rc == 0
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert Path(rows[0]["proteome_path"]).parent == custom.resolve()


def test_suffix_probing_finds_a_non_faa_protein_file(tmp_path):
    leaf = "100268_Paragonimus_heterotremus"
    d = tmp_path / "species_tree_data" / "ncbi_proteomes"
    _write_fasta(d / f"{leaf}.aa")
    mf = tmp_path / "m.tsv"
    _write_manifest(mf, REAL_ROWS[:1])
    out = tmp_path / "out.tsv"

    rc = bpm.main([
        "--manifest", str(mf), "--base-dir", str(tmp_path), "--out", str(out),
    ])
    assert rc == 0
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert rows[0]["proteome_path"].endswith(".aa")


# --------------------------------------------------------------------------
# the zero-resolution assertion
# --------------------------------------------------------------------------
def test_zero_resolution_is_loud_not_a_header_only_manifest(tmp_path, capsys):
    """THE silent-success guard.

    A manifest with species and a directory that resolves NONE of them is a
    path/convention failure, not N independently missing species. Before this
    guard the run exited 0 with a header-only manifest and the downstream
    SLURM scan array simply had nothing to do -- also succeeding.
    """
    mf = tmp_path / "m.tsv"
    _write_manifest(mf, REAL_ROWS)          # species present
    out = tmp_path / "out.tsv"              # ...but no proteomes anywhere

    rc = bpm.main([
        "--manifest", str(mf), "--base-dir", str(tmp_path), "--out", str(out),
    ])

    assert rc != 0, (
        "0 of 2 species resolved and the run still succeeded -- exactly the "
        "silent success this guard exists to prevent"
    )
    err = capsys.readouterr().err
    assert "species_tree_data/ncbi_proteomes" in err, (
        "the error must name the directories actually probed so the operator "
        "chases one bad path, not N missing species"
    )
    assert "0/2" in err or "0 of 2" in err


def test_partial_resolution_stays_a_warning(tmp_path, capsys):
    """A genuine per-species gap is not a path failure. Only a TOTAL
    non-resolution is."""
    d = tmp_path / "species_tree_data" / "ncbi_proteomes"
    _write_fasta(d / "100268_Paragonimus_heterotremus.faa")   # 1 of 2
    mf = tmp_path / "m.tsv"
    _write_manifest(mf, REAL_ROWS)
    out = tmp_path / "out.tsv"

    rc = bpm.main([
        "--manifest", str(mf), "--base-dir", str(tmp_path), "--out", str(out),
    ])
    assert rc == 0
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert len(rows) == 1
    assert "100452" in capsys.readouterr().err       # the gap is still named


def test_empty_manifest_does_not_trip_the_guard(tmp_path):
    """No species in, no species out -- that is a real measurement, not a
    path failure, so it must not be reported as one."""
    mf = tmp_path / "m.tsv"
    _write_manifest(mf, [])
    out = tmp_path / "out.tsv"

    rc = bpm.main([
        "--manifest", str(mf), "--base-dir", str(tmp_path), "--out", str(out),
    ])
    assert rc == 0


def test_allow_empty_escape_hatch(tmp_path):
    """Legitimate 'the download has not run yet' case stays expressible --
    but it must be requested explicitly, never be the default."""
    mf = tmp_path / "m.tsv"
    _write_manifest(mf, REAL_ROWS)
    out = tmp_path / "out.tsv"

    rc = bpm.main([
        "--manifest", str(mf), "--base-dir", str(tmp_path), "--out", str(out),
        "--allow-empty",
    ])
    assert rc == 0
    assert out.read_text().startswith("taxid\t")


def test_guard_does_not_fire_when_output_is_skipped_for_idempotency(tmp_path):
    """The idempotent early return happens before any resolution, so there is
    nothing measured and nothing to assert about."""
    out = tmp_path / "out.tsv"
    out.write_text("taxid\tbinomial\tproteome_path\n999\tFoo bar\t/p.faa\n")
    rc = bpm.main([
        "--manifest", str(tmp_path / "nonexistent.tsv"),
        "--base-dir", str(tmp_path),
        "--out", str(out),
    ])
    assert rc == 0
    assert "999" in out.read_text()
