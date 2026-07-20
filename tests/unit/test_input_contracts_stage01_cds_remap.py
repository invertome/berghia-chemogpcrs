"""Stage 01's CDS header remap moved a file nothing ever wrote (audit finding F9).

``scripts/update_headers.py`` always writes ``f"{fasta_file}_updated.fa"`` --
it appends to the WHOLE input path, extension included. For the CDS input
``reference_sequences/cds/all_references_cds.fna`` that is::

    reference_sequences/cds/all_references_cds.fna_updated.fa

Stage 01 moved ``reference_sequences/cds/all_references_cds_updated.fna``
instead: different stem (``_cds_updated`` vs ``_cds.fna_updated``) AND
different extension (``.fna`` vs ``.fa``). The mismatch was masked by
``2>/dev/null || true``, so the stage reported success while:

  * ``all_references_cds.fna`` kept its ORIGINAL headers -- CDS IDs were never
    remapped onto the protein ``id_map``'s short IDs, so every downstream
    CDS-to-protein join ran on unmapped IDs; and
  * a stray ``all_references_cds.fna_updated.fa`` accumulated in cds/.

The adjacent protein remap (01:465) was always correct -- ``all_references.fa``
-> ``all_references.fa_updated.fa`` -- which is what makes the CDS line's
divergence a copy-paste defect rather than a convention.

These tests run the REAL ``update_headers.py`` to derive the filename (rather
than trusting a description of it) and then execute the REAL ``mv`` lines
lifted out of the shipping stage 01 script.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
STAGE01 = PROJECT_ROOT / "01_reference_processing.sh"
UPDATE_HEADERS = PROJECT_ROOT / "scripts" / "update_headers.py"

CDS_FASTA = ">ref_seq_one taxid=6500\nATGAAATTTCTGGTGGCA\n>ref_seq_two taxid=6500\nATGAAATTTCTGGTGGCT\n"


def _run_update_headers(fasta: Path, id_map: Path, *extra: str) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        [sys.executable, str(UPDATE_HEADERS), str(fasta), str(id_map), "--source-type", "reference", *extra],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc


def _mv_lines_for(marker: str) -> list[str]:
    """Lift the `mv` statement whose FIRST operand mentions `marker` out of stage 01."""
    lines = STAGE01.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.lstrip().startswith("mv ") and marker in line:
            # the statement may be wrapped over one or more continuation lines
            block = [line.lstrip()]
            while block[-1].rstrip().endswith("\\"):
                block.append(lines[i + len(block)].strip())
            return block
    pytest.fail(f"stage 01 has no `mv` statement mentioning {marker!r}")


# --------------------------------------------------------------------------
# what update_headers.py ACTUALLY writes
# --------------------------------------------------------------------------
def test_update_headers_appends_updated_fa_to_the_whole_input_path(tmp_path):
    """The contract stage 01 has to match, measured rather than assumed."""
    cds = tmp_path / "all_references_cds.fna"
    cds.write_text(CDS_FASTA)

    _run_update_headers(cds, tmp_path / "id_map.csv")

    assert (tmp_path / "all_references_cds.fna_updated.fa").exists()
    assert not (tmp_path / "all_references_cds_updated.fna").exists(), (
        "the path stage 01 used to move is never produced"
    )


def test_update_headers_rewrites_the_cds_headers(tmp_path):
    """Confirms the remap is real work, i.e. that losing it actually loses data."""
    cds = tmp_path / "all_references_cds.fna"
    cds.write_text(CDS_FASTA)

    _run_update_headers(cds, tmp_path / "id_map.csv")

    updated = (tmp_path / "all_references_cds.fna_updated.fa").read_text()
    assert ">ref_seq_one" not in updated, "headers should have been replaced by short IDs"
    assert updated.startswith(">")


# --------------------------------------------------------------------------
# the regression: stage 01 must move the file that exists
# --------------------------------------------------------------------------
def test_stage01_cds_mv_relocates_the_real_update_headers_output(tmp_path):
    """End to end: run update_headers.py, then stage 01's REAL mv line."""
    results = tmp_path / "results"
    cds_dir = results / "reference_sequences" / "cds"
    cds_dir.mkdir(parents=True)
    cds = cds_dir / "all_references_cds.fna"
    cds.write_text(CDS_FASTA)
    _run_update_headers(cds, tmp_path / "id_map.csv")
    remapped = (cds_dir / "all_references_cds.fna_updated.fa").read_text()

    snippet = "\n".join(
        ["log() { :; }", f'RESULTS_DIR="{results}"', *_mv_lines_for("all_references_cds")]
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)

    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert cds.read_text() == remapped, "the remapped CDS must land at all_references_cds.fna"
    assert not (cds_dir / "all_references_cds.fna_updated.fa").exists(), (
        "the stray _updated.fa file must not accumulate"
    )


def test_stage01_cds_mv_fails_loudly_when_the_remap_produced_nothing(tmp_path):
    """`2>/dev/null || true` is what hid this for the whole life of the bug."""
    results = tmp_path / "results"
    (results / "reference_sequences" / "cds").mkdir(parents=True)

    snippet = "\n".join(
        ["log() { :; }", f'RESULTS_DIR="{results}"', *_mv_lines_for("all_references_cds")]
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)

    assert proc.returncode != 0, (
        "a missing update_headers.py output must not be silently swallowed:\n"
        + proc.stdout
        + proc.stderr
    )


def test_stage01_cds_mv_does_not_suppress_errors():
    """No `2>/dev/null` and no unconditional `|| true` on the CDS mv."""
    block = "\n".join(_mv_lines_for("all_references_cds"))
    assert "2>/dev/null" not in block, "stderr suppression re-hides this defect"
    assert "|| true" not in block, "`|| true` re-hides this defect"


def test_stage01_cds_mv_targets_the_documented_filename():
    """Pin stage 01 to update_headers.py's naming rule, not to a hand-typed path."""
    block = "\n".join(_mv_lines_for("all_references_cds"))
    assert "all_references_cds.fna_updated.fa" in block
    assert "all_references_cds_updated.fna" not in block, (
        "this path has no writer; moving it was always a no-op"
    )


def test_update_headers_naming_rule_is_unchanged():
    """If update_headers.py's rule ever changes, both stage 01 mv lines must too."""
    assert 'f"{fasta_file}_updated.fa"' in UPDATE_HEADERS.read_text(), (
        "update_headers.py's output-naming rule moved; re-derive stage 01's mv targets"
    )


# --------------------------------------------------------------------------
# the adjacent protein remap must stay correct (it always was)
# --------------------------------------------------------------------------
def test_stage01_protein_mv_still_matches_the_same_rule(tmp_path):
    """Regression guard on the line that was already right."""
    block = "\n".join(_mv_lines_for("all_references.fa_updated.fa"))
    assert "all_references.fa_updated.fa" in block

    results = tmp_path / "results"
    ref_dir = results / "reference_sequences"
    ref_dir.mkdir(parents=True)
    prot = ref_dir / "all_references.fa"
    prot.write_text(">ref_one taxid=6500\nMKFLVAAALL\n")
    _run_update_headers(prot, tmp_path / "id_map.csv")
    remapped = (ref_dir / "all_references.fa_updated.fa").read_text()

    snippet = "\n".join(["log() { :; }", f'RESULTS_DIR="{results}"', block])
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)

    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert prot.read_text() == remapped
