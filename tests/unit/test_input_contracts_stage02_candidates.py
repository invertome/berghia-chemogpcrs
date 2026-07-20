"""The 02 -> 02a -> 02b candidate-FASTA input contract (audit finding F8/F12.5).

Stage 02 writes its chemoreceptor candidates to::

    ${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa      (Berghia)
    ${RESULTS_DIR}/chemogpcrs/chemogpcrs_<taxid>.fa      (every other sample)

Stage 02a's input-discovery loop probed four OTHER paths
(``candidates/chemogpcr_candidates.fa``, ``.fasta``,
``candidates/all_candidates.fa``, ``chemogpcrs/candidates.fa``). Grepping the
repository turns up exactly one writer for any of them --
``tests/validate_pipeline.sh`` (fixture generation) -- so on real data the loop
matched nothing, 02a exited 1, and 02b then failed its hard
``step_completed_02a.txt`` check. The contract was satisfied only by the test
harness, never by the pipeline.

Berghia-only is deliberate. 02a's single output is the one canonical
``candidates_clustered.fa`` consumed by 02b, which drives classification and
the candidate ranking / HCR shortlist -- all Berghia-scoped. The per-species
``chemogpcrs_<taxid>.fa`` files are read directly out of
``results/chemogpcrs/`` by stage 06 (06:181, 06:203) and never flow through
02a, and CD-HIT at 0.98 over a pooled multi-species set would collapse
cross-species orthologs onto a single representative, silently deleting
species from the per-species sets stage 06 depends on.

F12.5, same contract: 02b hardcoded ``candidates_nr098.fa`` while 02a computes
that basename from ``${CDHIT_IDENTITY}``. They agree only at the default 0.98;
any other threshold makes 02b miss 02a's output.

These tests execute the REAL discovery loops lifted out of the shipping stage
scripts, so they break if the probe lists drift.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
STAGE02 = PROJECT_ROOT / "02_chemogpcrs_identification.sh"
STAGE02A = PROJECT_ROOT / "02a_cluster_sequences.sh"
STAGE02B = PROJECT_ROOT / "02b_classify_gpcrs.sh"


# --------------------------------------------------------------------------
# harness
# --------------------------------------------------------------------------
def _extract_for_loop(script: Path, var: str) -> str:
    """Lift a top-level `for <var> in ... ; do ... done` block out of a stage script."""
    lines = script.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith(f"for {var} in"):
            start = i
            break
    if start is None:
        pytest.fail(f"{script.name} has no top-level `for {var} in ...` discovery loop")
    for j in range(start + 1, len(lines)):
        if lines[j] == "done":
            return "\n".join(lines[start : j + 1])
    pytest.fail(f"could not find the closing `done` of `for {var}` in {script.name}")


def _run_discovery(script: Path, tmp_path: Path, **env_extra) -> str:
    """Run a stage's real candidate-discovery loop against a synthetic ${RESULTS_DIR}.

    Returns the resolved INPUT_FASTA (empty string when the loop matched nothing).
    """
    results = tmp_path / "results"
    snippet = "\n".join(
        [
            f'RESULTS_DIR="{results}"',
            'CANDIDATES_DIR="${RESULTS_DIR}/candidates"',
            'INPUT_FASTA=""',
            *(f'{k}="{v}"' for k, v in env_extra.items()),
            _extract_for_loop(script, "candidate_file"),
            'echo "INPUT_FASTA=${INPUT_FASTA}"',
        ]
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    for line in proc.stdout.splitlines():
        if line.startswith("INPUT_FASTA="):
            return line.split("=", 1)[1]
    pytest.fail("discovery loop produced no INPUT_FASTA line: " + proc.stdout + proc.stderr)


def _write_fasta(path: Path, seq_id: str = "s1") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{seq_id}\nMKFLVAAALLAVSSA\n")
    return path


# --------------------------------------------------------------------------
# F8: the real stage-02 output must satisfy 02a's discovery loop
# --------------------------------------------------------------------------
def test_02a_finds_the_berghia_fasta_stage02_actually_writes(tmp_path):
    """The core regression: only stage 02's real output exists, and 02a finds it."""
    real = _write_fasta(tmp_path / "results" / "chemogpcrs" / "chemogpcrs_berghia.fa")

    assert _run_discovery(STAGE02A, tmp_path) == str(real)


def test_02a_probes_the_path_stage02_documents_as_its_output():
    """Pin the two sides of the contract to each other, not to my description of it."""
    stage02 = STAGE02.read_text()
    assert '"${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"' in stage02, (
        "stage 02 no longer writes chemogpcrs_berghia.fa; re-derive 02a's probe list"
    )
    assert "chemogpcrs/chemogpcrs_berghia.fa" in STAGE02A.read_text(), (
        "02a must probe the path stage 02 actually writes, not only fixture paths"
    )


def test_02a_still_honours_the_pre_existing_probe_paths(tmp_path):
    """Additive fix: an explicitly-staged candidates/ file keeps working."""
    legacy = _write_fasta(tmp_path / "results" / "candidates" / "chemogpcr_candidates.fa")

    assert _run_discovery(STAGE02A, tmp_path) == str(legacy)


def test_02a_prefers_an_explicit_candidates_file_over_the_stage02_default(tmp_path):
    """Priority is stable: a deliberately-staged override wins over the default."""
    legacy = _write_fasta(tmp_path / "results" / "candidates" / "chemogpcr_candidates.fa")
    _write_fasta(tmp_path / "results" / "chemogpcrs" / "chemogpcrs_berghia.fa")

    assert _run_discovery(STAGE02A, tmp_path) == str(legacy)


def test_02a_finds_nothing_when_stage02_produced_nothing(tmp_path):
    """No candidates anywhere -> empty, so 02a's `[ -z ]` branch still exits 1."""
    (tmp_path / "results" / "chemogpcrs").mkdir(parents=True)

    assert _run_discovery(STAGE02A, tmp_path) == ""


# --------------------------------------------------------------------------
# F8, the multi-species question: Berghia-only, deliberately
# --------------------------------------------------------------------------
def test_02a_ignores_the_per_species_chemogpcrs_files(tmp_path):
    """Per-species stage-02 outputs must NOT be picked up as clustering input.

    They are stage 06's inputs (06:181, 06:203) and are read straight from
    results/chemogpcrs/. Feeding them to CD-HIT here would pool species into
    one non-redundant set and drop species from those per-species files.
    """
    _write_fasta(tmp_path / "results" / "chemogpcrs" / "chemogpcrs_9606.fa")
    _write_fasta(tmp_path / "results" / "chemogpcrs" / "chemogpcrs_6500.fa")

    assert _run_discovery(STAGE02A, tmp_path) == ""


def test_02a_picks_berghia_even_when_per_species_files_sit_beside_it(tmp_path):
    """The Berghia file is selected by name, not by glob order."""
    _write_fasta(tmp_path / "results" / "chemogpcrs" / "chemogpcrs_6500.fa")
    berghia = _write_fasta(tmp_path / "results" / "chemogpcrs" / "chemogpcrs_berghia.fa")
    _write_fasta(tmp_path / "results" / "chemogpcrs" / "chemogpcrs_9606.fa")

    assert _run_discovery(STAGE02A, tmp_path) == str(berghia)


def test_02a_does_not_glob_the_chemogpcrs_directory():
    """A wildcard probe would silently re-introduce the multi-species pooling."""
    src = STAGE02A.read_text()
    loop = _extract_for_loop(STAGE02A, "candidate_file")
    assert "chemogpcrs_*" not in loop, "02a must not glob per-species stage-02 outputs"
    assert "chemogpcrs_${" not in loop, "02a must not interpolate a taxid into its probe list"
    assert "cat " not in loop, "02a must not concatenate per-species files"
    assert src.count("chemogpcrs_berghia.fa") >= 1


# --------------------------------------------------------------------------
# F8 related: the symlink target directory is never created
# --------------------------------------------------------------------------
def test_02a_creates_the_candidates_directory_it_symlinks_into():
    """02a only ever mkdir -p'd clustering/; the ln -sf target dir was assumed."""
    lines = STAGE02A.read_text().splitlines()
    mkdir_at = next(
        (i for i, l in enumerate(lines) if 'mkdir -p "${CANDIDATES_DIR}"' in l and not l.lstrip().startswith("#")),
        None,
    )
    assert mkdir_at is not None, "02a symlinks into ${CANDIDATES_DIR} but never creates it"
    ln_at = next(i for i, l in enumerate(lines) if l.startswith("ln -sf"))
    assert mkdir_at < ln_at, "the candidates dir must be created BEFORE the symlink"


def test_02a_symlink_failure_is_fatal():
    """`ln -sf` had no `|| exit`, so the next line logged success unconditionally."""
    src = STAGE02A.read_text()
    ln_line_idx = next(i for i, l in enumerate(src.splitlines()) if l.startswith("ln -sf"))
    ln_block = "\n".join(src.splitlines()[ln_line_idx : ln_line_idx + 2])
    assert "||" in ln_block, "a failed `ln -sf` must not be reported as success"
    assert "exit 1" in ln_block


def test_02a_symlink_block_succeeds_on_a_clean_tree(tmp_path):
    """Behavioural: the real mkdir+ln pair produces a resolvable symlink."""
    results = tmp_path / "results"
    target = _write_fasta(results / "clustering" / "candidates_nr098.fa")
    src_lines = STAGE02A.read_text().splitlines()
    mkdir_line = next(l for l in src_lines if 'mkdir -p "${CANDIDATES_DIR}"' in l)
    ln_idx = next(i for i, l in enumerate(src_lines) if l.startswith("ln -sf"))

    snippet = "\n".join(
        [
            "log() { :; }",
            f'RESULTS_DIR="{results}"',
            'CANDIDATES_DIR="${RESULTS_DIR}/candidates"',
            f'OUTPUT_FASTA="{target}"',
            mkdir_line,
            'CANONICAL_OUTPUT="${CANDIDATES_DIR}/candidates_clustered.fa"',
            *src_lines[ln_idx : ln_idx + 2],
        ]
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)

    assert proc.returncode == 0, proc.stdout + proc.stderr
    link = results / "candidates" / "candidates_clustered.fa"
    assert link.is_symlink()
    assert link.resolve() == target.resolve()


# --------------------------------------------------------------------------
# F12.5: 02b must derive the clustered basename from CDHIT_IDENTITY
# --------------------------------------------------------------------------
def test_02b_finds_02a_output_at_the_default_identity(tmp_path):
    """Regression guard: the 0.98 default must keep resolving."""
    out = _write_fasta(tmp_path / "results" / "clustering" / "candidates_nr098.fa")

    resolved = _run_discovery(STAGE02B, tmp_path, CDHIT_IDENTITY="0.98")

    assert resolved == str(out)


def test_02b_finds_02a_output_at_a_non_default_identity(tmp_path):
    """The bug: 02b hardcoded nr098 and missed 02a's output at any other threshold."""
    out = _write_fasta(tmp_path / "results" / "clustering" / "candidates_nr095.fa")

    resolved = _run_discovery(STAGE02B, tmp_path, CDHIT_IDENTITY="0.95")

    assert resolved == str(out), "02b must compute the basename from ${CDHIT_IDENTITY}"


def test_02b_does_not_pick_up_a_stale_default_named_file(tmp_path):
    """At CDHIT_IDENTITY=0.90, a leftover nr098 file from a previous run is NOT ours."""
    _write_fasta(tmp_path / "results" / "clustering" / "candidates_nr098.fa")

    assert _run_discovery(STAGE02B, tmp_path, CDHIT_IDENTITY="0.90") == ""


def test_02b_uses_the_same_basename_expression_as_02a():
    """Guards against the two stages drifting apart again."""
    expr = 'candidates_nr${CDHIT_IDENTITY/./}.fa'
    assert expr in STAGE02A.read_text(), "02a's output basename expression moved"
    assert expr in STAGE02B.read_text(), "02b must reuse 02a's basename expression"
    assert "candidates_nr098.fa" not in STAGE02B.read_text(), (
        "02b must not hardcode the 0.98 default"
    )


def test_02b_still_prefers_the_canonical_symlink(tmp_path):
    """02a's canonical symlink stays the first choice, identity-independent."""
    link_target = _write_fasta(tmp_path / "results" / "clustering" / "candidates_nr095.fa")
    canonical = tmp_path / "results" / "candidates" / "candidates_clustered.fa"
    canonical.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(link_target, canonical)

    resolved = _run_discovery(STAGE02B, tmp_path, CDHIT_IDENTITY="0.95")

    assert resolved == str(canonical)
