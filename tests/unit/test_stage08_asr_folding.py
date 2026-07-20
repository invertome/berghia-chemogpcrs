"""Tests for the ASR-sequence folding path in 08_structural_analysis.sh (bead 444).

The defect
----------
Stage 08 intends to fold ancestral-sequence-reconstruction (ASR) sequences
alongside the extant top-ranked candidates, but they were silently dropped::

    cat "${RESULTS_DIR}/asr/"*_asr.fa > .../asr_seqs.fa 2>/dev/null     # line 59
    grep "^>" .../asr_seqs.fa | sed 's/>//' >> .../top_ids.txt          # line 60
    ...
    seqtk subseq "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" top_ids.txt

``chemogpcrs_berghia.fa`` holds only EXTANT candidates. Ancestral
reconstructions are, by definition, absent from it, so ``seqtk subseq``
silently dropped every ASR id. The per-candidate loop's
``grep -A1 "^>${id}$"`` then wrote a zero-byte input FASTA and hit
``continue``. Net effect: ASR ids inflated the reported ``final_count``
while never being folded, and ``asr_seqs.fa`` was written and never read
again.

The fix
-------
Build a combined extraction source (extant candidates + ``asr_seqs.fa``) and
run ``seqtk subseq`` against that, so ids from both sources resolve. Because a
duplicate id in the seqtk source is a silent correctness problem downstream
(seqtk would emit two records for one id, and the per-candidate
``grep -A1`` would fold whichever came first), the combined source is guarded
with the repo's existing ``assert_no_duplicate_fasta_ids`` fail-loud helper.

Testing approach
----------------
``08_structural_analysis.sh`` is a stage script, not a library: sourcing it
executes the whole pipeline stage. So -- following the extraction pattern
already used by tests/unit/test_channel_merge_wiring.py (which pulls a real
function out of rank_candidates.py rather than reimplementing it) -- these
tests pull the real shell functions out of the stage's source text and source
them in isolation alongside the real functions.sh. The logic under test is
therefore the shipping code, not a copy of it.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE08 = PROJECT_ROOT / "08_structural_analysis.sh"

# The stage-08 shell functions under test.
COLLECT_FN = "collect_asr_sequences"
BUILD_FN = "build_structural_seq_source"
COUNT_FN = "count_foldable_sequences"


def _extract_shell_function(name: str, src: str) -> str:
    """Exact source text of one top-level ``name() { ... }`` in ``src``.

    Anchored on a column-0 ``}`` terminator, which is the stage's style for
    every function it defines. Raises if the function is absent so a missing
    implementation fails as a clear error rather than a confusing shell syntax
    error further down.
    """
    match = re.search(
        rf"^{re.escape(name)}\(\)\s*\{{$.*?^\}}$",
        src,
        flags=re.MULTILINE | re.DOTALL,
    )
    if not match:
        raise AssertionError(
            f"{STAGE08.name} does not define a top-level shell function '{name}()'"
        )
    return match.group(0)


def _run_shell(body: str, tmp_path: Path, *fn_names: str):
    """Source functions.sh + the named stage-08 functions, then run ``body``.

    Sourcing functions.sh installs an EXIT trap (finalize_pipeline); it is
    disarmed with ``trap - EXIT`` so it cannot pollute the test's exit status,
    matching tests/unit/test_assert_fasta_non_empty.py and
    tests/unit/test_candidate_source_selection.sh.
    """
    src = STAGE08.read_text()
    fns = "\n".join(_extract_shell_function(n, src) for n in (fn_names or (COLLECT_FN, BUILD_FN, COUNT_FN)))
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    script = f'source "{FUNCTIONS_SH}"\ntrap - EXIT\n{fns}\n{body}\n'
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _ids_in(fasta: Path) -> list[str]:
    """Bare record ids (header up to first whitespace) in file order."""
    return [
        line[1:].split()[0]
        for line in fasta.read_text().splitlines()
        if line.startswith(">")
    ]


def _seq_of(fasta: Path, wanted: str) -> str | None:
    """Concatenated sequence for ``wanted``, or None if the id is absent.

    A tiny stand-in for ``seqtk subseq`` so these tests stay hermetic (seqtk is
    a cluster binary and is not assumed present locally).
    """
    seq, capturing = [], False
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            capturing = line[1:].split()[0] == wanted
        elif capturing:
            seq.append(line.strip())
    return "".join(seq) if seq else None


# --------------------------------------------------------------------------
# collect_asr_sequences: gathering stage 05's *_asr.fa robustly
# --------------------------------------------------------------------------

def test_collect_asr_gathers_every_asr_file(tmp_path):
    asr_dir = tmp_path / "asr"
    asr_dir.mkdir()
    (asr_dir / "OG0000001_Node5_asr.fa").write_text(">OG0000001_Node5\nMKFAA\n")
    (asr_dir / "OG0000002_Node9_asr.fa").write_text(">OG0000002_Node9\nMKVCC\n")
    (asr_dir / "OG0000001_asr_plot.png").write_text("not a fasta\n")

    out = tmp_path / "asr_seqs.fa"
    result = _run_shell(f'{COLLECT_FN} "{asr_dir}" "{out}"', tmp_path, COLLECT_FN)
    assert result.returncode == 0, result.stderr

    assert sorted(_ids_in(out)) == ["OG0000001_Node5", "OG0000002_Node9"]


def test_collect_asr_with_no_asr_directory_yields_empty_file(tmp_path):
    """Stage 05 ASR is optional -- an absent results/asr/ must not be fatal."""
    out = tmp_path / "asr_seqs.fa"
    result = _run_shell(f'{COLLECT_FN} "{tmp_path}/asr" "{out}"', tmp_path, COLLECT_FN)
    assert result.returncode == 0, result.stderr
    assert out.exists()
    assert out.read_text() == ""


def test_collect_asr_never_writes_a_literal_unexpanded_glob(tmp_path):
    """The old `cat .../*_asr.fa 2>/dev/null` idiom left a literal glob path.

    An empty asr/ directory leaves the glob unmatched. A bare glob would be
    passed through literally; the collector must yield an empty file with no
    '*' in it.
    """
    asr_dir = tmp_path / "asr"
    asr_dir.mkdir()
    out = tmp_path / "asr_seqs.fa"
    result = _run_shell(f'{COLLECT_FN} "{asr_dir}" "{out}"', tmp_path, COLLECT_FN)
    assert result.returncode == 0, result.stderr
    assert "*" not in out.read_text()
    assert out.read_text() == ""


def test_collect_asr_truncates_a_stale_previous_run(tmp_path):
    """Re-running with no ASR results must not leave the prior run's records."""
    asr_dir = tmp_path / "asr"
    asr_dir.mkdir()
    out = tmp_path / "asr_seqs.fa"
    out.write_text(">StaleNode\nMKFAA\n")
    result = _run_shell(f'{COLLECT_FN} "{asr_dir}" "{out}"', tmp_path, COLLECT_FN)
    assert result.returncode == 0, result.stderr
    assert out.read_text() == ""


# --------------------------------------------------------------------------
# build_structural_seq_source: the actual bug fix
# --------------------------------------------------------------------------

def test_asr_ids_resolve_in_the_combined_extraction_source(tmp_path):
    """THE regression test for the defect: ASR ids must be extractable."""
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA\n>cand2\nMKVCC\n")
    asr = tmp_path / "asr_seqs.fa"
    asr.write_text(">OG0000001_Node5\nMWWWW\n")
    combined = tmp_path / "extraction_source.fa"

    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{asr}" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 0, result.stderr

    # Every id stage 08 puts in top_ids.txt -- extant AND ancestral -- resolves.
    assert _ids_in(combined) == ["cand1", "cand2", "OG0000001_Node5"]
    # ...and resolves to the correct sequence, not a truncated/merged one.
    assert _seq_of(combined, "OG0000001_Node5") == "MWWWW"
    assert _seq_of(combined, "cand2") == "MKVCC"


def test_extant_fasta_without_trailing_newline_does_not_merge_records(tmp_path):
    """Naive `cat extant asr > combined` would fuse the last extant sequence
    line onto the first ASR header, corrupting BOTH records silently."""
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA")  # no trailing newline
    asr = tmp_path / "asr_seqs.fa"
    asr.write_text(">Node5\nMWWWW\n")
    combined = tmp_path / "extraction_source.fa"

    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{asr}" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 0, result.stderr
    assert _ids_in(combined) == ["cand1", "Node5"]
    assert _seq_of(combined, "cand1") == "MKFAA"
    assert _seq_of(combined, "Node5") == "MWWWW"


def test_no_asr_case_yields_just_the_extant_candidates(tmp_path):
    """ASR is optional: stage 08 must work when stage 05 produced no ASR."""
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA\n>cand2\nMKVCC\n")
    asr = tmp_path / "asr_seqs.fa"
    asr.write_text("")  # what collect_asr_sequences leaves when ASR was not run
    combined = tmp_path / "extraction_source.fa"

    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{asr}" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 0, result.stderr
    assert _ids_in(combined) == ["cand1", "cand2"]


def test_missing_asr_file_is_not_a_hard_dependency(tmp_path):
    """No hard dependency on stage 05: an absent asr_seqs.fa is fine."""
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA\n")
    combined = tmp_path / "extraction_source.fa"

    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{tmp_path}/nope.fa" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 0, result.stderr
    assert _ids_in(combined) == ["cand1"]


def test_headerless_asr_file_is_ignored(tmp_path):
    """Defensive: a non-empty but record-less asr_seqs.fa (e.g. a stray literal
    path written by the old redirection) must not corrupt the source."""
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA\n")
    asr = tmp_path / "asr_seqs.fa"
    asr.write_text("/results/asr/*_asr.fa\n")
    combined = tmp_path / "extraction_source.fa"

    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{asr}" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 0, result.stderr
    assert _ids_in(combined) == ["cand1"]
    assert "*_asr.fa" not in combined.read_text()


def test_duplicate_id_across_extant_and_asr_fails_loud(tmp_path):
    """A duplicate in the seqtk source is a silent correctness problem: seqtk
    emits both records and the per-candidate grep folds whichever came first."""
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA\n>Node5\nMKVCC\n")
    asr = tmp_path / "asr_seqs.fa"
    asr.write_text(">Node5\nMWWWW\n")
    combined = tmp_path / "extraction_source.fa"

    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{asr}" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 1
    assert "duplicate" in result.stderr.lower()
    assert "Node5" in result.stderr


def test_duplicate_id_between_two_asr_files_fails_loud(tmp_path):
    """Stage 05 writes ancestral headers as the bare node name, so two
    orthogroups can each emit a '>Node5'. That must fail loudly, not fold an
    arbitrary one of them."""
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.write_text(">cand1\nMKFAA\n")
    asr = tmp_path / "asr_seqs.fa"
    asr.write_text(">Node5\nMWWWW\n>Node5\nMYYYY\n")
    combined = tmp_path / "extraction_source.fa"

    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{asr}" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 1
    assert "duplicate" in result.stderr.lower()


def test_empty_extant_candidate_fasta_fails_loud(tmp_path):
    extant = tmp_path / "chemogpcrs_berghia.fa"
    extant.touch()
    combined = tmp_path / "extraction_source.fa"
    result = _run_shell(
        f'{BUILD_FN} "{extant}" "{tmp_path}/nope.fa" "{combined}"', tmp_path, BUILD_FN
    )
    assert result.returncode == 1
    assert result.stderr.strip()


def test_missing_extant_candidate_fasta_fails_loud(tmp_path):
    combined = tmp_path / "extraction_source.fa"
    result = _run_shell(
        f'{BUILD_FN} "{tmp_path}/gone.fa" "{tmp_path}/nope.fa" "{combined}"',
        tmp_path,
        BUILD_FN,
    )
    assert result.returncode == 1
    assert result.stderr.strip()


# --------------------------------------------------------------------------
# count_foldable_sequences: an honest final_count
# --------------------------------------------------------------------------

def test_reported_count_matches_what_was_actually_extracted(tmp_path):
    """final_count must reflect folded sequences, not requested ids."""
    top_ids = tmp_path / "top_ids.txt"
    top_ids.write_text("cand1\ncand2\nNode5\n")
    top_seqs = tmp_path / "top_seqs.fa"
    top_seqs.write_text(">cand1\nMKFAA\n>cand2\nMKVCC\n>Node5\nMWWWW\n")

    result = _run_shell(
        f'{COUNT_FN} "{top_seqs}" "{top_ids}"; echo "COUNT=${{FOLDABLE_COUNT}}"',
        tmp_path,
        COUNT_FN,
    )
    assert result.returncode == 0, result.stderr
    assert "COUNT=3" in result.stdout


def test_unresolved_ids_are_reported_not_silently_counted(tmp_path):
    """The original bug's symptom: ids inflating the count without folding."""
    top_ids = tmp_path / "top_ids.txt"
    top_ids.write_text("cand1\ncand2\nNode5\n")
    top_seqs = tmp_path / "top_seqs.fa"
    top_seqs.write_text(">cand1\nMKFAA\n>cand2\nMKVCC\n")  # Node5 did not resolve

    result = _run_shell(
        f'{COUNT_FN} "{top_seqs}" "{top_ids}"; echo "COUNT=${{FOLDABLE_COUNT}}"',
        tmp_path,
        COUNT_FN,
    )
    assert result.returncode == 0, result.stderr
    assert "COUNT=2" in result.stdout, result.stdout
    combined_output = result.stdout + result.stderr
    assert "WARN" in combined_output
    assert "3" in combined_output and "2" in combined_output


def test_count_handles_an_empty_extraction(tmp_path):
    top_ids = tmp_path / "top_ids.txt"
    top_ids.write_text("")
    top_seqs = tmp_path / "top_seqs.fa"
    top_seqs.touch()
    result = _run_shell(
        f'{COUNT_FN} "{top_seqs}" "{top_ids}"; echo "COUNT=${{FOLDABLE_COUNT}}"',
        tmp_path,
        COUNT_FN,
    )
    assert result.returncode == 0, result.stderr
    assert "COUNT=0" in result.stdout


# --------------------------------------------------------------------------
# Stage wiring: the functions must actually be used by the stage
# --------------------------------------------------------------------------

def test_stage08_extracts_from_the_combined_source_not_the_extant_fasta():
    """Regression guard for the defect itself."""
    text = STAGE08.read_text()
    seqtk_lines = [ln for ln in text.splitlines() if "subseq" in ln]
    assert seqtk_lines, "stage 08 must still extract sequences with seqtk subseq"
    for line in seqtk_lines:
        assert "chemogpcrs/chemogpcrs_berghia.fa" not in line, (
            "seqtk subseq must read the combined extant+ASR source, not the "
            f"extant-only candidate FASTA: {line.strip()}"
        )


def test_stage08_calls_each_helper():
    text = STAGE08.read_text()
    for fn in (COLLECT_FN, BUILD_FN, COUNT_FN):
        # Once for the definition, at least once for the call.
        assert text.count(fn) >= 2, f"stage 08 defines but never calls {fn}"


def test_stage08_guards_the_source_against_duplicate_ids():
    text = STAGE08.read_text()
    assert "assert_no_duplicate_fasta_ids" in text, (
        "the combined extraction source must be guarded with the repo's "
        "fail-loud duplicate-id helper"
    )


def test_stage08_feeds_asr_seqs_into_the_extraction_source():
    """asr_seqs.fa was written and never read again. It must now be an argument
    to the extraction-source builder, not merely mentioned somewhere."""
    text = STAGE08.read_text()
    # The builder invocation (the call, not the definition), including any
    # backslash-continued argument lines.
    call = re.search(
        rf"^{BUILD_FN} +(?:\\\n|.)*",
        text,
        flags=re.MULTILINE,
    )
    assert call, f"stage 08 never calls {BUILD_FN}"
    assert "ASR_SEQS" in call.group(0), (
        f"the ASR FASTA must be passed to {BUILD_FN}; got: {call.group(0)!r}"
    )


def test_stage08_still_parses():
    result = subprocess.run(["bash", "-n", str(STAGE08)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
