"""Repo-wide audit of the bead-444 `grep -c ... || echo 0` idiom.

``grep -c`` prints ``0`` to stdout AND exits non-zero when nothing matches, so::

    n=$(grep -c '^>' "$file" 2>/dev/null || echo 0)

appends a SECOND ``0`` on the no-match path and ``$n`` becomes the two-line
string ``"0\\n0"``. Three distinct downstream failure modes, all reproduced::

    [ "$n" -gt 0 ]          -> bash: [: 0\n0: integer expected
    $(( total + n ))        -> bash: arithmetic syntax error (error token "0")
    printf '%5d' "$n"       -> bash: printf: 0\n0: invalid number
    printf '%s\\t%s\\n' ...  -> embeds a raw newline, splitting a TSV row in two

The no-match path is ordinary operation, not a corner case: any zero-length or
header-less FASTA produced upstream (a filter that dropped every sequence, an
empty orthogroup, a failed extraction) reaches it.

The fix is the established in-repo pattern from ``count_taxid_occurrences`` in
functions.sh: capture with ``|| true`` so the count is kept without the exit
status, then read it back through ``${VAR:-0}`` so a missing file still yields
exactly one integer.

These tests execute the REAL assignment lines lifted out of the production
scripts, so they fail for the actual reason rather than a re-typed copy.
"""
from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Sites whose captured value feeds an integer test, arithmetic, or a %d
# conversion -- i.e. where the doubled zero is a hard failure, not a cosmetic
# blemish. `guard` is the real downstream expression from that same script.
HARMFUL_SITES = [
    pytest.param(
        "01_reference_processing.sh", "n_cds", '[ "${n_cds}" -lt 5000 ]',
        id="01_reference_processing-n_cds",
    ),
    pytest.param(
        "01_reference_processing.sh", "EXISTING_N", '[ "${EXISTING_N:-0}" -lt 5000 ]',
        id="01_reference_processing-EXISTING_N",
    ),
    pytest.param(
        "scripts/build_curated_chemo_hmms.sh", "n_seqs", '[ "$n_seqs" -lt 3 ]',
        id="build_curated_chemo_hmms-n_seqs",
    ),
    pytest.param(
        "scripts/run_aligner.sh", "N", '[ "$N" -lt 200 ]',
        id="run_aligner-N",
    ),
    # Structured machine-readable output: a raw newline splits the row.
    pytest.param(
        "scripts/unity/build_classification_reference_trees.sh", "n",
        '[ "$(echo -e "name\\t${n}" | wc -l)" -eq 1 ]',
        id="build_classification_reference_trees-manifest-row",
    ),
    pytest.param(
        "scripts/unity/merge_recover_proteome.sh", "ng",
        '[ "$(printf \'%s\\t%s\\n\' es "$ng" | wc -l)" -eq 1 ]',
        id="merge_recover_proteome-report-row",
    ),
    pytest.param(
        "scripts/unity/smoke_test_filter_stack_large.sh", "n",
        'printf "%5d\\n" "$n" >/dev/null',
        id="smoke_test_filter_stack_large-printf-percent-d",
    ),
]


def _extract_assignment(script: Path, var: str) -> str:
    """Return the production line that assigns *var* from a `grep -c` capture.

    Matches the *last* such assignment so that an unrelated earlier `grep -c`
    binding of the same (often single-letter) variable cannot shadow the site
    under test.
    """
    pattern = re.compile(rf'^\s*(?:local\s+)?{re.escape(var)}="?\$\(\s*grep -c\b.*\)"?\s*$')
    matches = [line.strip() for line in script.read_text().splitlines() if pattern.match(line)]
    if not matches:
        raise AssertionError(f"no `grep -c` assignment to {var!r} found in {script}")
    return matches[-1]


def _retarget(line: str, target: Path) -> str:
    """Point the line's `grep -c` at *target* instead of its production path."""
    new, count = re.subn(
        r'(grep -c\s+\S+\s+)"[^"]*"',
        lambda m: m.group(1) + shlex.quote(str(target)),
        line,
        count=1,
    )
    assert count == 1, f"could not retarget the grep argument in: {line}"
    return new


@pytest.mark.parametrize("rel_path,var,guard", HARMFUL_SITES)
def test_no_match_value_survives_its_own_downstream_guard(rel_path, var, guard, tmp_path):
    """Run the real assignment against a header-less FASTA, then its own guard."""
    script = PROJECT_ROOT / rel_path
    fasta = tmp_path / "headerless.fa"
    fasta.write_text("ACGT\nTTTT\n")  # exists, but `grep -c '^>'` matches nothing

    # Every `guard` above is written so that it evaluates TRUE exactly when the
    # capture is a clean single "0" -- so GUARD_FAILED means the doubled zero
    # corrupted the value even where it did not raise a shell error.
    body = _retarget(_extract_assignment(script, var), fasta)
    result = subprocess.run(
        ["bash", "-c", f"{body}\nif {guard}; then echo GUARD_OK; else echo GUARD_FAILED; fi"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"{rel_path}:{var} -> {result.stderr}"
    assert "integer expected" not in result.stderr, f"{rel_path}:{var} -> {result.stderr}"
    assert "integer expression expected" not in result.stderr, f"{rel_path}:{var} -> {result.stderr}"
    assert "arithmetic syntax error" not in result.stderr, f"{rel_path}:{var} -> {result.stderr}"
    assert "invalid number" not in result.stderr, f"{rel_path}:{var} -> {result.stderr}"
    assert result.stdout.strip() == "GUARD_OK", result.stdout


@pytest.mark.parametrize("rel_path,var,guard", HARMFUL_SITES)
def test_no_match_capture_is_exactly_one_zero(rel_path, var, guard, tmp_path):
    script = PROJECT_ROOT / rel_path
    fasta = tmp_path / "headerless.fa"
    fasta.write_text("ACGT\nTTTT\n")

    body = _retarget(_extract_assignment(script, var), fasta)
    result = subprocess.run(
        ["bash", "-c", f'{body}\nprintf "%s" "${{{var}:-0}}"'],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == "0", (
        f"{rel_path}:{var} captured {result.stdout!r}; expected a single '0'"
    )


# --- repo-wide regression guard ----------------------------------------------
#
# One shared assertion for every remaining (cosmetic) site rather than twenty
# near-identical tests. Cosmetic today is harmful the first time someone
# compares the value, so the idiom is banned outright.

IDIOM = re.compile(r"grep -c\b[^\n]*\|\|\s*echo\s+0")


# `preliminary/` is a gitignored legacy tree (a stale fork of the pipeline kept
# on the workstation only); it is not part of the tracked codebase and is out of
# scope for this guard.
EXCLUDED_DIRS = {".git", "reports", "preliminary"}


def _shell_files() -> list[Path]:
    return sorted(
        p for p in PROJECT_ROOT.rglob("*.sh") if EXCLUDED_DIRS.isdisjoint(p.parts)
    )


def test_no_shell_file_uses_the_grep_c_echo_zero_idiom():
    offenders = []
    for path in _shell_files():
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue  # the explanatory comments describing the defect
            if IDIOM.search(line):
                offenders.append(f"{path.relative_to(PROJECT_ROOT)}:{lineno}: {stripped}")
    assert not offenders, (
        "the `grep -c ... || echo 0` idiom doubles the zero on the no-match path; "
        "use `|| true` plus a `${VAR:-0}` read-back instead:\n  "
        + "\n  ".join(offenders)
    )


# `grep -c ... | tail -1` is the symptom-patch for the same defect: piping the
# count through `tail` discards the second "0", but a pipeline's exit status is
# the LAST command's, so `tail` also masks grep's failure. The `|| ...` fallback
# becomes dead code and a MISSING file now yields an empty string instead of 0
# -- which dies in an integer test exactly like the doubled zero did.
SYMPTOM_PATCH = re.compile(r"grep -c\b[^\n|]*\|\s*tail\b")


def test_no_shell_file_swallows_the_doubled_zero_with_tail():
    offenders = []
    for path in _shell_files():
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            if SYMPTOM_PATCH.search(line):
                offenders.append(f"{path.relative_to(PROJECT_ROOT)}:{lineno}: {stripped}")
    assert not offenders, (
        "`grep -c ... | tail -1` masks grep's exit status, so the `||` fallback "
        "never fires and a missing file captures an empty string; use "
        "`|| true` plus a `${VAR:-0}` read-back instead:\n  "
        + "\n  ".join(offenders)
    )


def test_every_touched_shell_file_still_parses():
    failures = []
    for path in _shell_files():
        result = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
        if result.returncode != 0:
            failures.append(f"{path.relative_to(PROJECT_ROOT)}: {result.stderr.strip()}")
    assert not failures, "\n".join(failures)
