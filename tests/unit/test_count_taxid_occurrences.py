"""Tests for count_taxid_occurrences in functions.sh (bead 444, item 1).

Stage 05 decided whether an orthogroup contains Berghia with::

    berghia_count=$(echo "$taxids" | grep -c "${BERGHIA_TAXID}" || echo 0)
    ...
    if [ "$berghia_count" -gt 0 ] && [ "$seq_count" -gt 2 ]; then

Two independent defects, both confirmed empirically:

(a) ``grep -c`` prints ``0`` *and* exits 1 when nothing matches, so the
    ``|| echo 0`` fallback appends a SECOND ``0``. The captured value is the
    two-line string ``"0\\n0"``, and ``[ "0\\n0" -gt 0 ]`` aborts with
    ``bash: [: 0\\n0: integer expected`` instead of evaluating cleanly.

(b) ``grep`` matches substrings, so BERGHIA_TAXID=1287507 is "found" inside an
    unrelated taxid such as 12875070 — a silent false positive that would run
    Berghia-LSE ASR on an orthogroup with no Berghia in it.

``count_taxid_occurrences`` replaces the idiom: it splits the whitespace-
separated taxid list into one field per line and counts only whole-field
(``-x``, literal ``-F``) matches, always emitting exactly one integer line.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE05 = PROJECT_ROOT / "05_selective_pressure_and_asr.sh"


def _run_bash(body: str, tmp_path: Path) -> subprocess.CompletedProcess:
    """Source functions.sh and run *body*.

    Sourcing functions.sh installs an EXIT trap (finalize_pipeline) that logs to
    stdout; `trap - EXIT` disarms it so only the helper's output is captured
    (same technique as tests/unit/test_candidate_source_selection.sh).
    """
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    script = f'source "{FUNCTIONS_SH}"; trap - EXIT; {body}'
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _count(taxid_list: str, taxid: str, tmp_path: Path) -> subprocess.CompletedProcess:
    return _run_bash(f'count_taxid_occurrences "{taxid_list}" "{taxid}"', tmp_path)


def test_counts_exact_match(tmp_path):
    result = _count("9606 1287507 10090", "1287507", tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


def test_no_match_emits_single_zero_line(tmp_path):
    # Defect (a): the old idiom emitted "0\n0" here.
    result = _count("9606 10090", "1287507", tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["0"], f"expected one '0' line, got {result.stdout!r}"


def test_empty_taxid_list_emits_single_zero_line(tmp_path):
    result = _count("", "1287507", tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["0"], f"expected one '0' line, got {result.stdout!r}"


def test_substring_taxid_is_not_a_match(tmp_path):
    # Defect (b): 12875070 must NOT count as 1287507.
    result = _count("9606 12875070 10090", "1287507", tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "0"


def test_taxid_embedded_as_prefix_is_not_a_match(tmp_path):
    # The reverse direction: 1287507 must not be found inside 51287507.
    result = _count("51287507", "1287507", tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "0"


def test_newline_separated_list_also_works(tmp_path):
    # get_taxids_from_fasta's fallback path can emit newline-separated ids.
    result = _count("9606\n1287507\n10090", "1287507", tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


def test_counts_repeats(tmp_path):
    result = _count("1287507 9606 1287507", "1287507", tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "2"


def test_result_is_usable_in_an_integer_test(tmp_path):
    """The whole point: the value must survive `[ "$n" -gt 0 ]` without error."""
    result = _run_bash(
        'n=$(count_taxid_occurrences "9606 10090" "1287507"); '
        'if [ "$n" -gt 0 ]; then echo YES; else echo NO; fi',
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "NO"
    assert "integer expression expected" not in result.stderr
    assert "integer expected" not in result.stderr


def test_stage05_uses_the_helper_not_the_buggy_grep():
    """Regression guard: the raw `grep -c "$BERGHIA_TAXID"` idiom is gone."""
    text = STAGE05.read_text()
    assert "count_taxid_occurrences" in text, \
        "stage 05 must use count_taxid_occurrences for the Berghia membership test"
    assert 'grep -c "${BERGHIA_TAXID}"' not in text, \
        "stage 05 still uses the substring-matching `grep -c` idiom"


def test_stage05_still_parses():
    result = subprocess.run(["bash", "-n", str(STAGE05)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
