"""Tests for stage 09 displaying ``rank_tier`` instead of ``confidence_tier``.

The pipeline carries TWO tier systems in the ranked CSV:

* ``confidence_tier`` --- a legacy heuristic in ``rank_candidates.py``
  (``assign_confidence_tier``) that thresholds ``evidence_completeness`` and a
  ``rank_score``/max ratio. Its 0.5/0.3 cutoffs are uncalibrated against the
  current score scale, and it derives from the hand-weighted composite that the
  label-free RRA ordering (``final_rank``) superseded.
* ``rank_tier`` --- produced by ``scripts/rank_confidence.py``, which resamples
  the SIGNAL SET with replacement, re-runs RRA on each of 1000 draws, and bins
  ``p_top_k`` (the bootstrap probability of being in the top k) into
  ``high`` / ``plausible`` / ``tail``. Lowercase, unlike ``confidence_tier``.

``rank_tier`` was computed and joined but read by nothing. These tests pin the
switch: stage 09's summary counts, its per-row badge, and its prose all key on
``rank_tier`` and describe what it MEANS (a bootstrap top-k probability), not
evidence completeness.

The critical case is the column being ABSENT: ``rank_confidence.py`` is invoked
non-fatally by stage 07 (``|| log --level=WARN``), so a run where it failed
yields a ranked CSV with no ``rank_tier`` column. The report must then say so,
rather than reporting three zeros (which reads as a real result) or painting
every candidate with the worst tier.

Both the shell summary block and the python row renderer are embedded in the
stage, so these tests EXTRACT and EXECUTE them against synthetic CSVs -- the
same pattern as ``test_stage09_poscontrol_status.py``.
"""
from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
STAGE_09 = REPO_ROOT / "09_report_generation.sh"
SUMMARY_LIB = REPO_ROOT / "scripts" / "_report_summary_lib.sh"

# The tier vocabulary rank_confidence.py actually emits. Asserted against the
# producer below so this fixture cannot drift from it.
TIER_HIGH = "high"
TIER_PLAUSIBLE = "plausible"
TIER_TAIL = "tail"
TIERS = (TIER_HIGH, TIER_PLAUSIBLE, TIER_TAIL)

# A minimal ranked-CSV schema: every column the two renderers read.
RANKED_COLUMNS = [
    "id", "rank_score", "confidence_tier", "rank_tier", "p_top_k",
    "evidence_completeness", "selection_significant", "busted_mh_significant",
    "hcr_probe_friendly", "tandem_cluster_size",
]


# --- Schema guard: the fixture vocabulary must match the producer ------------

def test_tier_vocabulary_matches_rank_confidence() -> None:
    """The tier strings and cutoffs under test must be the ones assign_tiers emits."""
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    try:
        import rank_confidence as rc
    finally:
        sys.path.pop(0)

    # Column name the augmenter joins.
    assert "rank_tier" in rc._JOIN_COLUMNS

    # Exact values + cutoffs, read off the real function rather than asserted
    # from a summary: 0.8 -> high, 0.2 -> plausible, below -> tail.
    intervals = {
        "a": {"p_top_k": 0.95}, "b": {"p_top_k": 0.80},
        "c": {"p_top_k": 0.79}, "d": {"p_top_k": 0.20},
        "e": {"p_top_k": 0.19}, "f": {"p_top_k": 0.0},
    }
    tiers = rc.assign_tiers(intervals)
    assert tiers == {
        "a": TIER_HIGH, "b": TIER_HIGH,
        "c": TIER_PLAUSIBLE, "d": TIER_PLAUSIBLE,
        "e": TIER_TAIL, "f": TIER_TAIL,
    }
    assert set(tiers.values()) <= set(TIERS)


# --- Harness: the python row renderer ---------------------------------------

def _extract_ranked_renderer() -> str:
    """Return the body of the ``$RANKED_FILE`` python heredoc in stage 09."""
    text = STAGE_09.read_text()
    m = re.search(
        r"^\s*python3 - \"\$RANKED_FILE\".*?<<'PY'\n(.*?)^PY$",
        text, flags=re.S | re.M,
    )
    assert m, "could not locate the ranked-candidates renderer heredoc in stage 09"
    return m.group(1)


def _ranked_row(**overrides: str) -> dict:
    row = {c: "" for c in RANKED_COLUMNS}
    row.update({"id": "cand1", "rank_score": "0.5", "tandem_cluster_size": "1"})
    row.update(overrides)
    return row


def _write_ranked_csv(path: Path, rows: list[dict], columns=None) -> Path:
    cols = list(columns if columns is not None else RANKED_COLUMNS)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})
    return path


def _render_rows(tmp_path: Path, rows: list[dict], columns=None) -> list[str]:
    """Run the extracted row renderer over a synthetic ranked CSV."""
    csv_path = _write_ranked_csv(tmp_path / "ranked.csv", rows, columns)
    script = tmp_path / "renderer.py"
    script.write_text(_extract_ranked_renderer())
    proc = subprocess.run([sys.executable, str(script), str(csv_path)],
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, f"renderer failed:\n{proc.stderr}"
    return [ln for ln in proc.stdout.splitlines() if ln.strip()]


def _badge(line: str) -> str:
    m = re.search(r"\\textcolor\{[^}]*\}\{[^}]*\}", line)
    assert m, f"row has no colour badge: {line}"
    return m.group(0)


# --- Harness: the shell summary block ---------------------------------------

def _extract_summary_block() -> str:
    """Return the stats block guarded by ``if [ -f "$RANKED_FILE" ]``.

    Two such guards exist in the stage; the leftmost match is the stats block
    (the other wraps the top-20 table), asserted by its use of the header-aware
    helpers from _report_summary_lib.sh.
    """
    text = STAGE_09.read_text()
    # Start at the variable initialisation so the harness inherits the stage's
    # own defaults (an absent-column run must be flagged by the STAGE, not by a
    # default the test supplied).
    m = re.search(r'^SIG_SELECTION="".*?^if \[ -f "\$RANKED_FILE" \]; then\n.*?^fi$',
                  text, flags=re.S | re.M)
    assert m, "could not locate the stats block in stage 09"
    block = m.group(0)
    assert "count_by_value" in block, "matched the wrong $RANKED_FILE block"
    return block


def _extract_tier_rows_block() -> str:
    """Return the conditional LaTeX emission for the summary tier rows."""
    text = STAGE_09.read_text()
    m = re.search(r'^if \[ "\$\{?RANK_TIER_AVAILABLE.*?^fi$',
                  text, flags=re.S | re.M)
    assert m, "could not locate the RANK_TIER_AVAILABLE emission block in stage 09"
    return m.group(0)


def _run_summary(tmp_path: Path, rows: list[dict], columns=None) -> dict:
    """Execute the stats block against a synthetic CSV; return the tier vars."""
    csv_path = _write_ranked_csv(tmp_path / "ranked.csv", rows, columns)
    script = tmp_path / "summary.sh"
    script.write_text(
        "#!/usr/bin/env bash\n"
        f'source "{SUMMARY_LIB}"\n'
        f'RANKED_FILE="{csv_path}"\n'
        + _extract_summary_block() +
        "\n"
        'echo "RANK_TIER_AVAILABLE=${RANK_TIER_AVAILABLE}"\n'
        'echo "RANK_TIER_HIGH=${RANK_TIER_HIGH}"\n'
        'echo "RANK_TIER_PLAUSIBLE=${RANK_TIER_PLAUSIBLE}"\n'
        'echo "RANK_TIER_TAIL=${RANK_TIER_TAIL}"\n'
        'echo "RANK_TIER_NA=${RANK_TIER_NA}"\n'
        'echo "TOP_N_RANKED=${TOP_N_RANKED}"\n'
    )
    proc = subprocess.run(["bash", str(script)], capture_output=True,
                          text=True, timeout=60)
    assert proc.returncode == 0, f"summary block failed:\n{proc.stderr}"
    out = {}
    for line in proc.stdout.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            out[k] = v
    return out


def _emit_tier_rows(tmp_path: Path, rows: list[dict], columns=None) -> list[str]:
    """Run stats block + the LaTeX tier-row emission; return the emitted lines."""
    csv_path = _write_ranked_csv(tmp_path / "ranked.csv", rows, columns)
    (tmp_path / "report").mkdir(exist_ok=True)
    script = tmp_path / "emit.sh"
    script.write_text(
        "#!/usr/bin/env bash\n"
        f'source "{SUMMARY_LIB}"\n'
        f'RANKED_FILE="{csv_path}"\n'
        f'RESULTS_DIR="{tmp_path}"\n'
        + _extract_summary_block()
        + "\n"
        + _extract_tier_rows_block()
        + "\n"
    )
    proc = subprocess.run(["bash", str(script)], capture_output=True,
                          text=True, timeout=60)
    assert proc.returncode == 0, f"tier-row emission failed:\n{proc.stderr}"
    tex = (tmp_path / "report" / "report.tex")
    assert tex.exists(), "emission block wrote no report.tex"
    return [ln for ln in tex.read_text().splitlines() if ln.strip()]


# --- Prose helpers ----------------------------------------------------------

def _tier_prose() -> str:
    """Return the tier-interpretation subsection as the stage emits it."""
    text = STAGE_09.read_text()
    # Ends at the heredoc terminator, so this is exactly the block the stage
    # emits -- no trailing shell scaffolding leaks into the LaTeX under test.
    m = re.search(r"(\\subsection\{Rank Confidence Tiers\}.*?)^EOF$",
                  text, flags=re.S | re.M)
    assert m, "could not locate the rank-tier interpretation prose in stage 09"
    return m.group(1)


def _rendered_tier_prose() -> str:
    """The tier prose AS LATEX SEES IT.

    The stage emits this from an UNQUOTED heredoc so the configured top-k is
    interpolated, which means every literal LaTeX ``$`` is written ``\\$``.
    Resolve the shell layer -- expand ``${VAR:-default}`` / ``${VAR}`` and
    unescape ``\\$`` -- so the remaining ``$`` are exactly the math delimiters.
    """
    prose = _tier_prose()
    prose = re.sub(r"\$\{[A-Z_]+:-([^}]*)\}", r"\1", prose)
    prose = re.sub(r"\$\{[A-Z_]+\}", "20", prose)
    assert not re.search(r"(?<!\\)\$\{", prose), (
        f"unresolved shell expansion in tier prose: {prose}")
    return prose.replace(r"\$", "$")


# ============================================================================
# The three tiers render distinctly
# ============================================================================

def test_high_tier_renders_green(tmp_path: Path) -> None:
    out = _render_rows(tmp_path, [_ranked_row(rank_tier=TIER_HIGH)])
    assert len(out) == 1
    assert r"\textcolor{highconf}{high}" in out[0]


def test_plausible_tier_renders_amber(tmp_path: Path) -> None:
    out = _render_rows(tmp_path, [_ranked_row(rank_tier=TIER_PLAUSIBLE)])
    assert r"\textcolor{medconf}{plausible}" in out[0]


def test_tail_tier_renders_neutral_not_red(tmp_path: Path) -> None:
    """``tail`` is a statistical statement, not a data-quality failure.

    Only k candidates can occupy the top k, so most of the list is necessarily
    tail. Painting the bulk of the table red (lowconf) would be alarmist and
    would misrepresent a correct, expected outcome as a problem.
    """
    out = _render_rows(tmp_path, [_ranked_row(rank_tier=TIER_TAIL)])
    assert r"\textcolor{neutral}{tail}" in out[0]
    assert "lowconf" not in out[0]


def test_three_tiers_are_visually_distinct(tmp_path: Path) -> None:
    out = _render_rows(tmp_path, [
        _ranked_row(id="a", rank_tier=TIER_HIGH),
        _ranked_row(id="b", rank_tier=TIER_PLAUSIBLE),
        _ranked_row(id="c", rank_tier=TIER_TAIL),
    ])
    assert len(out) == 3
    assert len({_badge(ln) for ln in out}) == 3


def test_tier_badge_colours_are_defined_in_preamble(tmp_path: Path) -> None:
    """Every colour a badge uses must actually be \\definecolor'd."""
    text = STAGE_09.read_text()
    out = _render_rows(tmp_path, [
        _ranked_row(id="a", rank_tier=TIER_HIGH),
        _ranked_row(id="b", rank_tier=TIER_PLAUSIBLE),
        _ranked_row(id="c", rank_tier=TIER_TAIL),
        _ranked_row(id="d", rank_tier=""),
    ])
    for line in out:
        colour = re.search(r"\\textcolor\{([^}]*)\}", line).group(1)
        assert re.search(rf"\\definecolor\{{{re.escape(colour)}}}", text), (
            f"badge colour '{colour}' is never defined in the preamble")


# ============================================================================
# The absent / blank column
# ============================================================================

def test_blank_tier_renders_not_available(tmp_path: Path) -> None:
    """A blank rank_tier means "no bootstrap statement", not "worst tier".

    rank_confidence.py leaves the cell empty for candidates covered by no
    signal (annotate_ranked_csv fillna("")), so this state is reachable even on
    a successful run.
    """
    out = _render_rows(tmp_path, [_ranked_row(rank_tier="")])
    assert "n/a" in out[0]
    assert "tail" not in out[0]
    assert "lowconf" not in out[0]


def test_absent_column_renders_not_available(tmp_path: Path) -> None:
    """No rank_tier column at all (rank_confidence.py failed) must not crash."""
    cols = [c for c in RANKED_COLUMNS if c != "rank_tier"]
    out = _render_rows(tmp_path, [_ranked_row(id="a")], columns=cols)
    assert len(out) == 1
    assert "n/a" in out[0]
    assert "tail" not in out[0]


def test_unavailable_is_distinct_from_all_three_tiers(tmp_path: Path) -> None:
    """The unavailable badge must not collide with any real tier badge."""
    out = _render_rows(tmp_path, [
        _ranked_row(id="a", rank_tier=TIER_HIGH),
        _ranked_row(id="b", rank_tier=TIER_PLAUSIBLE),
        _ranked_row(id="c", rank_tier=TIER_TAIL),
        _ranked_row(id="d", rank_tier=""),
    ])
    assert len({_badge(ln) for ln in out}) == 4


def test_unrecognised_tier_does_not_claim_a_tier(tmp_path: Path) -> None:
    """A future/unknown tier string must not be silently binned."""
    out = _render_rows(tmp_path, [_ranked_row(rank_tier="some_future_tier")])
    assert "n/a" in out[0]
    assert "highconf" not in out[0]


def test_summary_reports_unavailable_when_column_absent(tmp_path: Path) -> None:
    """Absent column -> flagged unavailable, NOT three zero counts."""
    cols = [c for c in RANKED_COLUMNS if c != "rank_tier"]
    got = _run_summary(tmp_path, [_ranked_row(id=f"c{i}") for i in range(5)],
                       columns=cols)
    assert got["RANK_TIER_AVAILABLE"] == "0"


def test_summary_table_says_not_computed_when_column_absent(tmp_path: Path) -> None:
    """The emitted LaTeX must state the statistic was not computed.

    Three ``0`` counts would read as a genuine finding ("no candidate is high
    confidence") rather than "this run never computed rank confidence".
    """
    cols = [c for c in RANKED_COLUMNS if c != "rank_tier"]
    lines = _emit_tier_rows(tmp_path, [_ranked_row(id=f"c{i}") for i in range(5)],
                            columns=cols)
    blob = "\n".join(lines)
    assert "not computed" in blob.lower()
    # No tier is claimed at all, in either direction.
    for tier in TIERS:
        assert not re.search(rf"\b{tier}\b.*\{{0\}}", blob), (
            f"absent column reported a '{tier}' count: {blob}")


# ============================================================================
# Counts
# ============================================================================

def test_summary_counts_each_tier(tmp_path: Path) -> None:
    rows = (
        [_ranked_row(id=f"h{i}", rank_tier=TIER_HIGH) for i in range(3)]
        + [_ranked_row(id=f"p{i}", rank_tier=TIER_PLAUSIBLE) for i in range(4)]
        + [_ranked_row(id=f"t{i}", rank_tier=TIER_TAIL) for i in range(5)]
        + [_ranked_row(id=f"n{i}", rank_tier="") for i in range(2)]
    )
    got = _run_summary(tmp_path, rows)
    assert got["RANK_TIER_AVAILABLE"] == "1"
    assert got["RANK_TIER_HIGH"] == "3"
    assert got["RANK_TIER_PLAUSIBLE"] == "4"
    assert got["RANK_TIER_TAIL"] == "5"
    assert got["RANK_TIER_NA"] == "2"


def test_summary_counts_partition_the_ranked_set(tmp_path: Path) -> None:
    """high + plausible + tail + not-scored must equal the ranked row count."""
    rows = (
        [_ranked_row(id=f"h{i}", rank_tier=TIER_HIGH) for i in range(2)]
        + [_ranked_row(id=f"p{i}", rank_tier=TIER_PLAUSIBLE) for i in range(3)]
        + [_ranked_row(id=f"t{i}", rank_tier=TIER_TAIL) for i in range(7)]
        + [_ranked_row(id=f"n{i}", rank_tier="") for i in range(1)]
    )
    got = _run_summary(tmp_path, rows)
    total = sum(int(got[f"RANK_TIER_{k}"])
                for k in ("HIGH", "PLAUSIBLE", "TAIL", "NA"))
    assert total == int(got["TOP_N_RANKED"]) == len(rows)


def test_summary_does_not_count_confidence_tier(tmp_path: Path) -> None:
    """The legacy tier must no longer drive any displayed count."""
    block = _extract_summary_block()
    assert "confidence_tier" not in block, (
        "stats block still counts the legacy confidence_tier")


def test_row_renderer_does_not_read_confidence_tier() -> None:
    """The per-row badge must no longer key on the legacy tier.

    Checks for an actual column ACCESS rather than any mention of the name, so
    a comment explaining why the legacy tier was retired stays allowed.
    """
    renderer = _extract_ranked_renderer()
    assert not re.search(r"""row\s*\.\s*get\(\s*["']confidence_tier""", renderer)
    assert not re.search(r"""row\s*\[\s*["']confidence_tier""", renderer)


def test_summary_table_rows_have_two_columns(tmp_path: Path) -> None:
    """The summary table is ``{lr}`` -- one '&' and EXACTLY one ``\\\\`` per row.

    The row terminator is checked exactly, not with ``endswith``: these rows are
    emitted from two different heredocs, and a quoted heredoc passes ``\\\\\\\\``
    through verbatim where an unquoted one collapses it to ``\\\\``. Four
    backslashes still ``endswith`` two, so a loose check silently accepts a
    doubled terminator -- which renders as a spurious blank row in the table.
    """
    for rows, cols in (
        ([_ranked_row(id="a", rank_tier=TIER_HIGH)], None),
        ([_ranked_row(id="a")], [c for c in RANKED_COLUMNS if c != "rank_tier"]),
    ):
        for line in _emit_tier_rows(tmp_path, rows, columns=cols):
            trailing = re.search(r"\\*$", line.rstrip()).group(0)
            assert trailing == "\\\\", (
                f"row terminator is {len(trailing)} backslashes, expected 2: {line!r}")
            assert line.count("&") == 1, f"expected 2 columns: {line}"


# ============================================================================
# Table structure
# ============================================================================

@pytest.mark.parametrize("tier", [TIER_HIGH, TIER_PLAUSIBLE, TIER_TAIL, ""])
def test_ranked_row_has_eight_columns(tmp_path: Path, tier: str) -> None:
    """The longtable is declared ``{rlrlcccr}`` -- 7 '&' and a ``\\\\`` per row."""
    out = _render_rows(tmp_path, [_ranked_row(rank_tier=tier)])
    line = out[0]
    assert line.rstrip().endswith("\\\\"), line
    assert line.count("&") == 7, f"expected 8 columns, got {line.count('&') + 1}"


def test_ranked_table_header_matches_column_spec() -> None:
    """Header row column count must match ``{rlrlcccr}`` and name the tier."""
    text = STAGE_09.read_text()
    block = re.search(r"\\begin\{longtable\}\{rlrlcccr\}.*?\\endhead",
                      text, flags=re.S)
    assert block, "top-20 longtable spec changed unexpectedly"
    headers = re.findall(r"^\\textbf\{\\#\}.*?\\\\$", block.group(0), flags=re.M)
    assert headers, "no header row found"
    for header in headers:
        assert header.count("&") == 7
        assert "Tier" in header, f"header does not name the tier column: {header}"
        assert "Conf" not in header, "header still labels the legacy Conf column"


def test_ranked_table_environment_is_balanced() -> None:
    text = STAGE_09.read_text()
    block = re.search(r"\\begin\{longtable\}\{rlrlcccr\}.*?\\end\{longtable\}",
                      text, flags=re.S)
    assert block, "top-20 longtable is not balanced"
    body = block.group(0)
    assert body.count(r"\begin{longtable}") == body.count(r"\end{longtable}") == 1


def test_summary_table_environment_is_balanced() -> None:
    """Splitting the summary heredoc must not orphan the table/tabular envs."""
    text = STAGE_09.read_text()
    block = re.search(
        r"\\begin\{table\}\[H\][^\\]*\\centering\s*\n"
        r"\\caption\{Pipeline Summary Statistics.*?\\end\{table\}",
        text, flags=re.S)
    assert block, "summary table block not found"
    body = block.group(0)
    for env in ("table", "tabular"):
        assert body.count(rf"\begin{{{env}}}") == body.count(rf"\end{{{env}}}") == 1
    assert body.count(r"\toprule") == 1
    assert body.count(r"\bottomrule") == 1


# ============================================================================
# Prose
# ============================================================================

def test_old_confidence_tier_section_is_gone() -> None:
    """The evidence-completeness tier explainer is now wrong and must go."""
    text = STAGE_09.read_text()
    assert r"\subsection{Confidence Tier Interpretation}" not in text


@pytest.mark.parametrize("stale", [
    "Strong evidence across multiple data sources",
    "Good evidence but incomplete data coverage",
    "Limited evidence or low scores",
])
def test_stale_evidence_completeness_prose_is_gone(stale: str) -> None:
    """Prose describing the OLD semantics must not survive anywhere in the stage."""
    assert stale not in STAGE_09.read_text(), (
        f"stage 09 still describes the old tier semantics: {stale!r}")


def test_prose_states_the_bootstrap_topk_semantics() -> None:
    prose = _tier_prose().lower()
    assert "bootstrap" in prose or "resampl" in prose
    assert "top" in prose
    assert "probability" in prose


def test_prose_names_all_three_tiers_and_their_cutoffs() -> None:
    prose = _tier_prose()
    for tier in TIERS:
        assert tier in prose, f"prose never names the '{tier}' tier"
    assert "0.8" in prose
    assert "0.2" in prose


def test_prose_explains_tail_is_not_a_failure() -> None:
    """The reader must be told ``tail`` is a rank statement, not bad data."""
    prose = _tier_prose().lower()
    assert "tail" in prose
    assert ("not a" in prose or "neither" in prose or "rather than" in prose), (
        "prose does not distinguish tail from a data-quality failure")


def test_prose_points_evidence_sufficiency_at_evidence_completeness() -> None:
    """Evidence sufficiency is a separate question with its own column."""
    prose = _tier_prose()
    assert "evidence" in prose.lower()
    assert re.search(r"evidence\\?_completeness", prose), (
        "prose does not direct the reader to evidence_completeness")


def test_prose_explains_the_unavailable_state() -> None:
    prose = _tier_prose().lower()
    assert "n/a" in prose or "not scored" in prose or "not computed" in prose


def test_recommendations_use_rank_tier() -> None:
    """The wet-lab shortlist criteria must key on the displayed tier."""
    text = STAGE_09.read_text()
    block = re.search(r"\\subsection\{High-Priority Wet-Lab Shortlist\}(.*?)"
                      r"\\subsection\{Follow-up Analyses\}", text, flags=re.S)
    assert block, "shortlist recommendations section not found"
    body = block.group(1)
    assert re.search(r"rank\\?_tier", body), (
        "shortlist criteria do not reference rank_tier")
    assert not re.search(r"confidence\\?_tier\s*=", body), (
        "shortlist criteria still gate on confidence_tier")
    assert "Medium confidence" not in body, (
        "shortlist still refers to the legacy Medium tier")


def test_ranking_figure_caption_describes_rank_tier() -> None:
    """Panel C now plots the rank-confidence tiers, so its caption must say so."""
    text = STAGE_09.read_text()
    m = re.search(r"\\caption\{Candidate ranking summary\.(.*?)\}\n", text, flags=re.S)
    assert m, "ranking figure caption not found"
    caption = m.group(1)
    assert "Confidence tier distribution" not in caption
    assert re.search(r"rank[- ]tier|rank confidence", caption, flags=re.I), (
        f"caption does not describe the rank-tier panel: {caption}")


# ============================================================================
# LaTeX validity
# ============================================================================

def test_amsmath_still_loaded() -> None:
    """A missing \\text macro previously made the stage exit non-zero."""
    preamble = STAGE_09.read_text().split(r"\begin{document}")[0]
    assert r"\usepackage{amsmath}" in preamble


def test_tier_prose_has_balanced_math_delimiters() -> None:
    """An odd number of unescaped ``$`` runs the rest of the doc into math mode."""
    dollars = _rendered_tier_prose().count("$")
    assert dollars % 2 == 0, f"unbalanced math delimiters in tier prose ({dollars})"


@pytest.mark.skipif(not shutil.which("pdflatex"), reason="pdflatex not installed")
def test_tier_display_compiles(tmp_path: Path) -> None:
    """Preamble + tier prose + summary rows + all four badge states must compile."""
    text = STAGE_09.read_text()
    preamble = re.search(r"(\\documentclass.*?)\n\\title\{", text, flags=re.S).group(1)

    prose = _rendered_tier_prose()

    rows = _render_rows(tmp_path, [
        _ranked_row(id="cand_a#1", rank_tier=TIER_HIGH, rank_score="0.91"),
        _ranked_row(id="cand_b", rank_tier=TIER_PLAUSIBLE, rank_score="0.55"),
        _ranked_row(id="cand_c", rank_tier=TIER_TAIL, rank_score="0.12"),
        _ranked_row(id="cand_d", rank_tier="", rank_score="0.05"),
    ])
    summary_rows = _emit_tier_rows(
        tmp_path, [_ranked_row(id="a", rank_tier=TIER_HIGH)])

    doc = "\n".join([
        preamble,
        r"\begin{document}",
        prose,
        r"\begin{tabular}{lr}", r"\toprule",
        "\n".join(summary_rows),
        r"\bottomrule", r"\end{tabular}",
        r"\begin{longtable}{rlrlcccr}", r"\toprule",
        r"\textbf{\#} & \textbf{ID} & \textbf{Score} & \textbf{Tier} & "
        r"\textbf{Sel} & \textbf{B-MH} & \textbf{HCR} & \textbf{TC size} \\",
        r"\midrule",
        "\n".join(rows),
        r"\bottomrule", r"\end{longtable}",
        r"\end{document}", "",
    ])
    tex = tmp_path / "doc.tex"
    tex.write_text(doc)
    proc = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "doc.tex"],
        cwd=tmp_path, capture_output=True, text=True, timeout=180)
    assert (tmp_path / "doc.pdf").exists(), (
        f"tier display did not compile:\n{proc.stdout[-4000:]}")
