"""Tests for the stage-09 positive-control table's three-state Status column.

``scripts/check_positive_controls.py`` (bead 444) stopped treating an ABSENT
control as an alert: only a control that is FOUND BUT BELOW the alert
percentile sets ``alert=True``. It gained a ``status`` column with exactly
three values — ``found_healthy`` / ``found_below_percentile`` / ``not_found``.

``09_report_generation.sh`` derived its Status column solely from ``alert``,
so after that change an absent control rendered as a green ``OK`` next to
``Found=N`` — a pass badge on a control that was never scored at all. The
prose above the table still claimed "ALERT rows indicate either a missing
control or one below the alert percentile", which is no longer true.

These tests pin the fix: the renderer reads ``status`` (not ``alert``) and
emits a THIRD, neutral state for absent controls, visually distinct from both
the healthy and the alert states; and the prose describes all three states and
says why absence is informational rather than a failure.

The renderer under test is a python heredoc embedded in the shell stage, so
these tests EXTRACT that block and execute it against synthetic TSVs — the
same behaviour the stage gets, rather than a grep over the source.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
STAGE_09 = REPO_ROOT / "09_report_generation.sh"

# The TSV schema check_positive_controls.py pins in OUTPUT_COLUMNS. Order and
# spelling are asserted against that module below so this fixture cannot drift.
POSCTRL_COLUMNS = [
    "gene_name", "found", "candidate_id", "rank", "percentile", "rank_score",
    "matched_strategy", "status", "alert", "notes",
]

STATUS_HEALTHY = "found_healthy"
STATUS_BELOW = "found_below_percentile"
STATUS_NOT_FOUND = "not_found"


# --- Harness ----------------------------------------------------------------

def _extract_poscontrol_renderer() -> str:
    """Return the body of the ``$POSCTRL_TSV`` python heredoc in stage 09.

    There are two ``<<'PY'`` blocks in the stage (the ranked-candidates table
    and this one); the opening line is matched on ``$POSCTRL_TSV`` so the
    right one is selected even if more are added later.
    """
    text = STAGE_09.read_text()
    m = re.search(
        r"^\s*python3 - \"\$POSCTRL_TSV\".*?<<'PY'\n(.*?)^PY$",
        text, flags=re.S | re.M,
    )
    assert m, "could not locate the positive-control renderer heredoc in stage 09"
    return m.group(1)


def _poscontrol_prose() -> str:
    """Return the explanatory prose paragraph rendered above the table."""
    text = STAGE_09.read_text()
    m = re.search(r"\\subsection\{Positive-Control Sanity Check\}(.*?)\\begin\{table\}",
                  text, flags=re.S)
    assert m, "could not locate the positive-control prose in stage 09"
    return m.group(1)


def _tsv_row(**overrides: str) -> dict:
    row = {c: "" for c in POSCTRL_COLUMNS}
    row["notes"] = "synthetic"
    row.update(overrides)
    return row


def _render(tmp_path: Path, rows: list[dict]) -> list[str]:
    """Run the extracted renderer over a synthetic TSV; return its stdout lines."""
    tsv = tmp_path / "positive_controls_check.tsv"
    lines = ["\t".join(POSCTRL_COLUMNS)]
    for row in rows:
        lines.append("\t".join(str(row.get(c, "")) for c in POSCTRL_COLUMNS))
    tsv.write_text("\n".join(lines) + "\n")

    script = tmp_path / "renderer.py"
    script.write_text(_extract_poscontrol_renderer())

    proc = subprocess.run([sys.executable, str(script), str(tsv)],
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, f"renderer failed:\n{proc.stderr}"
    return [ln for ln in proc.stdout.splitlines() if ln.strip()]


# --- Schema guard: the fixture must match the producer ----------------------

def test_fixture_schema_matches_producer() -> None:
    """The TSV columns these tests render must be the ones the producer writes."""
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    try:
        import check_positive_controls as cpc
    finally:
        sys.path.pop(0)
    assert cpc.OUTPUT_COLUMNS == POSCTRL_COLUMNS
    assert (cpc.STATUS_HEALTHY, cpc.STATUS_BELOW, cpc.STATUS_NOT_FOUND) == (
        STATUS_HEALTHY, STATUS_BELOW, STATUS_NOT_FOUND)


# --- The three rendered states ---------------------------------------------

def test_healthy_control_renders_ok(tmp_path: Path) -> None:
    """found_healthy -> green OK."""
    out = _render(tmp_path, [_tsv_row(
        gene_name="Ctrl_healthy", found="True", candidate_id="c1", rank="3",
        percentile="97.5", rank_score="0.9", matched_strategy="exact_id",
        status=STATUS_HEALTHY, alert="False")])
    assert len(out) == 1
    assert r"\textcolor{highconf}{OK}" in out[0]
    assert "ALERT" not in out[0]


def test_below_percentile_control_renders_alert(tmp_path: Path) -> None:
    """found_below_percentile -> red ALERT. This is the real drift signal."""
    out = _render(tmp_path, [_tsv_row(
        gene_name="Ctrl_sunk", found="True", candidate_id="c2", rank="700",
        percentile="11.4", rank_score="0.1", matched_strategy="exact_id",
        status=STATUS_BELOW, alert="True")])
    assert len(out) == 1
    assert r"\textcolor{lowconf}{ALERT}" in out[0]
    assert "OK" not in out[0]


def test_absent_control_renders_neutral_not_ok(tmp_path: Path) -> None:
    """not_found -> a THIRD neutral state; never the green OK pass badge.

    This is the regression: absence used to fall through to the `else` branch
    and render as OK beside `Found=N`.
    """
    out = _render(tmp_path, [_tsv_row(
        gene_name="Ga_olf", found="False", matched_strategy="none",
        status=STATUS_NOT_FOUND, alert="False")])
    assert len(out) == 1
    line = out[0]
    assert "{OK}" not in line, f"absent control rendered as a pass badge: {line}"
    assert "ALERT" not in line, f"absent control rendered as a failure: {line}"
    assert r"\textcolor{" in line, f"absent control has no colour treatment: {line}"


def test_three_states_are_visually_distinct(tmp_path: Path) -> None:
    """Each status maps to its own colour+label pair; no two collide."""
    out = _render(tmp_path, [
        _tsv_row(gene_name="A", found="True", rank="1", percentile="99.0",
                 status=STATUS_HEALTHY, alert="False"),
        _tsv_row(gene_name="B", found="True", rank="900", percentile="5.0",
                 status=STATUS_BELOW, alert="True"),
        _tsv_row(gene_name="C", found="False", status=STATUS_NOT_FOUND,
                 alert="False"),
    ])
    assert len(out) == 3
    badges = [re.search(r"\\textcolor\{[^}]*\}\{[^}]*\}", ln).group(0) for ln in out]
    assert len(set(badges)) == 3, f"states are not distinct: {badges}"


def test_absent_state_colour_is_defined_in_preamble(tmp_path: Path) -> None:
    """The colour the neutral badge uses must actually be \\definecolor'd."""
    out = _render(tmp_path, [_tsv_row(
        gene_name="Ga_olf", found="False", status=STATUS_NOT_FOUND,
        alert="False")])
    colour = re.search(r"\\textcolor\{([^}]*)\}", out[0]).group(1)
    text = STAGE_09.read_text()
    assert re.search(rf"\\definecolor\{{{re.escape(colour)}}}", text), (
        f"badge colour '{colour}' is never defined in the preamble")


def test_unknown_status_does_not_render_ok(tmp_path: Path) -> None:
    """A future/unrecognised status must not silently claim a pass."""
    out = _render(tmp_path, [_tsv_row(
        gene_name="Weird", found="True", status="some_future_status",
        alert="False")])
    assert "{OK}" not in out[0]


# --- Table structure --------------------------------------------------------

@pytest.mark.parametrize("status,alert,found", [
    (STATUS_HEALTHY, "False", "True"),
    (STATUS_BELOW, "True", "True"),
    (STATUS_NOT_FOUND, "False", "False"),
])
def test_row_has_five_columns_and_row_terminator(
        tmp_path: Path, status: str, alert: str, found: str) -> None:
    """The table is declared ``{llrrl}`` — every row needs 4 '&' and a ``\\\\``."""
    out = _render(tmp_path, [_tsv_row(
        gene_name="Ctrl", found=found, rank="5", percentile="80.0",
        status=status, alert=alert)])
    line = out[0]
    assert line.rstrip().endswith(r"\\"), line
    # '&' only ever appears as a column separator here (gene names are escaped).
    assert line.count("&") == 4, f"expected 5 columns, got {line.count('&') + 1}: {line}"


def test_gene_name_with_latex_specials_is_escaped(tmp_path: Path) -> None:
    """Escaping must survive the status rework (an unescaped _ aborts pdflatex)."""
    out = _render(tmp_path, [_tsv_row(
        gene_name="Ga_olf#1", found="False", status=STATUS_NOT_FOUND,
        alert="False")])
    assert r"Ga\_olf\#1" in out[0]
    assert out[0].count("&") == 4


def test_table_environment_is_balanced() -> None:
    """The surrounding LaTeX block stays well-formed."""
    text = STAGE_09.read_text()
    block = re.search(
        r"\\subsection\{Positive-Control Sanity Check\}.*?\\end\{table\}",
        text, flags=re.S).group(0)
    for env in ("table", "tabular"):
        assert block.count(rf"\begin{{{env}}}") == block.count(rf"\end{{{env}}}") == 1
    assert block.count(r"\toprule") == 1
    assert block.count(r"\midrule") == 1
    # Header row column count must match the {llrrl} column spec.
    header = re.search(r"\\textbf\{Control\}.*?\\\\", block, flags=re.S).group(0)
    assert header.count("&") == 4


# --- Prose ------------------------------------------------------------------

def test_prose_no_longer_calls_absence_an_alert() -> None:
    """The stale claim that ALERT covers 'a missing control' must be gone."""
    prose = _poscontrol_prose()
    assert "missing control" not in prose.lower(), (
        "prose still claims a missing control is an ALERT")


def test_prose_describes_all_three_states() -> None:
    """Prose must name the neutral state and explain why absence is not a failure."""
    prose = _poscontrol_prose()
    assert "ALERT" in prose
    assert "OK" in prose
    # The neutral badge label the renderer actually emits must appear in prose.
    renderer = _extract_poscontrol_renderer()
    labels = set(re.findall(r"\\textcolor\{[^}]*\}\{([A-Z ]+)\}", renderer))
    neutral = labels - {"OK", "ALERT"}
    assert neutral, "renderer emits no third badge label"
    for label in neutral:
        assert label in prose, f"prose never explains the '{label}' state"


def test_prose_explains_absence_is_informational() -> None:
    """The reader must be told absence is informational, not a pipeline failure."""
    prose = _poscontrol_prose().lower()
    assert "informational" in prose or "not a failure" in prose


# --- LaTeX validity ---------------------------------------------------------

def test_math_mode_text_macro_has_its_package() -> None:
    """The prose uses \\text{} in math mode, which requires amsmath.

    Found by compiling the section: without amsmath, pdflatex aborts on
    ``G$\\alpha_{\\text{olf}}$`` with a cascading "Undefined control sequence"
    and produces NO pdf, so stage 09's second latex pass exits 1.
    """
    text = STAGE_09.read_text()
    preamble = text.split(r"\begin{document}")[0]
    if re.search(r"\\text\{", text):
        assert r"\usepackage{amsmath}" in preamble, (
            r"prose uses \text{} but amsmath is not loaded")


@pytest.mark.skipif(not shutil.which("pdflatex"), reason="pdflatex not installed")
def test_positive_control_section_compiles(tmp_path: Path) -> None:
    """The real preamble + prose + all rendered states must compile to a PDF."""
    text = STAGE_09.read_text()
    preamble = re.search(r"(\\documentclass.*?)\n\\title\{", text, flags=re.S).group(1)
    section = re.search(
        r"(\\subsection\{Positive-Control Sanity Check\}.*?\\end\{table\})",
        text, flags=re.S).group(1)
    # The section spans TWO heredocs (prose+table head, then closing rules)
    # with the renderer between them. Drop the shell scaffolding so only the
    # LaTeX the stage actually emits reaches pdflatex.
    kept, skipping = [], False
    for line in section.splitlines():
        if line.rstrip() == "EOF":
            skipping = True
            continue
        if skipping:
            if line.rstrip().endswith("<<'EOF'"):
                skipping = False
            continue
        kept.append(line)
    section = "\n".join(kept)

    rows = _render(tmp_path, [
        _tsv_row(gene_name="Ga_olf#1", found="False", matched_strategy="none",
                 status=STATUS_NOT_FOUND, alert="False"),
        _tsv_row(gene_name="Ctrl_ok", found="True", rank="3", percentile="97.5",
                 status=STATUS_HEALTHY, alert="False"),
        _tsv_row(gene_name="Ctrl_sunk", found="True", rank="702",
                 percentile="11.4", status=STATUS_BELOW, alert="True"),
        _tsv_row(gene_name="Legacy", found="True", rank="10", percentile="88.0",
                 status="", alert="False"),
    ])
    section = section.replace(r"\midrule", "\\midrule\n" + "\n".join(rows), 1)

    tex = tmp_path / "doc.tex"
    tex.write_text(f"{preamble}\n\\begin{{document}}\n{section}\n\\end{{document}}\n")
    proc = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "doc.tex"],
        cwd=tmp_path, capture_output=True, text=True, timeout=180)
    assert (tmp_path / "doc.pdf").exists(), (
        f"section did not compile:\n{proc.stdout[-3000:]}")
