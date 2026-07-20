"""Stage 05 must read stage 04's PER-CLASS per-orthogroup outputs (bead F1).

Stage 04's per-class refactor writes every per-orthogroup alignment and tree
two levels deeper than they used to live::

    ${RESULTS_DIR}/phylogenies/protein/class_<CLASS>/<base>/<base>_trimmed.fa
    ${RESULTS_DIR}/phylogenies/protein/class_<CLASS>/<base>/<base>.treefile

Stage 05 still read the pre-refactor FLAT pair
(``protein/<base>_trimmed.fa`` / ``protein/<base>.treefile``). Nothing in the
repository writes those flat per-OG paths -- grepping turns up readers only --
so the ``[ -f "$protein_align" ] && [ -f "$tree" ]`` guard was never true and
NO aBSREL / GARD / BUSTED-S / BUSTED-MH / MEME / ASR output was ever produced.
Stage 05 nonetheless touched ``step_completed_05.txt``, so stage 07 ran on and
merely logged "aBSREL results not found": a silent, whole-stage no-op.

The fix makes stage 05 resolve the orthogroup's class exactly the way stage 04
does (``results/classification/og_class_majority.tsv``, majority class, with
stage 04's two fallbacks) and read from the resulting per-class directory.

These tests execute the REAL ``resolve_perog_paths`` shell function lifted out
of the shipping stage 05 script, so they break if the resolution logic drifts
away from stage 04's.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE04 = PROJECT_ROOT / "04_phylogenetic_analysis.sh"
STAGE05 = PROJECT_ROOT / "05_selective_pressure_and_asr.sh"


# --------------------------------------------------------------------------
# harness
# --------------------------------------------------------------------------
def _extract_shell_function(script: Path, name: str) -> str:
    """Lift a top-level `name() {` ... `}` definition out of a real stage script."""
    lines = script.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if re.match(rf"^{re.escape(name)}\(\)\s*\{{\s*$", line):
            start = i
            break
    if start is None:
        pytest.fail(
            f"{script.name} defines no `{name}()` function; stage 05 must resolve "
            f"per-OG paths through stage 04's per-class layout"
        )
    for j in range(start + 1, len(lines)):
        if lines[j] == "}":
            return "\n".join(lines[start : j + 1])
    pytest.fail(f"could not find the closing brace of `{name}()` in {script.name}")


def _run_resolver(tmp_path: Path, base: str) -> subprocess.CompletedProcess:
    """Run the real resolver against a synthetic ${RESULTS_DIR} tree."""
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path / "logs")
    (tmp_path / "logs").mkdir(exist_ok=True)
    script = f"""
source "{FUNCTIONS_SH}"; trap - EXIT
RESULTS_DIR="{tmp_path / 'results'}"
{_extract_shell_function(STAGE05, "resolve_perog_paths")}
if resolve_perog_paths "{base}"; then
    echo "RESOLVED_OK"
else
    echo "RESOLVED_FAIL"
fi
echo "CLASS=${{OG_CLASS}}"
echo "ALIGN=${{protein_align}}"
echo "TREE=${{tree}}"
"""
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _fields(proc: subprocess.CompletedProcess) -> dict[str, str]:
    out = {}
    for line in proc.stdout.splitlines():
        if "=" in line and line.split("=", 1)[0] in {"CLASS", "ALIGN", "TREE"}:
            k, v = line.split("=", 1)
            out[k] = v
    return out


def _make_og_outputs(tmp_path: Path, cls: str, base: str, *, align=True, tree=True) -> Path:
    """Create stage 04's per-class per-OG output directory and its products."""
    d = tmp_path / "results" / "phylogenies" / "protein" / f"class_{cls}" / base
    d.mkdir(parents=True, exist_ok=True)
    if align:
        (d / f"{base}_trimmed.fa").write_text(">s1\nMKFLV\n>s2\nMKFLA\n")
    if tree:
        (d / f"{base}.treefile").write_text("(s1:0.1,s2:0.2);\n")
    return d


def _write_class_majority(tmp_path: Path, rows: dict[str, str]) -> Path:
    p = tmp_path / "results" / "classification" / "og_class_majority.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("orthogroup\tclass\n" + "".join(f"{og}\t{c}\n" for og, c in rows.items()))
    return p


# --------------------------------------------------------------------------
# the core regression: the per-class per-OG pair is FOUND
# --------------------------------------------------------------------------
def test_resolves_class_a_perog_alignment_and_tree(tmp_path):
    """The exact layout stage 04 writes must now resolve."""
    _write_class_majority(tmp_path, {"OG0000123": "A"})
    d = _make_og_outputs(tmp_path, "A", "OG0000123")

    proc = _run_resolver(tmp_path, "OG0000123")

    assert "RESOLVED_OK" in proc.stdout, proc.stdout + proc.stderr
    f = _fields(proc)
    assert f["CLASS"] == "A"
    assert f["ALIGN"] == str(d / "OG0000123_trimmed.fa")
    assert f["TREE"] == str(d / "OG0000123.treefile")


def test_resolves_non_class_a_orthogroup(tmp_path):
    """Class routing is read from the TSV, not hard-coded to A."""
    _write_class_majority(tmp_path, {"OG0000123": "A", "OG0000900": "C"})
    d = _make_og_outputs(tmp_path, "C", "OG0000900")

    proc = _run_resolver(tmp_path, "OG0000900")

    assert "RESOLVED_OK" in proc.stdout, proc.stdout + proc.stderr
    f = _fields(proc)
    assert f["CLASS"] == "C"
    assert f["TREE"] == str(d / "OG0000900.treefile")


def test_flat_legacy_path_is_not_used(tmp_path):
    """A file at the old flat location must NOT satisfy the resolver.

    Pins the direction of the fix: the flat path has no writer, so honouring it
    would just re-open the hole this bug came from.
    """
    _write_class_majority(tmp_path, {"OG0000123": "A"})
    flat = tmp_path / "results" / "phylogenies" / "protein"
    flat.mkdir(parents=True, exist_ok=True)
    (flat / "OG0000123_trimmed.fa").write_text(">s1\nMKFLV\n")
    (flat / "OG0000123.treefile").write_text("(s1:0.1,s2:0.2);\n")

    proc = _run_resolver(tmp_path, "OG0000123")

    assert "RESOLVED_FAIL" in proc.stdout, proc.stdout + proc.stderr


# --------------------------------------------------------------------------
# stage 04's two class-resolution fallbacks, mirrored exactly
# --------------------------------------------------------------------------
def test_missing_class_tsv_falls_back_to_class_a(tmp_path):
    """Stage 04 (04:531-533): no og_class_majority.tsv -> class A, back-compat."""
    d = _make_og_outputs(tmp_path, "A", "OG0000123")

    proc = _run_resolver(tmp_path, "OG0000123")

    assert "RESOLVED_OK" in proc.stdout, proc.stdout + proc.stderr
    f = _fields(proc)
    assert f["CLASS"] == "A"
    assert f["TREE"] == str(d / "OG0000123.treefile")


def test_og_absent_from_tsv_routes_to_unclassified(tmp_path):
    """Stage 04 (04:527-530): OG missing from the TSV -> class_unclassified."""
    _write_class_majority(tmp_path, {"OG0000999": "A"})
    d = _make_og_outputs(tmp_path, "unclassified", "OG0000123")

    proc = _run_resolver(tmp_path, "OG0000123")

    assert "RESOLVED_OK" in proc.stdout, proc.stdout + proc.stderr
    f = _fields(proc)
    assert f["CLASS"] == "unclassified"
    assert f["TREE"] == str(d / "OG0000123.treefile")


def test_class_resolution_matches_stage04_source(tmp_path):
    """The awk lookup and both fallbacks must stay byte-identical to stage 04.

    Guards against the two stages drifting apart again, which is precisely the
    failure mode that produced this bug.
    """
    stage04 = STAGE04.read_text()
    stage05 = STAGE05.read_text()
    awk_expr = """awk -F'\\t' -v og="${base}" 'NR>1 && $1==og {print $2; exit}'"""
    assert awk_expr in stage04, "stage 04's class lookup moved; re-derive stage 05's"
    # stage 05 uses a local variable name, so compare the invariant core
    core = """'NR>1 && $1==og {print $2; exit}'"""
    assert core in stage05
    assert "og_class_majority.tsv" in stage05
    assert 'OG_CLASS="unclassified"' in stage04 and 'unclassified' in stage05
    assert '"${RESULTS_DIR}/phylogenies/protein/class_' in stage05


# --------------------------------------------------------------------------
# the "missing" branch must DISCRIMINATE, not cry wolf
# --------------------------------------------------------------------------
def test_orthogroup_with_no_tree_anywhere_is_handled_gracefully(tmp_path):
    """Stage 04 legitimately builds no tree for some OGs: expected, not a defect.

    This must NOT be logged as an error -- a message that fires for every
    orthogroup is exactly what hid the broken path for so long.
    """
    _write_class_majority(tmp_path, {"OG0000123": "A"})

    proc = _run_resolver(tmp_path, "OG0000123")

    assert proc.returncode == 0, proc.stderr
    assert "RESOLVED_FAIL" in proc.stdout
    assert _fields(proc)["TREE"] == ""
    assert "[ERROR]" not in proc.stdout + proc.stderr


def test_tree_under_a_different_class_is_reported_as_an_error(tmp_path):
    """Routing inconsistency between stages 04 and 05 must be LOUD.

    Happens when og_class_majority.tsv is (re)generated after stage 04 already
    ran under the class_A back-compat fallback. Distinguishable from "this OG
    has no tree", and the whole point of splitting the old single warning.
    """
    _write_class_majority(tmp_path, {"OG0000123": "C"})
    d = _make_og_outputs(tmp_path, "A", "OG0000123")  # stage 04 wrote it under A

    proc = _run_resolver(tmp_path, "OG0000123")

    combined = proc.stdout + proc.stderr
    assert "[ERROR]" in combined, combined
    assert "class" in combined.lower()
    # the single unambiguous match is adopted so the stage still produces output
    assert "RESOLVED_OK" in proc.stdout, combined
    assert _fields(proc)["TREE"] == str(d / "OG0000123.treefile")


def test_ambiguous_class_placement_refuses_to_guess(tmp_path):
    """Same base under two classes: error and skip rather than pick one."""
    _write_class_majority(tmp_path, {"OG0000123": "F"})
    _make_og_outputs(tmp_path, "A", "OG0000123")
    _make_og_outputs(tmp_path, "C", "OG0000123")

    proc = _run_resolver(tmp_path, "OG0000123")

    combined = proc.stdout + proc.stderr
    assert "[ERROR]" in combined, combined
    assert "RESOLVED_FAIL" in proc.stdout
    assert _fields(proc)["TREE"] == ""


def test_tree_without_alignment_is_flagged(tmp_path):
    """Stage 04 got as far as IQ-TREE but the trimmed alignment is gone."""
    _write_class_majority(tmp_path, {"OG0000123": "A"})
    _make_og_outputs(tmp_path, "A", "OG0000123", align=False)

    proc = _run_resolver(tmp_path, "OG0000123")

    combined = proc.stdout + proc.stderr
    assert "RESOLVED_FAIL" in proc.stdout
    assert "[WARN]" in combined or "[ERROR]" in combined


# --------------------------------------------------------------------------
# both call sites must go through the resolver
# --------------------------------------------------------------------------
def test_both_selection_and_asr_branches_use_the_resolver(tmp_path):
    """The dN/dS branch AND the ASR branch were both broken; both must be fixed."""
    src = STAGE05.read_text()
    # count invocations (`... resolve_perog_paths "$..."`), not the definition
    calls = re.findall(r'resolve_perog_paths\s+"', src)
    assert len(calls) >= 2, (
        f"expected the selection branch and the ASR branch to both call "
        f"resolve_perog_paths, found {len(calls)}"
    )


def test_no_flat_perog_paths_remain_in_stage05(tmp_path):
    """No writer exists for protein/<base>_trimmed.fa -- no reader may remain."""
    src = STAGE05.read_text()
    for dead in (
        '"${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"',
        '"${RESULTS_DIR}/phylogenies/protein/${base}.treefile"',
    ):
        assert dead not in src, f"stage 05 still reads the unwritten flat path {dead}"
