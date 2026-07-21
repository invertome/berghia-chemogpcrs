"""Stage 03d caps reconciliation at 100 gene trees and must SAY so.

The defect
----------
Gene-tree discovery and the processing cap were fused into one statement::

    GENE_TREES=$(find "$PHYLO_DIR" ... | head -100)
    TREE_COUNT=$(echo "$GENE_TREES" | wc -l)
    log "Found ${TREE_COUNT} gene trees to reconcile"

``TREE_COUNT`` therefore counts the TRUNCATED list. A project with 100 gene
trees and a project with 4,000 both log "Found 100 gene trees to reconcile",
so a run that reconciled 2.5% of the orthogroups is textually identical to one
that reconciled all of them. The truncation is invisible in the log AND in
``events_summary.json``, whose ``total_trees_analyzed`` / ``total_duplications``
/ ``total_losses`` keys are lower bounds presented as totals.

This is the cross-cutting invariant: a partial measurement must not be
reported as a complete one. The cap itself is deliberately NOT changed here --
raising or removing it is a resource decision for the user -- but the real
total must be visible so that decision can be made on evidence.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE03D = PROJECT_ROOT / "03d_notung_reconciliation.sh"

# The cap the stage currently applies. Pinned so that a change to the cap has
# to come with a deliberate edit here rather than sliding through unnoticed.
EXPECTED_CAP = 100


def _extract_discovery_block() -> str:
    """Lift the live gene-tree discovery + count + report block out of 03d.

    Bounded by the '# Find gene trees' banner and the GeneRax branch that
    follows it, so the test exercises the shipped code rather than a copy.
    """
    lines = STAGE03D.read_text().splitlines()
    start = end = None
    for i, line in enumerate(lines):
        if start is None and line.strip() == "# Find gene trees":
            start = i
        elif start is not None and line.startswith("# --- GeneRax branch"):
            end = i
            break
    if start is None or end is None:
        pytest.fail("could not locate the gene-tree discovery block in stage 03d")
    block = "\n".join(lines[start:end])
    # The empty-set guard exits the whole stage; harmless here because every
    # test below writes at least one tree.
    return block


def _run_discovery(tmp_path: Path, n_trees: int) -> subprocess.CompletedProcess:
    """Create ``n_trees`` real IQ-TREE-style .treefile records and run the block."""
    results = tmp_path / "results"
    phylo = results / "phylogenies" / "protein" / "class_A"
    phylo.mkdir(parents=True, exist_ok=True)
    for i in range(n_trees):
        base = f"OG{i:07d}"
        d = phylo / base
        d.mkdir(exist_ok=True)
        # Real Newick, as stage 04's IQ-TREE --prefix emits.
        (d / f"{base}.treefile").write_text("(s1:0.1,s2:0.2);\n")

    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path / "logs")
    (tmp_path / "logs").mkdir(exist_ok=True)
    script = f"""
source "{FUNCTIONS_SH}"; trap - EXIT
RESULTS_DIR="{results}"
create_checkpoint() {{ :; }}
{_extract_discovery_block()}
echo "TREE_COUNT=${{TREE_COUNT}}"
echo "TREE_TOTAL=${{TREE_TOTAL}}"
echo "N_SELECTED=$(printf '%s\\n' "$GENE_TREES" | grep -c . )"
"""
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _kv(proc: subprocess.CompletedProcess, key: str) -> str:
    for line in proc.stdout.splitlines():
        if line.startswith(f"{key}="):
            return line.split("=", 1)[1].strip()
    pytest.fail(f"{key} not emitted:\n{proc.stdout}\n{proc.stderr}")


# --------------------------------------------------------------------------
# the truncated case: the real total must survive the cap
# --------------------------------------------------------------------------
def test_truncated_run_reports_the_true_total_not_the_cap(tmp_path):
    """With more trees than the cap, the discovered total must be recoverable.

    This is the defect: TREE_COUNT alone cannot distinguish "exactly 100 trees
    exist" from "4,000 exist and we silently took 100".
    """
    n = EXPECTED_CAP + 37
    proc = _run_discovery(tmp_path, n)
    assert proc.returncode == 0, proc.stdout + proc.stderr

    assert _kv(proc, "TREE_TOTAL") == str(n), (
        "the number of gene trees DISCOVERED must be reported; without it a "
        "truncated reconciliation is indistinguishable from a complete one"
    )
    assert _kv(proc, "N_SELECTED") == str(EXPECTED_CAP), (
        "the cap must still be applied -- this test pins reporting, not the cap"
    )


def test_truncation_is_logged_loudly_with_both_numbers(tmp_path):
    """A truncated run must emit a WARN naming processed AND discovered."""
    n = EXPECTED_CAP + 37
    proc = _run_discovery(tmp_path, n)
    combined = proc.stdout + proc.stderr

    assert "WARN" in combined, f"truncation was not raised as a warning:\n{combined}"
    assert str(n) in combined, (
        f"the true total {n} never appears in the log:\n{combined}"
    )
    # The message must not present the capped figure as the population.
    assert f"Found {EXPECTED_CAP} gene trees to reconcile" not in combined, (
        "the capped count is still being announced as though it were the total"
    )


def test_truncated_run_is_not_worded_like_a_complete_one(tmp_path):
    """The complete and truncated cases must produce distinguishable text."""
    full = _run_discovery(tmp_path / "a", EXPECTED_CAP + 5)
    part = _run_discovery(tmp_path / "b", 7)
    assert (full.stdout + full.stderr) != (part.stdout + part.stderr)


# --------------------------------------------------------------------------
# the complete case: no false alarm
# --------------------------------------------------------------------------
def test_complete_run_does_not_warn(tmp_path):
    """Below the cap nothing is truncated, so no truncation warning."""
    proc = _run_discovery(tmp_path, 7)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    combined = proc.stdout + proc.stderr

    assert _kv(proc, "TREE_TOTAL") == "7"
    assert _kv(proc, "TREE_COUNT") == "7"
    assert "WARN" not in combined, f"spurious truncation warning:\n{combined}"


def test_exactly_at_the_cap_is_complete_not_truncated(tmp_path):
    """The boundary: exactly ``cap`` trees is a COMPLETE run.

    This is the case the old wording could never distinguish, so it is the one
    most worth pinning.
    """
    proc = _run_discovery(tmp_path, EXPECTED_CAP)
    combined = proc.stdout + proc.stderr

    assert _kv(proc, "TREE_TOTAL") == str(EXPECTED_CAP)
    assert _kv(proc, "TREE_COUNT") == str(EXPECTED_CAP)
    assert "WARN" not in combined, (
        "exactly-at-the-cap is a complete run and must not warn:\n" + combined
    )


# --------------------------------------------------------------------------
# the output artifact, not just the log
# --------------------------------------------------------------------------
def test_events_summary_records_discovery_and_truncation():
    """``events_summary.json`` must carry the coverage denominator.

    ``total_trees_analyzed``/``total_duplications``/``total_losses`` are lower
    bounds under truncation; a reader of the JSON alone cannot currently tell.
    """
    text = STAGE03D.read_text()
    assert "gene_trees_discovered" in text, (
        "events_summary.json does not record how many gene trees existed, so "
        "its 'total_*' keys cannot be read as partial"
    )
    assert "reconciliation_truncated" in text, (
        "events_summary.json does not flag that the tree set was truncated"
    )
