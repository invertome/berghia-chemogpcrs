"""Bead -0mhv: run_selection_stack.sh must not report success it did not
earn, must not advertise a branch pre-filter it does not implement, and must
never fabricate a MEME strict/lenient concordance signal.

(a) Every HyPhy call was suffixed `|| log --level=WARN`, so all five methods
    could fail while the script exited 0 and logged "Selection stack
    complete". The caller's `|| log WARN` guard at 05:396 was unreachable and
    the dN/dS axis silently contributed nothing.

(b) `--branches FOREGROUND` was passed against a tree that is never labelled
    (the code's own comment carried a TODO), so the header's claimed
    statistical-power advantage and BH-FDR q=0.10 did not exist.

(c) When ClipKit failed, both the strict and lenient codon alignments became
    `cp` copies of the SAME raw alignment. The dual-pass MEME concordance
    check then compared an alignment against itself, scored trivially 100%
    concordant, and fed that fabricated high-confidence signal into
    rank_candidates' positive_score boost.

The harness (stub hyphy/clipkit/parsers around the REAL script) is shared
with test_ih5u_selection_stack_concurrency.py.
"""
from __future__ import annotations

import sys
from pathlib import Path

try:  # pytest may import these as a package or as top-level modules
    from .test_ih5u_selection_stack_concurrency import (
        REAL_SCRIPT, make_project, run_stack,
    )
except ImportError:  # pragma: no cover - import-mode fallback
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_ih5u_selection_stack_concurrency import (  # type: ignore
        REAL_SCRIPT, make_project, run_stack,
    )

OG = "OGintegrity_main"

# Method names as the stub hyphy sees them ($1 of the hyphy command line).
ALL_HYPHY_METHODS = "gard busted absrel meme"


def setup(tmp_path: Path) -> tuple[Path, Path]:
    root = tmp_path / "proj"
    make_project(root)
    return root, root / "results" / "selective_pressure"


def read_status(out_dir: Path, og: str = OG) -> dict[str, str]:
    path = out_dir / f"{og}_selection_status.tsv"
    if not path.exists():
        return {}
    rows = [ln.split("\t") for ln in path.read_text().splitlines() if ln.strip()]
    assert rows[0] == ["og_name", "method", "status"], rows[0]
    return {r[1]: r[2] for r in rows[1:]}


def code_lines() -> list[str]:
    """Script lines with comments stripped."""
    return [ln for ln in REAL_SCRIPT.read_text().splitlines()
            if not ln.lstrip().startswith("#")]


def header_text() -> str:
    """The leading comment block, before `set -euo pipefail`."""
    out = []
    for ln in REAL_SCRIPT.read_text().splitlines():
        if ln.startswith("set -"):
            break
        out.append(ln)
    return "\n".join(out)


# ===========================================================================
# (a) failure must be detectable
# ===========================================================================

def test_all_methods_failing_exits_nonzero(tmp_path):
    """The whole point: total failure must be visible to the caller."""
    root, out_dir = setup(tmp_path)
    res = run_stack(root, OG, {"STUB_HYPHY_FAIL": ALL_HYPHY_METHODS})
    assert res.returncode != 0, (
        "every HyPhy method failed but the script still exited 0 — "
        "the caller's `|| log WARN` guard can never fire"
    )
    assert "NO usable inference output" in res.stderr


def test_all_methods_failing_does_not_claim_completion(tmp_path):
    root, _ = setup(tmp_path)
    res = run_stack(root, OG, {"STUB_HYPHY_FAIL": ALL_HYPHY_METHODS})
    assert "Selection stack complete" not in res.stderr


def test_status_tsv_reports_every_method(tmp_path):
    """A caller must be able to act on per-method outcomes, not just a
    single exit code."""
    root, out_dir = setup(tmp_path)
    run_stack(root, OG)
    status = read_status(out_dir)
    for method in ("gard", "busted_s", "busted_mh", "absrel",
                   "meme_strict", "meme_lenient"):
        assert method in status, f"{method} missing from status TSV: {status}"
        assert status[method] == "ok", f"{method}={status[method]}"


def test_status_tsv_marks_failures(tmp_path):
    root, out_dir = setup(tmp_path)
    run_stack(root, OG, {"STUB_HYPHY_FAIL": "busted absrel"})
    status = read_status(out_dir)
    assert status["busted_s"] == "failed"
    assert status["busted_mh"] == "failed"
    assert status["absrel"] == "failed"
    assert status["gard"] == "ok"


def test_exit_zero_without_output_is_not_success(tmp_path):
    """HyPhy can exit 0 having written nothing; a zero-length JSON would
    otherwise sail through the parse guards and yield an empty parse."""
    root, out_dir = setup(tmp_path)
    res = run_stack(root, OG, {"STUB_HYPHY_EMPTY": ALL_HYPHY_METHODS})
    status = read_status(out_dir)
    assert status["busted_s"] == "no_output"
    assert status["absrel"] == "no_output"
    assert res.returncode != 0, "all methods produced empty output but exit was 0"


def test_partial_failure_is_tolerated(tmp_path):
    """Policy: the stack is useful if ANY inference method succeeded, so a
    partial failure must not abort the orthogroup."""
    root, out_dir = setup(tmp_path)
    res = run_stack(root, OG, {"STUB_HYPHY_FAIL": "gard absrel meme"})
    assert res.returncode == 0, res.stderr[-2000:]
    assert "Selection stack complete" in res.stderr
    status = read_status(out_dir)
    assert status["busted_s"] == "ok"
    assert status["absrel"] == "failed"


def test_gard_failure_alone_does_not_gate_exit(tmp_path):
    """GARD is an advisory recombination screen; nothing downstream reads
    gard.json, so its failure must not fail the orthogroup."""
    root, _ = setup(tmp_path)
    res = run_stack(root, OG, {"STUB_HYPHY_FAIL": "gard"})
    assert res.returncode == 0, res.stderr[-2000:]


# ===========================================================================
# (b) the branch pre-filter must not be advertised or silently no-op'd
# ===========================================================================

def test_branch_list_argument_is_rejected(tmp_path):
    """Accepting the argument and ignoring it is what made this a no-op."""
    root, _ = setup(tmp_path)
    branch_file = root / "branches.txt"
    branch_file.write_text("Node1\nNode2\n")
    script = root / "scripts" / "hpc" / "run_selection_stack.sh"
    align = root / "aln.phy"
    align.write_text(">s1\nAAA\n")
    tree = root / "t.nwk"
    tree.write_text("(s1,s2);\n")

    import os
    import subprocess
    env = os.environ.copy()
    env["PATH"] = f"{root / 'bin'}{os.pathsep}{env['PATH']}"
    res = subprocess.run(
        ["bash", str(script), str(align), str(tree), OG, str(branch_file)],
        capture_output=True, text=True, env=env,
    )
    assert res.returncode == 2, (
        f"a branch list was supplied but not rejected (rc={res.returncode})"
    )
    assert "not implemented" in res.stderr.lower()


def test_foreground_flag_not_passed_to_absrel():
    """`--branches FOREGROUND` against an unlabelled tree is a no-op."""
    lines = code_lines()
    assert not any("ABSREL_BRANCH_ARG" in ln for ln in lines), (
        "the ABSREL_BRANCH_ARG indirection is still present"
    )
    # Extract the aBSREL invocation (continuation lines end with a backslash).
    block: list[str] = []
    for i, ln in enumerate(lines):
        if ln.startswith("hyphy absrel"):
            block.append(ln)
            j = i
            while lines[j].rstrip().endswith("\\"):
                j += 1
                block.append(lines[j])
            break
    assert block, "aBSREL invocation not found"
    joined = " ".join(block)
    assert "--branches" not in joined, f"aBSREL still passes --branches: {joined}"


def test_header_does_not_claim_unimplemented_bh_fdr():
    """The header asserted BH-FDR q=0.10 over a pre-filtered branch set.
    Neither exists in this script."""
    header = header_text()
    assert "q=0.10" not in header, (
        "header still claims a BH-FDR q=0.10 correction the script never applies"
    )


def test_header_has_no_dangling_todo():
    """A header must not document a property as pending-and-assumed."""
    assert "TODO" not in REAL_SCRIPT.read_text(), (
        "a TODO still stands in for unimplemented branch labelling"
    )


# ===========================================================================
# (c) the concordance signal must never be fabricated
# ===========================================================================

def test_no_raw_alignment_fallback_in_code():
    """`cp "$CODON_ALIGN" "$CODON_STRICT"` is what collapsed the two passes."""
    joined = "\n".join(code_lines())
    assert 'cp "$CODON_ALIGN"' not in joined, (
        "ClipKit failure still falls back to copying the raw alignment, "
        "which makes the two MEME passes the same file"
    )


def test_clipkit_strict_failure_produces_no_concordance(tmp_path):
    root, out_dir = setup(tmp_path)
    run_stack(root, OG, {"STUB_CLIPKIT_FAIL": "kpic-smart-gap"})

    assert not (out_dir / f"{OG}_codon_strict.fa").exists(), (
        "a strict alignment exists despite ClipKit failing — "
        "the raw alignment was substituted"
    )
    assert not (out_dir / f"{OG}_meme_concordance.csv").exists(), (
        "concordance was computed from a collapsed comparison"
    )
    assert not (out_dir / "meme_concordance.csv").exists(), (
        "a fabricated concordance row reached the cumulative CSV that "
        "rank_candidates.py reads"
    )


def test_clipkit_failure_is_loud(tmp_path):
    root, out_dir = setup(tmp_path)
    res = run_stack(root, OG, {"STUB_CLIPKIT_FAIL": "kpic-smart-gap"})
    assert "[ERROR]" in res.stderr, res.stderr[-2000:]
    assert "ClipKit strict-trim failed" in res.stderr
    assert read_status(out_dir)["clipkit_strict"] == "failed"


def test_clipkit_lenient_failure_produces_no_concordance(tmp_path):
    """The sibling path must be covered too, not just the one named first."""
    root, out_dir = setup(tmp_path)
    res = run_stack(root, OG, {"STUB_CLIPKIT_FAIL": "smart-gap"})
    assert not (out_dir / f"{OG}_codon_lenient.fa").exists()
    assert not (out_dir / f"{OG}_meme_concordance.csv").exists()
    assert read_status(out_dir)["clipkit_lenient"] == "failed"
    # The strict MEME pass still runs and still contributes.
    assert (out_dir / f"{OG}_meme.csv").exists()


def test_identical_trims_do_not_score_as_concordant(tmp_path):
    """Both regimes can legitimately remove the same columns. A
    self-comparison is trivially 100% concordant, so it must not be scored."""
    root, out_dir = setup(tmp_path)
    run_stack(root, OG, {"STUB_CLIPKIT_IDENTICAL": "1"})

    strict = out_dir / f"{OG}_codon_strict.fa"
    lenient = out_dir / f"{OG}_codon_lenient.fa"
    assert strict.read_bytes() == lenient.read_bytes(), "stub setup wrong"

    assert not (out_dir / f"{OG}_meme_concordance.csv").exists()
    assert not (out_dir / "meme_concordance.csv").exists()
    assert read_status(out_dir)["meme_concordance"] == "degenerate_identical_trims"


def test_stale_concordance_is_removed_not_reconcatenated(tmp_path):
    """A concordance CSV left by an earlier run must not be swept back into
    the cumulative once the comparison is judged degenerate."""
    root, out_dir = setup(tmp_path)
    stale = out_dir / f"{OG}_meme_concordance.csv"
    stale.write_text("og_name,high_confidence_sites_n\n" + OG + ",99\n")

    run_stack(root, OG, {"STUB_CLIPKIT_IDENTICAL": "1"})

    assert not stale.exists(), "stale concordance output survived"
    assert not (out_dir / "meme_concordance.csv").exists()


def test_distinct_trims_do_produce_concordance(tmp_path):
    """The happy path must still work — the fix must not disable the feature."""
    root, out_dir = setup(tmp_path)
    res = run_stack(root, OG)
    assert res.returncode == 0, res.stderr[-2000:]

    strict = out_dir / f"{OG}_codon_strict.fa"
    lenient = out_dir / f"{OG}_codon_lenient.fa"
    assert strict.read_bytes() != lenient.read_bytes()

    assert (out_dir / f"{OG}_meme_concordance.csv").exists()
    assert (out_dir / "meme_concordance.csv").exists()
    assert read_status(out_dir)["meme_concordance"] == "ok"


def test_concordance_never_emitted_from_identical_inputs(tmp_path):
    """Invariant across all ClipKit outcomes: if a concordance row was
    published, the two passes were genuinely different alignments."""
    scenarios = [
        {},
        {"STUB_CLIPKIT_IDENTICAL": "1"},
        {"STUB_CLIPKIT_FAIL": "kpic-smart-gap"},
        {"STUB_CLIPKIT_FAIL": "smart-gap"},
        {"STUB_CLIPKIT_FAIL": "kpic-smart-gap smart-gap"},
    ]
    for env_extra in scenarios:
        root, out_dir = setup(tmp_path / f"s{len(env_extra)}{'-'.join(env_extra)}")
        run_stack(root, OG, env_extra)
        if not (out_dir / f"{OG}_meme_concordance.csv").exists():
            continue
        strict = out_dir / f"{OG}_codon_strict.fa"
        lenient = out_dir / f"{OG}_codon_lenient.fa"
        assert strict.exists() and lenient.exists(), env_extra
        assert strict.read_bytes() != lenient.read_bytes(), (
            f"concordance published from identical alignments under {env_extra}"
        )
