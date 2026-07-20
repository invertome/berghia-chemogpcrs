"""A completion flag must assert the stage PRODUCED something, not that it ran.

Of 18 stages, exactly one verified a declared output before writing its
completion flag. ``step_completed_05.txt`` meant "the script reached its last
line" -- so stage 05, reading a path a refactor had moved two levels deeper,
produced nothing for any orthogroup, wrote its flag anyway, and stage 07
proceeded and merely logged "aBSREL results not found".

The mechanism here is the one the repo already proved out in
``06c_classify_non_chemoreceptors.sh`` (named sub-step tracking through
``is_step_completed "<substep>"``) -- propagated rather than reinvented, and
extended so that the ONLY thing that writes a flag is a verifier that checked a
declaration.

Two design points the tests pin:

  * Cardinality, not "non-empty". Some stages legitimately produce nothing for
    a given dataset, so a stage declares ``0..N`` or ``1..N`` (or an exact
    ``k..k``) rather than being force-fed a non-emptiness assumption that would
    turn a legitimate empty result into a false alarm.
  * No literal paths the resolver can derive. A declaration may name an output
    by its RELATIONSHIP (``og_product:tree``) so the stage script never repeats
    a path string that can drift.

Everything is exercised with ``functions.sh`` sourced ALONE, because stages are
submitted with a bare ``sbatch`` and never see ``run_pipeline.sh``.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"


def _run(tmp_path: Path, body: str) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    logs = tmp_path / "logs"
    logs.mkdir(exist_ok=True)
    env["LOGS_DIR"] = str(logs)
    results = tmp_path / "results"
    results.mkdir(exist_ok=True)
    script = f"""
source "{FUNCTIONS_SH}"; trap - EXIT
RESULTS_DIR="{results}"
{body}
"""
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _touch(tmp_path: Path, rel: str, content: str = "x\n") -> Path:
    p = tmp_path / "results" / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content)
    return p


def _flag(tmp_path: Path, stage: str) -> Path:
    return tmp_path / "results" / f"step_completed_{stage}.txt"


# --------------------------------------------------------------------------
# cardinality: 1..N
# --------------------------------------------------------------------------
def test_one_or_more_passes_when_output_present(tmp_path: Path) -> None:
    _touch(tmp_path, "selective_pressure/OG1_absrel.json")
    r = _run(tmp_path, """
declare_stage_output 05 'aBSREL results' 1..N 'glob:${RESULTS_DIR}/selective_pressure/*_absrel.json'
verify_stage_outputs 05 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_OK" in r.stdout, r.stdout + r.stderr


def test_one_or_more_fails_loudly_when_nothing_produced(tmp_path: Path) -> None:
    """This is the exact stage-05 failure: zero outputs, flag written anyway."""
    r = _run(tmp_path, """
declare_stage_output 05 'aBSREL results' 1..N 'glob:${RESULTS_DIR}/selective_pressure/*_absrel.json'
verify_stage_outputs 05 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_FAIL" in r.stdout, r.stdout + r.stderr
    combined = r.stdout + r.stderr
    assert "ERROR" in combined
    assert "05" in combined, "the failure must name the stage"
    assert "aBSREL results" in combined, "the failure must name the unmet expectation"


# --------------------------------------------------------------------------
# cardinality: 0..N tolerates legitimate emptiness
# --------------------------------------------------------------------------
def test_zero_or_more_tolerates_zero(tmp_path: Path) -> None:
    r = _run(tmp_path, """
declare_stage_output 06 'tandem clusters' 0..N 'glob:${RESULTS_DIR}/synteny/*_tandem.tsv'
verify_stage_outputs 06 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_OK" in r.stdout, r.stdout + r.stderr
    assert "ERROR" not in r.stdout + r.stderr


def test_zero_or_more_and_one_or_more_differ_on_the_same_empty_tree(tmp_path: Path) -> None:
    """The discriminator between 'legitimately empty' and 'silently broken'."""
    body = """
declare_stage_output S 'out' {card} 'glob:${{RESULTS_DIR}}/nothing/*.tsv'
verify_stage_outputs S && echo VERIFY_OK || echo VERIFY_FAIL
"""
    lenient = _run(tmp_path, body.format(card="0..N"))
    strict = _run(tmp_path, body.format(card="1..N"))
    assert "VERIFY_OK" in lenient.stdout
    assert "VERIFY_FAIL" in strict.stdout


# --------------------------------------------------------------------------
# cardinality: floors above one, and exact counts
# --------------------------------------------------------------------------
def test_floor_above_one_fails_below_it(tmp_path: Path) -> None:
    _touch(tmp_path, "phylogenies/a.treefile")
    r = _run(tmp_path, """
declare_stage_output 04 'class trees' 2..N 'glob:${RESULTS_DIR}/phylogenies/*.treefile'
verify_stage_outputs 04 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_FAIL" in r.stdout
    assert "2" in r.stdout + r.stderr


def test_exact_cardinality_passes_on_the_nose(tmp_path: Path) -> None:
    """`k..k` lets a stage assert one-output-per-input from a RUNTIME count."""
    for i in range(3):
        _touch(tmp_path, f"phylogenies/og{i}.treefile")
    r = _run(tmp_path, """
n=3
declare_stage_output 04 'one tree per OG' "${n}..${n}" 'glob:${RESULTS_DIR}/phylogenies/*.treefile'
verify_stage_outputs 04 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_OK" in r.stdout, r.stdout + r.stderr


def test_exact_cardinality_fails_when_over(tmp_path: Path) -> None:
    for i in range(4):
        _touch(tmp_path, f"phylogenies/og{i}.treefile")
    r = _run(tmp_path, """
declare_stage_output 04 'one tree per OG' 3..3 'glob:${RESULTS_DIR}/phylogenies/*.treefile'
verify_stage_outputs 04 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_FAIL" in r.stdout


def test_malformed_cardinality_is_rejected_loudly(tmp_path: Path) -> None:
    r = _run(tmp_path, """
declare_stage_output 04 'x' lots 'glob:${RESULTS_DIR}/*' && echo DECL_OK || echo DECL_FAIL
""")
    assert "DECL_FAIL" in r.stdout
    assert "ERROR" in r.stdout + r.stderr


# --------------------------------------------------------------------------
# spec kinds
# --------------------------------------------------------------------------
def test_path_spec_counts_a_single_named_file(tmp_path: Path) -> None:
    _touch(tmp_path, "ranking/ranked.csv")
    r = _run(tmp_path, """
declare_stage_output 07 'ranked csv' 1..N 'path:${RESULTS_DIR}/ranking/ranked.csv'
verify_stage_outputs 07 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_OK" in r.stdout, r.stdout + r.stderr


def test_empty_file_does_not_satisfy_a_declaration(tmp_path: Path) -> None:
    """A zero-byte output is the failure mode, not evidence of success."""
    _touch(tmp_path, "ranking/ranked.csv", content="")
    r = _run(tmp_path, """
declare_stage_output 07 'ranked csv' 1..N 'path:${RESULTS_DIR}/ranking/ranked.csv'
verify_stage_outputs 07 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_FAIL" in r.stdout


def test_og_product_spec_needs_no_literal_path(tmp_path: Path) -> None:
    """The declaration names the RELATIONSHIP; the resolver supplies location."""
    d = tmp_path / "results" / "phylogenies" / "protein" / "class_A" / "OG0000123"
    d.mkdir(parents=True)
    (d / "OG0000123.treefile").write_text("(s1,s2);\n")
    r = _run(tmp_path, """
declare_stage_output 04 'per-OG trees' 1..N 'og_product:tree'
verify_stage_outputs 04 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_OK" in r.stdout, r.stdout + r.stderr


def test_og_product_spec_counts_across_every_class(tmp_path: Path) -> None:
    for cls, base in (("A", "OG1"), ("C", "OG2"), ("unclassified", "OG3")):
        d = tmp_path / "results" / "phylogenies" / "protein" / f"class_{cls}" / base
        d.mkdir(parents=True)
        (d / f"{base}.treefile").write_text("(s1,s2);\n")
    r = _run(tmp_path, """
declare_stage_output 04 'per-OG trees' 3..3 'og_product:tree'
verify_stage_outputs 04 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_OK" in r.stdout, r.stdout + r.stderr


def test_og_product_spec_ignores_original_treefiles(tmp_path: Path) -> None:
    d = tmp_path / "results" / "phylogenies" / "protein" / "class_A" / "OG1"
    d.mkdir(parents=True)
    (d / "OG1.original.treefile").write_text("(s1,s2);\n")
    r = _run(tmp_path, """
declare_stage_output 04 'per-OG trees' 1..N 'og_product:tree'
verify_stage_outputs 04 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_FAIL" in r.stdout


def test_unknown_spec_kind_is_rejected(tmp_path: Path) -> None:
    r = _run(tmp_path, """
declare_stage_output 04 'x' 1..N 'telepathy:tree' && echo DECL_OK || echo DECL_FAIL
""")
    assert "DECL_FAIL" in r.stdout
    assert "telepathy" in r.stdout + r.stderr


# --------------------------------------------------------------------------
# multiple declarations
# --------------------------------------------------------------------------
def test_all_declarations_must_pass(tmp_path: Path) -> None:
    _touch(tmp_path, "a/one.txt")
    r = _run(tmp_path, """
declare_stage_output 05 'first'  1..N 'glob:${RESULTS_DIR}/a/*.txt'
declare_stage_output 05 'second' 1..N 'glob:${RESULTS_DIR}/b/*.txt'
verify_stage_outputs 05 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    assert "VERIFY_FAIL" in r.stdout
    assert "second" in r.stdout + r.stderr


def test_every_unmet_declaration_is_reported_not_just_the_first(tmp_path: Path) -> None:
    """One run should surface the whole truth, not one failure per re-run."""
    r = _run(tmp_path, """
declare_stage_output 05 'alpha' 1..N 'glob:${RESULTS_DIR}/a/*.txt'
declare_stage_output 05 'beta'  1..N 'glob:${RESULTS_DIR}/b/*.txt'
verify_stage_outputs 05 && echo VERIFY_OK || echo VERIFY_FAIL
""")
    combined = r.stdout + r.stderr
    assert "alpha" in combined and "beta" in combined


def test_declarations_are_scoped_per_stage(tmp_path: Path) -> None:
    """Stage 06's unmet declaration must not sink stage 05's verification."""
    _touch(tmp_path, "a/one.txt")
    r = _run(tmp_path, """
declare_stage_output 05 'mine'    1..N 'glob:${RESULTS_DIR}/a/*.txt'
declare_stage_output 06 'theirs'  1..N 'glob:${RESULTS_DIR}/b/*.txt'
verify_stage_outputs 05 && echo FIVE_OK || echo FIVE_FAIL
verify_stage_outputs 06 && echo SIX_OK || echo SIX_FAIL
""")
    assert "FIVE_OK" in r.stdout, r.stdout + r.stderr
    assert "SIX_FAIL" in r.stdout


# --------------------------------------------------------------------------
# complete_stage -- the ONLY writer of the completion flag
# --------------------------------------------------------------------------
def test_complete_stage_writes_the_flag_when_verified(tmp_path: Path) -> None:
    _touch(tmp_path, "a/one.txt")
    r = _run(tmp_path, """
declare_stage_output 05 'out' 1..N 'glob:${RESULTS_DIR}/a/*.txt'
complete_stage 05 && echo DONE_OK || echo DONE_FAIL
""")
    assert "DONE_OK" in r.stdout, r.stdout + r.stderr
    assert _flag(tmp_path, "05").exists(), "verified stage must be flagged completed"


def test_complete_stage_refuses_the_flag_when_unverified(tmp_path: Path) -> None:
    """The regression that started all of this: flag written despite no output."""
    r = _run(tmp_path, """
declare_stage_output 05 'out' 1..N 'glob:${RESULTS_DIR}/a/*.txt'
complete_stage 05 && echo DONE_OK || echo DONE_FAIL
""")
    assert "DONE_FAIL" in r.stdout, r.stdout + r.stderr
    assert not _flag(tmp_path, "05").exists(), (
        "an unverified stage must NOT be flagged completed -- this is the bug"
    )


def test_complete_stage_refuses_a_stage_that_declared_nothing(tmp_path: Path) -> None:
    """No declaration means no assertion; silently flagging is what broke 05."""
    r = _run(tmp_path, 'complete_stage 05 && echo DONE_OK || echo DONE_FAIL')
    assert "DONE_FAIL" in r.stdout, r.stdout + r.stderr
    assert not _flag(tmp_path, "05").exists()
    assert "ERROR" in r.stdout + r.stderr


def test_completed_stage_is_visible_to_the_existing_is_step_completed(tmp_path: Path) -> None:
    """Propagates 06c's contract rather than inventing a parallel one."""
    _touch(tmp_path, "a/one.txt")
    r = _run(tmp_path, """
declare_stage_output 05 'out' 1..N 'glob:${RESULTS_DIR}/a/*.txt'
complete_stage 05
is_step_completed 05 && echo SEEN_YES || echo SEEN_NO
""")
    assert "SEEN_YES" in r.stdout, r.stdout + r.stderr


def test_complete_stage_works_with_a_named_substep(tmp_path: Path) -> None:
    """06c's `is_step_completed "06c_classify_hmm"` granularity must survive."""
    _touch(tmp_path, "classification/hmm_hits.tsv")
    r = _run(tmp_path, """
declare_stage_output 06c_classify_hmm 'hmm hits' 1..N 'glob:${RESULTS_DIR}/classification/*.tsv'
complete_stage 06c_classify_hmm
is_step_completed 06c_classify_hmm && echo SEEN_YES || echo SEEN_NO
""")
    assert "SEEN_YES" in r.stdout, r.stdout + r.stderr


def test_whole_contract_runs_with_only_functions_sh(tmp_path: Path) -> None:
    """A bare `sbatch 05_...sh` sources functions.sh -- no orchestrator exists."""
    _touch(tmp_path, "a/one.txt")
    r = _run(tmp_path, """
declare_stage_output 05 'out' 1..N 'glob:${RESULTS_DIR}/a/*.txt'
complete_stage 05 && echo DONE_OK || echo DONE_FAIL
""")
    assert "DONE_OK" in r.stdout
    assert "command not found" not in (r.stdout + r.stderr).lower()
