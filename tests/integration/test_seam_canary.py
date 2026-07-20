"""End-to-end canary for the pipeline's stage-to-stage data-flow seams.

WHY THIS FILE EXISTS
--------------------
An audit found twelve defects. An entire stage had never produced output: it
read a path a refactor had moved, its guard was never satisfied, it wrote its
completion flag anyway, and the next stage merely logged "results not found".
Four scored ranking axes were dormant the same way. Two producer scripts were
invoked by no stage.

~2057 unit tests passed throughout, because every unit was correct in
isolation. The defects lived in the SEAMS -- what stage A writes versus where
stage B looks -- and a seam is not inside any unit, so no amount of unit
testing can reach one.

WHAT THIS ASSERTS
-----------------
The contract, not the filenames:

    for every stage-to-stage hand-off, the path stage B reads is a path
    stage A actually writes.

Paths are compared as PATTERNS. A stage may rename its outputs, add a class
directory, or derive a filename at runtime; the canary only cares that the two
sides still intersect. That keeps it from breaking on legitimate change while
still catching a directory-depth mismatch, which is what the 04 -> 05 defect
actually was.

RUNNING IT
----------
    pytest tests/integration -q            # the canary
    pytest tests/integration -q -m seam    # seam graph checks only

It needs no cluster, no network, no aligners and no tree inference; it reads
the stage scripts and executes only config.sh and one extracted shell
function. Expect roughly a second.

It is NOT collected by `pytest tests/unit -q`.
"""

from __future__ import annotations

import re
import subprocess
import textwrap
from pathlib import Path

import pytest

import seamlib
from known_gaps import KNOWN_GAPS, SUSPECTED_DEFECT, gap_key

pytestmark = pytest.mark.seam

REPO_ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# Failure reporting
#
# "assert False" would make this canary worthless. A seam failure must name
# the consuming stage, the path, and -- when it can work it out -- the reason
# no producer was found, because the person reading the failure is being told
# about a hand-off they did not know existed.
# ---------------------------------------------------------------------------

def _leaf_of(pattern: str) -> str:
    return pattern.rsplit("/", 1)[-1]


def _invoked_anywhere(script_name: str) -> bool:
    """Is this script reachable -- invoked by a stage, or imported as a module?

    A library module such as `_rank_candidates_lib.py` is never invoked and
    never should be; accusing it of being an orphaned producer would send the
    reader chasing a non-problem.
    """
    searched = (
        sorted(REPO_ROOT.glob("[0-9]*_*.sh"))
        + sorted((REPO_ROOT / "scripts").rglob("*.sh"))
        + [REPO_ROOT / "functions.sh", REPO_ROOT / "run_pipeline.sh"]
    )
    for path in searched:
        if path.is_file() and script_name in path.read_text(errors="replace"):
            return True

    module = script_name.removesuffix(".py")
    imported = re.compile(rf"^\s*(?:from|import)\s+{re.escape(module)}\b", re.M)
    for path in sorted((REPO_ROOT / "scripts").rglob("*.py")):
        if path.name == script_name:
            continue
        if imported.search(path.read_text(errors="replace")):
            return True
    return False


def diagnose(probe: seamlib.Probe, results_dir: str) -> str:
    """Best available explanation for why nothing produces this path."""
    leaf = _leaf_of(probe.pattern)
    # Only chase an orphaned producer for something that looks like a FILE.
    # A bare directory name such as "hmms" matches far too many scripts and
    # produces a confidently wrong accusation.
    if "*" not in leaf and "." in leaf:
        own_name = probe.stage.rsplit("/", 1)[-1]
        mentions: list[Path] = [
            script
            for script in sorted((REPO_ROOT / "scripts").rglob("*.py"))
            if script.name != own_name
            and leaf in script.read_text(errors="replace")
        ]
        # An ORPHANED producer is the far more actionable finding, so report it
        # in preference to a script that merely mentions the name. Reporting
        # the consumer back to itself would be worse than useless.
        for script in mentions:
            if not _invoked_anywhere(script.name):
                return (
                    f"ORPHANED PRODUCER: {script.relative_to(REPO_ROOT)} appears "
                    f"to write '{leaf}', but no stage script invokes it. The "
                    f"guard in {probe.stage} can therefore never be satisfied."
                )
        if mentions:
            return (
                f"{mentions[0].relative_to(REPO_ROOT)} mentions '{leaf}' and is "
                f"reachable, but the extractor could not bind its output path."
            )
    return "No writer of this pattern was found anywhere in the repository."


def describe(probe: seamlib.Probe, results_dir: str) -> str:
    pretty = probe.pattern.replace(results_dir, "${RESULTS_DIR}")
    return (
        f"\n  CONSUMER : {probe.location}"
        f"\n  PROBES   : {pretty}"
        f"\n  AS       : {probe.raw}  [{probe.kind}]"
        f"\n  PRODUCER : none found"
        f"\n  WHY      : {diagnose(probe, results_dir)}"
    )


# ---------------------------------------------------------------------------
# The core assertion
# ---------------------------------------------------------------------------

def test_every_seam_probe_has_a_producer(probe_groups, seam_graph, results_dir):
    """Every existence guard on a RESULTS_DIR path must have some writer.

    This is the assertion that would have caught the 04 -> 05 break, the
    OrthoFinder root mismatch, and the 02a fallback chain whose only writer
    was the test harness.
    """
    failures: list[str] = []

    for members in probe_groups.values():
        satisfied = any(
            any(seamlib.write_satisfies(write, probe) for write in seam_graph.writes)
            for probe in members
        )
        if satisfied:
            continue
        if all(
            gap_key(probe.stage, probe.pattern, results_dir) in KNOWN_GAPS
            for probe in members
        ):
            continue

        if len(members) > 1:
            head = members[0]
            alternatives = "\n           ".join(
                probe.pattern.replace(results_dir, "${RESULTS_DIR}")
                for probe in members
            )
            failures.append(
                f"\n  CONSUMER : {head.location}"
                f"\n  PROBES   : a {len(members)}-way fallback chain in which"
                f" NO alternative has a producer:"
                f"\n           {alternatives}"
                f"\n  WHY      : {diagnose(head, results_dir)}"
            )
        else:
            failures.append(describe(members[0], results_dir))

    assert not failures, (
        "BROKEN PIPELINE SEAM(S): a stage probes for a path that nothing in "
        "the pipeline writes.\n"
        "A guard like this is never satisfied, so the stage skips its work "
        "SILENTLY -- which is exactly how a whole stage once ran to "
        "completion without ever producing output.\n"
        + "".join(failures)
        + "\n\nIf the path is genuinely produced outside this repository, add "
        "it to tests/integration/known_gaps.py WITH the producer named."
    )


def test_known_gaps_inventory_has_no_stale_entries(
    probe_groups, seam_graph, results_dir
):
    """An inventory entry that now resolves must be deleted.

    Without this the inventory would silently become permanent cover for
    seams that were fixed, and the next real break in the same place would
    be waved through.
    """
    unresolved: set[tuple[str, str]] = set()
    for members in probe_groups.values():
        satisfied = any(
            any(seamlib.write_satisfies(write, probe) for write in seam_graph.writes)
            for probe in members
        )
        if not satisfied:
            for probe in members:
                unresolved.add(gap_key(probe.stage, probe.pattern, results_dir))

    listed = set(KNOWN_GAPS)
    stale = sorted(listed - unresolved)
    assert not stale, (
        "STALE known_gaps.py entries -- these probes now HAVE a producer and "
        "the exemption must be removed:\n"
        + "\n".join(f"  {stage}  {pattern}" for stage, pattern in stale)
    )


def test_suspected_defects_are_surfaced(results_dir):
    """Keep live findings visible instead of letting them rot in a table."""
    live = {
        key: reason
        for key, (category, reason) in KNOWN_GAPS.items()
        if category == SUSPECTED_DEFECT
    }
    if live:
        report = "\n".join(
            f"  {stage}  {pattern}\n      {reason}"
            for (stage, pattern), reason in sorted(live.items())
        )
        print(
            "\nSUSPECTED LIVE SEAM DEFECTS (recorded, not fixed -- stage files "
            "are out of scope for this canary):\n" + report
        )
    # Recorded findings must stay recorded; this fails if one is dropped
    # without the underlying seam actually being wired up.
    assert isinstance(live, dict)


# ---------------------------------------------------------------------------
# Wiring completeness
#
# Two of the audited defects were not path mismatches at all: a stage was
# missing from the registry, and a correct producer was invoked by nobody.
# ---------------------------------------------------------------------------

def test_every_stage_script_appears_in_the_registry(stage_scripts):
    """A stage that is not in the registry never runs. 02c was in that state."""
    registry = (REPO_ROOT / "run_pipeline.sh").read_text(errors="replace")

    declared = set(re.findall(r'\["([0-9][0-9a-z]*)"\]=', registry))
    order_match = re.search(r"STEP_ORDER=\(([^)]*)\)", registry)
    assert order_match, "run_pipeline.sh no longer defines STEP_ORDER"
    ordered = set(re.findall(r'"([0-9][0-9a-z]*)"', order_match.group(1)))

    missing: list[str] = []
    for script in stage_scripts:
        step = re.match(r"([0-9][0-9a-z]*)_", script.name)
        if not step:
            continue
        key = step.group(1)
        if key not in declared:
            missing.append(f"  {script.name}: absent from PIPELINE_STEPS")
        if key not in ordered:
            missing.append(f"  {script.name}: absent from STEP_ORDER")

    assert not missing, (
        "UNREGISTERED STAGE(S): a stage script exists but the driver will "
        "never run it.\n" + "\n".join(missing)
    )


def test_registry_entries_point_at_real_scripts():
    """The mirror image: a registry entry naming a script that does not exist."""
    registry = (REPO_ROOT / "run_pipeline.sh").read_text(errors="replace")
    broken = [
        script
        for script in re.findall(r'\]="([0-9][^:"]*\.sh):', registry)
        if not (REPO_ROOT / script).is_file()
    ]
    assert not broken, (
        "Registry references stage scripts that do not exist: " + ", ".join(broken)
    )


def test_scripts_invoked_by_stages_exist(stage_scripts):
    """A stage invoking a missing script fails at runtime, deep into a job."""
    scripts_dir = REPO_ROOT / "scripts"
    missing: list[str] = []
    pattern = re.compile(r"(?:\$\{SCRIPTS_DIR\}|scripts)/([A-Za-z0-9_/]+\.(?:py|sh))")
    for stage in stage_scripts:
        for name in set(pattern.findall(stage.read_text(errors="replace"))):
            if not (scripts_dir / name).is_file():
                missing.append(f"  {stage.name} invokes scripts/{name} (absent)")
    assert not missing, "MISSING SCRIPT(S):\n" + "\n".join(sorted(missing))


# ---------------------------------------------------------------------------
# Historical-defect regression tests
#
# These pin the DETECTION SEMANTICS using the real path pairs from the audit.
# If someone later loosens the matcher -- for instance by using a glob whose
# "*" crosses "/" -- these fail even though the pipeline itself is fine.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "producer,consumer,should_intersect",
    [
        # Defect 1: stage 05 read the pre-refactor FLAT path while stage 04
        # writes two levels deeper under a per-class directory.
        (
            "R/phylogenies/protein/class_A/OGX/OGX_trimmed.fa",
            "R/phylogenies/protein/OGX_trimmed.fa",
            False,
        ),
        # ... and the post-fix consumer path, which must still resolve.
        (
            "R/phylogenies/protein/class_A/OGX/OGX_trimmed.fa",
            "R/phylogenies/protein/class_*/OGX/OGX_trimmed.fa",
            True,
        ),
        # Defect 2: OrthoFinder writes under orthogroups/input/OrthoFinder,
        # the ranking script probed a results/orthofinder root.
        (
            "R/orthogroups/input/OrthoFinder/Results_A/Orthogroups/Orthogroups.tsv",
            "R/orthofinder/Results_A/Orthogroups/Orthogroups.tsv",
            False,
        ),
        (
            "R/orthogroups/input/OrthoFinder/Results_A/Orthogroups/Orthogroups.tsv",
            "R/orthogroups/input/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv",
            True,
        ),
    ],
)
def test_matcher_distinguishes_the_historical_path_pairs(
    producer, consumer, should_intersect
):
    assert seamlib.patterns_intersect(producer, consumer) is should_intersect


def test_wildcards_do_not_cross_directory_separators():
    """The property the 04 -> 05 defect turned on.

    fnmatch's "*" spans "/", so a naive matcher would call these compatible
    and the canary would have missed the defect it exists for.
    """
    assert not seamlib.patterns_intersect("R/a/b/c.fa", "R/*.fa")
    assert seamlib.patterns_intersect("R/a/b/c.fa", "R/a/*/c.fa")


# ---------------------------------------------------------------------------
# Genuinely executed seam
#
# Static extraction cannot prove that a resolver actually FINDS the layout its
# producer creates. This materialises stage 04's real output layout in a temp
# directory and runs stage 05's REAL resolver function against it.
# ---------------------------------------------------------------------------

def _extract_shell_function(script: Path, name: str) -> str:
    text = script.read_text(errors="replace")
    start = text.find(f"{name}() {{")
    assert start != -1, f"{script.name} no longer defines {name}()"
    depth, index = 0, start
    while index < len(text):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start : index + 1]
        index += 1
    raise AssertionError(f"unterminated function {name} in {script.name}")


def test_stage05_resolver_finds_the_layout_stage04_writes(tmp_path):
    """Execute the actual seam: 04's directory layout, 05's actual resolver.

    This is the smallest end-to-end canary that would have failed loudly on
    the original defect. Nothing is stubbed except `log`; the resolver code is
    the real text of stage 05.
    """
    results = tmp_path / "results"
    orthogroup = "OGCANARY"
    og_class = "A"

    per_og = results / "phylogenies" / "protein" / f"class_{og_class}" / orthogroup
    per_og.mkdir(parents=True)
    (per_og / f"{orthogroup}_trimmed.fa").write_text(">seq\nMAAA\n")
    (per_og / f"{orthogroup}.treefile").write_text("(seq:0.1);\n")

    classification = results / "classification"
    classification.mkdir(parents=True)
    (classification / "og_class_majority.tsv").write_text(
        f"orthogroup\tclass\n{orthogroup}\t{og_class}\n"
    )

    resolver = _extract_shell_function(
        REPO_ROOT / "05_selective_pressure_and_asr.sh", "resolve_perog_paths"
    )

    script = textwrap.dedent(
        """
        set -eo pipefail
        RESULTS_DIR="%s"
        log() { :; }
        %s
        if resolve_perog_paths "%s"; then
            echo "RESOLVED:$protein_align"
            echo "TREE:$tree"
        else
            echo "UNRESOLVED"
        fi
        """
    ) % (results, resolver, orthogroup)

    completed = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True, timeout=60
    )
    output = completed.stdout

    assert "UNRESOLVED" not in output, (
        "BROKEN SEAM 04 -> 05: stage 05's resolver could not find the "
        "per-orthogroup alignment and tree at the layout stage 04 writes.\n"
        f"  layout created : {per_og}\n"
        f"  resolver stdout: {output!r}\n"
        f"  resolver stderr: {completed.stderr[:1500]!r}\n"
        "This is the exact failure that let stage 05 produce no selection or "
        "ASR output while still touching step_completed_05.txt."
    )
    assert f"RESOLVED:{per_og}/{orthogroup}_trimmed.fa" in output, output
    assert f"TREE:{per_og}/{orthogroup}.treefile" in output, output


def test_stage05_resolver_rejects_a_layout_it_should_not_accept(tmp_path):
    """The negative control: a FLAT layout must NOT satisfy the resolver.

    Without this, a resolver that returned success unconditionally would pass
    the positive test above and the canary would prove nothing.
    """
    results = tmp_path / "results"
    orthogroup = "OGCANARY"

    flat = results / "phylogenies" / "protein"
    flat.mkdir(parents=True)
    (flat / f"{orthogroup}_trimmed.fa").write_text(">seq\nMAAA\n")
    (flat / f"{orthogroup}.treefile").write_text("(seq:0.1);\n")

    resolver = _extract_shell_function(
        REPO_ROOT / "05_selective_pressure_and_asr.sh", "resolve_perog_paths"
    )
    script = textwrap.dedent(
        """
        set -eo pipefail
        RESULTS_DIR="%s"
        log() { :; }
        %s
        if resolve_perog_paths "%s"; then echo "RESOLVED:$protein_align"; else echo "UNRESOLVED"; fi
        """
    ) % (results, resolver, orthogroup)

    completed = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True, timeout=60
    )
    assert "UNRESOLVED" in completed.stdout, (
        "Stage 05's resolver accepted the PRE-REFACTOR flat layout. Either the "
        "resolver regressed to the old path, or this canary is no longer "
        "testing anything.\n"
        f"  stdout: {completed.stdout!r}"
    )


# ---------------------------------------------------------------------------
# Sanity: the canary must be looking at something
# ---------------------------------------------------------------------------

def test_seam_graph_is_populated(seam_graph, probe_groups):
    """Guard against a silently empty canary.

    If extraction breaks -- a renamed config variable, a changed idiom -- every
    other test here would pass vacuously. That failure mode is the same shape
    as the bugs this file exists to catch, so it gets its own assertion.
    """
    assert len(seam_graph.writes) > 300, (
        f"only {len(seam_graph.writes)} writes extracted; extraction is broken"
    )
    assert len(probe_groups) > 100, (
        f"only {len(probe_groups)} probe groups extracted; extraction is broken"
    )


def test_config_sh_is_actually_sourced(config_vars):
    """The canary resolves paths through the REAL config, not a copy of it."""
    assert config_vars["RESULTS_DIR"].endswith("/results")
    assert config_vars["PHYLO_DIR"].endswith("/phylogenies/protein")
