"""ONE shared resolver for stage 04's per-class / per-orthogroup output layout.

Stage 04 writes two different shapes under ``${RESULTS_DIR}/phylogenies/protein``::

    per-class (global)  class_<C>/class_<C>.treefile      + class_<C>_trimmed.fa
    per-orthogroup      class_<C>/<base>/<base>.treefile  + <base>_trimmed.fa

Three consumers re-derived that location independently and ALL THREE got it
wrong against the post-refactor layout:

  * ``05_selective_pressure_and_asr.sh`` read the pre-refactor FLAT pair, so
    the ``[ -f ... ]`` guard was never true and the whole stage silently
    produced nothing while still writing its completion flag (bead F1).
  * ``03d_notung_reconciliation.sh`` looked up alignments at the same dead flat
    paths, so ``families.txt`` was left as a bare ``[FAMILIES]`` header and
    GeneRax reconciled nothing (bead F6).
  * ``scripts/rank_candidates.py`` hardcodes ``class_A/class_A.treefile``.

The council's discriminator is that machinery storing literal path VALUES
drifts exactly the way those three copies drifted, whereas a resolver encoding
the RELATIONSHIP (class comes from the class map; product name is a function of
the base) cannot. These tests pin the relationship.

The resolver lives in ``functions.sh`` because stages are submitted with a bare
``sbatch``, not only through the orchestrator, so every test here sources
``functions.sh`` ALONE -- no ``config.sh``, no ``run_pipeline.sh``.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"


# --------------------------------------------------------------------------
# harness -- functions.sh sourced alone, against a synthetic ${RESULTS_DIR}
# --------------------------------------------------------------------------
def _run(tmp_path: Path, body: str) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    logs = tmp_path / "logs"
    logs.mkdir(exist_ok=True)
    env["LOGS_DIR"] = str(logs)
    script = f"""
source "{FUNCTIONS_SH}"; trap - EXIT
RESULTS_DIR="{tmp_path / 'results'}"
{body}
"""
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _make_perog(tmp_path: Path, cls: str, base: str, *, align=True, tree=True) -> Path:
    d = tmp_path / "results" / "phylogenies" / "protein" / f"class_{cls}" / base
    d.mkdir(parents=True, exist_ok=True)
    if align:
        (d / f"{base}_trimmed.fa").write_text(">s1\nMKFLV\n>s2\nMKFLA\n")
    if tree:
        (d / f"{base}.treefile").write_text("(s1:0.1,s2:0.2);\n")
    return d


def _make_perclass(tmp_path: Path, cls: str) -> Path:
    d = tmp_path / "results" / "phylogenies" / "protein" / f"class_{cls}"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"class_{cls}.treefile").write_text("(s1:0.1,s2:0.2);\n")
    (d / f"class_{cls}_trimmed.fa").write_text(">s1\nMKFLV\n")
    return d


def _write_class_map(tmp_path: Path, rows: dict[str, str]) -> Path:
    p = tmp_path / "results" / "classification" / "og_class_majority.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("orthogroup\tclass\n" + "".join(f"{og}\t{c}\n" for og, c in rows.items()))
    return p


# --------------------------------------------------------------------------
# resolve_og_class -- the class is READ, never hardcoded
# --------------------------------------------------------------------------
def test_class_comes_from_the_class_map(tmp_path: Path) -> None:
    _write_class_map(tmp_path, {"OG0000123": "C"})
    r = _run(tmp_path, 'resolve_og_class OG0000123')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().splitlines()[-1] == "C"


def test_class_falls_back_to_unclassified_when_og_absent_from_map(tmp_path: Path) -> None:
    """Stage 04 routes OGs missing from the map to class_unclassified."""
    _write_class_map(tmp_path, {"OG0000999": "A"})
    r = _run(tmp_path, 'resolve_og_class OG0000123')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().splitlines()[-1] == "unclassified"


def test_class_falls_back_to_A_when_map_is_absent(tmp_path: Path) -> None:
    """Back-compat: runs predating per-OG class assignment routed to class_A."""
    r = _run(tmp_path, 'resolve_og_class OG0000123')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().splitlines()[-1] == "A"


# --------------------------------------------------------------------------
# resolve_og_dir -- derives location without requiring the products to exist
# --------------------------------------------------------------------------
def test_dir_is_derived_for_a_perog_base(tmp_path: Path) -> None:
    _write_class_map(tmp_path, {"OG0000123": "C"})
    r = _run(tmp_path, 'resolve_og_dir OG0000123')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("phylogenies/protein/class_C/OG0000123")


def test_dir_for_a_perclass_base_is_not_nested(tmp_path: Path) -> None:
    """`class_A` is stage 04's GLOBAL output; it must not become class_A/class_A."""
    r = _run(tmp_path, 'resolve_og_dir class_A')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("phylogenies/protein/class_A")


# --------------------------------------------------------------------------
# resolve_og_product -- the drift-prone lookup, in one place
# --------------------------------------------------------------------------
def test_product_found_at_the_class_the_map_names(tmp_path: Path) -> None:
    _write_class_map(tmp_path, {"OG0000123": "C"})
    _make_perog(tmp_path, "C", "OG0000123")
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("class_C/OG0000123/OG0000123.treefile")


def test_alignment_product_uses_the_trimmed_suffix(tmp_path: Path) -> None:
    _write_class_map(tmp_path, {"OG0000123": "A"})
    _make_perog(tmp_path, "A", "OG0000123")
    r = _run(tmp_path, 'resolve_og_product OG0000123 alignment')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("class_A/OG0000123/OG0000123_trimmed.fa")


def test_product_found_under_the_unclassified_fallback(tmp_path: Path) -> None:
    _write_class_map(tmp_path, {"OG0000999": "A"})
    _make_perog(tmp_path, "unclassified", "OG0000123")
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("class_unclassified/OG0000123/OG0000123.treefile")


def test_perclass_global_product_resolves(tmp_path: Path) -> None:
    _make_perclass(tmp_path, "A")
    r = _run(tmp_path, 'resolve_og_product class_A tree')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("class_A/class_A.treefile")


def test_product_elsewhere_is_found_but_warns(tmp_path: Path) -> None:
    """Map says C, stage 04 actually wrote under A: use it, but say so LOUDLY.

    This is the stage-04/stage-05 disagreement that happens when the class map
    is regenerated after stage 04 already ran under the back-compat fallback.
    """
    _write_class_map(tmp_path, {"OG0000123": "C"})
    _make_perog(tmp_path, "A", "OG0000123")
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("class_A/OG0000123/OG0000123.treefile")
    assert "WARN" in r.stdout + r.stderr


def test_same_og_under_two_classes_refuses_to_guess(tmp_path: Path) -> None:
    """Choosing one would silently attach results to the wrong tree."""
    _write_class_map(tmp_path, {"OG0000123": "C"})
    _make_perog(tmp_path, "A", "OG0000123")
    _make_perog(tmp_path, "B", "OG0000123")
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 2, f"expected ambiguity rc=2, got {r.returncode}"
    combined = r.stdout + r.stderr
    assert "ERROR" in combined
    assert "OG0000123" in combined


def test_missing_product_is_quiet_and_not_an_error(tmp_path: Path) -> None:
    """Stage 04 legitimately builds no tree for some OGs -- must not cry wolf."""
    _write_class_map(tmp_path, {"OG0000123": "A"})
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 1, f"expected not-found rc=1, got {r.returncode}"
    assert "ERROR" not in r.stdout + r.stderr


def test_original_treefile_is_not_mistaken_for_a_product(tmp_path: Path) -> None:
    """TreeShrink's `<base>.original.treefile` is an INPUT copy, not an output.

    03d had to special-case this; the shared resolver must handle it so the
    next caller does not rediscover it.
    """
    _write_class_map(tmp_path, {"OG0000123": "A"})
    d = tmp_path / "results" / "phylogenies" / "protein" / "class_A" / "OG0000123"
    d.mkdir(parents=True)
    (d / "OG0000123.original.treefile").write_text("(s1,s2);\n")
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 1, "the .original.treefile copy must not satisfy `tree`"


def test_empty_product_file_does_not_count_as_produced(tmp_path: Path) -> None:
    """A zero-byte treefile IS the 'produced nothing' failure, not a success."""
    _write_class_map(tmp_path, {"OG0000123": "A"})
    d = tmp_path / "results" / "phylogenies" / "protein" / "class_A" / "OG0000123"
    d.mkdir(parents=True)
    (d / "OG0000123.treefile").write_text("")
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 1


def test_unknown_product_kind_is_a_usage_error(tmp_path: Path) -> None:
    r = _run(tmp_path, 'resolve_og_product OG0000123 sandwich')
    assert r.returncode == 3, f"expected usage rc=3, got {r.returncode}"
    assert "sandwich" in r.stdout + r.stderr


def test_resolver_needs_no_orchestrator(tmp_path: Path) -> None:
    """A bare `sbatch 05_...sh` sources functions.sh and nothing else."""
    _write_class_map(tmp_path, {"OG0000123": "A"})
    _make_perog(tmp_path, "A", "OG0000123")
    r = _run(tmp_path, 'resolve_og_product OG0000123 tree')
    assert r.returncode == 0
    assert "run_pipeline" not in r.stderr


# --- stdout is the resolver's RETURN CHANNEL, not a log sink ---------------
#
# `log()` routes INFO and WARN to STDOUT (only ERROR goes to stderr). The
# resolver returns its answer by `echo`-ing a path, so any diagnostic it logs
# lands in the caller's command substitution and corrupts the value.
#
# The class-mismatch branch is the dangerous one: it returns 0 (SUCCESS) and
# echoes the path, but logs a WARN first -- so `p=$(resolve_og_product ...)`
# captures "<warn line>\n<path>", `[ -f "$p" ]` fails, and the caller silently
# skips the orthogroup. That is precisely the silent-attrition bug this whole
# mechanism exists to prevent, occurring on the branch meant to RESCUE a
# class-map mismatch.


def _resolve_capture(tmp_path, base, kind, layout, class_map):
    """Call resolve_og_product exactly as a real caller would: $( )."""
    results = tmp_path / "results"
    (results / "classification").mkdir(parents=True, exist_ok=True)
    rows = ["orthogroup\tclass"] + [f"{og}\t{cls}" for og, cls in class_map.items()]
    (results / "classification" / "og_class_majority.tsv").write_text("\n".join(rows) + "\n")
    for rel, content in layout.items():
        p = results / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content)
    script = (
        f'source "{FUNCTIONS_SH}"; trap - EXIT; '
        f'p=$(resolve_og_product "{base}" "{kind}"); rc=$?; '
        f'printf "%s" "$p"; exit $rc'
    )
    env = os.environ.copy()
    env["RESULTS_DIR"] = str(results)
    env["LOGS_DIR"] = str(tmp_path)
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def test_class_mismatch_returns_a_usable_path_not_a_log_line(tmp_path):
    """Success path: the captured value must be openable by the caller."""
    res = _resolve_capture(
        tmp_path, "OGY", "tree",
        {"phylogenies/protein/class_C/OGY/OGY.treefile": "(A:0.1,B:0.1);\n"},
        {"OGY": "A"},
    )
    assert res.returncode == 0, res.stderr
    captured = res.stdout
    assert "\n" not in captured.strip(), (
        "resolver emitted more than one line on stdout; a diagnostic leaked into "
        f"its return channel: {captured!r}"
    )
    assert "WARN" not in captured and "[20" not in captured, (
        f"a log line was captured as the path: {captured!r}"
    )
    assert captured.strip().endswith("OGY.treefile"), captured


def test_not_found_returns_an_empty_value_not_a_log_line(tmp_path):
    """Failure path: a caller ignoring rc must get '', never a plausible string."""
    res = _resolve_capture(tmp_path, "OGX", "tree", {}, {"OGX": "A"})
    assert res.returncode == 1
    assert res.stdout.strip() == "", (
        "not-found leaked a log line into the return channel, so a caller that "
        f"does not check rc gets a non-empty non-path: {res.stdout!r}"
    )


def test_resolver_diagnostics_still_reach_the_operator_on_stderr(tmp_path):
    """Redirecting to stderr must not silence the diagnostic."""
    res = _resolve_capture(
        tmp_path, "OGY", "tree",
        {"phylogenies/protein/class_C/OGY/OGY.treefile": "(A:0.1,B:0.1);\n"},
        {"OGY": "A"},
    )
    assert "disagree" in res.stderr, (
        f"the class-mismatch warning must still be visible on stderr: {res.stderr!r}"
    )
