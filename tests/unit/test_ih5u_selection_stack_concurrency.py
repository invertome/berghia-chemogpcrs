"""Bead -ih5u: run_selection_stack.sh must not corrupt cumulative CSVs
when many SLURM array tasks publish concurrently.

Stage 05 runs this script under `#SBATCH --array=0-999%50`, so up to 50
copies rebuild the cumulative results CSVs in the same shared directory at
the same time. The original code staged every rebuild through a FIXED path
``"${cum}.tmp"``: task A's `: >` truncated the staging file while task B was
mid-append, and B then `mv`-published the truncated result over the
cumulative CSV that rank_candidates.py reads.

These tests exercise the REAL script against stub `hyphy` / `clipkit` /
parser binaries, so they cover the actual production control flow rather
than a re-implementation of it. They are fast (seconds) and require no
external tools.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import textwrap
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
REAL_SCRIPT = PROJECT_ROOT / "scripts" / "hpc" / "run_selection_stack.sh"

# Fixture orthogroup labels are an explicit literal set, never minted by a
# counter: this repo treats sequential ID generation as a hard-blocked pattern
# because a re-run would silently re-point the same label at a different record.
OG_FIXTURES = [
    "OGfix_alpha", "OGfix_bravo", "OGfix_charlie", "OGfix_delta",
    "OGfix_echo", "OGfix_foxtrot", "OGfix_golf", "OGfix_hotel",
    "OGfix_india", "OGfix_juliet", "OGfix_kilo", "OGfix_lima",
]

# Header shared by the stub parser output and the pre-seeded per-OG CSVs.
BUSTED_HEADER = "og_name,p_value,is_significant"


def _w(path: Path, body: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(body).lstrip())
    if executable:
        path.chmod(0o755)


def make_project(root: Path) -> Path:
    """Build a minimal fake project root wrapping the REAL script.

    config.sh derives every path from its own location, so dropping a stub
    config.sh at ``root`` makes the script write into ``root/results``.
    Returns the path to the copied script.
    """
    _w(root / "config.sh", """
        export BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        export RESULTS_DIR="${BASE_DIR}/results"
        export SCRIPTS_DIR="${BASE_DIR}/scripts"
        export LOGS_DIR="${RESULTS_DIR}/logs"
        mkdir -p "$LOGS_DIR"
    """)
    # Stub log(): stderr only. The real log() also appends to a shared
    # pipeline.log, which is irrelevant to what these tests measure.
    _w(root / "functions.sh", """
        log() {
            local level="INFO"
            if [[ "$1" == --level=* ]]; then level="${1#--level=}"; shift; fi
            echo "[$level] $1" >&2
        }
    """)

    script = root / "scripts" / "hpc" / "run_selection_stack.sh"
    script.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(REAL_SCRIPT, script)

    # --- stub hyphy -------------------------------------------------------
    # STUB_HYPHY_FAIL / STUB_HYPHY_EMPTY hold space-delimited method names.
    _w(root / "bin" / "hyphy", """
        #!/bin/bash
        method="$1"
        out=""
        prev=""
        for a in "$@"; do
            [ "$prev" = "--output" ] && out="$a"
            prev="$a"
        done
        case " ${STUB_HYPHY_FAIL:-} " in *" $method "*) exit 1 ;; esac
        case " ${STUB_HYPHY_EMPTY:-} " in *" $method "*) : > "$out"; exit 0 ;; esac
        printf '{"method":"%s"}\\n' "$method" > "$out"
    """, executable=True)

    # --- stub clipkit -----------------------------------------------------
    # STUB_CLIPKIT_FAIL holds space-delimited -m modes to fail on.
    # STUB_CLIPKIT_IDENTICAL=1 makes both trim regimes emit the same bytes.
    _w(root / "bin" / "clipkit", """
        #!/bin/bash
        shift                       # positional input alignment
        mode=""; out=""
        while [ $# -gt 0 ]; do
            case "$1" in
                -m) mode="$2"; shift 2 ;;
                -o) out="$2"; shift 2 ;;
                *)  shift ;;
            esac
        done
        case " ${STUB_CLIPKIT_FAIL:-} " in *" $mode "*) exit 1 ;; esac
        if [ "${STUB_CLIPKIT_IDENTICAL:-0}" = "1" ]; then
            printf '>s1\\nAAACCC\\n' > "$out"
        elif [ "$mode" = "kpic-smart-gap" ]; then
            printf '>s1\\nAAA\\n' > "$out"
        else
            printf '>s1\\nAAACCC\\n' > "$out"
        fi
    """, executable=True)

    # --- stub parsers -----------------------------------------------------
    # Written atomically so the concurrency test measures the script's own
    # staging behaviour, not the parsers'.
    _w(root / "scripts" / "parse_busted.py", f"""
        import argparse, os
        p = argparse.ArgumentParser()
        p.add_argument("--json"); p.add_argument("--og-name")
        p.add_argument("--variant"); p.add_argument("--out")
        a = p.parse_args()
        tmp = a.out + ".stubtmp"
        with open(tmp, "w") as f:
            f.write("{BUSTED_HEADER}\\n" + a.og_name + ",0.01,1\\n")
        os.replace(tmp, a.out)
    """)
    _w(root / "scripts" / "parse_meme.py", """
        import argparse, os
        p = argparse.ArgumentParser()
        for flag in ("--json", "--lenient-json", "--strict-fa", "--lenient-fa",
                     "--og-name", "--out-og", "--out-sites",
                     "--out-og-concordance", "--out-sites-concordance"):
            p.add_argument(flag)
        p.add_argument("--codon", action="store_true")
        a = p.parse_args()

        def put(path, header, row):
            if not path:
                return
            tmp = path + ".stubtmp"
            with open(tmp, "w") as f:
                f.write(header + "\\n" + row + "\\n")
            os.replace(tmp, path)

        put(a.out_og, "og_name,n_sites_total,n_sites_episodic,fraction_episodic",
            a.og_name + ",10,2,0.2")
        put(a.out_sites, "og_name,site,p_value", a.og_name + ",1,0.01")
        put(a.out_og_concordance,
            "og_name,high_confidence_sites_n,alignment_robustness_index",
            a.og_name + ",2,1.0")
        put(a.out_sites_concordance, "og_name,lenient_site,concordance_tier",
            a.og_name + ",1,high_confidence")
    """)

    (root / "results" / "selective_pressure").mkdir(parents=True, exist_ok=True)
    return script


def run_stack(root: Path, og: str, env_extra: dict | None = None,
              array_task_id: str = "0") -> subprocess.CompletedProcess:
    """Invoke the script for one orthogroup, as stage 05 does."""
    script = root / "scripts" / "hpc" / "run_selection_stack.sh"
    align = root / f"{og}_codon.phy"
    tree = root / f"{og}.treefile"
    align.write_text(">s1\nAAACCCGGG\n")
    tree.write_text("(s1,s2);\n")

    env = os.environ.copy()
    env["PATH"] = f"{root / 'bin'}{os.pathsep}{env['PATH']}"
    env["SLURM_JOB_ID"] = "999001"
    env["SLURM_ARRAY_TASK_ID"] = array_task_id
    if env_extra:
        env.update(env_extra)
    return subprocess.run(
        ["bash", str(script), str(align), str(tree), og],
        capture_output=True, text=True, env=env,
    )


def seed_per_og_csvs(out_dir: Path, n: int) -> list[str]:
    """Pre-populate per-OG CSVs so the concat loop has real work to do,
    widening the window in which a shared staging file could be truncated.

    Names are random, not counter-derived, for the reason noted on
    OG_FIXTURES above.
    """
    names = [f"OGseed_{uuid.uuid4().hex[:10]}" for _ in range(n)]
    for name in names:
        (out_dir / f"{name}_busted_s.csv").write_text(
            f"{BUSTED_HEADER}\n{name},0.02,1\n"
        )
    return names


def read_cumulative(path: Path) -> list[str]:
    return [ln for ln in path.read_text().splitlines() if ln.strip()]


# ---------------------------------------------------------------------------
# Structural guards
# ---------------------------------------------------------------------------

def _code_lines() -> list[str]:
    """Script lines with comments stripped (the fix rationale is documented
    in comments that legitimately quote the old, broken path)."""
    return [ln for ln in REAL_SCRIPT.read_text().splitlines()
            if not ln.lstrip().startswith("#")]


def test_no_shared_fixed_temp_path_remains():
    """The fixed `"${cum}.tmp"` staging path must be gone from the code."""
    code = "\n".join(_code_lines())
    assert '"${cum}.tmp"' not in code, (
        "a fixed shared staging path is still present; "
        "concurrent array tasks will truncate each other"
    )
    assert '${cum}.tmp.${TASK_TAG}' in code


def test_task_tag_is_unique_across_nodes():
    """$$ alone is not enough: PIDs are per-node and OUT_DIR is shared."""
    tag_line = next(ln for ln in _code_lines() if ln.startswith("TASK_TAG="))
    for component in ("SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID", "$$"):
        assert component in tag_line, (
            f"TASK_TAG lacks {component}; it is not unique cluster-wide: {tag_line}"
        )


def test_status_file_also_staged_uniquely(tmp_path):
    """Every staged publish in the script, not just the cumulative loop,
    must use the per-task tag."""
    staged = [ln.strip() for ln in _code_lines() if ".tmp" in ln and "=" in ln]
    assert staged, "expected at least one staging assignment"
    for line in staged:
        assert "TASK_TAG" in line, f"unstaged/shared temp path: {line}"


# ---------------------------------------------------------------------------
# Behavioural: real concurrency against the real script
# ---------------------------------------------------------------------------

def _run_concurrently(root: Path, ogs: list[str]) -> list[subprocess.CompletedProcess]:
    with ThreadPoolExecutor(max_workers=len(ogs)) as pool:
        return list(pool.map(
            lambda t: run_stack(root, t[1], array_task_id=str(t[0])),
            enumerate(ogs),
        ))


def test_concurrent_tasks_publish_complete_cumulative(tmp_path):
    """Concurrent invocations must leave every cumulative CSV complete:
    one header plus exactly one row per contributing orthogroup."""
    root = tmp_path / "proj"
    make_project(root)
    out_dir = root / "results" / "selective_pressure"
    n_seed = 80
    seed_per_og_csvs(out_dir, n_seed)

    ogs = list(OG_FIXTURES)
    n_tasks = len(ogs)
    results = _run_concurrently(root, ogs)

    for og, res in zip(ogs, results):
        assert res.returncode == 0, f"{og} failed: {res.stderr[-2000:]}"

    cumulative = out_dir / "busted_s_results.csv"
    assert cumulative.exists()
    lines = read_cumulative(cumulative)

    # 1 header + every seeded OG + every task's own OG.
    assert lines[0] == BUSTED_HEADER
    assert len(lines) == 1 + n_seed + n_tasks, (
        f"cumulative CSV is truncated/corrupt: {len(lines)} lines, "
        f"expected {1 + n_seed + n_tasks}"
    )
    names = [ln.split(",")[0] for ln in lines[1:]]
    assert len(names) == len(set(names)), "duplicate rows in cumulative CSV"
    assert set(ogs).issubset(set(names)), "a task's own row is missing"
    # No interleaved header rows — the signature of a partial republish.
    assert BUSTED_HEADER not in names


def test_no_staging_files_left_behind(tmp_path):
    """Staging files must always be published or removed, never orphaned."""
    root = tmp_path / "proj"
    make_project(root)
    out_dir = root / "results" / "selective_pressure"

    _run_concurrently(root, OG_FIXTURES[:6])

    leftovers = sorted(p.name for p in out_dir.glob("*.tmp*"))
    assert not leftovers, f"orphaned staging files: {leftovers}"


def test_meme_and_concordance_cumulatives_also_complete(tmp_path):
    """The race was in a loop over four variants; all four must be safe,
    not only the two the original bug report named."""
    root = tmp_path / "proj"
    make_project(root)
    out_dir = root / "results" / "selective_pressure"

    ogs = OG_FIXTURES[:8]
    n = len(ogs)
    _run_concurrently(root, ogs)

    for name in ("busted_s_results.csv", "busted_mh_results.csv",
                 "meme_results.csv", "meme_concordance.csv"):
        path = out_dir / name
        assert path.exists(), f"{name} was never published"
        lines = read_cumulative(path)
        assert len(lines) == 1 + n, f"{name}: {len(lines)} lines, expected {1 + n}"
        rows = [ln.split(",")[0] for ln in lines[1:]]
        assert sorted(rows) == sorted(ogs), f"{name} rows: {rows}"
