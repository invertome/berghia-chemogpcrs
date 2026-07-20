"""scripts/hpc/foldseek_against_gpcrdb.sh must refuse a zero-structure query dir.

The script takes QUERY_DIR as $1 and hands it straight to `foldseek createdb`
with no check that anything is in it. `foldseek createdb` on an empty (or
simply wrong) directory does not fail: it builds an empty DB, the subsequent
search returns nothing, convertalis writes a zero-row hits.tsv, and the whole
job exits 0. Every downstream consumer then reads an empty hits file and
treats "no structural evidence" as a legitimate result -- the same silent
dormancy that bead 5ubd fixed in scripts/unity/run_foldseek_candidates.sh,
which DOES guard its query staging (N_MODELS == 0 -> exit 3).

These tests drive the REAL script with foldseek + conda stubbed, so they
exercise the guard as written rather than a reimplementation of it.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REAL_SCRIPT = REPO_ROOT / "scripts" / "hpc" / "foldseek_against_gpcrdb.sh"

# The staged query filename convention this repo uses: one structure per
# candidate, named <candidate_id>.<ext> (scripts/unity/run_foldseek_candidates.sh
# symlinks stage 08's nested AF3 models in under exactly this scheme).
CANDIDATES = ("BersteEVm001t1", "BersteEVm042t3")


@pytest.fixture
def sandbox(tmp_path):
    """A self-contained repo copy so the script's own $0-relative PROJECT_DIR
    resolution lands inside tmp_path rather than the real checkout."""
    repo = tmp_path / "repo"
    (repo / "scripts" / "hpc").mkdir(parents=True)
    shutil.copy2(REAL_SCRIPT, repo / "scripts" / "hpc" / REAL_SCRIPT.name)

    results = repo / "results"
    results.mkdir()
    (repo / "config.sh").write_text(f'RESULTS_DIR="{results}"\n')

    # A present-but-stub GPCRdb foldseek DB dir, so the script gets past its
    # existing DB check and reaches the query-dir handling under test.
    db = repo / "references" / "foldseek" / "gpcrdb_2025"
    db.mkdir(parents=True)
    (db / "stub_db").write_text("stub")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    calls = tmp_path / "foldseek_calls.txt"
    stub = bin_dir / "foldseek"
    stub.write_text(
        "#!/bin/bash\n"
        f'echo "$@" >> "{calls}"\n'
        # convertalis writes the hits TSV ($4); everything else is a no-op that
        # succeeds, mirroring foldseek's real behaviour on an empty query DB.
        'if [ "$1" = "convertalis" ]; then : > "$4"; fi\n'
        "exit 0\n"
    )
    stub.chmod(0o755)

    return {"tmp": tmp_path, "repo": repo, "results": results,
            "bin": bin_dir, "db": db, "calls": calls}


def run_script(sandbox, query_dir, extra_env=None):
    env = dict(os.environ)
    env["PATH"] = f"{sandbox['bin']}{os.pathsep}{env['PATH']}"
    env["GPCRDB_FOLDSEEK_DB"] = str(sandbox["db"])
    env.update(extra_env or {})
    return subprocess.run(
        ["bash", str(sandbox["repo"] / "scripts" / "hpc" / REAL_SCRIPT.name),
         str(query_dir)],
        cwd=sandbox["repo"], env=env, capture_output=True, text=True,
    )


def _populate(query_dir: Path, ext="cif"):
    query_dir.mkdir(parents=True, exist_ok=True)
    for cand in CANDIDATES:
        (query_dir / f"{cand}.{ext}").write_text("data_stub\n")
    return query_dir


# --------------------------------------------------------------------------- #
# The guard
# --------------------------------------------------------------------------- #

def test_empty_query_dir_is_refused_not_silently_searched(sandbox):
    """An existing but EMPTY query dir must abort before `foldseek createdb`."""
    query_dir = sandbox["tmp"] / "empty_queries"
    query_dir.mkdir()

    proc = run_script(sandbox, query_dir)

    assert proc.returncode != 0, (
        "an empty query dir exited 0 -- a zero-structure run that looks like "
        f"success.\nstdout={proc.stdout}\nstderr={proc.stderr}"
    )
    # And it must not have reached foldseek at all.
    assert not sandbox["calls"].exists() or "createdb" not in sandbox["calls"].read_text()


def test_missing_query_dir_is_refused(sandbox):
    """A wrong/nonexistent path must abort too, with the path in the message."""
    query_dir = sandbox["tmp"] / "does_not_exist"

    proc = run_script(sandbox, query_dir)

    assert proc.returncode != 0
    assert str(query_dir) in proc.stderr


def test_dir_of_non_structure_files_is_refused(sandbox):
    """The 'wrong directory' case: files present, but none foldseek can ingest.

    Pointing at e.g. a FASTA dir, or at stage 08's alphafold/ root (whose depth-1
    entries are all DIRECTORIES -- the bead 5ubd failure), must not read as a
    valid zero-hit search.
    """
    query_dir = sandbox["tmp"] / "wrong_dir"
    query_dir.mkdir()
    (query_dir / "candidates.fasta").write_text(">x\nMEEP\n")
    (query_dir / "BersteEVm001t1").mkdir()  # a nested per-candidate dir

    proc = run_script(sandbox, query_dir)

    assert proc.returncode != 0, (
        "a directory containing no structure files exited 0.\n"
        f"stdout={proc.stdout}\nstderr={proc.stderr}"
    )


def test_error_message_names_the_dir_and_the_count(sandbox):
    query_dir = sandbox["tmp"] / "empty_queries"
    query_dir.mkdir()

    proc = run_script(sandbox, query_dir)

    assert str(query_dir) in proc.stderr, proc.stderr
    assert "structure" in proc.stderr.lower(), proc.stderr


# --------------------------------------------------------------------------- #
# The guard must not block the legitimate path
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("ext", ["cif", "pdb", "mmcif", "ent"])
def test_populated_query_dir_still_runs(sandbox, ext):
    """Every structure extension foldseek ingests must pass the guard."""
    query_dir = _populate(sandbox["tmp"] / f"queries_{ext}", ext=ext)

    proc = run_script(sandbox, query_dir)

    assert proc.returncode == 0, f"stdout={proc.stdout}\nstderr={proc.stderr}"
    calls = sandbox["calls"].read_text()
    assert "createdb" in calls
    assert "search" in calls
    assert "convertalis" in calls


def test_gzipped_structures_pass_the_guard(sandbox):
    """foldseek reads .cif.gz/.pdb.gz natively; the guard must not reject them."""
    query_dir = sandbox["tmp"] / "queries_gz"
    query_dir.mkdir()
    (query_dir / f"{CANDIDATES[0]}.cif.gz").write_bytes(b"\x1f\x8b\x08")

    proc = run_script(sandbox, query_dir)

    assert proc.returncode == 0, f"stdout={proc.stdout}\nstderr={proc.stderr}"


def test_symlinked_structures_pass_the_guard(sandbox):
    """run_foldseek_candidates.sh stages queries as SYMLINKS, not regular files,
    so a `-type f`-only count would reject the pipeline's own query dir."""
    real = sandbox["tmp"] / "real_models"
    real.mkdir()
    (real / "af3_job_model.cif").write_text("data_stub\n")

    query_dir = sandbox["tmp"] / "queries_links"
    query_dir.mkdir()
    for cand in CANDIDATES:
        (query_dir / f"{cand}.cif").symlink_to(real / "af3_job_model.cif")

    proc = run_script(sandbox, query_dir)

    assert proc.returncode == 0, f"stdout={proc.stdout}\nstderr={proc.stderr}"


def test_reports_how_many_structures_it_will_search(sandbox):
    """Positive confirmation of the count, so a 1-of-790 run is visible in the
    log rather than indistinguishable from a full one."""
    query_dir = _populate(sandbox["tmp"] / "queries_count")

    proc = run_script(sandbox, query_dir)

    assert proc.returncode == 0, proc.stderr
    combined = proc.stdout + proc.stderr
    assert str(len(CANDIDATES)) in combined, combined
