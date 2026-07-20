"""RESID REPAIR 1 — parse_absrel.py must publish its per-OG CSV atomically.

``scripts/parse_absrel.py`` produces ``*_absrel.csv``, the per-branch dN/dS
table that is the PRIMARY selection axis of the ranking. Its sibling
parsers (parse_busted.py, parse_meme.py) were fixed to publish through
``write_csv_atomic`` (unique temp + fsync + ``os.replace``); parse_absrel
was left on a bare ``open(output_csv, mode)``, which truncates the
published file at open() and refills it incrementally.

Why that is a live corruption path, verified in this repo rather than
assumed:

  * ``05_selective_pressure_and_asr.sh:14`` -> ``#SBATCH --array=0-999%50``
    so up to 50 tasks run concurrently on DIFFERENT nodes.
  * ``05_selective_pressure_and_asr.sh:453`` invokes parse_absrel.py to
    write ``<base>_absrel.csv`` into ``${RESULTS_DIR}/selective_pressure``.
  * ``05_selective_pressure_and_asr.sh:454`` then runs
    ``concat_per_og_csv "${RESULTS_DIR}/selective_pressure" "_absrel.csv" ...``
    which globs that same shared directory.

So one task's truncate-and-refill window is another task's glob input, and
a half-written CSV gets concatenated into ``absrel_results.csv`` with no
error anywhere. Stage 05's own comment block (lines 292-297) already
recorded parse_absrel as "NOT ATOMIC" while its module docstring claimed
to support "atomic per-OG CSV writes" -- that phrase meant per-orthogroup
isolation, not write atomicity.

These tests encode BEHAVIOUR (a failed write must not damage the published
file; a successful write must leave no debris in the globbed directory),
not the spelling of the implementation.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"

# A previously-published per-OG aBSREL CSV, in the real schema.
SENTINEL = (
    "branch_id,omega,omega_max,omega_mean,weight_at_max,n_rate_classes,"
    "p_value,corrected_p_value,is_already_corrected,is_selected,cds_source\n"
    "BRANCH_PREVIOUS,0.21,0.21,0.21,1.0,1,0.9,0.9,1,0,native\n"
)


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def absrel():
    return _load("parse_absrel")


def _tmp_debris(directory: Path) -> list[str]:
    """Any staging file left behind in the directory stage 05 globs."""
    return sorted(
        p.name for p in directory.iterdir()
        if ".tmp" in p.name or p.name.endswith("~")
    )


def _absrel_json(tmp_path: Path, name: str = "absrel.json") -> Path:
    """A minimal but real-shaped aBSREL JSON (partitioned branch attributes)."""
    data = {
        "branch attributes": {
            "0": {
                "NODE1": {
                    "original name": "NODE1",
                    "Rate Distributions": [[0.4, 0.9], [3.1, 0.1]],
                    "Uncorrected P-value": 0.01,
                    "Corrected P-value": 0.03,
                }
            }
        }
    }
    p = tmp_path / name
    p.write_text(json.dumps(data))
    return p


# --------------------------------------------------------------------------
# the shared helper
# --------------------------------------------------------------------------
def test_parse_absrel_exposes_an_atomic_csv_writer(absrel):
    """parse_absrel must route its CSV publish through the same helper the
    other two stage-05 parsers use -- one mechanism, not a second one."""
    assert hasattr(absrel, "write_csv_atomic"), (
        "parse_absrel.py must expose write_csv_atomic(); the *_absrel.csv it "
        "publishes is globbed by concurrent stage-05 array tasks"
    )


def test_failed_write_leaves_the_published_file_untouched(absrel, tmp_path, monkeypatch):
    """A crash mid-write must not truncate the previously published CSV.

    This is the discriminating test. With ``open(path, "w")`` the file is
    already 0 bytes by the time the writer raises, so a concurrent globbing
    task concatenates an empty CSV into absrel_results.csv. With tmp +
    os.replace the published file is never opened at all.
    """
    out = tmp_path / "OG0000007_absrel.csv"
    out.write_text(SENTINEL)

    class Boom(Exception):
        pass

    def exploding_writer(*a, **kw):
        raise Boom("simulated crash partway through serialisation")

    monkeypatch.setattr(absrel.csv, "DictWriter", exploding_writer)

    with pytest.raises(Boom):
        absrel.write_csv_atomic(out, absrel.FIELDNAMES, [{"branch_id": "NODE1"}])

    assert out.read_text() == SENTINEL, (
        "the published per-OG aBSREL CSV was damaged by a failed write; a "
        "concurrent array task would concatenate the damaged file into "
        "absrel_results.csv -- the primary dN/dS ranking axis"
    )


def test_successful_write_publishes_and_leaves_no_debris(absrel, tmp_path):
    out = tmp_path / "OG0000007_absrel.csv"
    out.write_text(SENTINEL)

    absrel.write_csv_atomic(
        out, absrel.FIELDNAMES,
        [{"branch_id": "NODE1", "omega": 0.5}],
    )

    text = out.read_text()
    assert "BRANCH_PREVIOUS" not in text, "the new content did not replace the old"
    assert "NODE1" in text
    assert _tmp_debris(tmp_path) == [], (
        "staging files were left in the published directory; "
        "05_selective_pressure_and_asr.sh:454 globs this directory"
    )


def test_staging_name_does_not_collide_across_array_tasks(absrel, tmp_path, monkeypatch):
    """Two tasks on DIFFERENT nodes can hold the same PID, and this is
    shared storage -- so the staging name must fold in the SLURM identity."""
    out = tmp_path / "OG0000007_absrel.csv"
    seen: list[str] = []

    real_replace = absrel.os.replace

    def spy_replace(src, dst):
        seen.append(str(src))
        return real_replace(src, dst)

    monkeypatch.setattr(absrel.os, "replace", spy_replace)
    monkeypatch.setattr(absrel.os, "getpid", lambda: 4242)  # same PID, two nodes

    for task_id in ("11", "12"):
        monkeypatch.setenv("SLURM_JOB_ID", "61999999")
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", task_id)
        absrel.write_csv_atomic(out, absrel.FIELDNAMES, [{"branch_id": "NODE1"}])

    assert len(seen) == 2
    assert seen[0] != seen[1], (
        f"both array tasks staged through the same name {seen[0]!r}; on shared "
        f"storage that is a lost-update race"
    )


def test_staging_file_is_invisible_to_the_stage05_concat_glob(absrel, tmp_path, monkeypatch):
    """concat_per_og_csv globs '*_absrel.csv'. The staging name must not
    match it, or the concat would swallow a half-written temp."""
    out = tmp_path / "OG0000007_absrel.csv"
    staged: list[Path] = []

    real_replace = absrel.os.replace

    def spy_replace(src, dst):
        staged.append(Path(src))
        return real_replace(src, dst)

    monkeypatch.setattr(absrel.os, "replace", spy_replace)
    absrel.write_csv_atomic(out, absrel.FIELDNAMES, [{"branch_id": "NODE1"}])

    assert staged, "write_csv_atomic did not stage through os.replace"
    for p in staged:
        assert not p.name.endswith("_absrel.csv"), (
            f"staging name {p.name!r} matches the stage-05 concat glob "
            f"'*_absrel.csv'"
        )
        assert p.parent == out.parent, (
            "the temp must live in the destination directory, or os.replace "
            "can cross a filesystem boundary and lose its atomicity"
        )


# --------------------------------------------------------------------------
# end-to-end through the real parse entry point
# --------------------------------------------------------------------------
def test_parse_absrel_json_publishes_atomically(absrel, tmp_path, monkeypatch):
    """The public entry point -- not just the helper -- must go through
    os.replace. A fix applied only to the helper would be dead code."""
    j = _absrel_json(tmp_path)
    out = tmp_path / "OG0000007_absrel.csv"

    replaced: list[str] = []
    real_replace = absrel.os.replace

    def spy_replace(src, dst):
        replaced.append(str(dst))
        return real_replace(src, dst)

    monkeypatch.setattr(absrel.os, "replace", spy_replace)
    absrel.parse_absrel_json(str(j), str(out))

    assert str(out) in replaced, (
        "parse_absrel_json published without os.replace; the write is still "
        "a truncate-and-refill on the file stage 05 globs"
    )
    rows = list(csv.DictReader(out.open()))
    assert [r["branch_id"] for r in rows] == ["NODE1"]
    assert _tmp_debris(tmp_path) == []


def test_append_mode_is_also_atomic_and_preserves_prior_rows(absrel, tmp_path, monkeypatch):
    """Legacy append mode must not regress to a non-atomic 'a' open either.

    Appending in place has the same visibility problem: a globbing task can
    read the file between the two writes. The published file must go from
    'old rows' to 'old + new rows' in one atomic step.
    """
    j1 = _absrel_json(tmp_path, "a.json")
    out = tmp_path / "OG0000007_absrel.csv"
    absrel.parse_absrel_json(str(j1), str(out))

    data = {
        "branch attributes": {
            "0": {
                "NODE2": {
                    "original name": "NODE2",
                    "Rate Distributions": [[0.2, 1.0]],
                    "Uncorrected P-value": 0.4,
                }
            }
        }
    }
    j2 = tmp_path / "b.json"
    j2.write_text(json.dumps(data))

    replaced: list[str] = []
    real_replace = absrel.os.replace

    def spy_replace(src, dst):
        replaced.append(str(dst))
        return real_replace(src, dst)

    monkeypatch.setattr(absrel.os, "replace", spy_replace)
    absrel.parse_absrel_json(str(j2), str(out), append_mode="append")

    assert str(out) in replaced, "append mode still publishes non-atomically"
    rows = list(csv.DictReader(out.open()))
    assert [r["branch_id"] for r in rows] == ["NODE1", "NODE2"], (
        "append mode must preserve the previously published rows"
    )
    assert _tmp_debris(tmp_path) == []


def test_append_mode_still_rejects_a_schema_mismatch(absrel, tmp_path):
    """The bead -mqt guard must survive the atomicity change."""
    out = tmp_path / "OG0000007_absrel.csv"
    out.write_text("branch_id,omega,p_value\nfoo,1.0,0.5\n")
    j = _absrel_json(tmp_path)

    with pytest.raises(ValueError, match="schema"):
        absrel.parse_absrel_json(str(j), str(out), append_mode="append")

    assert out.read_text() == "branch_id,omega,p_value\nfoo,1.0,0.5\n", (
        "a rejected append must leave the published file byte-identical"
    )


def test_docstring_no_longer_claims_write_atomicity_it_does_not_have(absrel):
    """The module docstring said 'Supports atomic per-OG CSV writes' while the
    write was a truncate-and-refill. Whatever it says now must be true."""
    doc = absrel.__doc__ or ""
    if "atomic" in doc.lower():
        assert "os.replace" in doc or "temp" in doc.lower(), (
            "the docstring claims atomicity; it must name the mechanism that "
            "delivers it so the claim stays checkable"
        )
