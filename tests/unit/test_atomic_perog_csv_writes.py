"""REPAIR 2 — the per-OG CSV writers must publish atomically.

Stage 05 runs as ``#SBATCH --array=0-999%50``: up to 50 tasks run on
DIFFERENT nodes against the same shared ``results/selective_pressure/``.
Each task writes its own per-OG CSV and then re-runs a concatenation loop
that globs *every* per-OG CSV in that directory.

``05_selective_pressure_and_asr.sh``'s staging comment asserted the per-OG
CSVs written by ``parse_busted.py`` / ``parse_meme.py`` were "fast, atomic".
They were not: both used a bare ``open(path, "w")``, which truncates the
published file at open() and refills it incrementally. The window between
truncate and flush is a silent data-corruption path that the *staging* fix
does not close, because the corruption is in the glob's INPUT:

  task A: open("OG7_busted_s.csv", "w")   <- file is now 0 bytes
  task B: globs the directory, cats OG7's EMPTY/partial csv into the
          cumulative, and publishes it. No error anywhere.

The fix is the project standard: write a per-task-unique temp in the same
directory, then ``os.replace`` it into place. ``os.replace`` is atomic
within a filesystem, so a concurrent reader observes either the whole old
file or the whole new file, never a partial one.

These tests encode the *behaviour* (a failed write must not damage the
published file, and a successful write must leave no debris), not the
spelling of the implementation.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"

SENTINEL = "og_name,model_variant,p_value\nOG_PREVIOUS,S,0.001\n"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def busted():
    return _load("parse_busted")


@pytest.fixture
def meme():
    return _load("parse_meme")


def _tmp_debris(directory: Path) -> list[str]:
    """Any staging file left behind in the published directory."""
    return sorted(
        p.name for p in directory.iterdir()
        if ".tmp" in p.name or p.name.endswith("~")
    )


# --------------------------------------------------------------------------
# the shared helper
# --------------------------------------------------------------------------
@pytest.mark.parametrize("modname", ["parse_busted", "parse_meme", "parse_absrel"])
def test_parser_exposes_an_atomic_csv_writer(modname):
    """Both parsers must route every CSV publish through one atomic helper."""
    mod = _load(modname)
    assert hasattr(mod, "write_csv_atomic"), (
        f"{modname}.py must expose write_csv_atomic(); every per-OG CSV it "
        f"publishes is globbed by concurrent stage-05 array tasks"
    )


@pytest.mark.parametrize("modname", ["parse_busted", "parse_meme", "parse_absrel"])
def test_failed_write_leaves_the_published_file_untouched(modname, tmp_path, monkeypatch):
    """A crash mid-write must not truncate the previously published CSV.

    This is the discriminating test: with ``open(path, "w")`` the file is
    already 0 bytes by the time the writer raises, so a concurrent reader
    sees an empty CSV. With tmp + os.replace the published file is never
    opened at all.
    """
    mod = _load(modname)
    out = tmp_path / "OG0000007_busted_s.csv"
    out.write_text(SENTINEL)

    class Boom(Exception):
        pass

    def exploding_writer(*a, **kw):
        raise Boom("simulated crash partway through serialisation")

    monkeypatch.setattr(mod.csv, "DictWriter", exploding_writer)

    with pytest.raises(Boom):
        mod.write_csv_atomic(out, ["og_name"], [{"og_name": "OG0000007"}])

    assert out.read_text() == SENTINEL, (
        "the published per-OG CSV was damaged by a failed write; a concurrent "
        "array task would concatenate the damaged file into the cumulative"
    )


@pytest.mark.parametrize("modname", ["parse_busted", "parse_meme", "parse_absrel"])
def test_successful_write_publishes_and_leaves_no_debris(modname, tmp_path):
    mod = _load(modname)
    out = tmp_path / "OG0000007_busted_s.csv"
    out.write_text(SENTINEL)

    mod.write_csv_atomic(out, ["og_name", "p_value"],
                         [{"og_name": "OG0000007", "p_value": 0.5}])

    text = out.read_text()
    assert "OG_PREVIOUS" not in text, "the new content did not replace the old"
    assert "OG0000007" in text
    assert _tmp_debris(tmp_path) == [], (
        "staging files were left in the published directory; the stage-05 "
        "concat globs this directory"
    )


def test_staging_name_is_unique_per_process(busted, tmp_path, monkeypatch):
    """Two concurrent writers must not share a staging path.

    Array tasks land on different nodes, so a staging name derived from the
    PID alone can collide. The name must fold in the SLURM array identity.
    """
    seen: list[str] = []
    real_replace = os.replace

    def spy(src, dst):
        seen.append(Path(src).name)
        real_replace(src, dst)

    monkeypatch.setattr(busted.os, "replace", spy)

    out = tmp_path / "OG0000007_busted_s.csv"
    for job, task, pid in (("1001", "7", 4242), ("1001", "8", 4242)):
        monkeypatch.setenv("SLURM_JOB_ID", job)
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", task)
        monkeypatch.setattr(busted.os, "getpid", lambda pid=pid: pid)
        busted.write_csv_atomic(out, ["og_name"], [{"og_name": "OG0000007"}])

    assert len(seen) == 2
    assert seen[0] != seen[1], (
        f"both array tasks staged through the same name {seen[0]!r}; on shared "
        f"storage that is a lost-update race"
    )


# --------------------------------------------------------------------------
# end-to-end through main()
# --------------------------------------------------------------------------
def test_busted_main_publishes_atomically(busted, tmp_path, monkeypatch):
    batch = tmp_path / "batch"
    batch.mkdir()
    out = tmp_path / "out" / "busted_s_results.csv"
    monkeypatch.setattr(
        sys, "argv",
        ["parse_busted.py", "--batch", str(batch), "--variant", "S",
         "--out", str(out)],
    )
    assert busted.main() == 0
    assert out.exists()
    assert _tmp_debris(out.parent) == []


def test_meme_main_publishes_both_csvs_atomically(meme, tmp_path, monkeypatch):
    batch = tmp_path / "batch"
    batch.mkdir()
    outdir = tmp_path / "out"
    out_og = outdir / "OG0000007_meme.csv"
    out_sites = outdir / "OG0000007_meme_sites.csv"
    monkeypatch.setattr(
        sys, "argv",
        ["parse_meme.py", "--batch", str(batch),
         "--out-og", str(out_og), "--out-sites", str(out_sites)],
    )
    assert meme.main() == 0
    assert out_og.exists() and out_sites.exists()
    assert _tmp_debris(outdir) == []


@pytest.mark.parametrize(
    "modname,attr",
    [("parse_busted", "out"), ("parse_meme", "out_og")],
)
def test_no_raw_truncating_open_remains(modname, attr):
    """Guard against a future edit reintroducing the non-atomic publish."""
    text = (SCRIPTS / f"{modname}.py").read_text()
    assert f'open(args.{attr}, "w"' not in text, (
        f"{modname}.py publishes args.{attr} with a truncating open(); it must "
        f"go through write_csv_atomic()"
    )
