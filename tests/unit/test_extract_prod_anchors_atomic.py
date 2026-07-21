"""The anchor-subset FASTA must be published by rename, never written in place.

This output is a production input: the reference npz is built from it, and the family
prototypes candidates are scored against are built from that. Writing the destination
with mode "w" truncates it at open and refills incrementally, so an interruption
leaves a SHORT BUT VALID FASTA. A short FASTA parses, so no downstream consumer
rejects it -- the corruption is silent all the way to the prototypes.

That is not a hypothetical. The anchor TSV this script reads was itself left
half-written by an interrupted run (55 of 63 intended row removals applied), and the
only signal was a failing equality test in a neighbouring artifact. Had the same
interruption hit this script instead, the result would have been a reference set
silently missing sequences with nothing to compare it against.
"""
from __future__ import annotations

import importlib.util
import pathlib
import sys

import pytest

SCRIPT = pathlib.Path(__file__).resolve().parents[1].parent / "scripts" / "extract_prod_anchors.py"


def _load():
    spec = importlib.util.spec_from_file_location("extract_prod_anchors", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["extract_prod_anchors"] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_source(path: pathlib.Path, ids) -> None:
    with open(path, "w") as fh:
        for i in ids:
            fh.write(f">{i}\nMKTAYIAKQRQISFVKSHFSRQ\n")


def test_source_never_opened_for_writing_at_destination(tmp_path, monkeypatch) -> None:
    """The destination path must never be passed to open() in a write mode.

    This is the discriminating check: a fix that merely adds an fsync, or writes to
    a temp in the WRONG directory, would still leave a truncation window or a
    cross-filesystem rename that is not atomic.
    """
    mod = _load()
    src = tmp_path / "all.fasta"
    out = tmp_path / "subset.fasta"
    _write_source(src, [f"ANCHOR_A_2_ACC{i:03d}" for i in range(6)])
    wanted = {f"ANCHOR_A_2_ACC{i:03d}" for i in range(4)}

    real_open = open
    offending = []

    def watched_open(file, mode="r", *a, **k):
        if str(file) == str(out) and any(m in mode for m in ("w", "a", "+")):
            offending.append((str(file), mode))
        return real_open(file, mode, *a, **k)

    monkeypatch.setattr("builtins.open", watched_open)
    mod.subset_fasta_by_ids(str(src), wanted, str(out))

    assert not offending, (
        f"destination opened for writing directly: {offending}; it must be staged "
        f"to a temp and renamed"
    )
    assert out.read_text().count(">") == 4


def test_interrupted_write_leaves_the_previous_file_intact(tmp_path, monkeypatch) -> None:
    """A crash mid-write must not shorten an existing destination.

    Under the pre-fix code the destination was already truncated to zero bytes by
    the time the failure occurred, so a consumer globbing the directory would find
    an empty-but-present reference set.
    """
    mod = _load()
    src = tmp_path / "all.fasta"
    out = tmp_path / "subset.fasta"
    _write_source(src, [f"ANCHOR_A_2_ACC{i:03d}" for i in range(6)])
    out.write_text(">PREVIOUS_GOOD\nMKTAYIAKQRQ\n")
    before = out.read_text()

    real_replace = mod.os.replace

    def boom(*a, **k):
        raise RuntimeError("interrupted before publish")

    monkeypatch.setattr(mod.os, "replace", boom)
    with pytest.raises(RuntimeError):
        mod.subset_fasta_by_ids(str(src), {"ANCHOR_A_2_ACC000"}, str(out))

    assert out.read_text() == before, "the prior reference set was damaged"
    monkeypatch.setattr(mod.os, "replace", real_replace)


def test_failed_write_leaves_no_temp_debris(tmp_path, monkeypatch) -> None:
    """A partial file must not survive where a glob could pick it up.

    Several consumers in this pipeline locate inputs by globbing a directory, so
    debris is not merely untidy -- it is a candidate input.
    """
    mod = _load()
    src = tmp_path / "all.fasta"
    out = tmp_path / "subset.fasta"
    _write_source(src, [f"ANCHOR_A_2_ACC{i:03d}" for i in range(3)])

    monkeypatch.setattr(mod.os, "replace", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")))
    with pytest.raises(RuntimeError):
        mod.subset_fasta_by_ids(str(src), {"ANCHOR_A_2_ACC000"}, str(out))

    leftovers = [p.name for p in tmp_path.iterdir() if p.name.startswith("subset.fasta.")]
    assert not leftovers, f"temp debris left behind: {leftovers}"


def test_write_once_violation_still_raises_before_any_write(tmp_path) -> None:
    """Validation must precede publication, not follow it.

    A requested anchor absent from the source is a write-once violation. If that
    check ran after the destination was opened, the failure would still have
    truncated the previous good file.
    """
    mod = _load()
    src = tmp_path / "all.fasta"
    out = tmp_path / "subset.fasta"
    _write_source(src, ["ANCHOR_A_2_ACC000"])
    out.write_text(">PREVIOUS_GOOD\nMKTAYIAKQRQ\n")

    with pytest.raises(ValueError, match="missing"):
        mod.subset_fasta_by_ids(str(src), {"ANCHOR_A_2_ACC000", "ANCHOR_A_2_ABSENT"}, str(out))

    assert out.read_text().startswith(">PREVIOUS_GOOD"), (
        "a write-once violation damaged the existing file"
    )


def test_output_is_deterministic_across_runs(tmp_path) -> None:
    """Byte-identical output for identical input.

    The reference npz is keyed on these headers; a reordering would look like a
    changed reference set to anything comparing checksums.
    """
    mod = _load()
    src = tmp_path / "all.fasta"
    _write_source(src, [f"ANCHOR_A_2_ACC{i:03d}" for i in range(8)])
    wanted = {f"ANCHOR_A_2_ACC{i:03d}" for i in (5, 1, 7, 2)}

    a, b = tmp_path / "a.fasta", tmp_path / "b.fasta"
    mod.subset_fasta_by_ids(str(src), wanted, str(a))
    mod.subset_fasta_by_ids(str(src), wanted, str(b))
    assert a.read_bytes() == b.read_bytes()


def test_the_function_under_test_actually_exists() -> None:
    """Pin the real API name.

    These tests were first written against an assumed function name and passed
    nothing, which is the same defect they exist to prevent: a fixture encoding
    what the author expected rather than what the module provides.
    """
    mod = _load()
    assert hasattr(mod, "subset_fasta_by_ids"), sorted(
        n for n in dir(mod) if not n.startswith("_")
    )
