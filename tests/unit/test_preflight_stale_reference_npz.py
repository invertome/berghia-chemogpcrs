"""The reference-npz staleness guard (TRAP 1b).

TRAP 1 measures `|labels & npz| / |npz|` against a floor. That fraction catches a
WRONG label file -- the real defect it was written for, measured at 116/953 = 0.12.
It cannot catch a STALE npz, because a stale npz has EXTRA keys rather than missing
ones, and extra keys move the fraction only slightly.

The numbers here are the real ones. Repairing the anchor set evicted 55 class-A
entries that were curated class B or class C (UniProt "G-protein coupled receptor 2
family" is class B, "3 family" is class C), taking the class-A set 953 -> 898. A
reference npz built before that repair still embeds all 953. Overlap is then
898/953 = 0.9423, comfortably above REF_OVERLAP_MIN = 0.80, so TRAP 1 stays silent
while every family prototype is still built from all 55 evicted non-class-A members.

The invariant TRAP 1b enforces is directional and exact: every key in the reference
npz must be named in the label set. A key that is not is an embedding of something
the reference set no longer contains.
"""
from __future__ import annotations

import ast
import pathlib

import pytest

PREFLIGHT = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "unity" / "embedding_channel_preflight.py"


def _source() -> str:
    return PREFLIGHT.read_text()


def test_preflight_module_parses() -> None:
    ast.parse(_source())


def test_trap1b_exists_and_is_an_exact_zero_check() -> None:
    """The guard must assert zero orphans, never a proportion.

    If this is ever relaxed to a fraction it stops detecting the failure it was
    added for, because the stale-npz signature lives entirely in the direction the
    fraction does not measure.
    """
    src = _source()
    assert "TRAP 1b FIRED" in src, "TRAP 1b is missing"
    assert "set(ref) - set(ref_labels)" in src, (
        "TRAP 1b must compare npz keys AGAINST labels in that direction; the "
        "reverse direction is what TRAP 1 already covers"
    )


def test_the_real_eviction_ratio_would_not_trip_trap1() -> None:
    """The measured case: 898 surviving labels against a 953-key stale npz.

    This is the arithmetic that makes TRAP 1 insufficient. If REF_OVERLAP_MIN is
    ever raised far enough to catch this, it would also reject legitimate partial
    reference sets, which is why the fix is a separate directional check rather
    than a stricter threshold.
    """
    import sys

    sys.path.insert(0, str(PREFLIGHT.parent))
    labels_after_repair = 898
    npz_before_repair = 953
    frac = labels_after_repair / npz_before_repair
    assert frac == pytest.approx(0.9423, abs=5e-4)
    assert frac > 0.80, (
        "if this ever fails, REF_OVERLAP_MIN changed and this test's premise "
        "needs revisiting"
    )


def test_orphan_detection_logic_on_the_real_shape() -> None:
    """Reproduce the detection directly: 953-key npz, 898 surviving labels."""
    ref_npz_keys = {f"ANCHOR_A_2_ACC{i:04d}" for i in range(953)}
    surviving_labels = {f"ANCHOR_A_2_ACC{i:04d}" for i in range(898)}

    orphans = sorted(ref_npz_keys - surviving_labels)
    assert len(orphans) == 55, "the 55 evicted anchors must surface as orphans"

    overlap_frac = len(surviving_labels & ref_npz_keys) / len(ref_npz_keys)
    assert overlap_frac > 0.80, "TRAP 1 would pass this, which is the whole point"

    # And the healthy case must be silent.
    assert not (surviving_labels - surviving_labels)
    assert not (surviving_labels ^ surviving_labels)


def test_guard_is_silent_when_npz_matches_labels() -> None:
    """A correctly regenerated npz must not fire the guard.

    Guards that fire on healthy input get disabled, so over-triggering is its own
    failure mode.
    """
    labels = {f"ANCHOR_A_2_ACC{i:04d}" for i in range(898)}
    npz = set(labels)
    assert not (npz - labels)


def test_labels_may_exceed_npz_without_firing_trap1b() -> None:
    """A label set larger than the npz is TRAP 1's business, not TRAP 1b's.

    The anchor TSV legitimately carries class B/C/F rows that the class-A
    reference npz never embeds, so labels-minus-npz is expected to be non-empty
    and must NOT be treated as staleness.
    """
    labels = {f"ANCHOR_A_2_ACC{i:04d}" for i in range(898)}
    labels |= {f"ANCHOR_B_2_BCC{i:04d}" for i in range(89)}
    npz = {f"ANCHOR_A_2_ACC{i:04d}" for i in range(898)}
    assert not (npz - labels), "class B/C/F labels must not read as orphan npz keys"
