"""A per-item warning that fires for 100% of items is indistinguishable from noise.

That is precisely what hid the stage-05 outage. Its loop logged

    "Missing alignment or tree for OG0000123"

once per orthogroup, for every orthogroup, forever. Each individual message was
true and plausible -- stage 04 really does skip orthogroups it cannot build a
tree for -- so nobody could tell the difference between "a few OGs were
skipped" and "this stage has produced nothing since the refactor".

The guard here tracks per-item outcomes across a loop and escalates to a
STAGE-LEVEL error when the failure RATE crosses a threshold. Rate, not count:
a count threshold tuned for 3000 orthogroups never fires on a 3-orthogroup test
run, and one tuned for 3 screams on every large run.

Threshold rationale pinned by these tests: the default is 0.9, not 1.0.
100% is the unambiguous case, but a near-total failure is the same bug wearing
a disguise -- one stale directory that happens to resolve is enough to pull a
broken stage off 100% and back into silence. Below N=10 a rate above 0.9
requires 100% anyway, so the default cannot fire on a handful of items unless
literally everything failed. Legitimate per-item skips in this pipeline are a
minority of items by construction, so 0.9 sits far above the healthy band.
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


def _loop(label: str, n_ok: int, n_fail: int, threshold: str = "") -> str:
    """Emit a loop that records n_ok successes and n_fail failures."""
    thr = f" {threshold}" if threshold else ""
    return f"""
saturation_begin {label}{thr}
for i in $(seq 1 {n_ok}); do saturation_record {label} ok "item_ok_$i"; done
for i in $(seq 1 {n_fail}); do saturation_record {label} fail "item_bad_$i"; done
saturation_report {label} && echo GUARD_OK || echo GUARD_SATURATED
"""


# --------------------------------------------------------------------------
# the unambiguous case
# --------------------------------------------------------------------------
def test_total_failure_escalates(tmp_path: Path) -> None:
    """The stage-05 signature: every single item failed."""
    r = _run(tmp_path, _loop("perog", n_ok=0, n_fail=20))
    assert "GUARD_SATURATED" in r.stdout, r.stdout + r.stderr
    assert "ERROR" in r.stdout + r.stderr


def test_total_failure_escalates_at_small_n(tmp_path: Path) -> None:
    """Rate-based, so a 3-item run is as protected as a 3000-item run."""
    r = _run(tmp_path, _loop("perog", n_ok=0, n_fail=3))
    assert "GUARD_SATURATED" in r.stdout, r.stdout + r.stderr


def test_total_failure_escalates_at_large_n(tmp_path: Path) -> None:
    r = _run(tmp_path, _loop("perog", n_ok=0, n_fail=300))
    assert "GUARD_SATURATED" in r.stdout, r.stdout + r.stderr


# --------------------------------------------------------------------------
# the healthy band stays quiet -- 'legitimately empty' is not a false alarm
# --------------------------------------------------------------------------
def test_healthy_failure_rate_stays_quiet(tmp_path: Path) -> None:
    """Stage 04 skipping a few small orthogroups is normal and must not shout."""
    r = _run(tmp_path, _loop("perog", n_ok=90, n_fail=10))
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr
    assert "ERROR" not in r.stdout + r.stderr


def test_substantial_but_sub_threshold_failure_stays_quiet(tmp_path: Path) -> None:
    """Half failing is bad news but it is NOT the silent-no-op signature."""
    r = _run(tmp_path, _loop("perog", n_ok=50, n_fail=50))
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr
    assert "ERROR" not in r.stdout + r.stderr


def test_no_failures_at_all_stays_quiet(tmp_path: Path) -> None:
    r = _run(tmp_path, _loop("perog", n_ok=10, n_fail=0))
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr
    assert "ERROR" not in r.stdout + r.stderr


# --------------------------------------------------------------------------
# the near-total case is why the default is 0.9 rather than 1.0
# --------------------------------------------------------------------------
def test_near_total_failure_escalates(tmp_path: Path) -> None:
    """95/100 failing is the same bug as 100/100, with one stale leftover."""
    r = _run(tmp_path, _loop("perog", n_ok=5, n_fail=95))
    assert "GUARD_SATURATED" in r.stdout, r.stdout + r.stderr
    assert "ERROR" in r.stdout + r.stderr


def test_just_under_the_default_threshold_stays_quiet(tmp_path: Path) -> None:
    """89% must not fire, or the guard becomes the noise it exists to remove."""
    r = _run(tmp_path, _loop("perog", n_ok=11, n_fail=89))
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr


def test_small_n_cannot_trip_the_default_without_total_failure(tmp_path: Path) -> None:
    """Below N=10, >0.9 requires 100% -- so 2-of-3 failing stays quiet."""
    r = _run(tmp_path, _loop("perog", n_ok=1, n_fail=2))
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr


# --------------------------------------------------------------------------
# reporting quality -- the message must be actionable
# --------------------------------------------------------------------------
def test_escalation_names_the_label_and_the_rate(tmp_path: Path) -> None:
    r = _run(tmp_path, _loop("perog_trees", n_ok=0, n_fail=8))
    combined = r.stdout + r.stderr
    assert "perog_trees" in combined
    assert "8" in combined, "the message must carry the observed counts"


def test_escalation_shows_example_failed_items(tmp_path: Path) -> None:
    """Naming a few offenders is what turns an alarm into a diagnosis."""
    r = _run(tmp_path, _loop("perog", n_ok=0, n_fail=5))
    assert "item_bad_1" in r.stdout + r.stderr


# --------------------------------------------------------------------------
# configuration + isolation
# --------------------------------------------------------------------------
def test_custom_threshold_is_honoured(tmp_path: Path) -> None:
    r = _run(tmp_path, _loop("perog", n_ok=70, n_fail=30, threshold="0.25"))
    assert "GUARD_SATURATED" in r.stdout, r.stdout + r.stderr


def test_custom_threshold_can_be_made_stricter_than_default(tmp_path: Path) -> None:
    """Opting into 1.0 keeps the guard to the strictly unambiguous case."""
    r = _run(tmp_path, _loop("perog", n_ok=5, n_fail=95, threshold="1.0"))
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr


def test_zero_items_is_not_a_saturation_error(tmp_path: Path) -> None:
    """An empty loop has no rate; emptiness is the cardinality check's job."""
    r = _run(tmp_path, """
saturation_begin perog
saturation_report perog && echo GUARD_OK || echo GUARD_SATURATED
""")
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr
    assert "ERROR" not in r.stdout + r.stderr


def test_two_loops_do_not_contaminate_each_other(tmp_path: Path) -> None:
    r = _run(tmp_path, """
saturation_begin trees
saturation_begin aligns
for i in $(seq 1 10); do saturation_record trees fail "t$i"; done
for i in $(seq 1 10); do saturation_record aligns ok "a$i"; done
saturation_report trees  && echo TREES_OK  || echo TREES_SATURATED
saturation_report aligns && echo ALIGNS_OK || echo ALIGNS_SATURATED
""")
    assert "TREES_SATURATED" in r.stdout, r.stdout + r.stderr
    assert "ALIGNS_OK" in r.stdout


def test_begin_resets_a_reused_label(tmp_path: Path) -> None:
    """A second loop under the same label must not inherit the first's tally."""
    r = _run(tmp_path, """
saturation_begin perog
for i in $(seq 1 10); do saturation_record perog fail "x$i"; done
saturation_report perog >/dev/null 2>&1
saturation_begin perog
for i in $(seq 1 10); do saturation_record perog ok "y$i"; done
saturation_report perog && echo GUARD_OK || echo GUARD_SATURATED
""")
    assert "GUARD_OK" in r.stdout, r.stdout + r.stderr


def test_recording_without_begin_is_a_usage_error(tmp_path: Path) -> None:
    """Silently inventing a counter would recreate the invisible-failure class."""
    r = _run(tmp_path, 'saturation_record perog fail item1 && echo REC_OK || echo REC_FAIL')
    assert "REC_FAIL" in r.stdout, r.stdout + r.stderr


def test_invalid_outcome_is_a_usage_error(tmp_path: Path) -> None:
    r = _run(tmp_path, """
saturation_begin perog
saturation_record perog maybe item1 && echo REC_OK || echo REC_FAIL
""")
    assert "REC_FAIL" in r.stdout, r.stdout + r.stderr


def test_guard_works_with_only_functions_sh(tmp_path: Path) -> None:
    r = _run(tmp_path, _loop("perog", n_ok=0, n_fail=4))
    assert "command not found" not in (r.stdout + r.stderr).lower()
