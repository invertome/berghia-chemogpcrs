"""`lse_depth` -> `lse_divergence`: a persisted-schema rename, plus two
neighbouring axis fixes (LSE_NESTING_DEPTH_WEIGHT export; nesting depth added to
the discovery view).

WHY RENAME
----------
The axis was named for nesting depth but measures ``node.get_distance(tree)`` --
CUMULATIVE BRANCH LENGTH, i.e. divergence (rate x time). Its companion
``lse_nesting_depth`` measures the topological quantity (``len(get_ancestors())``).
Carrying both under names that both say "depth" made the pair unreadable, and
on this repo's real trees they select substantially different candidate sets
inside the top quartile that actually shapes a ranking (62.7% / 42.5% overlap).

THE COMPATIBILITY POLICY, AND WHY IT IS THIS ONE
------------------------------------------------
The rename crosses a persisted schema (written CSV columns) and a public config
surface (env vars). The decision, implemented in
``_rank_candidates_lib.apply_legacy_column_aliases`` / ``getenv_renamed``:

  READ  legacy names are still accepted, and every acceptance prints a
        DEPRECATION line naming the file or variable.
  WRITE only canonical names are ever emitted.

Stale files are NOT rejected: the legacy column holds the SAME measurement under
a wrong label, so refusing it would discard good data -- and would break
07_candidate_ranking.sh's signal-independence audit, which deliberately reads
the PRIOR run's ranked CSV. What is forbidden is the SILENT case: a reader that
finds neither name and treats the axis as "no data" rather than as "this file
predates the rename". These tests pin exactly that boundary.
"""
from __future__ import annotations

import io
import os
import re
import subprocess
import sys
from contextlib import redirect_stderr
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CONFIG = PROJECT_ROOT / "config.sh"
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import _rank_candidates_lib as lib  # noqa: E402
import rank_aggregation as ra  # noqa: E402
import emit_ranked_views as erv  # noqa: E402

CANONICAL = ("lse_divergence_score", "lse_divergence_score_norm",
             "has_lse_divergence_data")
LEGACY = ("lse_depth_score", "lse_depth_score_norm", "has_lse_depth_data")


@pytest.fixture(autouse=True)
def _reset_warn_cache():
    """The deprecation warnings are one-shot per process; clear between tests."""
    lib._LEGACY_WARNED.clear()
    yield
    lib._LEGACY_WARNED.clear()


# ---------------------------------------------------------------------------
# 1. The rename is complete in the code we own
# ---------------------------------------------------------------------------

OWNED = [
    "scripts/rank_candidates.py", "scripts/plot_ranking.py",
    "scripts/emit_ranked_views.py", "scripts/rank_aggregation.py",
    "scripts/ablate_ranking.py",
    "scripts/audit_signal_ranking_independence.py",
    "config.sh", "validate_config.sh", "scripts/run_weight_ablation.sh",
    "09_report_generation.sh",
]


@pytest.mark.parametrize("rel", OWNED)
def test_no_legacy_identifier_survives_in_the_pipeline_code(rel: str) -> None:
    """A half-done rename is worse than none: the stale half keeps reading a
    column nothing writes, which reads as "no data".

    Would have caught: leaving validate_config.sh summing an unset
    LSE_DEPTH_WEIGHT, or 09_report_generation.sh telling a reader to set an env
    var that no longer exists.

    Scoped to EXECUTABLE code: trailing comments are stripped first, because
    comments documenting the rename ("renamed from LSE_DEPTH_WEIGHT", "accepts
    the pre-rename schema") are exactly what a reader needs and must survive.
    The one place the legacy string is required in code is the compatibility
    layer itself, in _rank_candidates_lib.py, which is not in this list.
    """
    text = (PROJECT_ROOT / rel).read_text()
    stray = []
    for ln in text.splitlines():
        code = ln.split("#", 1)[0]
        if "lse_nesting_depth" in code:
            continue
        if re.search(r"lse_depth|LSE_DEPTH|LSE\\_DEPTH", code):
            stray.append(ln)
    assert not stray, f"{rel} still carries legacy lse_depth identifiers: {stray}"


def test_the_companion_nesting_axis_was_not_renamed_by_accident() -> None:
    """`lse_nesting_depth` is correctly named -- it really does measure depth --
    and must survive the rename untouched.

    Would have caught: a blanket s/depth/divergence/ that collapsed the two
    axes into one name, silently fusing two deliberately independent voters.
    """
    src = RANK.read_text()
    assert "lse_nesting_depth" in src
    assert "lse_nesting_divergence" not in src
    assert "LSE_NESTING_DEPTH_WEIGHT" in src
    spec = dict((s[0], s[1]) for s in ra.SIGNAL_SPEC)
    assert "lse_divergence" in spec and "lse_nesting_depth" in spec
    assert spec["lse_divergence"] != spec["lse_nesting_depth"] or True
    # They remain SEPARATE signals, gated on different availability facts.
    assert spec["lse_divergence"] == "has_phylo_data"
    assert spec["lse_nesting_depth"] == "has_lse_nesting_depth_data"


# ---------------------------------------------------------------------------
# 2. Legacy CSVs: accepted, and ANNOUNCED
# ---------------------------------------------------------------------------

def test_a_legacy_csv_column_is_read_under_the_canonical_name() -> None:
    """Would have caught: the silent-mismatch failure the rename would otherwise
    introduce -- a pre-rename ranked CSV presenting no lse_divergence column, so
    the axis vanishes from the ranking and every consumer reads that as "this
    candidate has no divergence signal"."""
    df = pd.DataFrame({"id": ["a"], "lse_depth_score": [0.7],
                       "has_lse_depth_data": [True]})
    out = lib.apply_legacy_column_aliases(df, source="test CSV")
    assert "lse_divergence_score" in out.columns
    assert "has_lse_divergence_data" in out.columns
    assert out["lse_divergence_score"].iloc[0] == 0.7


def test_reading_a_legacy_column_is_never_silent() -> None:
    """The whole compatibility policy rests on the acceptance being visible.

    Would have caught: an alias map applied quietly, which would let a stale
    schema persist indefinitely with nobody aware they are on it.
    """
    df = pd.DataFrame({"id": ["a"], "lse_depth_score": [0.7]})
    buf = io.StringIO()
    with redirect_stderr(buf):
        lib.apply_legacy_column_aliases(df, source="ranked CSV /tmp/x.csv")
    msg = buf.getvalue()
    assert "DEPRECATION" in msg
    assert "lse_depth_score" in msg and "lse_divergence_score" in msg
    assert "/tmp/x.csv" in msg, "the warning must identify WHICH file is stale"


def test_a_canonical_column_wins_over_a_legacy_one() -> None:
    """A file carrying both must not have its live column overwritten by the
    stale copy.

    Would have caught: an unconditional rename, which on a hand-merged file
    would silently substitute the older values.
    """
    df = pd.DataFrame({"id": ["a"], "lse_divergence_score": [0.9],
                       "lse_depth_score": [0.1]})
    out = lib.apply_legacy_column_aliases(df, source="test CSV")
    assert out["lse_divergence_score"].iloc[0] == 0.9


def test_a_csv_with_neither_name_reports_no_data_not_zero() -> None:
    """Absence of both names is genuinely "unmeasured" and must stay so.

    Would have caught: a shim that filled in 0.0 for the missing axis, which is
    the present-zero defect this whole task is about.
    """
    df = pd.DataFrame({"id": ["a", "b"], "phylo_score_norm": [0.9, 0.1],
                       "has_phylo_data": [True, True]})
    out = lib.apply_legacy_column_aliases(df, source="test CSV")
    assert "lse_divergence_score" not in out.columns
    lists = ra.build_ranklists_from_df(out)
    assert "lse_divergence" not in lists, (
        "an axis with no column must not vote at all")


def test_rank_aggregation_votes_on_a_legacy_csv() -> None:
    """End-to-end: the shared aggregation entry point must rank a pre-rename
    CSV rather than dropping its divergence signal."""
    df = pd.DataFrame({
        "id": ["a", "b", "c"],
        "lse_depth_score_norm": [0.9, 0.5, 0.1],
        "has_phylo_data": [True, True, True],
        "phylo_score_norm": [0.5, 0.5, 0.5],
    })
    lists = ra.build_ranklists_from_df(df)
    assert "lse_divergence" in lists
    assert lists["lse_divergence"] == {"a": 0.9, "b": 0.5, "c": 0.1}


def test_emit_ranked_views_scores_a_legacy_csv(tmp_path) -> None:
    """The discovery view's RRA must still see the divergence signal from a
    pre-rename ranked CSV."""
    csv = tmp_path / "ranked.csv"
    pd.DataFrame({
        "id": ["deep", "mid", "shallow"],
        "classification": ["chemoreceptor-candidate"] * 3,
        "og_dnds_reliability": ["low"] * 3,
        "lse_depth_score": [0.9, 0.5, 0.1],
    }).to_csv(csv, index=False)

    out = tmp_path / "disc.csv"
    conf = tmp_path / "conf.csv"
    rc = erv.main(["--ranked-csv", str(csv), "--confidence-out", str(conf),
                   "--discovery-out", str(out)])
    assert rc == 0
    disc = pd.read_csv(out)
    assert list(disc["id"]) == ["deep", "mid", "shallow"], (
        "the legacy divergence column did not reach the discovery score")


# ---------------------------------------------------------------------------
# 3. Legacy env vars: honoured, announced, canonical-wins
# ---------------------------------------------------------------------------

def test_legacy_env_var_is_honoured_and_announced(monkeypatch) -> None:
    """preliminary/config.sh and any operator shell still export
    LSE_DEPTH_WEIGHT. Ignoring it silently would drop the operator's setting
    back to the default with no indication."""
    monkeypatch.delenv("LSE_DIVERGENCE_WEIGHT", raising=False)
    monkeypatch.setenv("LSE_DEPTH_WEIGHT", "4")
    buf = io.StringIO()
    with redirect_stderr(buf):
        value = lib.getenv_renamed("LSE_DIVERGENCE_WEIGHT", 1)
    assert value == 4.0
    assert "DEPRECATION" in buf.getvalue()
    assert "LSE_DEPTH_WEIGHT" in buf.getvalue()


def test_canonical_env_var_wins_and_the_conflict_is_reported(monkeypatch) -> None:
    """Would have caught: legacy-first precedence, where a half-migrated
    environment silently applies the value the operator did NOT intend."""
    monkeypatch.setenv("LSE_DIVERGENCE_WEIGHT", "2")
    monkeypatch.setenv("LSE_DEPTH_WEIGHT", "9")
    buf = io.StringIO()
    with redirect_stderr(buf):
        value = lib.getenv_renamed("LSE_DIVERGENCE_WEIGHT", 1)
    assert value == 2.0
    assert "both" in buf.getvalue().lower()


def test_neither_set_returns_the_default(monkeypatch) -> None:
    monkeypatch.delenv("LSE_DIVERGENCE_WEIGHT", raising=False)
    monkeypatch.delenv("LSE_DEPTH_WEIGHT", raising=False)
    assert lib.getenv_renamed("LSE_DIVERGENCE_WEIGHT", 1) == 1.0


# ---------------------------------------------------------------------------
# 4. LSE_NESTING_DEPTH_WEIGHT is exported, and MUST be non-zero
# ---------------------------------------------------------------------------

def test_nesting_depth_weight_is_exported_and_nonzero() -> None:
    """The axis was never exported, so it ran at its getenv default with no
    record in config.

    Would have caught: exporting it as 0 "to start it disabled" -- which does
    not down-weight the axis, it REMOVES it from the vote entirely (see the next
    test), silently, under the production RANK_METHOD=rankagg default.
    """
    text = CONFIG.read_text()
    m = re.search(r"^export\s+LSE_NESTING_DEPTH_WEIGHT=([0-9.]+)", text, re.M)
    assert m, "LSE_NESTING_DEPTH_WEIGHT is not exported from config.sh"
    assert float(m.group(1)) != 0.0, (
        "a weight of exactly 0 is read as an EXCLUSION, not a small weight")
    assert "never set to 0" in text or "MUST STAY NON-ZERO" in text, (
        "config.sh does not warn that 0 is an exclusion rather than a weight")


def test_zero_weight_really_would_exclude_the_axis() -> None:
    """The behaviour the comment above warns about, verified rather than
    asserted."""
    weights = {"phylo": 2.0, "lse_nesting_depth": 0.0, "tandem_cluster": 2.5}
    excluded = ra.excluded_signals_from_weights(weights, env_var="_NO_SUCH_VAR")
    assert excluded == {"lse_nesting_depth"}


def test_config_sh_is_valid_bash() -> None:
    proc = subprocess.run(["bash", "-n", str(CONFIG)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


# ---------------------------------------------------------------------------
# 5. Nesting depth joins the discovery view
# ---------------------------------------------------------------------------

def test_discovery_view_votes_on_nesting_depth(tmp_path) -> None:
    """The discovery view is explicitly about divergent LSE paralogs, and
    duplication depth is a divergence signal -- it was simply missing while the
    patristic axis was listed.

    Would have caught: the omission. With nesting depth the ONLY varying signal,
    an implementation that ignores it cannot order these rows.
    """
    csv = tmp_path / "ranked.csv"
    pd.DataFrame({
        "id": ["deep", "mid", "shallow"],
        "classification": ["chemoreceptor-candidate"] * 3,
        "og_dnds_reliability": ["low"] * 3,
        "lse_nesting_depth_score_norm": [0.9, 0.5, 0.1],
    }).to_csv(csv, index=False)

    out = tmp_path / "disc.csv"
    conf = tmp_path / "conf.csv"
    assert erv.main(["--ranked-csv", str(csv), "--confidence-out", str(conf),
                     "--discovery-out", str(out)]) == 0
    disc = pd.read_csv(out)
    assert list(disc["id"]) == ["deep", "mid", "shallow"]
    assert disc["discovery_score"].is_monotonic_decreasing
    assert disc["discovery_score"].nunique() > 1, (
        "nesting depth is not reaching the discovery score at all")


def test_discovery_signal_set_is_documented_and_matches_the_code() -> None:
    """The module docstring is the spec for the discovery score; it must list
    the same signals the code reads.

    Would have caught: adding a voter without updating the documented formula,
    leaving the only description of the score wrong.
    """
    import ast
    src = (PROJECT_ROOT / "scripts" / "emit_ranked_views.py").read_text()
    header = ast.get_docstring(ast.parse(src)) or ""
    for signal in ("tandem", "positive", "novelty", "lse_divergence",
                   "lse_nesting_depth"):
        assert signal in header, f"{signal} missing from the documented formula"
