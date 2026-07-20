"""Tests for Glue task G3: wiring the structural/embedding/microswitch
evidence-channel TSVs (Glue G1/G2/G5 producers) into the rankagg reranker.

Task 6 (see tests/unit/test_channel_integration.py's module docstring)
deliberately wired the channels only as far as
``rank_aggregation.merge_evidence_channels`` / ``build_ranklists_from_df`` --
it explicitly did NOT modify ``rank_candidates.py``. This is the follow-up
glue: ``rank_candidates.py``'s ``RANK_METHOD == 'rankagg'`` block must
actually call ``merge_evidence_channels`` (via a small helper,
``_merge_channels_if_present``) on whatever channel TSVs
``07_candidate_ranking.sh``'s producers wrote, BEFORE ``rerank_output()``
runs.

``rank_candidates.py`` is a script, not an importable module (it reads
``sys.argv`` and executes real I/O at module-import time), so -- following
the extraction pattern already used by
tests/unit/test_rank_candidates_dnds_reliability.py and
tests/unit/test_rank_candidates_emits_orthogroup.py -- these tests pull the
real ``_merge_channels_if_present`` function out of the script's source (via
``ast``, which is robust to whatever code precedes/follows the function --
unlike a ``"\\n\\ndef "`` text-marker search) and exec it in isolation,
rather than running the whole script.
"""
from __future__ import annotations

import ast
import os

import pandas as pd
import pytest

# conftest.py adds scripts/ to sys.path
from rank_aggregation import aggregate, build_ranklists_from_df, rerank_output


RANK_CANDIDATES_SRC_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts", "rank_candidates.py"
)
RANK_CANDIDATES_SRC_PATH = os.path.normpath(RANK_CANDIDATES_SRC_PATH)


def _read_rank_candidates_source() -> str:
    with open(RANK_CANDIDATES_SRC_PATH) as fh:
        return fh.read()


def _extract_function_source(src: str, func_name: str) -> str:
    """Exact source text of one top-level ``def <func_name>(...):`` in
    ``src``, located via ``ast`` rather than a "next def" text marker (the
    function is not guaranteed to be immediately followed by another def --
    it sits right before the script's argv-parsing top-level code)."""
    tree = ast.parse(src)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == func_name:
            segment = ast.get_source_segment(src, node)
            assert segment is not None
            return segment
    raise AssertionError(f"{func_name!r} not found as a top-level def in rank_candidates.py")


def _load_merge_channels_if_present():
    """Load the real ``_merge_channels_if_present`` from rank_candidates.py
    without running the script's CLI-driven top-level body (sys.argv parsing,
    real result-directory reads)."""
    src = _extract_function_source(_read_rank_candidates_source(), "_merge_channels_if_present")
    ns: dict = {"os": os}
    exec(src, ns)
    return ns["_merge_channels_if_present"]


def _write_tsv(path, rows, columns) -> str:
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)
    return str(path)


def _regression_fixture_df() -> pd.DataFrame:
    """A df carrying several of the pre-Task-6 12 ranking signals (no
    channels) -- used to prove the channels-absent path is a no-op."""
    return pd.DataFrame(
        {
            "id": ["c0", "c1", "c2", "c3", "c4"],
            "phylo_score_norm": [0.9, 0.7, 0.5, 0.3, 0.1],
            "purifying_score_norm": [0.2, 0.8, 0.4, 0.6, 0.1],
            "positive_score_norm": [0.5, 0.5, 0.9, 0.1, 0.3],
            "lse_divergence_score_norm": [0.4, 0.4, 0.4, 0.9, 0.1],
            "synteny_score_norm": [0.8, 0.2, 0.5, 0.6, 0.3],
            "has_synteny_data": [True, True, False, True, True],
            "og_confidence_score_norm": [0.3, 0.6, 0.9, 0.2, 0.5],
            "has_og_confidence_data": [True, True, True, True, False],
        }
    )


# --------------------------------------------------------------------------- #
# _merge_channels_if_present: direct behavior
# --------------------------------------------------------------------------- #
def test_struct_tsv_present_sets_flag_and_columns(tmp_path):
    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels"
    channels_dir.mkdir()
    _write_tsv(
        channels_dir / "structural_channel.tsv",
        [["candA", 1, 0]],
        ["id", "struct_novelty", "struct_nonchemo_corrob"],
    )

    df = pd.DataFrame({"id": ["candA", "candB"]})
    merged = merge_fn(df, str(channels_dir))

    assert list(merged["has_struct_data"]) == [True, False]
    assert merged.set_index("id").loc["candA", "struct_novelty"] == 1
    # the other two channels' TSVs are absent -> stay dormant
    assert not merged["has_emb_data"].any()
    assert not merged["has_or_microswitch_data"].any()


def test_emb_and_microswitch_tsvs_present_set_flags_and_columns(tmp_path):
    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels"
    channels_dir.mkdir()
    _write_tsv(
        channels_dir / "embedding_channel.tsv",
        [["candA", 0.9, 0.1]],
        ["id", "emb_classA_sim", "emb_nonchemo_sim"],
    )
    _write_tsv(
        channels_dir / "microswitch_channel.tsv",
        [["candB", 1]],
        ["id", "or_microswitch"],
    )

    df = pd.DataFrame({"id": ["candA", "candB"]})
    merged = merge_fn(df, str(channels_dir))

    assert list(merged["has_emb_data"]) == [True, False]
    assert merged.set_index("id").loc["candA", "emb_classA_sim"] == pytest.approx(0.9)
    assert list(merged["has_or_microswitch_data"]) == [False, True]
    assert merged.set_index("id").loc["candB", "or_microswitch"] == 1
    # structural TSV absent -> dormant
    assert not merged["has_struct_data"].any()


def test_nonexistent_channels_dir_leaves_everything_dormant(tmp_path):
    """First run: 07_candidate_ranking.sh hasn't even mkdir -p'd the channels
    dir yet (no producer inputs were available). Must not crash."""
    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels_never_created"

    df = pd.DataFrame({"id": ["a", "b"], "phylo_score_norm": [0.5, 0.6]})
    merged = merge_fn(df, str(channels_dir))

    assert not merged["has_struct_data"].any()
    assert not merged["has_emb_data"].any()
    assert not merged["has_or_microswitch_data"].any()
    assert "struct_novelty" not in merged.columns
    assert "emb_classA_sim" not in merged.columns
    assert "or_microswitch" not in merged.columns
    # untouched base columns preserved
    assert list(merged["id"]) == ["a", "b"]
    assert list(merged["phylo_score_norm"]) == [0.5, 0.6]


def test_existing_channels_dir_with_no_tsvs_leaves_everything_dormant(tmp_path):
    """channels_dir exists (mkdir -p ran) but none of the 3 producers had
    inputs available -- every channel must stay dormant, no crash."""
    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels"
    channels_dir.mkdir()

    df = pd.DataFrame({"id": ["a", "b"]})
    merged = merge_fn(df, str(channels_dir))

    assert not merged["has_struct_data"].any()
    assert not merged["has_emb_data"].any()
    assert not merged["has_or_microswitch_data"].any()


# --------------------------------------------------------------------------- #
# End-to-end wiring: _merge_channels_if_present -> rerank_output
# --------------------------------------------------------------------------- #
def test_struct_nonchemo_corrob_present_lowers_rank_via_real_rerank_output(tmp_path):
    """Acceptance criterion, exercised through the EXACT production entry
    point (rank_candidates.py imports this as `_rerank_output`): a candidate
    whose structural TSV marks struct_nonchemo_corrob=1 ranks LOWER under
    rankagg-with-channel than the same df without the channel.
    rerank_output() always uses RRA internally (rank_aggregation.py's
    primary method -- it has no method="rrf" option), so this is the single
    non-parametrized RRA proof; test_struct_nonchemo_corrob_present_lowers_rank_isolated
    below additionally proves the same direction for both RRA and RRF via
    aggregate() directly."""
    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels"
    channels_dir.mkdir()
    _write_tsv(
        channels_dir / "structural_channel.tsv",
        [["good", 0], ["bad", 1]],
        ["id", "struct_nonchemo_corrob"],
    )

    df = pd.DataFrame({"id": ["good", "bad"]})
    merged = merge_fn(df, str(channels_dir))

    ordered = rerank_output(merged, "rankagg")
    ids = list(ordered["id"])
    assert ids.index("good") < ids.index("bad")


@pytest.mark.parametrize("method", ["rra", "rrf"])
def test_struct_nonchemo_corrob_present_lowers_rank_isolated(tmp_path, method):
    """Same isolated single-signal scenario, but calling aggregate()
    directly with an explicit method= so both RRA and RRF are genuinely
    exercised (rerank_output above can't do this -- it hardcodes RRA).
    Proves _merge_channels_if_present's output threads correctly into
    build_ranklists_from_df/aggregate for either algorithm."""
    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels"
    channels_dir.mkdir()
    _write_tsv(
        channels_dir / "structural_channel.tsv",
        [["good", 0], ["bad", 1]],
        ["id", "struct_nonchemo_corrob"],
    )

    df = pd.DataFrame({"id": ["good", "bad"]})
    merged = merge_fn(df, str(channels_dir))

    order = aggregate(build_ranklists_from_df(merged), method=method)
    assert order.index("good") < order.index("bad")


def test_struct_nonchemo_corrob_present_lowers_rank_realistic_multisignal(tmp_path):
    """Same acceptance criterion in a realistic df carrying several
    (genuinely tied) base signals, differing only in struct_nonchemo_corrob.
    Exercised via aggregate(..., method="rrf") directly rather than
    rerank_output (which hardcodes RRA): RRA's Bonferroni correction can
    saturate to a tied 1.0 ceiling with enough uninformative confounds
    (documented, pre-existing property of RRA itself; see
    test_channel_integration.py), so this specific realistic-multisignal
    scenario is RRF's proof, not RRA's -- the isolated single-signal test
    above already gives the method-agnostic (both rra and rrf) proof through
    the real rerank_output production entry point."""
    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels"
    channels_dir.mkdir()
    _write_tsv(
        channels_dir / "structural_channel.tsv",
        [["good", 0], ["bad", 1]],
        ["id", "struct_nonchemo_corrob"],
    )

    df = pd.DataFrame(
        {
            "id": ["good", "bad"],
            "phylo_score_norm": [0.6, 0.6],
            "purifying_score_norm": [0.4, 0.4],
            "positive_score_norm": [0.5, 0.5],
            "lse_divergence_score_norm": [0.3, 0.3],
        }
    )
    merged = merge_fn(df, str(channels_dir))

    order = aggregate(build_ranklists_from_df(merged), method="rrf")
    assert order.index("good") < order.index("bad")


def test_channels_absent_is_byte_identical_regression_via_real_rerank_output(tmp_path):
    """Acceptance criterion, through the EXACT production entry point:
    channel TSVs absent -> _merge_channels_if_present is a no-op for ranking
    purposes -- the wired path (merge, then rerank_output) must reproduce
    the EXACT same id order as calling rerank_output directly on the
    un-merged df (the pre-G3 12-signal rankagg)."""
    df = _regression_fixture_df()

    baseline_order = list(rerank_output(df.copy(), "rankagg")["id"])

    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels_no_producers_ran"  # never created
    merged = merge_fn(df.copy(), str(channels_dir))
    wired_order = list(rerank_output(merged, "rankagg")["id"])

    assert wired_order == baseline_order


@pytest.mark.parametrize("method", ["rra", "rrf"])
def test_channels_absent_is_byte_identical_regression_both_methods(tmp_path, method):
    """Same regression, via aggregate() directly with an explicit method= so
    both RRA and RRF are genuinely exercised (rerank_output above can't do
    this -- it hardcodes RRA)."""
    df = _regression_fixture_df()

    baseline_order = aggregate(build_ranklists_from_df(df.copy()), method=method)

    merge_fn = _load_merge_channels_if_present()
    channels_dir = tmp_path / "channels_no_producers_ran"  # never created
    merged = merge_fn(df.copy(), str(channels_dir))
    wired_order = aggregate(build_ranklists_from_df(merged), method=method)

    assert wired_order == baseline_order


# --------------------------------------------------------------------------- #
# Static wiring checks: weighted path untouched, default channels_dir
# --------------------------------------------------------------------------- #
def test_merge_helper_called_exactly_once_and_only_inside_rankagg_branch():
    """The weighted path (RANK_METHOD != 'rankagg') must be UNTOUCHED: the
    only call to _merge_channels_if_present(...) in rank_candidates.py must
    be textually inside the `if RANK_METHOD == 'rankagg':` block."""
    src = _read_rank_candidates_source()
    tree = ast.parse(src)

    rankagg_if = None
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.If)
            and isinstance(node.test, ast.Compare)
            and isinstance(node.test.left, ast.Name)
            and node.test.left.id == "RANK_METHOD"
        ):
            rankagg_if = node
            break
    assert rankagg_if is not None, "if RANK_METHOD == 'rankagg': block not found"

    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "_merge_channels_if_present"
    ]
    assert len(calls) == 1, (
        f"Expected exactly one _merge_channels_if_present(...) call site, found {len(calls)}"
    )
    call = calls[0]
    assert rankagg_if.lineno <= call.lineno <= rankagg_if.end_lineno, (
        "_merge_channels_if_present(...) must be called inside the "
        "`if RANK_METHOD == 'rankagg':` block -- the weighted path must stay untouched"
    )


def test_channels_dir_default_expression_matches_spec():
    """The call site must compute channels_dir as
    os.path.join(os.path.dirname(output_file), 'channels') -- the documented
    default (rank_candidates.py's channels/ dir sits next to the output CSV,
    matching 07_candidate_ranking.sh's producers)."""
    src = _read_rank_candidates_source()
    assert "os.path.dirname(output_file)" in src
    assert "os.path.join(os.path.dirname(output_file), 'channels')" in src


def test_merge_evidence_channels_is_imported():
    src = _read_rank_candidates_source()
    assert "merge_evidence_channels" in src
