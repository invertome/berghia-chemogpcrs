"""Production row order must follow `final_rank`, not the demoted `rank_score`.

Two coupled defects with one consequence: the emitted CSV's ROW ORDER silently
reverts to the weighted composite even when RANK_METHOD=rankagg is the
production default.

1. THE COLUMN NEVER SURVIVES. rank_candidates.py builds `output_cols` at
   module level twice (the "Update output columns if sensitivity was run"
   block re-assigns a hand-maintained SECOND literal list). The second list
   drifted: it omits `final_rank`, `dnds_reliability_weight` and
   `n_signals_present`. Because the second write
   (``df_sorted[output_cols].to_csv(output_file)``) is unconditional and runs
   LAST, the shipped CSV has no `final_rank` column at all. This is the same
   hand-maintained-second-copy drift documented in
   tests/unit/test_tier_normalization_denominator.py.

   Downstream, emit_ranked_views.py's "honor the PRODUCTION order" branch
   (`if "final_rank" in out.columns`) is therefore DEAD in production and
   falls through to its rank_score fallback.

2. THE SORT IS RE-DERIVED. add_classification_columns.py then re-sorts the
   whole CSV by `rank_score` -- exactly what rank_candidates.py's invariant
   comment forbids, because under rankagg `rank_score` stays the DEMOTED
   weighted value while row order carries the RRA order. With nothing
   reclassified (every factor 1.0) the sort is still destructive: it is a
   pure weighted-order restore.

08_structural_analysis.sh selects the AlphaFold set by ROW POSITION
(``head -n N | cut -d, -f1``), so the structural stage runs on the top-N by
the demoted composite rather than the production ranking.

The fix keeps f2e's classification suppression (suppress, don't drop) but
enacts it on the PRODUCTION order key instead of the weighted score, so it
composes with either RANK_METHOD.

rank_candidates.py is not import-safe (top-level ``sys.argv``), so the
structural guards read the source via ast, following the convention in
tests/unit/test_tier_normalization_denominator.py.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd

import add_classification_columns as acc

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"


# --------------------------------------------------------------------------
# 1. the column must survive to the LAST write
# --------------------------------------------------------------------------

def _output_cols_assignments() -> list[ast.Assign]:
    tree = ast.parse(RANK.read_text())
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and node.col_offset == 0:
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == "output_cols":
                    found.append(node)
    return found


def _literal_cols(node: ast.Assign) -> list[str] | None:
    """Column names if the assignment is a plain list literal, else None."""
    if not isinstance(node.value, ast.List):
        return None
    return [e.value for e in node.value.elts if isinstance(e, ast.Constant)]


def test_output_cols_defined_exactly_once_as_a_literal() -> None:
    """A second hand-maintained literal list is how `final_rank` got dropped.

    Exactly one module-level assignment may spell the column list out; any
    later assignment must DERIVE from it (e.g. `output_cols + [...]`) so the
    two can never drift apart again.
    """
    literals = [n for n in _output_cols_assignments() if _literal_cols(n) is not None]
    assert len(literals) == 1, (
        f"{len(literals)} module-level `output_cols` list literals at lines "
        f"{[n.lineno for n in literals]}; a second copy drifts (it already "
        f"dropped final_rank). Derive later assignments from the first."
    )


def test_final_write_carries_the_production_order_columns() -> None:
    """Whatever `output_cols` holds at the LAST write must include final_rank.

    The last `df_sorted[output_cols].to_csv(output_file)` is the shipped
    artifact. Union the literals so this fails while a drifted second copy
    omits the column, and passes once there is a single source of truth.
    """
    literals = _output_cols_assignments()
    last = literals[-1]
    cols = _literal_cols(last)
    if cols is None:  # derived form: check the base literal it extends
        cols = _literal_cols(literals[0]) or []
    for required in ("final_rank", "n_signals_present", "dnds_reliability_weight"):
        assert required in cols, (
            f"{required!r} missing from the output_cols in effect at the final "
            f"CSV write (line {last.lineno}); downstream consumers key on it."
        )


def test_purifying_score_norm_is_emitted_alongside_positive() -> None:
    """`positive_score_norm` is emitted but `purifying_score_norm` is not.

    Both are shrunk by the same dN/dS reliability multiplier upstream, so a
    consumer re-reading the CSV compares a shrunk `positive_score_norm`
    against an UNSHRUNK `purifying_score` -- which is why the signal
    independence audit cannot see that one confirmed shared-multiplier
    coupling.
    """
    cols = _literal_cols(_output_cols_assignments()[0]) or []
    assert "positive_score_norm" in cols  # guard: the asymmetry's other half
    assert "purifying_score_norm" in cols, (
        "purifying_score_norm missing while positive_score_norm is present; "
        "the shared dN/dS reliability shrink is then invisible to consumers."
    )


# --------------------------------------------------------------------------
# 2. add_classification_columns must not re-derive order from rank_score
# --------------------------------------------------------------------------

def _ranked_csv(tmp_path: Path) -> Path:
    """A CSV whose production order (final_rank) DISAGREES with rank_score.

    This is exactly the rankagg case: rank_score stays the weighted composite
    while final_rank carries the RRA order.
    """
    csv = tmp_path / "ranked.csv"
    csv.write_text(
        "id,final_rank,rank_score,rank\n"
        "cand_A,1,0.10,1\n"
        "cand_B,2,0.90,2\n"
        "cand_C,3,0.50,3\n"
    )
    return csv


def _empty_consensus(tmp_path: Path) -> Path:
    tsv = tmp_path / "classifications.tsv"
    tsv.write_text(
        "candidate_id\tclassification\tclassification_confidence\t"
        "classification_family\tclassification_subfamily\t"
        "classification_evidence\n"
    )
    return tsv


def test_nothing_reclassified_preserves_production_order(tmp_path: Path) -> None:
    """With every suppression factor 1.0 the pass must be order-preserving.

    Reproduces the audit: row order comes out [cand_B, cand_C, cand_A]
    (weighted order) instead of [cand_A, cand_B, cand_C] (final_rank order).
    """
    out = tmp_path / "out.csv"
    acc.add_classification_columns(
        str(_ranked_csv(tmp_path)), str(_empty_consensus(tmp_path)), str(out))

    df = pd.read_csv(out)
    assert list(df["id"]) == ["cand_A", "cand_B", "cand_C"], (
        "row order was re-derived from the demoted rank_score; "
        "08_structural_analysis.sh selects AlphaFold targets by row position."
    )


def test_suppression_demotes_but_keeps_production_order_within_a_class(
        tmp_path: Path) -> None:
    """f2e suppression must still work -- and still be a suppression.

    A classified non-chemoreceptor drops below the unclassified candidates,
    but the surviving candidates keep their PRODUCTION relative order rather
    than being re-sorted into weighted order.
    """
    csv = tmp_path / "ranked.csv"
    csv.write_text(
        "id,final_rank,rank_score,rank\n"
        "cand_A,1,0.10,1\n"
        "cand_B,2,0.90,2\n"
        "cand_C,3,0.50,3\n"
    )
    tsv = tmp_path / "classifications.tsv"
    tsv.write_text(
        "candidate_id\tclassification\tclassification_confidence\t"
        "classification_family\tclassification_subfamily\t"
        "classification_evidence\n"
        "cand_A\tnon-chemoreceptor\thigh\taminergic\t5HT\thmm:aminergic\n"
    )
    out = tmp_path / "out.csv"
    acc.add_classification_columns(str(csv), str(tsv), str(out))

    df = pd.read_csv(out)
    assert list(df["id"]) == ["cand_B", "cand_C", "cand_A"], (
        "the suppressed candidate must sink to the tail while B and C keep "
        "their production (final_rank) order"
    )
    # suppressed, not dropped -- still recoverable from the full list
    assert len(df) == 3


def test_final_rank_is_renumbered_and_the_prior_order_preserved(
        tmp_path: Path) -> None:
    """After re-ranking, final_rank must describe the NEW production order.

    emit_ranked_views.py sorts its confidence view by final_rank, so a stale
    final_rank would re-scatter the suppressed rows back through the view.
    The pre-suppression order is kept for provenance, mirroring the existing
    `rank_score_prefilter` column.
    """
    csv = tmp_path / "ranked.csv"
    csv.write_text(
        "id,final_rank,rank_score,rank\n"
        "cand_A,1,0.10,1\n"
        "cand_B,2,0.90,2\n"
        "cand_C,3,0.50,3\n"
    )
    tsv = tmp_path / "classifications.tsv"
    tsv.write_text(
        "candidate_id\tclassification\tclassification_confidence\t"
        "classification_family\tclassification_subfamily\t"
        "classification_evidence\n"
        "cand_A\tnon-chemoreceptor\thigh\taminergic\t5HT\thmm:aminergic\n"
    )
    out = tmp_path / "out.csv"
    acc.add_classification_columns(str(csv), str(tsv), str(out))

    df = pd.read_csv(out)
    assert list(df["final_rank"]) == [1, 2, 3]
    assert "final_rank_prefilter" in df.columns
    assert list(df["final_rank_prefilter"]) == [2, 3, 1]


def test_legacy_csv_without_final_rank_still_sorts(tmp_path: Path) -> None:
    """Pre-final_rank CSVs must keep working (rank_score fallback)."""
    csv = tmp_path / "ranked.csv"
    csv.write_text("id,rank_score\ncand_A,0.10\ncand_B,0.90\ncand_C,0.50\n")
    out = tmp_path / "out.csv"
    acc.add_classification_columns(
        str(csv), str(_empty_consensus(tmp_path)), str(out))
    df = pd.read_csv(out)
    assert list(df["id"]) == ["cand_B", "cand_C", "cand_A"]
