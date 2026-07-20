"""Regression test for bead -6hl: rank_candidates.py must emit an
``orthogroup`` column so that add_og_coverage_columns.py can join OG
coverage data.

Before the fix, ``rank_candidates.py`` computed ``cand_og`` internally but
never included it in either ``output_cols`` list, and never appended it to
each candidate's result dict.  That caused ``add_og_coverage_columns`` to
hit its ``if og_column not in df.columns`` branch and silently assign
og_n_ref_cds=0, og_n_total=0, og_dnds_reliability='low' for **every** row.

Strategy: approach (b) — test the downstream contract directly.
Given a ranked CSV that *already* has an ``orthogroup`` column (the
column rank_candidates.py should produce), confirm that
``add_og_coverage_columns.add_coverage_columns`` returns non-zero
``og_n_total`` / non-'low' ``og_dnds_reliability`` for a well-populated OG.

This proves:
    1. The column name expected by ``add_og_coverage_columns`` is exactly
       ``'orthogroup'`` (not ``'og'`` or ``'orthogroup_id'`` etc.)
    2. Once ``rank_candidates.py`` emits that column the entire contract
       works end-to-end.

A companion test at the bottom (test_rank_candidates_result_dict_has_orthogroup)
does a lighter-weight check: it reads the source of ``rank_candidates.py``,
extracts the fragment that builds the per-candidate result dict, executes it
in a minimal namespace, and asserts that ``'orthogroup'`` is a key in the
produced dict.  This is cheaper than running ``main()`` and directly targets
the line that was missing.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd
import pytest

# conftest.py adds scripts/ to sys.path
import add_og_coverage_columns as aoc


# ---------------------------------------------------------------------------
# Approach (b): downstream contract — add_og_coverage_columns must work when
# the ranked CSV contains an 'orthogroup' column.
# ---------------------------------------------------------------------------

def test_add_coverage_joins_when_orthogroup_column_present(tmp_path: Path) -> None:
    """When the ranked CSV has an 'orthogroup' column populated with real OG
    IDs, add_coverage_columns must produce non-zero og_n_total and a
    reliability label that is NOT always 'low'.

    This is the contract that rank_candidates.py's orthogroup column must
    satisfy.  Pre-fix: the column was absent and every row got og_n_total=0,
    og_dnds_reliability='low'.  Post-fix: with the column present, rows in
    well-covered OGs get the correct reliability.
    """
    # Build synthetic reference CDS FASTA with 10 sequences all in OG0000001
    # (so the OG will get HIGH reliability >= 10).
    cds_ids = [f"ref_species_{i}" for i in range(10)]
    cds_text = "".join(f">{gid}\nATGAAA\n" for gid in cds_ids)
    cds_fa = tmp_path / "all_references_cds.fna"
    cds_fa.write_text(cds_text)

    # OrthoGroups.tsv: OG0000001 has all 10 ref genes + 2 Berghia candidates.
    og_tsv = tmp_path / "Orthogroups.tsv"
    og_tsv.write_text(
        "Orthogroup\trefspecies\tberghia\n"
        "OG0000001\t" + ", ".join(cds_ids) + "\tbste_cand_a, bste_cand_b\n"
        "OG0000002\tbste_cand_c\t\n"
    )

    # Ranked CSV — this is what rank_candidates.py SHOULD produce after the fix.
    # It contains 'orthogroup' (which was missing before).
    ranked_csv = tmp_path / "ranked.csv"
    ranked_csv.write_text(
        "id,orthogroup,rank_score\n"
        "bste_cand_a,OG0000001,0.95\n"
        "bste_cand_b,OG0000001,0.92\n"
        "bste_cand_c,OG0000002,0.80\n"
    )
    out_csv = tmp_path / "ranked_with_coverage.csv"

    aoc.add_coverage_columns(
        ranked_csv_path=str(ranked_csv),
        cds_fasta_path=str(cds_fa),
        orthogroups_tsv_path=str(og_tsv),
        out_path=str(out_csv),
    )

    df = pd.read_csv(out_csv)

    # OG0000001 has 10 ref CDS → og_n_total >= 10, og_n_ref_cds = 10 → 'high'
    row_a = df[df["id"] == "bste_cand_a"].iloc[0]
    assert row_a["og_n_ref_cds"] == 10, (
        "Expected 10 reference CDS in OG0000001 but got "
        f"{row_a['og_n_ref_cds']!r} — 'orthogroup' column join failed"
    )
    assert row_a["og_n_total"] == 12, (
        "Expected 12 total OG members (10 ref + 2 Berghia) "
        f"but got {row_a['og_n_total']!r}"
    )
    assert row_a["og_dnds_reliability"] == "high", (
        f"Expected 'high' reliability for 10-member OG, got "
        f"{row_a['og_dnds_reliability']!r}"
    )

    # OG0000002 has 0 ref CDS → low reliability (control row)
    row_c = df[df["id"] == "bste_cand_c"].iloc[0]
    assert row_c["og_dnds_reliability"] == "low"


# ---------------------------------------------------------------------------
# Lightweight source-level check: the result dict in rank_candidates.py must
# contain an 'orthogroup' key.
# ---------------------------------------------------------------------------

def _extract_results_append_dict() -> dict:
    """Read rank_candidates.py source, find the ``results.append({...})``
    block, exec just that literal dict in a minimal namespace, and return it.

    This approach is taken from test_rank_candidates_dnds_reliability.py which
    similarly exec()s fragments of the script to avoid running main().
    """
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.normpath(os.path.join(here, '..', '..'))
    src_path = os.path.join(repo, 'scripts', 'rank_candidates.py')
    with open(src_path) as f:
        src = f.read()

    # Locate the per-candidate `results.append({` — this is at 4-space indent
    # (inside the `for cand_id in candidates:` loop).  Searching for the
    # exact indented form avoids matching `stability_results.append({` or
    # `cv_results.append({` which are substrings of `results.append({` that
    # appear earlier in the file.
    start = src.find('\n    results.append({')
    assert start != -1, "Could not find '    results.append({' in rank_candidates.py"
    start += 1  # skip the leading newline so `start` points to the spaces
    # Find the opening '{' of the dict by scanning forward from start.
    # (Can't use fixed offset because indentation may vary.)
    dict_start = src.index('{', start)
    i = dict_start
    brace_depth = 0
    assert src[i] == '{', f"Expected '{{' at index {i}, got {src[i]!r}"
    end = None
    for j in range(i, len(src)):
        if src[j] == '{':
            brace_depth += 1
        elif src[j] == '}':
            brace_depth -= 1
            if brace_depth == 0:
                end = j + 1
                break
    assert end is not None, "Could not find matching '}' for results.append dict"
    dict_literal = src[dict_start:end]

    # Build a namespace with dummy values for all names the dict references.
    # We only care that 'orthogroup' is a key — not the actual values.
    ns: dict = {
        'cand_id': 'dummy_id',
        'leaf_lookup': {},  # o98: has_phylo = cand_id in leaf_lookup
        # Bound at rank_candidates.py:1918-1919, well before the dict literal at
        # :2040-2041 uses them. This harness execs only the literal, so the
        # loop-local names must be supplied here.
        'has_phylo': False,
        'has_dnds': False,
        'completeness': 0.0,
        'phylo_score': 0.0,
        'purifying_score': 0.0,
        'positive_score': 0.0,
        'selection_significant': False,
        'busted_s_p': float('nan'),
        'busted_s_sig': False,
        'busted_mh_p': float('nan'),
        'busted_mh_sig': False,
        'meme_n_episodic': 0,
        'meme_fraction_episodic': 0.0,
        'meme_high_confidence_sites': 0,
        'meme_lenient_only_sites': 0,
        'meme_strict_only_sites': 0,
        'meme_alignment_robustness_index': float('nan'),
        'synteny_score': 0.0,
        'has_synteny': False,
        'expr_score': 0.0,
        'has_expression': False,
        'lse_score': 0.0,
        'raw_depth': 0.0,
        'has_lse_depth': False,   # hf3u: lse_depth availability flag
        # hf3u: the topological nesting-depth axis, emitted alongside the
        # patristic one (see tests/unit/test_hf3u_nesting_depth_axis.py).
        'lse_nesting_score': 0.0,
        'raw_nesting_depth': 0,
        'has_lse_nesting_depth': False,
        'chemo_expr_score': 0.0,
        'has_chemo_expr': False,
        'gprotein_score': 0.0,
        'has_gprotein': False,
        'ecl_score': 0.0,
        'has_ecl': False,
        'expansion_score': 0.0,
        'has_expansion': False,
        'og_conf': 0.0,
        'has_og_data': False,
        'tandem_score': 0.0,
        'tandem_size': 0,
        'tandem_cluster_label': '',
        'has_tandem': False,
        'cds_source': 'native',
        'og_dnds_reliability': {},
        'cand_og': 'OG_TEST',
        'evidence_completeness': 0.0,
        '__builtins__': __builtins__,
        'float': float,
    }
    result_dict = eval(dict_literal, ns)
    return result_dict


def test_rank_candidates_result_dict_has_orthogroup() -> None:
    """The per-candidate result dict in rank_candidates.py must have an
    'orthogroup' key so that the output CSV contains an 'orthogroup' column.

    Pre-fix: this key was absent.  Post-fix: it must be present and its value
    must reflect ``cand_og`` (which holds the orthogroup name for each
    candidate).
    """
    d = _extract_results_append_dict()
    assert 'orthogroup' in d, (
        "results.append({...}) in rank_candidates.py is missing the "
        "'orthogroup' key.  Add \"'orthogroup': cand_og or '',\" to the dict."
    )
    # Also verify it takes the value from cand_og (not hardcoded empty string).
    assert d['orthogroup'] == 'OG_TEST', (
        f"Expected orthogroup value 'OG_TEST' (from cand_og) but got {d['orthogroup']!r}"
    )


def test_rank_candidates_output_cols_include_orthogroup() -> None:
    """Both output_cols lists in rank_candidates.py must include 'orthogroup'.

    rank_candidates.py has two separate output_cols definitions (one at ~line
    2044 for the main write path and one at ~line 2259 for the post-sensitivity
    re-write path).  Both must include 'orthogroup' so the column survives in
    both code paths.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.normpath(os.path.join(here, '..', '..'))
    src_path = os.path.join(repo, 'scripts', 'rank_candidates.py')
    with open(src_path) as f:
        src = f.read()

    # Find all output_cols = [ ... ] list definitions
    occurrences = [m.start() for m in re.finditer(r'output_cols\s*=\s*\[', src)]
    # There must be exactly ONE. There used to be two: the second was a drifted
    # hand-maintained literal that omitted 'final_rank', and because it was
    # assigned last it won -- so the production rank order never reached the
    # shipped CSV at all. A second copy is the defect, not the design.
    assert len(occurrences) == 1, (
        f"Expected exactly ONE 'output_cols = [' definition in rank_candidates.py "
        f"(a second, drifting copy is how 'final_rank' was lost), "
        f"found {len(occurrences)}"
    )

    for pos in occurrences:
        # Find the end of this list literal
        bracket_depth = 0
        list_start = src.index('[', pos)
        end = None
        for j in range(list_start, len(src)):
            if src[j] == '[':
                bracket_depth += 1
            elif src[j] == ']':
                bracket_depth -= 1
                if bracket_depth == 0:
                    end = j + 1
                    break
        assert end is not None
        list_src = src[list_start:end]
        assert "'orthogroup'" in list_src or '"orthogroup"' in list_src, (
            f"output_cols list at position {pos} in rank_candidates.py does not "
            f"include 'orthogroup'.  Both output_cols lists must include it so "
            f"the column is written in both code paths.\n"
            f"List fragment: {list_src[:200]!r}"
        )


def test_output_cols_include_norm_signals_for_discovery_view() -> None:
    """Bead 1nr: emit_ranked_views.py's discovery score needs
    positive_score_norm and tandem_cluster_score_norm in the ranked CSV. Both
    output_cols lists in rank_candidates.py must include them so the discovery
    view scores on BOTH novelty signals (not tandem-only via the fallback)."""
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.normpath(os.path.join(here, '..', '..'))
    with open(os.path.join(repo, 'scripts', 'rank_candidates.py')) as f:
        src = f.read()
    occurrences = [m.start() for m in re.finditer(r'output_cols\s*=\s*\[', src)]
    assert len(occurrences) == 1  # one definition; a second drifting copy is the defect
    for pos in occurrences:
        bracket_depth = 0
        list_start = src.index('[', pos)
        end = None
        for j in range(list_start, len(src)):
            if src[j] == '[':
                bracket_depth += 1
            elif src[j] == ']':
                bracket_depth -= 1
                if bracket_depth == 0:
                    end = j + 1
                    break
        list_src = src[list_start:end]
        assert "'positive_score_norm'" in list_src, "output_cols missing positive_score_norm"
        assert "'tandem_cluster_score_norm'" in list_src, "output_cols missing tandem_cluster_score_norm"
