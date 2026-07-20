"""Nath quarantine: orthogroup-derived signals must go UNAVAILABLE, not zero.

WHY THE QUARANTINE EXISTS
-------------------------
Measured 2026-07-20: 426 of the 427 FASTAs fed to the stage-03 OrthoFinder run
are Nath et al. reference data (239 from ``one_to_one_ortholog/``, 187 from
``lse/``, verified by identical sequence counts and accessions against
``references/nath_et_al/``). The only non-Nath input is the Berghia file. That
run predates by 13 days the 2026-05-28 decision that Nath references are not
consumed by stages 03/04, and was never redone. So every orthogroup -- and
every signal derived from one -- encodes a dataset the project has decided does
not exist. ``ORTHOLOGY_SOURCE_TRUSTED`` (default 0) quarantines them pending a
re-run on scan-derived candidates; nothing is deleted.

WHAT THESE TESTS PIN, AND WHAT EACH WOULD HAVE CAUGHT
-----------------------------------------------------
The whole hazard of a quarantine is that the obvious implementation -- set the
score to 0 -- is WRONG, and wrong in a way that produces no error. This project
has repeatedly been bitten by exactly that (a present 0.0 read as a measurement:
the ete3 constant-support bug, the ``return 1.0, True`` on no bootstrap data,
the ``lse_threshold = 0.5`` fallback). A quarantined axis must be
indistinguishable from an axis that does not exist:

  * it must not contribute to the score NUMERATOR      (obvious)
  * it must not contribute to the weight DENOMINATOR   (NOT obvious, and this
    is the half that silently breaks things: a zeroed-but-present axis depresses
    every candidate's evidence_completeness by a constant, which cannot reorder
    anything -- so it looks harmless -- but IS decisive at the absolute
    CONFIDENCE_MIN_COMPLETENESS=0.7 gate, where it can quietly empty the
    confidence view.)

Each test below states the specific silent-wrong outcome it blocks.
"""
from __future__ import annotations

import ast
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RANK = PROJECT_ROOT / "scripts" / "rank_candidates.py"
CONFIG = PROJECT_ROOT / "config.sh"
STAGE07 = PROJECT_ROOT / "07_candidate_ranking.sh"
STAGE06C = PROJECT_ROOT / "06c_classify_non_chemoreceptors.sh"

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from _rank_candidates_lib import calculate_fair_rank_score  # noqa: E402
import rank_aggregation as ra  # noqa: E402


# The axes that are derived from the OrthoFinder orthogroups. Traced, not
# assumed: og_confidence reads Gene_Trees/ + Orthogroups.tsv; expansion is CAFE
# over orthogroup gene counts joined by Orthogroups.tsv; purifying/positive come
# from stage 05, which is a SLURM array over Orthogroup_Sequences/ (so every
# omega is estimated on a codon alignment whose membership IS an orthogroup).
QUARANTINED = {"purifying", "positive", "expansion", "og_confidence"}


def _rank_source() -> str:
    return RANK.read_text()


def _function_source(name: str) -> str:
    tree = ast.parse(_rank_source())
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.get_source_segment(_rank_source(), node)
    raise AssertionError(f"{name}() not found in rank_candidates.py")


def _evidence_completeness():
    """Exec just ``evidence_completeness`` out of rank_candidates.py."""
    ns: dict = {}
    exec(compile(_function_source("evidence_completeness"), "<ec>", "exec"), ns)
    return ns["evidence_completeness"]


ALL_AXES = dict(
    has_phylo_data=True, has_dnds_data=True, has_synteny=True,
    has_expression=True, has_lse_data=True, has_chemo_expr=True,
    has_gprotein=True, has_ecl=True, has_expansion=True, has_og_data=True,
    has_tandem=True,
)


# ---------------------------------------------------------------------------
# 1. The switch itself
# ---------------------------------------------------------------------------

def test_quarantine_is_off_by_default_and_documented() -> None:
    """config.sh must default ORTHOLOGY_SOURCE_TRUSTED to 0 and say what would
    make it trustworthy again.

    Would have caught: a quarantine that defaults to TRUSTED, i.e. no quarantine
    at all; and an undocumented flag, which turns lifting it into code
    archaeology instead of a one-line decision.
    """
    text = CONFIG.read_text()
    assert re.search(
        r'export\s+ORTHOLOGY_SOURCE_TRUSTED="\$\{ORTHOLOGY_SOURCE_TRUSTED:-0\}"',
        text), "ORTHOLOGY_SOURCE_TRUSTED must be exported with default 0"
    # The condition for flipping it must be stated, not implied.
    for token in ("Phase 1f", "re-run", "nath_et_al"):
        assert token.lower() in text.lower(), (
            f"config.sh does not state {token!r} as part of the condition that "
            f"makes the orthology source trustworthy again")


def test_rank_candidates_reads_the_switch_and_names_the_axes() -> None:
    """rank_candidates.py must derive QUARANTINED_AXES from the env switch.

    Would have caught: a hard-coded axis list that cannot be turned off, or a
    switch that is read but never applied.
    """
    src = _rank_source()
    assert "ORTHOLOGY_SOURCE_TRUSTED" in src
    m = re.search(r"QUARANTINED_AXES\s*=\s*frozenset\(\)\s*if\s+"
                  r"ORTHOLOGY_SOURCE_TRUSTED\s+else\s+frozenset\(\{(.*?)\}\)",
                  src, re.S)
    assert m, "QUARANTINED_AXES must be empty iff ORTHOLOGY_SOURCE_TRUSTED"
    named = set(re.findall(r"'([a-z_]+)'", m.group(1)))
    assert named == QUARANTINED, (
        f"quarantined axis set drifted: {named} != {QUARANTINED}. dN/dS is "
        f"included because stage 05 is an array over Orthogroup_Sequences/, so "
        f"omega is estimated on orthogroup membership.")


# ---------------------------------------------------------------------------
# 2. Unavailable != zero, in the SCORE
# ---------------------------------------------------------------------------

def test_passing_None_is_NOT_by_itself_enough_to_make_an_axis_absent() -> None:
    """MEASURED, not assumed: in the fair scorer a ``None`` score is numerically
    IDENTICAL to ``0.0`` whenever the completeness floor is not binding.

    The formula is ``(weighted_sum / avail_weight) * (avail_weight /
    total_weight)``, which algebraically collapses to ``weighted_sum /
    total_weight`` -- the avail_weight cancels. So dropping an axis from the
    numerator buys exactly nothing while the axis remains in ``total_weight``.

    This is the single most important fact about this change, and it is
    counter-intuitive: the pre-existing ``has_*_data`` gating, on its own, does
    NOT make an unavailable axis behave like an absent one in the weighted
    composite. Only removing the axis from the WEIGHTS does. Anyone
    re-implementing the quarantine by passing ``None`` and stopping there would
    ship a silent no-op.
    """
    weights = {"phylo": 2.0, "positive": 1.0, "synteny": 3.0}
    live = {"phylo": 0.8, "synteny": 0.6}

    measured_zero = calculate_fair_rank_score({**live, "positive": 0.0}, weights)
    gated_none = calculate_fair_rank_score({**live, "positive": None}, weights)
    removed = calculate_fair_rank_score(
        live, {k: v for k, v in weights.items() if k != "positive"})

    assert gated_none == pytest.approx(measured_zero), (
        "if this ever stops holding the scorer's formula changed; re-derive "
        "whether the weight-removal step is still required")
    assert removed > gated_none, (
        "removing the quarantined axis from the weights is what actually makes "
        "it absent -- this is the step the quarantine depends on")


def test_None_and_zero_do_diverge_once_the_completeness_floor_binds() -> None:
    """The cancellation above stops at the 0.4 floor: for a sparse candidate,
    ``max(floor, avail/total)`` no longer equals ``avail/total``, so an absent
    axis really is scored differently from a measured zero.

    Recorded so the boundary of the previous test's claim is explicit rather
    than something a future reader has to rediscover.
    """
    weights = {k: 1.0 for k in ("a", "b", "c", "d", "e", "positive")}
    sparse = {"a": 0.9, "b": None, "c": None, "d": None, "e": None}

    measured_zero = calculate_fair_rank_score({**sparse, "positive": 0.0}, weights)
    unavailable = calculate_fair_rank_score({**sparse, "positive": None}, weights)

    assert unavailable > measured_zero


def test_quarantined_axes_leave_the_weight_denominator_not_just_the_numerator() -> None:
    """Removing an axis from SCORING_WEIGHTS must be a removal, not a zeroing.

    With the axis merely zeroed it still sits in ``total_weight``, so
    ``evidence_completeness`` (avail/total) is depressed for EVERY candidate by
    a constant. That cannot reorder anything, so it looks harmless -- and then
    silently changes who clears the absolute CONFIDENCE_MIN_COMPLETENESS gate.

    Would have caught: `SCORING_WEIGHTS['positive'] = 0.0` as the quarantine.
    """
    full = {"phylo": 2.0, "positive": 1.0, "synteny": 3.0}
    scores_full = {"phylo": 1.0, "positive": None, "synteny": 1.0}

    zeroed = calculate_fair_rank_score(
        scores_full, {**full, "positive": 0.0},
        return_diagnostics=True)
    removed = calculate_fair_rank_score(
        {k: v for k, v in scores_full.items() if k != "positive"},
        {k: v for k, v in full.items() if k != "positive"},
        return_diagnostics=True)

    # A candidate perfect on every axis that still exists must read as fully
    # complete once the quarantined axis is genuinely gone.
    assert removed["evidence_completeness_raw"] == pytest.approx(1.0)
    assert zeroed["evidence_completeness_raw"] == pytest.approx(1.0)
    # ...and the source must actually delete the keys rather than zero them.
    src = _rank_source()
    assert "k not in QUARANTINED_AXES" in src, (
        "SCORING_WEIGHTS must have quarantined axes REMOVED, not set to 0")


def test_zero_weight_is_a_different_concept_from_unavailable() -> None:
    """A weight of exactly 0 is an explicit EXCLUSION (PURIFYING_WEIGHT=0),
    which is a policy statement about an axis that HAS data. Unavailable is a
    statement about data that does not exist. Conflating them would make the
    quarantine invisible in provenance.

    Would have caught: implementing the quarantine by setting weights to 0,
    which excluded_signals_from_weights would then report as a weighting
    decision rather than as a data-trust decision.
    """
    excluded = ra.excluded_signals_from_weights(
        {"phylo": 2.0, "purifying": 0.0, "positive": 1.0}, env_var="_NO_SUCH_VAR")
    assert excluded == {"purifying"}


# ---------------------------------------------------------------------------
# 3. Unavailable != zero, in evidence_completeness
# ---------------------------------------------------------------------------

def test_completeness_denominator_shrinks_when_axes_are_quarantined() -> None:
    """A candidate with every LIVE axis present must read as 100% complete.

    Would have caught: leaving the quarantined axes in the denominator, which
    caps completeness at 8/11 = 0.727 -- barely above the 0.7 confidence-view
    floor, so a single further missing axis silently drops a fully-evidenced
    candidate out of the shortlist for a reason that has nothing to do with it.
    """
    ec = _evidence_completeness()
    dark = {**ALL_AXES, "has_dnds_data": False,
            "has_expansion": False, "has_og_data": False}

    assert ec(**dark) == pytest.approx(8 / 11), (
        "without exclusion the quarantined axes still count against everyone")
    assert ec(**dark, excluded_axes=QUARANTINED) == pytest.approx(1.0), (
        "a candidate present on every live axis must be fully complete")


def test_completeness_still_penalizes_a_genuinely_missing_live_axis() -> None:
    """The exclusion must not become a blanket amnesty: an axis that is IN
    SCOPE and missing for this candidate must still cost completeness.

    Would have caught: excluding by "any False flag" rather than by the named
    quarantine set, which would make every sparse candidate look complete.
    """
    ec = _evidence_completeness()
    dark = {**ALL_AXES, "has_dnds_data": False,
            "has_expansion": False, "has_og_data": False}
    missing_synteny = {**dark, "has_synteny": False}

    assert ec(**missing_synteny, excluded_axes=QUARANTINED) == pytest.approx(7 / 8)


def test_naming_both_selection_axes_excludes_the_shared_dnds_source() -> None:
    """purifying and positive share ONE availability fact (did aBSREL report),
    counted once in completeness under the name ``dnds``. The exclusion must
    resolve that aliasing.

    Would have caught: an exclusion set naming 'purifying'/'positive' that
    silently fails to match the 'dnds' source key, leaving dN/dS in the
    denominator while it is out of the score -- a mismatch between the two
    halves of the same decision.
    """
    ec = _evidence_completeness()
    dark = {**ALL_AXES, "has_dnds_data": False,
            "has_expansion": False, "has_og_data": False}
    # 11 sources; naming both selection axes must drop exactly ONE ('dnds'),
    # leaving 10 with 8 present. If the aliasing failed, nothing would be
    # dropped and the answer would be 8/11.
    assert ec(**dark, excluded_axes={"purifying", "positive"}) == pytest.approx(8 / 10)
    assert ec(**dark) == pytest.approx(8 / 11)
    assert ec(**dark, excluded_axes=QUARANTINED) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# 4. The dN/dS axes must actually be gated
# ---------------------------------------------------------------------------

def test_dnds_axes_gate_on_has_dnds_data_in_the_production_scorer() -> None:
    """purifying/positive were UNGATED: a candidate aBSREL never reported on
    contributed a full-weight present 0.0.

    Would have caught: the quarantine emptying dnds_data (so has_dnds_data is
    False everywhere) while the scorer still read the score column -- pinning
    the entire cohort's selection axes at a measured 0 instead of dropping them.
    """
    body = _function_source("calculate_fair_rank_score")
    for axis in ("purifying", "positive"):
        assert re.search(
            rf"'{axis}':\s*\(?row\.get\('{axis}_score_norm'\)\s*"
            rf"if\s+row\.get\('has_dnds_data',\s*True\)\s+else\s+None", body), (
            f"the {axis} axis does not gate on has_dnds_data")


def test_legacy_row_without_the_flag_still_contributes_dnds() -> None:
    """Gating must not break a ranked CSV that predates the flag column.

    Would have caught: `row.get('has_dnds_data', False)`, which would silently
    drop dN/dS from every legacy row rather than treating an absent COLUMN as
    "this file predates the flag" (the established _OPTIONAL_FLAG_SIGNALS rule).
    """
    body = _function_source("calculate_fair_rank_score")
    assert "row.get('has_dnds_data', True)" in body, (
        "an absent has_dnds_data column must default to contributing")


# ---------------------------------------------------------------------------
# 5. Loaders are starved, so the existing has_*_data machinery does the work
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("loader", [
    "load_absrel_with_fdr", "load_busted_signals", "load_cafe_expansion",
    "load_gene_to_orthogroup", "load_og_dnds_reliability",
])
def test_every_orthogroup_derived_loader_is_gated(loader: str) -> None:
    """Each OG-derived loader must be conditioned on ORTHOLOGY_SOURCE_TRUSTED.

    Implementing the quarantine at the LOADERS (rather than as a second set of
    conditionals at each scoring site) is what makes every downstream
    has_*_data flag go False through the existing, already-tested path.

    Would have caught: a loader added later that quietly reintroduces
    Nath-derived data through a path the score-site guards do not cover.
    """
    src = _rank_source()
    assert re.search(
        rf"{loader}\([^)]*\)\s*if\s+ORTHOLOGY_SOURCE_TRUSTED\s+else\s*\{{\}}", src), (
        f"{loader}() is not gated on ORTHOLOGY_SOURCE_TRUSTED")


def test_dnds_reliability_weight_is_nan_not_zero_under_quarantine() -> None:
    """0.0 is a real instruction here ("this omega is unsupported, shrink it to
    the cohort median"); NaN is "not measured". They must not be the same value.

    Would have caught: writing 0.0, which asserts maximal UNRELIABILITY on no
    measurement and instructs reliability_shrink to rewrite every candidate's
    selection score toward the median.
    """
    src = _rank_source()
    assert re.search(r"float\('nan'\)\s*if\s+not\s+ORTHOLOGY_SOURCE_TRUSTED", src), (
        "dnds_reliability_weight must be NaN, not 0.0, when unmeasured")
    assert "if not ORTHOLOGY_SOURCE_TRUSTED:" in src
    assert "skipped the dN/dS reliability shrink" in src, (
        "the shrink step must be skipped, not run against a NaN weight")


# ---------------------------------------------------------------------------
# 6. Rank aggregation
# ---------------------------------------------------------------------------

def test_rankagg_excludes_the_quarantined_signals_explicitly() -> None:
    """They would fall out anyway (empty ranklists), but relying on that leaves
    the quarantine invisible in provenance and reactivates the axis if a stale
    CSV ever supplies a stale flag.

    Would have caught: dropping the union, so a ranked CSV carrying an old
    has_og_confidence_data=True column would let a Nath-derived axis vote.
    """
    src = _rank_source()
    assert "_excluded(SCORING_WEIGHTS) | set(QUARANTINED_AXES)" in src


def test_a_gated_signal_with_all_flags_false_produces_no_ranklist() -> None:
    """Behavioural check of the mechanism the quarantine leans on."""
    import pandas as pd
    df = pd.DataFrame({
        "id": ["a", "b", "c"],
        "og_confidence_score_norm": [0.9, 0.5, 0.1],
        "has_og_confidence_data": [False, False, False],
        "phylo_score_norm": [0.9, 0.5, 0.1],
        "has_phylo_data": [True, True, True],
    })
    lists = ra.build_ranklists_from_df(df)
    assert "og_confidence" not in lists, (
        "a signal with no available data must not vote at all")
    assert "phylo" in lists


# ---------------------------------------------------------------------------
# 7. Shell integration points
# ---------------------------------------------------------------------------

def test_stage07_skips_og_coverage_columns_under_quarantine() -> None:
    """og_dnds_reliability=='low' is a discovery-view MEMBERSHIP disjunct, so
    running the augmenter without a trusted OG mapping labels every candidate
    'low' and sweeps the whole cohort into the discovery shortlist.

    Would have caught: leaving the augmenter running and treating its output as
    merely cosmetic. The column must be ABSENT so emit_ranked_views takes its
    documented "disjunct skipped" path.
    """
    text = STAGE07.read_text()
    assert re.search(r'if\s+\[\s+"\$\{ORTHOLOGY_SOURCE_TRUSTED:-0\}"\s+!=\s+"1"\s+\];\s*then\s*\n\s*COV_OG_TSV=""',
                     text), "stage 07 still builds OG-coverage columns under quarantine"


def test_stage06c_disables_the_og_vote_and_states_the_consequence() -> None:
    """The OG vote is entirely Nath-fed. Disabling it is correct -- but
    classify_consensus reaches 'non-chemoreceptor' only on 3-of-3 and
    'likely-non-chemoreceptor' only when HMM *and* OG agree, so with OG dark
    NEITHER demotion tier is reachable and every candidate stays
    'chemoreceptor-candidate'.

    Would have caught: silently disabling the vote without recording that the
    non-chemoreceptor filter now provides zero suppression -- a reviewer would
    read an all-'chemoreceptor-candidate' column as a finding rather than as a
    disabled classifier.
    """
    text = STAGE06C.read_text()
    assert 'ORTHOLOGY_SOURCE_TRUSTED:-0' in text
    assert "chemoreceptor-candidate" in text and "demotion tier" in text, (
        "06c does not state that neither demotion tier remains reachable")


def test_consensus_really_cannot_demote_without_the_og_source() -> None:
    """Verify the consequence above against the real consensus function rather
    than trusting the comment that documents it."""
    from classify_consensus import _consensus
    # HMM + placement agree perfectly; OG is dark.
    family, label, n = _consensus("bioamine", "", "bioamine")
    assert label == "chemoreceptor-candidate", (
        "HMM+placement must NOT reach a demotion tier -- if this ever changes, "
        "the quarantine's stated cost in 06c is wrong and must be updated")
    assert n == 1


# ---------------------------------------------------------------------------
# 8. Nothing was deleted
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("symbol", [
    "load_og_trees", "get_og_confidence_score", "load_cafe_expansion",
    "get_expansion_score", "load_og_dnds_reliability", "load_absrel_with_fdr",
])
def test_quarantined_code_paths_still_exist(symbol: str) -> None:
    """This is a quarantine pending replacement data, not a removal. The code
    must be intact so lifting the switch restores the axes.

    Would have caught: 'cleaning up' the dormant paths, which turns a one-line
    reversal into a re-implementation.
    """
    assert re.search(rf"def\s+{symbol}\s*\(", _rank_source()), (
        f"{symbol}() was removed; the quarantine must not delete code paths")


def test_config_and_stage_scripts_are_syntactically_valid() -> None:
    """bash -n on every shell file this change touched."""
    for path in (CONFIG, STAGE07, STAGE06C):
        proc = subprocess.run(["bash", "-n", str(path)],
                              capture_output=True, text=True)
        assert proc.returncode == 0, f"{path.name}: {proc.stderr}"
