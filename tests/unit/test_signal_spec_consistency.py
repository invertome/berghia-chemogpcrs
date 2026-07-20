"""Drift guard: the 13 CORE ranking-signal keys are declared independently
in three modules, each in its own representation:

- ``audit_signal_ranking_independence.SIGNAL_COLUMNS`` -- a flat list of
  ranked-CSV column names: the 13 core signals as ``"<key>_score"`` (the
  rank_candidates.py CSV convention), plus 5 Task-6 evidence-channel
  signals (struct_novelty, struct_nonchemo_corrob, emb_classA_sim,
  emb_nonchemo_sim, or_microswitch) that do NOT follow the ``_score``
  suffix. Filtering on the suffix isolates exactly the core 13 by name,
  with no need to hard-code the 5 channel names in this test too.
- ``rank_aggregation.SIGNAL_SPEC`` -- a list of per-signal specs; the 13
  core signals are plain ``(key, has_flag_or_None)`` 2-tuples, the 5 Task-6
  channel signals are ``(key, has_flag, column, invert)`` 4-tuples.
  Filtering on tuple length isolates the core 13.
- ``ablate_ranking._SIGNAL_SPEC`` -- a list of 13
  ``(weight_key, score_norm_column, has_flag_or_None, env_var, default)``
  5-tuples. ALL of its entries are core signals -- ablate_ranking has no
  Task-6 evidence-channel entries at all (see its own module docstring:
  "signal axes of the production fair scorer").

This test imports the live constants (not copies), so an edit to any one
module's core-13 set -- e.g. adding/renaming/dropping a signal in
rank_aggregation.py's SIGNAL_SPEC without mirroring it into
audit_signal_ranking_independence.py's SIGNAL_COLUMNS or
ablate_ranking.py's _SIGNAL_SPEC -- fails this test immediately instead of
silently producing a rank-aggregation / audit / ablation harness that
disagree about which signals exist.
"""
import ablate_ranking as ablate
import audit_signal_ranking_independence as audit
import rank_aggregation as ra

# Was 12; bead hf3u added `lse_nesting_depth` as a 13th core signal --
# topological nesting depth (root-to-tip node count) alongside the existing
# `lse_divergence` (cumulative branch length). They are deliberately NOT fused: on
# this repo's real trees they correlate +0.897/+0.604 population-wide but only
# +0.078/+0.034 within the scored top quartile, and there is no positive
# control that could arbitrate them, so both vote and the independence audit
# decides. Bumping this number is the intended way to add a core signal -- it
# is what forces all three modules to be updated together.
CORE_SIGNAL_COUNT = 13


def _core_keys_from_audit():
    """SIGNAL_COLUMNS entries ending in "_score" are the 13 core signals;
    the 5 Task-6 evidence-channel signals in the same list do not follow
    that suffix (see module docstring)."""
    return {c[: -len("_score")] for c in audit.SIGNAL_COLUMNS if c.endswith("_score")}


def _core_keys_from_rank_aggregation():
    """SIGNAL_SPEC's core-13 entries are plain (key, flag) 2-tuples; the 5
    Task-6 channel entries are 4-tuples -- tuple length isolates the
    core-13 without a second hard-coded name list."""
    return {spec[0] for spec in ra.SIGNAL_SPEC if len(spec) == 2}


def _core_keys_from_ablate():
    """_SIGNAL_SPEC contains ONLY the 13 core signals."""
    return {spec[0] for spec in ablate._SIGNAL_SPEC}


def test_core_signal_key_sets_are_exactly_twelve():
    # Guards the filters above themselves: if a future edit changed what
    # counts as "core" (e.g. a 6th evidence channel without a "_score"
    # suffix, or a channel entry that happens to be a 2-tuple), the count
    # would drift away from CORE_SIGNAL_COUNT and this fails before the cross-module
    # comparison below could give a misleading pass.
    assert len(_core_keys_from_audit()) == CORE_SIGNAL_COUNT
    assert len(_core_keys_from_rank_aggregation()) == CORE_SIGNAL_COUNT
    assert len(_core_keys_from_ablate()) == CORE_SIGNAL_COUNT


def test_core_signal_keys_match_across_audit_rankagg_and_ablate():
    from_audit = _core_keys_from_audit()
    from_rankagg = _core_keys_from_rank_aggregation()
    from_ablate = _core_keys_from_ablate()

    assert from_audit == from_rankagg == from_ablate
