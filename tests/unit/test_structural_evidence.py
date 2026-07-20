"""Unit tests for scripts/structural_evidence.py.

Verifies the Foldseek structural-evidence channel (Task 4 of the ML/PLM
chemoreceptor ranking plan, docs/plans/2026-07-01-ml-plm-chemoreceptor-
ranking.md). Council rule: cross-species structural resemblance is an
honest RECALL (novelty) + EXCLUSION (non-chemoreceptor corroboration)
signal only — never a positive "looks like a known chemoreceptor" score.
There is no such state to test for, by design.
"""
from pathlib import Path

import pytest

from structural_evidence import (
    NON_CHEMORECEPTOR_FAMILIES,
    classify_hit,
    parse_foldseek,
    structural_channel,
)


def _write(tmp_path: Path, lines: list) -> Path:
    p = tmp_path / "foldseek_hits.tsv"
    p.write_text("\n".join(lines) + "\n")
    return p


# --- parse_foldseek ---------------------------------------------------------

def test_parse_keeps_higher_alntmscore_hit_per_query(tmp_path):
    """Two hits for the same query -> the one with higher alntmscore wins,
    even though it is listed second and has a lower fident."""
    path = _write(tmp_path, [
        "cand1\ttargetA\t0.30\t0.40\t1e-5",
        "cand1\ttargetB\t0.55\t0.72\t1e-9",
    ])
    result = parse_foldseek(str(path))
    assert set(result) == {"cand1"}
    assert result["cand1"]["target"] == "targetB"
    assert result["cand1"]["alntmscore"] == pytest.approx(0.72)
    assert result["cand1"]["fident"] == pytest.approx(0.55)
    assert result["cand1"]["evalue"] == pytest.approx(1e-9)


def test_parse_two_queries_get_independent_best_hits(tmp_path):
    path = _write(tmp_path, [
        "cand1\ttargetA\t0.30\t0.40\t1e-5",
        "cand2\ttargetC\t0.90\t0.95\t1e-20",
    ])
    result = parse_foldseek(str(path))
    assert set(result) == {"cand1", "cand2"}
    assert result["cand2"]["target"] == "targetC"


def test_parse_skips_blank_lines(tmp_path):
    path = _write(tmp_path, [
        "cand1\ttargetA\t0.30\t0.40\t1e-5",
        "",
        "   ",
        "cand2\ttargetC\t0.90\t0.95\t1e-20",
    ])
    result = parse_foldseek(str(path))
    assert set(result) == {"cand1", "cand2"}


def test_parse_skips_short_and_garbage_lines(tmp_path):
    path = _write(tmp_path, [
        "cand1\ttargetA\t0.30\t0.40\t1e-5",
        "not\tenough\tcolumns",
        "cand2\ttargetB\tNaN\toops\tnotanumber",
        "cand3\ttargetC\t0.90\t0.95\t1e-20",
    ])
    result = parse_foldseek(str(path))
    assert set(result) == {"cand1", "cand3"}


def test_parse_missing_file_returns_empty_dict(tmp_path):
    missing = tmp_path / "does_not_exist.tsv"
    assert parse_foldseek(str(missing)) == {}


# --- classify_hit ------------------------------------------------------------

def test_classify_hit_none_is_novel():
    assert classify_hit(None, family_map={}) == "novel"


def test_classify_hit_below_threshold_is_novel():
    best = {"target": "t1", "alntmscore": 0.49, "fident": 0.5, "evalue": 1e-5}
    assert classify_hit(best, family_map={"t1": "class-A-orphan"},
                         tm_threshold=0.5) == "novel"


def test_classify_hit_at_threshold_boundary_counts_as_hit():
    """alntmscore == tm_threshold is a confident hit (>= threshold wins),
    not novel."""
    best = {"target": "t1", "alntmscore": 0.5, "fident": 0.5, "evalue": 1e-5}
    result = classify_hit(best, family_map={"t1": "aminergic_5HT"},
                           tm_threshold=0.5)
    assert result != "novel"
    assert result == "known_non_chemoreceptor"


def test_classify_hit_non_chemoreceptor_family_corroborates_exclusion():
    best = {"target": "t1", "alntmscore": 0.9, "fident": 0.8, "evalue": 1e-30}
    assert classify_hit(best, family_map={"t1": "class-B-secretin"}) == \
        "known_non_chemoreceptor"


def test_classify_hit_known_other_when_family_not_nonchemoreceptor():
    best = {"target": "t1", "alntmscore": 0.9, "fident": 0.8, "evalue": 1e-30}
    assert classify_hit(
        best, family_map={"t1": "class-A-chemoreceptor-like"}) == "known_other"


def test_classify_hit_known_other_when_target_absent_from_family_map():
    best = {"target": "unmapped_target", "alntmscore": 0.9, "fident": 0.8,
            "evalue": 1e-30}
    assert classify_hit(best, family_map={}) == "known_other"


@pytest.mark.parametrize("label_form", ["{family}", "{family}_subfamilyX"])
@pytest.mark.parametrize("family", sorted(NON_CHEMORECEPTOR_FAMILIES))
def test_classify_hit_recognizes_every_nonchemoreceptor_tag(family, label_form):
    """NON_CHEMORECEPTOR_FAMILIES is a module constant so it is testable and
    overridable (per the spec) -- prove every tag in it actually triggers
    the exclusion-corroboration state.

    Both label forms the vocabulary actually uses are exercised: the bare
    coarse family (what references/anchors/anchor_set.tsv's `family` column
    really holds -- family and subfamily are separate columns there) and the
    documented "<coarse>_<subfamily>" composite. The previous
    "{family}-subfamily-x" form was invented -- nothing emits it -- and it
    only passed because the match was an unanchored substring test, the very
    defect that would have let "Class A (Rhodopsin)" match "opsin".
    """
    best = {"target": "t1", "alntmscore": 1.0, "fident": 1.0, "evalue": 0.0}
    family_map = {"t1": label_form.format(family=family)}
    assert classify_hit(best, family_map) == "known_non_chemoreceptor"


def test_nonchemoreceptor_families_matches_authoritative_taxonomy():
    """The exclusion set must equal the classifier's canonical coarse families
    (COARSE_FAMILIES in scripts/validate_classification_hmms.py). Guards the
    completeness gap the fix closed: no missing live families, no stray tags."""
    assert NON_CHEMORECEPTOR_FAMILIES == {
        "aminergic", "peptide", "opsin", "lipid", "nucleotide",
        "metabotropic-neurotransmitter", "glycoprotein-hormone",
        "class-B-secretin", "class-C", "class-F-frizzled",
    }
    # The vestigial colloquialism must be gone (taxonomy emits 'aminergic').
    assert "bioamine" not in NON_CHEMORECEPTOR_FAMILIES


def test_classify_hit_glycoprotein_hormone_is_nonchemoreceptor_regression():
    """Regression for the exact completeness gap: a confident structural hit to
    the glycoprotein-hormone family (TSHR/FSHR/LHCGR) MUST corroborate
    exclusion. Before the taxonomy fix, NON_CHEMORECEPTOR_FAMILIES omitted this
    LIVE family, so the hit was silently misclassified 'known_other' and the
    exclusion signal was lost. Tests both the bare coarse label and the
    realistic '<coarse>_<subfamily>' form."""
    best = {"target": "tshr_model", "alntmscore": 0.88, "fident": 0.70,
            "evalue": 1e-28}
    assert classify_hit(best, {"tshr_model": "glycoprotein-hormone"}) == \
        "known_non_chemoreceptor"
    assert classify_hit(best, {"tshr_model": "glycoprotein-hormone_TSH"}) == \
        "known_non_chemoreceptor"


def test_classify_hit_metabotropic_neurotransmitter_is_nonchemoreceptor_regression():
    """Companion gap: metabotropic-neurotransmitter (mGlu / GABA-B) was also
    absent from the pre-fix set and would have leaked to 'known_other'."""
    best = {"target": "mglur_model", "alntmscore": 0.82, "fident": 0.60,
            "evalue": 1e-20}
    assert classify_hit(
        best, {"mglur_model": "metabotropic-neurotransmitter"}) == \
        "known_non_chemoreceptor"


# --- structural_channel -------------------------------------------------------

def test_structural_channel_flags_novel_and_corroborated(tmp_path):
    path = _write(tmp_path, [
        "cand_novel\ttargetX\t0.10\t0.20\t1e-2",
        "cand_nonchemo\ttargetY\t0.85\t0.90\t1e-25",
        "cand_other\ttargetZ\t0.85\t0.90\t1e-25",
    ])
    family_map = {
        "targetY": "opsin",
        "targetZ": "class-A-chemoreceptor-like",
    }
    channel = structural_channel(str(path), family_map)

    assert channel["cand_novel"] == {
        "struct_state": "novel", "struct_novelty": 1,
        "struct_nonchemo_corrob": 0, "has_struct_data": True,
    }
    assert channel["cand_nonchemo"] == {
        "struct_state": "known_non_chemoreceptor", "struct_novelty": 0,
        "struct_nonchemo_corrob": 1, "has_struct_data": True,
    }
    assert channel["cand_other"] == {
        "struct_state": "known_other", "struct_novelty": 0,
        "struct_nonchemo_corrob": 0, "has_struct_data": True,
    }


def test_structural_channel_only_emits_queries_present_in_foldseek_output(tmp_path):
    """A candidate with no Foldseek row at all is simply absent from the
    returned dict (has_struct_data=False is a downstream-integration
    contract, not fabricated here)."""
    path = _write(tmp_path, [
        "cand_present\ttargetX\t0.90\t0.95\t1e-20",
    ])
    channel = structural_channel(str(path), family_map={})
    assert set(channel) == {"cand_present"}
    assert "cand_absent" not in channel


def test_structural_channel_missing_file_returns_empty_dict(tmp_path):
    missing = tmp_path / "nope.tsv"
    assert structural_channel(str(missing), family_map={}) == {}


def test_structural_channel_respects_custom_tm_threshold(tmp_path):
    path = _write(tmp_path, [
        "cand1\ttargetX\t0.60\t0.60\t1e-10",
    ])
    lenient = structural_channel(str(path), family_map={}, tm_threshold=0.5)
    strict = structural_channel(str(path), family_map={}, tm_threshold=0.9)
    assert lenient["cand1"]["struct_state"] != "novel"
    assert strict["cand1"]["struct_state"] == "novel"
