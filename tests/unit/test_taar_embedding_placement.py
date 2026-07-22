"""Tests for scripts/analysis/taar_embedding_placement.py (bead en2r / h0y0).

The module answers, from embedding evidence alone, whether the 36 vertebrate
TAAR2-9 should be ADDED to the reference anchors (as their own prototype).

These tests are SYNTHETIC and LOCAL: no GPU, no Unity, no real embeddings.
The two planted cases are the must-accept / must-reject guards for the decision
rule:

  * PLANTED-ADD  — TAARs drawn FROM the aminergic Gaussian. Q1 must find them
    statistically indistinguishable from aminergic, Q2 must find no suppression
    of candidate novelty beyond the random-36 null -> verdict ADD.
  * PLANTED-HOLD — TAARs sit ON the candidate cloud, far from aminergic. Q1 must
    reject membership, Q2 must show suppression AND candidates' nearest anchor
    flips to TAAR -> verdict HOLD.

Plus fail-loud guards (missing prefix / empty set / single-member family) and
unit checks of the exposed statistics (BH-FDR, two-sample permutation tests,
family mapping from the anchor TSV).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

# conftest puts scripts/ on sys.path; the module lives under scripts/analysis/.
SCRIPTS = Path(__file__).resolve().parent.parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from analysis.taar_embedding_placement import (  # noqa: E402
    Config,
    bh_fdr,
    map_anchor_families,
    run_and_write,
    run_consensus,
    run_per_model,
    split_by_prefix,
    two_sample_tests,
    validate_sets,
)

DIM = 12
FAST = Config(n_perm=200, n_random=80, knn=7, k=3, seed=0, split_frac=0.5,
              min_aminergic=4, topk=(5, 10))


def _cluster(rng, n, center, sd=1.0):
    return rng.randn(n, DIM) * sd + np.asarray(center, dtype=float)


def _center(axis, mag=6.0):
    c = np.zeros(DIM)
    c[axis] = mag
    return c


def _planted(kind, seed=0):
    """Build (embeddings, anchor_families) for a planted scenario.

    Anchors: aminergic / peptide / opsin Gaussian clusters on distinct axes.
    Probes: their own cluster. Candidates + TAARs placed per `kind`.
    """
    rng = np.random.RandomState(seed)
    emb, fam = {}, {}

    def add_anchors(name, axis, n):
        block = _cluster(rng, n, _center(axis))
        for i, v in enumerate(block):
            aid = f"ANCHOR_A_10_{name[:3].upper()}{i:03d}"
            emb[aid] = v
            fam[aid] = name

    # Aminergic needs enough calibration points for split-conformal to reach a
    # p below alpha (min achievable p = 1/(n_cal+1)); 60 -> ~30 calibration.
    add_anchors("aminergic", 0, 60)
    add_anchors("peptide", 1, 25)
    add_anchors("opsin", 2, 25)

    probes = _cluster(rng, 30, _center(3))
    for i, v in enumerate(probes):
        emb[f"PROBE_A_P{i:03d}"] = v

    novel_center = _center(6, 12.0)  # far region no anchor family covers
    if kind == "add":
        # TAARs live INSIDE the aminergic Gaussian; candidates sit on the peptide
        # family (+ a small novelty bump) so their nearest prototype is peptide,
        # never the (aminergic-co-located) TAAR -> Q1 indistinguishable, no flips.
        cand_center = _center(1) + _center(6, 4.0)
        cand = _cluster(rng, 40, cand_center)
        taar = _cluster(rng, 24, _center(0), sd=0.8)  # clearly inside aminergic
    elif kind == "hold":
        # TAARs sit ON the candidate cloud, far from aminergic.
        cand = _cluster(rng, 40, novel_center)
        taar = _cluster(rng, 24, novel_center)
    else:
        raise ValueError(kind)

    for i, v in enumerate(cand):
        emb[f"BersteEVm{i:06d}t1"] = v
    for i, v in enumerate(taar):
        emb[f"TAAR_ACC{i:03d}"] = v
    return emb, fam


# --------------------------------------------------------------------------- #
# split + coverage + fail-loud                                                 #
# --------------------------------------------------------------------------- #
def test_taar_embedding_split_by_prefix_partitions_the_four_sets():
    emb, _ = _planted("add")
    sets = split_by_prefix(emb)
    assert set(sets) == {"anchor", "candidate", "probe", "taar"}
    assert len(sets["anchor"]) == 110
    assert len(sets["candidate"]) == 40
    assert len(sets["probe"]) == 30
    assert len(sets["taar"]) == 24
    # nothing dropped, nothing double-counted
    assert sum(len(s) for s in sets.values()) == len(emb)


def test_taar_embedding_split_rejects_an_unknown_prefix():
    emb, _ = _planted("add")
    emb["WEIRD_1"] = np.zeros(DIM)
    with pytest.raises(ValueError, match="prefix"):
        split_by_prefix(emb)


def test_taar_embedding_failloud_on_empty_candidate_set():
    emb, fam = _planted("add")
    for k in [k for k in emb if k.startswith("BersteEVm")]:
        del emb[k]
    with pytest.raises(ValueError, match="empty|candidate"):
        run_per_model(emb, fam, "m", FAST)


def test_taar_embedding_failloud_on_missing_taar_prefix():
    emb, fam = _planted("add")
    for k in [k for k in emb if k.startswith("TAAR_")]:
        del emb[k]
    with pytest.raises(ValueError, match="empty|taar|TAAR"):
        run_per_model(emb, fam, "m", FAST)


def test_taar_embedding_failloud_on_single_member_aminergic_family():
    emb, fam = _planted("add")
    amin = [k for k in fam if fam[k] == "aminergic"]
    for k in amin[1:]:                    # leave exactly one aminergic anchor
        fam[k] = "peptide"
    with pytest.raises(ValueError, match="aminergic|too few|member"):
        run_per_model(emb, fam, "m", FAST)


# --------------------------------------------------------------------------- #
# exposed statistics                                                           #
# --------------------------------------------------------------------------- #
def test_taar_embedding_bh_fdr_matches_known_values():
    p = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    adj = bh_fdr(p)
    # classic BH: smallest stays *5, all others collapse to 0.05 here
    assert np.isclose(adj[0], 0.05)
    assert np.all(adj <= 1.0)
    # monotone non-decreasing in the original p-order
    assert np.all(np.diff(adj) >= -1e-12)


def test_taar_embedding_two_sample_separates_and_merges():
    rng = np.random.RandomState(3)
    A = rng.randn(40, DIM)
    B_same = rng.randn(40, DIM)
    B_far = rng.randn(40, DIM) + 12.0
    same = two_sample_tests(A, B_same, n_perm=300, seed=1)
    far = two_sample_tests(A, B_far, n_perm=300, seed=1)
    # identical distributions -> not significant; separated -> significant
    assert same["mmd_p"] > 0.05 and same["energy_p"] > 0.05
    assert far["mmd_p"] < 0.05 and far["energy_p"] < 0.05
    assert far["mmd_stat"] > same["mmd_stat"]


def test_taar_embedding_map_anchor_families_from_tsv(tmp_path):
    tsv = tmp_path / "anchors.tsv"
    tsv.write_text(
        "accession\ttier\ttaxid\tspecies\tfamily\tclass\tevidence\n"
        "Q1\t10\t9606\tHomo\taminergic\tA\treviewed\n"
        "Q2\t5\t9606\tHomo\tpeptide\tA\treviewed\n"
    )
    ids = ["ANCHOR_A_10_Q1", "ANCHOR_A_5_Q2"]
    fam = map_anchor_families(ids, str(tsv))
    assert fam["ANCHOR_A_10_Q1"] == "aminergic"
    assert fam["ANCHOR_A_5_Q2"] == "peptide"


def test_taar_embedding_map_anchor_families_failloud_on_unmapped(tmp_path):
    tsv = tmp_path / "anchors.tsv"
    tsv.write_text("accession\tfamily\nQ1\taminergic\n")
    with pytest.raises(ValueError, match="Q9|unmapped|overlap"):
        map_anchor_families(["ANCHOR_A_10_Q1", "ANCHOR_A_10_Q9"], str(tsv))


# --------------------------------------------------------------------------- #
# planted ADD (must-accept)                                                    #
# --------------------------------------------------------------------------- #
def test_taar_embedding_planted_add_verdict():
    emb, fam = _planted("add")
    res = run_per_model(emb, fam, "planted_add", FAST)
    q1, q2 = res["q1"], res["q2"]
    # Q1: TAARs indistinguishable from aminergic
    assert q1["indistinguishable"] is True
    assert q1["n_rejected"] == 0
    assert q1["mmd_taar_aminergic"]["p"] > 0.05
    # Q2: no suppression beyond the random-36 null; candidates do not localise to TAAR
    assert q2["no_excess_suppression"] is True
    assert q2["nearest_anchor_flips"]["count"] == 0
    assert res["verdict"] == "ADD"


def test_taar_embedding_planted_add_reports_coverage():
    emb, fam = _planted("add")
    res = run_per_model(emb, fam, "planted_add", FAST)
    cov = res["coverage"]
    assert cov["anchor"] == 110 and cov["candidate"] == 40
    assert cov["probe"] == 30 and cov["taar"] == 24
    assert cov["anchor_families"]["aminergic"] == 60


# --------------------------------------------------------------------------- #
# planted HOLD (must-reject)                                                   #
# --------------------------------------------------------------------------- #
def test_taar_embedding_planted_hold_verdict():
    emb, fam = _planted("hold")
    res = run_per_model(emb, fam, "planted_hold", FAST)
    q1, q2 = res["q1"], res["q2"]
    # Q1: TAARs are distinguishable from aminergic
    assert q1["indistinguishable"] is False
    assert q1["n_rejected"] > 0
    # Q1 support: TAAR-vs-aminergic separates
    assert q1["mmd_taar_aminergic"]["p"] < 0.05
    # Q2: suppression beyond null AND candidates localise onto the TAAR prototype
    assert q2["no_excess_suppression"] is False
    assert q2["nearest_anchor_flips"]["count"] > 0
    # adding TAAR lowers candidate novelty (paired deltas negative on median)
    deltas = np.array([r["delta"] for r in q2["novelty"]])
    assert np.median(deltas) < 0
    assert res["verdict"] == "HOLD"


# --------------------------------------------------------------------------- #
# robustness                                                                   #
# --------------------------------------------------------------------------- #
def test_taar_embedding_robustness_reports_three_methods():
    emb, fam = _planted("hold")
    res = run_per_model(emb, fam, "planted_hold", FAST)
    rob = res["robustness"]
    assert set(rob["methods"]) == {"mahalanobis", "knn", "lof"}
    for m in rob["methods"].values():
        assert "suppression_sign" in m and "kendall_tau" in m
    assert isinstance(rob["stable"], bool)
    assert "relative_maha" in rob


# --------------------------------------------------------------------------- #
# consensus                                                                    #
# --------------------------------------------------------------------------- #
def test_taar_embedding_consensus_agrees_when_both_add():
    emb, fam = _planted("add")
    r1 = run_per_model(emb, fam, "proteinclip3b", FAST)
    r2 = run_per_model(emb, fam, "protrek", FAST)
    con = run_consensus([r1, r2])
    assert con["agreement"] is True
    assert con["consensus_verdict"] == "ADD"


def test_taar_embedding_consensus_flags_disagreement_as_a_finding():
    add_emb, add_fam = _planted("add")
    hold_emb, hold_fam = _planted("hold")
    r1 = run_per_model(add_emb, add_fam, "proteinclip3b", FAST)
    r2 = run_per_model(hold_emb, hold_fam, "protrek", FAST)
    con = run_consensus([r1, r2])
    assert con["agreement"] is False
    assert con["consensus_verdict"] == "DISAGREEMENT"
    assert "proteinclip3b" in con["per_model"] and "protrek" in con["per_model"]


# --------------------------------------------------------------------------- #
# end-to-end file writing (exercises load_embeddings + TSV join + outputs)     #
# --------------------------------------------------------------------------- #
def test_taar_embedding_run_and_write_produces_all_outputs(tmp_path):
    emb, fam = _planted("add")
    # emit an npz keyed by id and a matching anchor TSV (accession join).
    npz = tmp_path / "emb.npz"
    np.savez(npz, **emb)
    rows = ["accession\tfamily"]
    seen = set()
    for aid, f in fam.items():
        acc = aid.split("_", 3)[3]
        if acc not in seen:
            rows.append(f"{acc}\t{f}")
            seen.add(acc)
    tsv = tmp_path / "anchors.tsv"
    tsv.write_text("\n".join(rows) + "\n")

    out = tmp_path / "out"
    res = run_and_write(str(npz), str(tsv), str(out), "proteinclip3b", FAST)
    for name in [
        "q1_conformal.tsv", "q1_twosample.tsv", "q1_probe_posterior.tsv",
        "q1_knn_purity.tsv", "q2_novelty_deltas.tsv", "q2_rank_agreement.tsv",
        "q2_nearest_anchor.tsv", "q2_random36_null.tsv", "robustness.tsv",
        "summary.json", "VERDICT.md",
    ]:
        assert (out / name).exists(), f"missing output {name}"
    summary = json.loads((out / "summary.json").read_text())
    assert summary["verdict"] in {"ADD", "HOLD"}
    assert res["verdict"] == summary["verdict"]
    assert res["verdict"] in (out / "VERDICT.md").read_text()
