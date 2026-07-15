"""Tests for the Mahalanobis (`maha`) producer path in build_embedding_channel.py.

Phase-0 wire-up (bead cw3): `build_embedding_channel_maha` glues the tested
`embedding_evidence.mahalanobis_channel` scorer to real .npz + anchor_set.tsv
inputs and emits a TSV with the novelty column (`emb_novelty`) that
rank_aggregation consumes as a positive ranking axis. The cosine path is
untouched.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from build_embedding_channel import (
    MAHA_TSV_COLUMNS,
    build_embedding_channel_maha,
    main,
    novelty_reference_labels,
    write_channel_tsv,
)


def _write_reference(tmp_path):
    """A 2-family reference (aminergic +x, peptide -x), keyed by the composite
    ANCHOR_<class>_<tier>_<accession> ids load_ref_labels reconstructs."""
    rng = np.random.RandomState(0)
    ref, rows = {}, []
    for i in range(24):
        acc = f"P{i:03d}"
        family = "aminergic" if i < 12 else "peptide"
        sign = 1.0 if i < 12 else -1.0
        ref[f"ANCHOR_A_1_{acc}"] = rng.randn(8) + np.eye(8)[0] * 5.0 * sign
        rows.append({"accession": acc, "tier": 1, "class": "A", "family": family,
                     "taxid": 1, "species": "x", "evidence": "reviewed"})
    ref_npz = tmp_path / "ref.npz"
    np.savez(ref_npz, **ref)
    labels_tsv = tmp_path / "labels.tsv"
    pd.DataFrame(rows).to_csv(labels_tsv, sep="\t", index=False)
    return str(ref_npz), str(labels_tsv)


def test_build_embedding_channel_maha_emits_novelty_and_family(tmp_path):
    ref_npz, labels_tsv = _write_reference(tmp_path)
    cand = {"like_A": np.eye(8)[0] * 5.0, "novel": np.eye(8)[1] * 50.0}
    cand_npz = tmp_path / "cand.npz"
    np.savez(cand_npz, **cand)

    channel = build_embedding_channel_maha(str(cand_npz), ref_npz, labels_tsv, k=3)

    assert set(channel) == {"like_A", "novel"}
    for entry in channel.values():
        assert "emb_novelty" in entry
        assert "emb_nonchemo_family" in entry
        assert entry["has_emb_data"] is True
        assert entry["emb_leakage_flag"] is True
        assert not any("chemoreceptor" in key.lower() for key in entry)
    # the divergent candidate is ranked more novel (the positive axis)
    assert channel["novel"]["emb_novelty"] > channel["like_A"]["emb_novelty"]


def test_write_channel_tsv_maha_columns(tmp_path):
    ref_npz, labels_tsv = _write_reference(tmp_path)
    cand_npz = tmp_path / "cand.npz"
    np.savez(cand_npz, c=np.eye(8)[0])
    channel = build_embedding_channel_maha(str(cand_npz), ref_npz, labels_tsv)
    out = tmp_path / "chan.tsv"
    write_channel_tsv(channel, str(out), columns=MAHA_TSV_COLUMNS)
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == MAHA_TSV_COLUMNS
    assert "emb_novelty" in df.columns


def test_main_scorer_maha_writes_maha_columns(tmp_path):
    ref_npz, labels_tsv = _write_reference(tmp_path)
    cand_npz = tmp_path / "cand.npz"
    np.savez(cand_npz, c=np.eye(8)[0] * 5.0)
    out = tmp_path / "chan.tsv"
    main(["--scorer", "maha", "--candidate-npz", str(cand_npz),
          "--ref-npz", ref_npz, "--ref-labels", labels_tsv, "--out", str(out)])
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == MAHA_TSV_COLUMNS


def test_main_defaults_to_cosine_columns(tmp_path):
    from build_embedding_channel import TSV_COLUMNS
    ref_npz, labels_tsv = _write_reference(tmp_path)
    cand_npz = tmp_path / "cand.npz"
    np.savez(cand_npz, c=np.eye(8)[0] * 5.0)
    out = tmp_path / "chan.tsv"
    main(["--candidate-npz", str(cand_npz), "--ref-npz", ref_npz,
          "--ref-labels", labels_tsv, "--out", str(out)])
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == TSV_COLUMNS   # unchanged default


def _write_mixed_class_reference(tmp_path):
    """Reference spanning the three cases the novelty restriction must tell apart:
    CHARACTERIZED class-A families (aminergic, peptide), an OUT-OF-CLASS family
    (class-C) and an UNCHARACTERIZED class-A family (orphan). Each family is a
    tight cluster around its own axis so distances are unambiguous."""
    rng = np.random.RandomState(0)
    ref, rows = {}, []

    def add(n, prefix, family, cls, center):
        for i in range(n):
            acc = f"{prefix}{i:03d}"
            ref[f"ANCHOR_{cls}_1_{acc}"] = rng.randn(8) * 0.1 + center
            rows.append({"accession": acc, "tier": 1, "class": cls, "family": family,
                         "taxid": 1, "species": "x", "evidence": "reviewed"})

    add(12, "A", "aminergic", "A", np.eye(8)[0] * 5.0)
    add(12, "P", "peptide", "A", -np.eye(8)[0] * 5.0)
    add(12, "C", "class-C", "C", np.eye(8)[1] * 50.0)      # out-of-class: gate only
    add(12, "O", "orphan", "A", np.eye(8)[2] * 50.0)       # class-A but uncharacterized
    ref_npz = tmp_path / "mixed_ref.npz"
    np.savez(ref_npz, **ref)
    labels_tsv = tmp_path / "mixed_labels.tsv"
    pd.DataFrame(rows).to_csv(labels_tsv, sep="\t", index=False)
    return str(ref_npz), str(labels_tsv)


def _maha(tmp_path, cand):
    ref_npz, labels_tsv = _write_mixed_class_reference(tmp_path)
    cand_npz = tmp_path / "mixed_cand.npz"
    np.savez(cand_npz, **cand)
    return build_embedding_channel_maha(str(cand_npz), ref_npz, labels_tsv)


def test_candidate_near_out_of_class_centroid_still_scores_novel(tmp_path):
    """A divergent class-A candidate that happens to land near a class-C centroid
    must stay HIGH novelty: out-of-class (B/C/F) references serve the upstream
    class-A gate, they must never act as novelty prototypes that explain a
    candidate away."""
    channel = _maha(tmp_path, {"on_aminergic": np.eye(8)[0] * 5.0,
                               "on_classC": np.eye(8)[1] * 50.0})
    assert channel["on_classC"]["emb_novelty"] > channel["on_aminergic"]["emb_novelty"]


def test_candidate_near_uncharacterized_orphan_still_scores_novel(tmp_path):
    """An orphan has no known ligand, so it is NOT a known non-chemoreceptor and
    must not define 'known' space: a candidate resembling one stays novel."""
    channel = _maha(tmp_path, {"on_aminergic": np.eye(8)[0] * 5.0,
                               "on_orphan": np.eye(8)[2] * 50.0})
    assert channel["on_orphan"]["emb_novelty"] > channel["on_aminergic"]["emb_novelty"]


def test_nearest_family_annotation_is_always_a_characterized_classA_family(tmp_path):
    """emb_nonchemo_family annotates which KNOWN non-chemoreceptor family a
    candidate most resembles -- it must never report an out-of-class or orphan
    family."""
    channel = _maha(tmp_path, {"on_classC": np.eye(8)[1] * 50.0,
                               "on_orphan": np.eye(8)[2] * 50.0})
    for cid in ("on_classC", "on_orphan"):
        assert channel[cid]["emb_nonchemo_family"] in {"aminergic", "peptide"}


def test_novelty_reference_labels_keeps_only_characterized_classA():
    """Contract of the novelty restriction: out-of-class (B/C/F) and
    uncharacterized (orphan) references are excluded from the novelty
    prototypes; characterized class-A families are kept."""
    labels = {
        "ANCHOR_A_1_P1": "peptide",
        "ANCHOR_A_1_A1": "aminergic",
        "ANCHOR_A_1_K1": "chemokine",
        "ANCHOR_A_1_O1": "orphan",
        "ANCHOR_B_1_B1": "class-B-secretin",
        "ANCHOR_C_1_C1": "class-C",
        "ANCHOR_F_1_F1": "class-F-frizzled",
    }
    kept = novelty_reference_labels(labels)
    assert sorted(kept.values()) == ["aminergic", "chemokine", "peptide"]


def test_novelty_reference_labels_empty_when_only_out_of_class():
    assert novelty_reference_labels({"ANCHOR_C_1_C1": "class-C"}) == {}
