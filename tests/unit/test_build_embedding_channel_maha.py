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
