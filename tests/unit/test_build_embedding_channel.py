"""Tests for scripts/build_embedding_channel.py.

Glue task G2 — the ESM-C embedding-channel PRODUCER. Wires the pure scorer
(embedding_evidence.py: load_embeddings/centroid/nearest_family/
embedding_channel) to real inputs: a candidate embeddings .npz, a reference
embeddings .npz, and a family/class-labeled reference TSV
(references/anchors/anchor_set.tsv), and writes the per-candidate channel out
as a TSV rank_aggregation.merge_evidence_channels can left-join by id.

HARD GUARDRAIL (embedding_evidence.py's module contract, inherited here):
this producer must never build a chemoreceptor centroid or emit any positive
similarity-to-known-chemoreceptor key. The reference set consumed here
(anchor_set.tsv / curate_gpcr_references.py's output) is the curated
NON-chemoreceptor GPCR set by construction -- every family it contains is
eligible for nonchemo_centroids -- but build_family_centroids defensively
drops any family literally naming "chemoreceptor" anyway (contamination
guard), and the emitted TSV's column set is fixed to exactly the 5 keys
embedding_channel() documents, so no chemoreceptor-similarity key can ever
appear in this producer's output.

Coverage:
    - family_to_class: mirrors build_anchor_set.py's documented
      "family -> class rule: class-B*->B, class-C*->C, class-F*->F, else A"
    - load_ref_labels: {accession: family} from a TSV; missing file -> {}
    - build_family_centroids: one centroid per non-chemoreceptor family +
      a classA centroid over every class-A-mapped ref; drops any
      "chemoreceptor"-named family defensively (both outputs)
    - build_embedding_channel: end-to-end candidate_npz + ref_npz +
      ref_labels_tsv -> embedding_evidence.embedding_channel() output
    - write_channel_tsv: output TSV columns EXACTLY match
      rank_aggregation.merge_evidence_channels' expected emb columns, and
      genuinely interop with it
    - guard: no *chemoreceptor_sim* key/column anywhere in this producer's
      output
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import rank_aggregation as ra
from build_embedding_channel import (
    build_embedding_channel,
    build_family_centroids,
    family_to_class,
    load_ref_labels,
    main,
    write_channel_tsv,
)
from embedding_evidence import centroid


# ---- family_to_class ---------------------------------------------------------

@pytest.mark.parametrize(
    "family,expected_class",
    [
        ("aminergic", "A"),
        ("class-A-other", "A"),
        ("class-B-secretin", "B"),
        ("class-C", "C"),
        ("class-F-frizzled", "F"),
        ("lipid", "A"),
        ("opsin", "A"),
        ("peptide", "A"),
    ],
)
def test_family_to_class_matches_build_anchor_set_rule(family, expected_class):
    """Drift guard: must reproduce build_anchor_set.py's documented rule
    ("family -> class rule: class-B*->B, class-C*->C, class-F*->F, else A")
    for every family literal actually present in anchor_set.tsv."""
    assert family_to_class(family) == expected_class


# ---- load_ref_labels ----------------------------------------------------------

def test_load_ref_labels_missing_file_returns_empty_dict(tmp_path):
    assert load_ref_labels(str(tmp_path / "missing.tsv")) == {}


def test_load_ref_labels_reads_accession_and_family_columns(tmp_path):
    tsv = tmp_path / "labels.tsv"
    pd.DataFrame(
        {
            "accession": ["P1", "P2", "P3"],
            "tier": [1, 1, 2],
            "family": ["opsin", "peptide", "class-B-secretin"],
            "class": ["A", "A", "B"],
        }
    ).to_csv(tsv, sep="\t", index=False)
    labels = load_ref_labels(str(tsv))
    assert labels == {"P1": "opsin", "P2": "peptide", "P3": "class-B-secretin"}


# ---- build_family_centroids ----------------------------------------------------

def test_build_family_centroids_one_centroid_per_nonchemo_family():
    ref_embeddings = {
        "o1": np.array([1.0, 0.0]),
        "o2": np.array([1.0, 0.0]),
        "p1": np.array([0.0, 1.0]),
    }
    ref_labels = {"o1": "opsin", "o2": "opsin", "p1": "peptide"}
    nonchemo, _ = build_family_centroids(ref_embeddings, ref_labels)
    assert set(nonchemo) == {"opsin", "peptide"}
    assert np.allclose(nonchemo["opsin"], centroid(ref_embeddings, ["o1", "o2"]))
    assert np.allclose(nonchemo["peptide"], centroid(ref_embeddings, ["p1"]))


def test_build_family_centroids_classA_centroid_includes_only_class_a_families():
    ref_embeddings = {
        "a1": np.array([1.0, 0.0, 0.0]),
        "b1": np.array([0.0, 1.0, 0.0]),
        "c1": np.array([0.0, 0.0, 1.0]),
    }
    ref_labels = {"a1": "opsin", "b1": "class-B-secretin", "c1": "class-C"}
    _, classA_centroid = build_family_centroids(ref_embeddings, ref_labels)
    # only a1 (opsin -> class A) should contribute
    assert np.allclose(classA_centroid, centroid(ref_embeddings, ["a1"]))


def test_build_family_centroids_drops_family_containing_chemoreceptor():
    """Defensive contamination guard: even if a future reference-set edit
    smuggled in a family literally naming 'chemoreceptor', it must be
    excluded from BOTH nonchemo_centroids AND classA_centroid -- never
    silently centroided as if it were a legitimate non-chemoreceptor
    exclusion family."""
    ref_embeddings = {
        "x1": np.array([1.0, 0.0]),
        "ok1": np.array([0.0, 1.0]),
    }
    ref_labels = {"x1": "candidate-chemoreceptor", "ok1": "opsin"}
    nonchemo, classA_centroid = build_family_centroids(ref_embeddings, ref_labels)
    assert "candidate-chemoreceptor" not in nonchemo
    assert set(nonchemo) == {"opsin"}
    # x1 must not leak into classA_centroid either, even though the
    # fallback family->class rule would otherwise call it class A
    assert np.allclose(classA_centroid, centroid(ref_embeddings, ["ok1"]))


def test_build_family_centroids_empty_ref_labels_returns_empty_and_zero_vector():
    ref_embeddings = {"a1": np.array([1.0, 0.0])}
    nonchemo, classA_centroid = build_family_centroids(ref_embeddings, {})
    assert nonchemo == {}
    assert np.array_equal(classA_centroid, np.zeros(2))


# ---- build_embedding_channel (end-to-end via .npz + TSV fixtures) --------------

def _write_npz(path, vectors):
    np.savez(path, **vectors)
    return str(path)


def _write_labels_tsv(path, rows):
    pd.DataFrame(rows, columns=["accession", "family", "class"]).to_csv(
        path, sep="\t", index=False
    )
    return str(path)


def test_build_embedding_channel_end_to_end(tmp_path):
    candidate_npz = _write_npz(
        tmp_path / "candidates.npz",
        {"cand_1": np.array([1.0, 0.0]), "cand_2": np.array([0.0, 1.0])},
    )
    ref_npz = _write_npz(
        tmp_path / "ref.npz",
        {
            "P1": np.array([1.0, 0.0]),  # opsin
            "P2": np.array([0.0, 1.0]),  # peptide
        },
    )
    ref_labels_tsv = _write_labels_tsv(
        tmp_path / "labels.tsv",
        [["P1", "opsin", "A"], ["P2", "peptide", "A"]],
    )

    channel = build_embedding_channel(candidate_npz, ref_npz, ref_labels_tsv)

    assert set(channel) == {"cand_1", "cand_2"}
    assert channel["cand_1"]["emb_nonchemo_family"] == "opsin"
    assert channel["cand_1"]["has_emb_data"] is True
    assert channel["cand_1"]["emb_leakage_flag"] is True


def test_build_embedding_channel_missing_candidate_npz_is_empty_not_a_crash(tmp_path):
    ref_npz = _write_npz(tmp_path / "ref.npz", {"P1": np.array([1.0, 0.0])})
    ref_labels_tsv = _write_labels_tsv(tmp_path / "labels.tsv", [["P1", "opsin", "A"]])
    channel = build_embedding_channel(
        str(tmp_path / "does_not_exist.npz"), ref_npz, ref_labels_tsv
    )
    assert channel == {}


# ---- write_channel_tsv: schema MUST match merge_evidence_channels -------------

_EXPECTED_TSV_COLUMNS = [
    "id",
    "emb_nonchemo_sim",
    "emb_nonchemo_family",
    "emb_classA_sim",
    "has_emb_data",
    "emb_leakage_flag",
]


def test_write_channel_tsv_columns_exactly_match_merge_evidence_channels(tmp_path):
    channel = {
        "cand_1": {
            "emb_nonchemo_sim": 0.4,
            "emb_nonchemo_family": "opsin",
            "emb_classA_sim": 0.8,
            "has_emb_data": True,
            "emb_leakage_flag": True,
        }
    }
    out = tmp_path / "channel.tsv"
    write_channel_tsv(channel, str(out))
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _EXPECTED_TSV_COLUMNS


def test_write_channel_tsv_empty_channel_writes_header_only(tmp_path):
    out = tmp_path / "channel.tsv"
    write_channel_tsv({}, str(out))
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _EXPECTED_TSV_COLUMNS
    assert len(df) == 0


def test_write_channel_tsv_rows_are_sorted_by_id(tmp_path):
    channel = {
        "z_cand": {"emb_nonchemo_sim": 0.1, "emb_nonchemo_family": "opsin",
                    "emb_classA_sim": 0.1, "has_emb_data": True, "emb_leakage_flag": True},
        "a_cand": {"emb_nonchemo_sim": 0.2, "emb_nonchemo_family": "peptide",
                    "emb_classA_sim": 0.2, "has_emb_data": True, "emb_leakage_flag": True},
    }
    out = tmp_path / "channel.tsv"
    write_channel_tsv(channel, str(out))
    df = pd.read_csv(out, sep="\t")
    assert list(df["id"]) == ["a_cand", "z_cand"]


def test_write_channel_tsv_interops_with_merge_evidence_channels(tmp_path):
    """Real proof of interop, not just column-name matching: feed the
    written TSV straight into rank_aggregation.merge_evidence_channels and
    confirm the left-join + has_emb_data flag behave as documented there."""
    candidate_npz = _write_npz(
        tmp_path / "candidates.npz",
        {"cand_1": np.array([1.0, 0.0]), "cand_2": np.array([0.0, 1.0])},
    )
    ref_npz = _write_npz(tmp_path / "ref.npz", {"P1": np.array([1.0, 0.0])})
    ref_labels_tsv = _write_labels_tsv(tmp_path / "labels.tsv", [["P1", "opsin", "A"]])
    channel = build_embedding_channel(candidate_npz, ref_npz, ref_labels_tsv)
    emb_tsv = tmp_path / "emb_channel.tsv"
    write_channel_tsv(channel, str(emb_tsv))

    df = pd.DataFrame({"id": ["cand_1", "cand_2", "cand_3_no_embedding"]})
    merged = ra.merge_evidence_channels(df, emb_tsv=str(emb_tsv))

    assert list(merged["has_emb_data"]) == [True, True, False]
    assert merged.set_index("id").loc["cand_1", "emb_classA_sim"] == pytest.approx(1.0)
    assert pd.isna(merged.set_index("id").loc["cand_3_no_embedding", "emb_classA_sim"])


# ---- HARD GUARDRAIL: no chemoreceptor-similarity key anywhere -----------------

def test_no_chemoreceptor_sim_key_anywhere_in_channel_or_tsv(tmp_path):
    candidate_npz = _write_npz(
        tmp_path / "candidates.npz",
        {"cand_1": np.array([1.0, 0.0]), "cand_2": np.array([0.0, 1.0])},
    )
    ref_npz = _write_npz(
        tmp_path / "ref.npz",
        {"P1": np.array([1.0, 0.0]), "P2": np.array([0.0, 1.0])},
    )
    ref_labels_tsv = _write_labels_tsv(
        tmp_path / "labels.tsv", [["P1", "opsin", "A"], ["P2", "peptide", "A"]]
    )
    channel = build_embedding_channel(candidate_npz, ref_npz, ref_labels_tsv)
    assert channel  # sanity: guard below must not vacuously pass

    for candidate_id, fields in channel.items():
        for key in fields:
            assert "chemoreceptor_sim" not in key.lower(), (
                f"banned positive-similarity key {key!r} found for {candidate_id}"
            )

    out = tmp_path / "channel.tsv"
    write_channel_tsv(channel, str(out))
    df = pd.read_csv(out, sep="\t")
    for col in df.columns:
        assert "chemoreceptor_sim" not in col.lower()
    assert list(df.columns) == _EXPECTED_TSV_COLUMNS


# ---- main() CLI -----------------------------------------------------------------

def test_main_cli_writes_expected_tsv(tmp_path):
    candidate_npz = _write_npz(
        tmp_path / "candidates.npz", {"cand_1": np.array([1.0, 0.0])}
    )
    ref_npz = _write_npz(tmp_path / "ref.npz", {"P1": np.array([1.0, 0.0])})
    ref_labels_tsv = _write_labels_tsv(tmp_path / "labels.tsv", [["P1", "opsin", "A"]])
    out = tmp_path / "out.tsv"

    main(
        [
            "--candidate-npz", candidate_npz,
            "--ref-npz", ref_npz,
            "--ref-labels", ref_labels_tsv,
            "--out", str(out),
        ]
    )

    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == _EXPECTED_TSV_COLUMNS
    assert list(df["id"]) == ["cand_1"]
