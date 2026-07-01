"""Unit tests for scripts/build_structural_channel.py.

The Foldseek structural-evidence SCORER (scripts/structural_evidence.py,
Rank Task 4) parses one Foldseek easy-search tab file and classifies a
candidate's best hit. This module is the PRODUCER that sits upstream of it:
it merges the three per-DB hit files scripts/unity/run_foldseek_candidates.sh
emits (PDB, AFDB50, GPCRdb) into one best-hit-per-query map, builds a
target-id -> family label map from the anchor set (+ optional GPCRdb
metadata), classifies every candidate, and writes the TSV that
scripts/rank_aggregation.py's merge_evidence_channels() reads to wire the
channel into the rank-aggregation reranker.

Same council rule as structural_evidence.py: structural resemblance is
honest RECALL (novelty) + EXCLUSION (non-chemoreceptor corroboration) only,
never a positive "looks like a known chemoreceptor" score.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

# conftest.py adds scripts/ to sys.path
import build_structural_channel as bsc
import rank_aggregation as ra

ANCHOR_FIELDS = ["accession", "tier", "taxid", "species", "family", "class", "evidence"]


def _write_foldseek(tmp_path: Path, name: str, lines: list) -> str:
    p = tmp_path / name
    p.write_text("\n".join(lines) + ("\n" if lines else ""))
    return str(p)


def _write_anchor_set(tmp_path: Path, rows: list, name: str = "anchor_set.tsv") -> str:
    """rows: list of dicts with (a subset of) ANCHOR_FIELDS keys."""
    p = tmp_path / name
    with open(p, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ANCHOR_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return str(p)


def _anchor_row(accession, family, **overrides):
    row = {
        "accession": accession, "tier": "1", "taxid": "6637",
        "species": "Todarodes pacificus (Japanese flying squid)",
        "family": family, "class": "A", "evidence": "reviewed",
    }
    row.update(overrides)
    return row


# --------------------------------------------------------------------------- #
# merge_best_hits
# --------------------------------------------------------------------------- #

def test_merge_best_hits_picks_max_alntmscore_across_dbs(tmp_path):
    pdb = _write_foldseek(tmp_path, "PDB.tsv", ["cand1\ttargetP\t0.40\t0.55\t1e-8"])
    afdb = _write_foldseek(tmp_path, "AFDB50.tsv", ["cand1\ttargetA\t0.70\t0.85\t1e-15"])
    gpcrdb = _write_foldseek(tmp_path, "GPCRdb.tsv", ["cand1\ttargetG\t0.60\t0.65\t1e-10"])

    merged = bsc.merge_best_hits([pdb, afdb, gpcrdb])

    assert set(merged) == {"cand1"}
    assert merged["cand1"]["target"] == "targetA"
    assert merged["cand1"]["alntmscore"] == pytest.approx(0.85)


def test_merge_best_hits_independent_per_query(tmp_path):
    pdb = _write_foldseek(tmp_path, "PDB.tsv", ["cand1\ttargetP\t0.40\t0.90\t1e-8"])
    afdb = _write_foldseek(tmp_path, "AFDB50.tsv", ["cand2\ttargetA\t0.70\t0.30\t1e-15"])

    merged = bsc.merge_best_hits([pdb, afdb])

    assert merged["cand1"]["target"] == "targetP"
    assert merged["cand2"]["target"] == "targetA"


def test_merge_best_hits_missing_files_return_empty_dict(tmp_path):
    missing1 = str(tmp_path / "nope1.tsv")
    missing2 = str(tmp_path / "nope2.tsv")
    assert bsc.merge_best_hits([missing1, missing2]) == {}


def test_merge_best_hits_empty_list_returns_empty_dict():
    assert bsc.merge_best_hits([]) == {}


def test_merge_best_hits_tolerates_some_dbs_missing(tmp_path):
    """A DB the wrapper skipped (never wrote an output file) shouldn't crash
    the merge -- the other DBs' hits still come through."""
    afdb = _write_foldseek(tmp_path, "AFDB50.tsv", ["cand1\ttargetA\t0.70\t0.85\t1e-15"])
    missing_pdb = str(tmp_path / "PDB.tsv")
    missing_gpcrdb = str(tmp_path / "GPCRdb.tsv")

    merged = bsc.merge_best_hits([missing_pdb, afdb, missing_gpcrdb])

    assert set(merged) == {"cand1"}
    assert merged["cand1"]["target"] == "targetA"


# --------------------------------------------------------------------------- #
# build_family_map
# --------------------------------------------------------------------------- #

def test_build_family_map_from_anchor_set(tmp_path):
    anchor_tsv = _write_anchor_set(tmp_path, [_anchor_row("P31356", "opsin")])

    fam = bsc.build_family_map(anchor_tsv)

    assert fam["P31356"] == "opsin"


def test_build_family_map_covers_multiple_accessions(tmp_path):
    anchor_tsv = _write_anchor_set(tmp_path, [
        _anchor_row("P31356", "opsin"),
        _anchor_row("Q5W9T5", "class-A-other"),
        _anchor_row("O15973", "aminergic_5HT"),
    ])

    fam = bsc.build_family_map(anchor_tsv)

    assert fam == {
        "P31356": "opsin",
        "Q5W9T5": "class-A-other",
        "O15973": "aminergic_5HT",
    }


def test_build_family_map_merges_optional_gpcrdb_meta(tmp_path):
    anchor_tsv = _write_anchor_set(tmp_path, [_anchor_row("P31356", "opsin")])
    gpcrdb_meta = tmp_path / "gpcrdb_meta.tsv"
    gpcrdb_meta.write_text("target\tfamily\ngpcrdb_101\taminergic_5HT\n")

    fam = bsc.build_family_map(anchor_tsv, gpcrdb_meta=str(gpcrdb_meta))

    assert fam["P31356"] == "opsin"
    assert fam["gpcrdb_101"] == "aminergic_5HT"


def test_build_family_map_gpcrdb_meta_none_is_fine(tmp_path):
    anchor_tsv = _write_anchor_set(tmp_path, [_anchor_row("Q5W9T5", "class-A-other")])

    fam = bsc.build_family_map(anchor_tsv, gpcrdb_meta=None)

    assert fam == {"Q5W9T5": "class-A-other"}


def test_build_family_map_missing_anchor_set_returns_empty_dict(tmp_path):
    missing = str(tmp_path / "nope.tsv")
    assert bsc.build_family_map(missing) == {}


def test_build_family_map_missing_gpcrdb_meta_is_graceful(tmp_path):
    anchor_tsv = _write_anchor_set(tmp_path, [_anchor_row("P31356", "opsin")])
    missing_meta = str(tmp_path / "nope_meta.tsv")

    fam = bsc.build_family_map(anchor_tsv, gpcrdb_meta=missing_meta)

    assert fam == {"P31356": "opsin"}


# --------------------------------------------------------------------------- #
# build_structural_channel
# --------------------------------------------------------------------------- #

def test_build_structural_channel_end_to_end_across_dbs(tmp_path):
    """cand_nonchemo's best hit is in AFDB50 (highest alntmscore of its three
    per-DB hits) and maps via the anchor family map to an opsin -> exclusion
    corroboration. cand_novel has no confident hit anywhere. cand_other hits
    a target absent from the family map -> known_other, graceful."""
    pdb = _write_foldseek(tmp_path, "PDB.tsv", [
        "cand_nonchemo\ttargetP\t0.30\t0.40\t1e-5",
        "cand_novel\ttargetQ\t0.10\t0.15\t1e-2",
    ])
    afdb = _write_foldseek(tmp_path, "AFDB50.tsv", [
        "cand_nonchemo\tP31356\t0.85\t0.90\t1e-25",
    ])
    gpcrdb = _write_foldseek(tmp_path, "GPCRdb.tsv", [
        "cand_other\tunmapped_target\t0.85\t0.90\t1e-25",
    ])
    anchor_tsv = _write_anchor_set(tmp_path, [_anchor_row("P31356", "opsin")])
    family_map = bsc.build_family_map(anchor_tsv)

    channel = bsc.build_structural_channel([pdb, afdb, gpcrdb], family_map)

    assert channel["cand_nonchemo"]["struct_state"] == "known_non_chemoreceptor"
    assert channel["cand_nonchemo"]["struct_nonchemo_corrob"] == 1
    assert channel["cand_nonchemo"]["struct_novelty"] == 0

    assert channel["cand_novel"]["struct_state"] == "novel"
    assert channel["cand_novel"]["struct_novelty"] == 1
    assert channel["cand_novel"]["struct_nonchemo_corrob"] == 0

    assert channel["cand_other"]["struct_state"] == "known_other"
    assert channel["cand_other"]["struct_novelty"] == 0
    assert channel["cand_other"]["struct_nonchemo_corrob"] == 0

    for row in channel.values():
        assert row["has_struct_data"] is True


def test_build_structural_channel_empty_inputs_returns_empty_dict():
    assert bsc.build_structural_channel([], {}) == {}


# --------------------------------------------------------------------------- #
# write_channel_tsv <-> rank_aggregation.merge_evidence_channels schema contract
#
# rank_aggregation.merge_evidence_channels() reads a struct TSV via a
# ("has_struct_data", ["struct_novelty", "struct_nonchemo_corrob",
# "struct_state"]) entry in its `channels` tuple, joined on id_col="id" (see
# scripts/rank_aggregation.py). Pinned here as a literal list (read directly
# from that source on 2026-07-01); the round-trip test below is the
# load-bearing proof since it exercises the real function, not this pin.
# --------------------------------------------------------------------------- #
EXPECTED_STRUCT_VALUE_COLS = ["struct_novelty", "struct_nonchemo_corrob", "struct_state"]


def test_write_channel_tsv_columns_match_merge_evidence_channels_contract(tmp_path):
    channel = {
        "cand1": {"struct_state": "novel", "struct_novelty": 1,
                  "struct_nonchemo_corrob": 0, "has_struct_data": True},
    }
    out = tmp_path / "struct_channel.tsv"

    bsc.write_channel_tsv(channel, str(out))

    with open(out) as fh:
        header = fh.readline().strip().split("\t")
    assert "id" in header
    for col in EXPECTED_STRUCT_VALUE_COLS:
        assert col in header, f"missing column {col!r} merge_evidence_channels expects"


def test_write_channel_tsv_round_trips_through_merge_evidence_channels(tmp_path):
    """Load-bearing integration proof: the TSV write_channel_tsv() produces is
    actually consumed correctly by the real rank_aggregation.merge_evidence_
    channels() -- values, not just column names, survive the round trip -- and
    the exclusion signal it carries actually lowers rank-aggregated standing."""
    channel = {
        "good": {"struct_state": "novel", "struct_novelty": 1,
                 "struct_nonchemo_corrob": 0, "has_struct_data": True},
        "bad": {"struct_state": "known_non_chemoreceptor", "struct_novelty": 0,
                "struct_nonchemo_corrob": 1, "has_struct_data": True},
    }
    out = tmp_path / "struct_channel.tsv"
    bsc.write_channel_tsv(channel, str(out))

    df = pd.DataFrame({"id": ["good", "bad", "absent"]})
    merged = ra.merge_evidence_channels(df, struct_tsv=str(out))

    assert list(merged["has_struct_data"]) == [True, True, False]
    indexed = merged.set_index("id")
    assert indexed.loc["good", "struct_novelty"] == 1
    assert indexed.loc["bad", "struct_nonchemo_corrob"] == 1
    assert indexed.loc["good", "struct_state"] == "novel"
    assert pd.isna(indexed.loc["absent", "struct_novelty"])

    order = ra.aggregate(ra.build_ranklists_from_df(merged), method="rrf")
    assert order.index("good") < order.index("bad")


def test_write_channel_tsv_empty_channel_writes_header_only(tmp_path):
    out = tmp_path / "empty.tsv"

    bsc.write_channel_tsv({}, str(out))

    df = pd.read_csv(out, sep="\t")
    assert len(df) == 0
    for col in ["id"] + EXPECTED_STRUCT_VALUE_COLS:
        assert col in df.columns


def test_write_channel_tsv_sorted_by_id_for_deterministic_output(tmp_path):
    channel = {
        "zeta": {"struct_state": "novel", "struct_novelty": 1,
                 "struct_nonchemo_corrob": 0, "has_struct_data": True},
        "alpha": {"struct_state": "novel", "struct_novelty": 1,
                  "struct_nonchemo_corrob": 0, "has_struct_data": True},
    }
    out = tmp_path / "struct_channel.tsv"

    bsc.write_channel_tsv(channel, str(out))

    df = pd.read_csv(out, sep="\t")
    assert list(df["id"]) == ["alpha", "zeta"]


# --------------------------------------------------------------------------- #
# CLI (main)
# --------------------------------------------------------------------------- #

def test_cli_end_to_end_writes_consumable_tsv(tmp_path):
    pdb = _write_foldseek(tmp_path, "PDB.tsv", ["cand1\tP31356\t0.85\t0.90\t1e-25"])
    missing_afdb = str(tmp_path / "AFDB50.tsv")
    missing_gpcrdb = str(tmp_path / "GPCRdb.tsv")
    anchor_tsv = _write_anchor_set(tmp_path, [_anchor_row("P31356", "opsin")])
    out = tmp_path / "struct_channel.tsv"

    rc = bsc.main([
        "--foldseek-tsvs", pdb, missing_afdb, missing_gpcrdb,
        "--anchor-set", anchor_tsv,
        "--out", str(out),
    ])

    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    row = df.set_index("id").loc["cand1"]
    assert row["struct_state"] == "known_non_chemoreceptor"
    assert int(row["struct_novelty"]) == 0
    assert int(row["struct_nonchemo_corrob"]) == 1

    # And the CLI's output is directly consumable by merge_evidence_channels.
    merged = ra.merge_evidence_channels(pd.DataFrame({"id": ["cand1"]}), struct_tsv=str(out))
    assert merged.loc[0, "has_struct_data"] is True or bool(merged.loc[0, "has_struct_data"]) is True


def test_cli_optional_gpcrdb_meta_flag(tmp_path):
    pdb = _write_foldseek(tmp_path, "PDB.tsv", ["cand1\tgpcrdb_101\t0.85\t0.90\t1e-25"])
    missing_afdb = str(tmp_path / "AFDB50.tsv")
    missing_gpcrdb = str(tmp_path / "GPCRdb.tsv")
    anchor_tsv = _write_anchor_set(tmp_path, [_anchor_row("P31356", "opsin")])
    gpcrdb_meta = tmp_path / "gpcrdb_meta.tsv"
    gpcrdb_meta.write_text("target\tfamily\ngpcrdb_101\taminergic_5HT\n")
    out = tmp_path / "struct_channel.tsv"

    rc = bsc.main([
        "--foldseek-tsvs", pdb, missing_afdb, missing_gpcrdb,
        "--anchor-set", anchor_tsv,
        "--gpcrdb-meta", str(gpcrdb_meta),
        "--out", str(out),
    ])

    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert df.set_index("id").loc["cand1", "struct_state"] == "known_non_chemoreceptor"
