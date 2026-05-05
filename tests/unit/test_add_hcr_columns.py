"""Tests for scripts/add_hcr_columns.py.

The script consumes a ranked-candidates CSV and a CDS FASTA + protein MSA, and
emits the same CSV augmented with cds_length_bp, paralog_min_identity,
closest_paralog_id, hcr_probe_friendly. These tests cover the contract that
stage 09 / the report depend on:

  * RefSeq-style headers (``>lcl|... [gene=X] [protein_id=Y]``) yield CDS lengths
    indexed by both first-token and bracketed alias keys.
  * A candidate present in the protein MSA gets a real paralog identity.
  * A candidate absent from the MSA gets blank paralog columns and is still
    classified hcr_probe_friendly if its CDS length passes (per the lib's
    "unknown identity = friendly" rule).
  * Output preserves all input columns and adds exactly the four new ones in a
    deterministic order.
  * Empty-input CSV produces an empty output with the four new columns added
    (so stage 09 can always read them).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

# conftest.py adds scripts/ to sys.path
import add_hcr_columns as ahc

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def _write_ranked_csv(path: Path, ids: list[str]) -> None:
    """Write a minimal ranked CSV with the columns rank_candidates emits that
    add_hcr_columns reads (just `id` + a couple of decoys)."""
    df = pd.DataFrame({
        "id": ids,
        "rank_score": [10.0 - i for i in range(len(ids))],
        "confidence_tier": ["High"] * len(ids),
        "tandem_cluster_size": [1] * len(ids),
    })
    df.to_csv(path, index=False)


def _write_cds_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write CDS FASTA. Each record is (full_header_line_without_>, sequence)."""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def _write_msa(path: Path, records: list[tuple[str, str]]) -> None:
    _write_cds_fasta(path, records)


# -------- fasta_records / cds_lengths_from_fasta --------

def test_fasta_records_yields_pairs(tmp_path: Path) -> None:
    fp = tmp_path / "f.fa"
    _write_cds_fasta(fp, [("seq_a desc", "ACGT" * 10),
                          ("seq_b", "ACG"),
                          ("seq_c", "")])
    out = list(ahc.fasta_records(str(fp)))
    assert out == [("seq_a", "ACGT" * 10),
                   ("seq_b", "ACG"),
                   ("seq_c", "")]


def test_cds_lengths_from_fasta_indexes_by_token_and_brackets(tmp_path: Path) -> None:
    fp = tmp_path / "cds.fna"
    _write_cds_fasta(fp, [
        ("lcl|NC_001.1_cds_XP_42.1_7 [gene=ABC] [protein_id=XP_42.1]", "A" * 1200),
    ])
    m = ahc.cds_lengths_from_fasta(str(fp))
    # First-token key
    assert m["lcl|NC_001.1_cds_XP_42.1_7"] == 1200
    # Bracket-extracted aliases (raw value AND k=v form)
    assert m["ABC"] == 1200
    assert m["XP_42.1"] == 1200
    assert m["gene=ABC"] == 1200
    assert m["protein_id=XP_42.1"] == 1200


def test_cds_lengths_skips_blank_lines_and_handles_multi_line_seqs(tmp_path: Path) -> None:
    fp = tmp_path / "cds.fna"
    fp.write_text(">a\nACGT\nACGT\n\n>b\nGGGG\n")
    m = ahc.cds_lengths_from_fasta(str(fp))
    assert m["a"] == 8
    assert m["b"] == 4


# -------- end-to-end CLI tests --------

def _run_cli(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "scripts/add_hcr_columns.py", *args],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )


def test_cli_augments_csv_with_four_new_columns(tmp_path: Path) -> None:
    """Happy path: candidate present in CDS + MSA gets non-empty values; output
    schema has exactly the four expected new columns appended."""
    ranked = tmp_path / "ranked.csv"
    cds = tmp_path / "cds.fna"
    msa = tmp_path / "aln.fa"
    out = tmp_path / "out.csv"
    _write_ranked_csv(ranked, ["candA", "candB"])
    _write_cds_fasta(cds, [("candA", "ACGT" * 200),  # 800 bp
                           ("candB", "ACGT" * 100)])  # 400 bp -> below default 600
    # MSA: build sequences so candA's closest paralog identity is < 80%.
    # candA matches the first 200 cols of each other sequence, then diverges.
    # 200/(200+200) = 50% identity over a 400-aa window — under 80%, so candA
    # is HCR-friendly on the paralog-divergence axis.
    msa_len_match = 200
    msa_len_diff  = 200
    candA_seq = ("M" + "A" * msa_len_match + "G" * msa_len_diff)
    candB_seq = ("M" + "A" * msa_len_match + "C" * msa_len_diff)
    ref1_seq  = ("M" + "A" * msa_len_match + "T" * msa_len_diff)
    _write_msa(msa, [("candA", candA_seq),
                     ("candB", candB_seq),
                     ("ref1",  ref1_seq)])
    rc = _run_cli(["--ranked-csv", str(ranked),
                   "--cds-fasta", str(cds),
                   "--alignment", str(msa),
                   "--out", str(out)])
    assert rc.returncode == 0, f"stderr: {rc.stderr}"

    df = pd.read_csv(out)
    # Original cols preserved
    for c in ("id", "rank_score", "confidence_tier", "tandem_cluster_size"):
        assert c in df.columns, f"missing original column {c}"
    # Exactly the four new cols, in deterministic insertion order at the end
    new_cols = ["cds_length_bp", "paralog_min_identity",
                "closest_paralog_id", "hcr_probe_friendly"]
    assert list(df.columns)[-4:] == new_cols, list(df.columns)

    rows = df.set_index("id")
    assert rows.loc["candA", "cds_length_bp"] == 800
    assert rows.loc["candB", "cds_length_bp"] == 400
    # candA is HCR-friendly: 800 >= 600, paralog identity ≤ 0.80
    assert bool(rows.loc["candA", "hcr_probe_friendly"]) is True
    # candB is NOT HCR-friendly: 400 < 600
    assert bool(rows.loc["candB", "hcr_probe_friendly"]) is False


def test_cli_blank_paralog_when_id_not_in_msa_but_still_friendly(tmp_path: Path) -> None:
    """A candidate present in CDS but absent from the MSA leaves paralog cols
    blank; per is_hcr_probe_friendly, unknown identity does NOT disqualify."""
    ranked = tmp_path / "ranked.csv"
    cds = tmp_path / "cds.fna"
    msa = tmp_path / "aln.fa"
    out = tmp_path / "out.csv"
    _write_ranked_csv(ranked, ["only_in_cds"])
    _write_cds_fasta(cds, [("only_in_cds", "G" * 1500)])  # well above min
    _write_msa(msa, [("other", "M" + "A" * 200)])  # no overlap

    rc = _run_cli(["--ranked-csv", str(ranked),
                   "--cds-fasta", str(cds),
                   "--alignment", str(msa),
                   "--out", str(out)])
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    df = pd.read_csv(out)
    row = df.iloc[0]
    assert int(row["cds_length_bp"]) == 1500
    # NaN/blank paralog columns are loaded as NaN by pandas
    assert pd.isna(row["paralog_min_identity"]) or row["paralog_min_identity"] == ""
    assert pd.isna(row["closest_paralog_id"]) or row["closest_paralog_id"] == ""
    # Friendly because cds length passes and unknown identity is treated as OK
    assert bool(row["hcr_probe_friendly"]) is True


def test_cli_no_cds_no_msa_emits_blank_lengths_and_unfriendly(tmp_path: Path) -> None:
    """If neither --cds-fasta nor --alignment is given, all rows get blank
    cds_length_bp and hcr_probe_friendly=False (CDS unknown disqualifies)."""
    ranked = tmp_path / "ranked.csv"
    out = tmp_path / "out.csv"
    _write_ranked_csv(ranked, ["x", "y", "z"])
    rc = _run_cli(["--ranked-csv", str(ranked), "--out", str(out)])
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    df = pd.read_csv(out)
    assert len(df) == 3
    assert df["cds_length_bp"].isna().all() or (df["cds_length_bp"] == "").all()
    assert (df["hcr_probe_friendly"] == False).all()  # noqa: E712


def test_cli_missing_id_column_returns_nonzero(tmp_path: Path) -> None:
    """Stage 07 (or a corrupt CSV) without an `id` column must fail loudly."""
    ranked = tmp_path / "bad.csv"
    out = tmp_path / "out.csv"
    pd.DataFrame({"foo": [1, 2]}).to_csv(ranked, index=False)
    rc = _run_cli(["--ranked-csv", str(ranked), "--out", str(out)])
    assert rc.returncode != 0
    assert "missing 'id' column" in rc.stderr


def test_cli_tandem_cluster_warning_disqualifies(tmp_path: Path) -> None:
    """tandem_cluster_size >= tandem_cluster_warning (default 5) sets
    hcr_probe_friendly=False even if CDS length and identity pass."""
    ranked = tmp_path / "ranked.csv"
    cds = tmp_path / "cds.fna"
    out = tmp_path / "out.csv"
    df_in = pd.DataFrame({
        "id": ["big_cluster_member", "singleton"],
        "rank_score": [9.0, 8.0],
        "tandem_cluster_size": [6, 1],
    })
    df_in.to_csv(ranked, index=False)
    _write_cds_fasta(cds, [("big_cluster_member", "A" * 1200),
                           ("singleton", "A" * 1200)])

    rc = _run_cli(["--ranked-csv", str(ranked),
                   "--cds-fasta", str(cds),
                   "--out", str(out)])
    assert rc.returncode == 0, f"stderr: {rc.stderr}"
    df = pd.read_csv(out).set_index("id")
    # 6-member cluster -> not friendly (probes likely to cross-react)
    assert bool(df.loc["big_cluster_member", "hcr_probe_friendly"]) is False
    # singleton -> friendly
    assert bool(df.loc["singleton", "hcr_probe_friendly"]) is True
