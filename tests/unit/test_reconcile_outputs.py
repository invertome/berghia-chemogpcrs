"""Tests for scripts/reconcile_candidates.py — Task 5.

Genome-track + candidate reconciliation feature (bead: genome-track
stage 02c). Task 5 covers only the output writers (``write_tsv`` /
``write_faa`` / ``write_report``) and the ``main`` argparse CLI. The
placement cascade, ``reconcile()``, loci grouping, and the aligner
parsers belong to earlier tasks and are exercised by their own suites.

Per the plan's verify these tests focus on the WRITERS + ``--help``;
the full file-parsing wiring of ``main`` is validated by the Task-6
stage smoke test and real runs, not here.

Serialization rules under test (design §8):
  * ``pct_identity`` / ``pct_coverage``: None -> ""; float -> ``f"{x:g}"``.
  * ``best_vs_second_margin``: None -> ""; ``inf`` -> "inf"; float -> ``f"{x:g}"``.
  * tuple fields (isoform id lists, qc_flags): ``";".join``; () -> "".
  * TSV header = the ``ReconciledGene`` field names in §8 order; LF endings.
"""
from __future__ import annotations

import dataclasses

import pytest

import reconcile_candidates as rc


# ---- builders -------------------------------------------------------------

def _gene(**overrides):
    """A ReconciledGene with sensible defaults, overridable per field."""
    base = dict(
        gene_id="chr1:100-200:+",
        representative_id="ref_1",
        representative_source="genome",
        provenance="both",
        genome_locus="chr1:100-200:+",
        placement_method="minimap2+gmap_concordant",
        placement_confidence="high",
        pct_identity=97.5,
        pct_coverage=99.0,
        best_vs_second_margin=5.0,
        transcriptome_isoform_ids=("tx_1", "tx_2"),
        genome_isoform_ids=("ref_1",),
        n_tm=7,
        completeness="complete",
        qc_flags=(),
    )
    base.update(overrides)
    return rc.ReconciledGene(**base)


_FIELD_NAMES = [f.name for f in dataclasses.fields(rc.ReconciledGene)]
_COL = {name: i for i, name in enumerate(_FIELD_NAMES)}


# ---- write_tsv ------------------------------------------------------------

def test_write_tsv_header_is_schema_order(tmp_path):
    path = tmp_path / "out.tsv"
    rc.write_tsv([_gene()], str(path))
    header = path.read_text().splitlines()[0]
    assert header == "\t".join(_FIELD_NAMES)
    # Pin the literal §8 column set/order so a schema drift is caught here.
    assert header.split("\t") == [
        "gene_id", "representative_id", "representative_source", "provenance",
        "genome_locus", "placement_method", "placement_confidence",
        "pct_identity", "pct_coverage", "best_vs_second_margin",
        "transcriptome_isoform_ids", "genome_isoform_ids", "n_tm",
        "completeness", "qc_flags"]


def test_write_tsv_serializes_margin_inf_and_none(tmp_path):
    genes = [
        _gene(gene_id="g_inf", best_vs_second_margin=float("inf")),
        _gene(gene_id="g_none", pct_identity=None, pct_coverage=None,
              best_vs_second_margin=None),
        _gene(gene_id="g_finite", best_vs_second_margin=2.5, pct_identity=100.0),
    ]
    path = tmp_path / "out.tsv"
    rc.write_tsv(genes, str(path))
    rows = {r.split("\t")[0]: r.split("\t")
            for r in path.read_text().splitlines()[1:]}

    # inf = single-locus unambiguous placement (Task-4 hand-off).
    assert rows["g_inf"][_COL["best_vs_second_margin"]] == "inf"
    # None = genome_native / unplaced -> empty cell, never "None".
    assert rows["g_none"][_COL["best_vs_second_margin"]] == ""
    assert rows["g_none"][_COL["pct_identity"]] == ""
    assert rows["g_none"][_COL["pct_coverage"]] == ""
    # finite floats via %g (trailing .0 stripped).
    assert rows["g_finite"][_COL["best_vs_second_margin"]] == "2.5"
    assert rows["g_finite"][_COL["pct_identity"]] == "100"


def test_write_tsv_joins_tuple_fields_with_semicolons(tmp_path):
    genes = [_gene(gene_id="g1",
                   transcriptome_isoform_ids=("tx_1", "tx_2"),
                   genome_isoform_ids=(),
                   qc_flags=("low_margin", "multi_mapping"),
                   n_tm=6)]
    path = tmp_path / "out.tsv"
    rc.write_tsv(genes, str(path))
    row = path.read_text().splitlines()[1].split("\t")
    assert row[_COL["transcriptome_isoform_ids"]] == "tx_1;tx_2"
    assert row[_COL["genome_isoform_ids"]] == ""            # empty tuple -> ""
    assert row[_COL["qc_flags"]] == "low_margin;multi_mapping"
    assert row[_COL["n_tm"]] == "6"


def test_write_tsv_row_order_follows_input(tmp_path):
    genes = [_gene(gene_id="zzz"), _gene(gene_id="aaa"), _gene(gene_id="mmm")]
    path = tmp_path / "out.tsv"
    rc.write_tsv(genes, str(path))
    ids = [r.split("\t")[0] for r in path.read_text().splitlines()[1:]]
    assert ids == ["zzz", "aaa", "mmm"]


def test_write_tsv_uses_lf_endings(tmp_path):
    path = tmp_path / "out.tsv"
    rc.write_tsv([_gene()], str(path))
    raw = path.read_bytes()
    assert b"\r\n" not in raw
    assert raw.endswith(b"\n")


# ---- write_faa ------------------------------------------------------------

def test_write_faa_one_record_per_gene_with_sequence(tmp_path):
    genes = [
        _gene(gene_id="g1", representative_id="ref_1"),
        _gene(gene_id="g2", representative_id="tx_9"),
        _gene(gene_id="g3", representative_id="missing_rep"),
    ]
    seqs = {"ref_1": "MAAAA", "tx_9": "MBBBB"}   # missing_rep intentionally absent
    path = tmp_path / "out.faa"
    rc.write_faa(genes, seqs, str(path))
    content = path.read_text()

    headers = [ln for ln in content.splitlines() if ln.startswith(">")]
    assert headers == [">ref_1", ">tx_9"]        # gene order, missing skipped
    assert content.count(">") == 2
    assert ">missing_rep" not in content
    assert "MAAAA" in content and "MBBBB" in content
    assert b"\r\n" not in path.read_bytes()


def test_write_faa_all_missing_writes_empty_file(tmp_path):
    path = tmp_path / "out.faa"
    rc.write_faa([_gene(representative_id="nope")], {}, str(path))
    assert path.read_text() == ""


# ---- write_report ---------------------------------------------------------

def test_write_report_case_counts_and_review_lists(tmp_path):
    genes = [
        _gene(gene_id="both_1", provenance="both", genome_locus="chr1:1-2:+"),
        _gene(gene_id="genome_1", provenance="genome_only",
              genome_locus="chr2:1-2:+", representative_source="genome"),
        _gene(gene_id="tx_placed_1", provenance="transcriptome_only",
              genome_locus="chr3:1-2:+"),
        _gene(gene_id="tx_unplaced_1", provenance="transcriptome_only",
              genome_locus="unplaced"),
        _gene(gene_id="chimeric_1", provenance="both",
              genome_locus="chr4:1-2:+", qc_flags=("chimeric",)),
        _gene(gene_id="lowmargin_1", provenance="transcriptome_only",
              genome_locus="chr5:1-2:+", qc_flags=("low_margin",)),
    ]
    path = tmp_path / "report.md"
    rc.write_report(genes, str(path))
    text = path.read_text()

    assert f"Total reconciled genes: {len(genes)}" in text
    # per-case provenance counts
    assert "- both (shared): 2" in text                    # both_1 + chimeric_1
    assert "- genome_only (recovered): 1" in text
    assert "- transcriptome_only, placed: 2" in text       # tx_placed_1 + lowmargin_1
    assert "- transcriptome_only, unplaced: 1" in text
    # flag counts
    assert "- chimeric: 1" in text
    assert "- low_margin: 1" in text
    assert "- source_disagreement: 0" in text
    # flagged review lists carry the right gene_ids ("- " => list line)
    assert "- chimeric_1" in text
    assert "- lowmargin_1" in text


def test_write_report_is_deterministic(tmp_path):
    genes = [_gene(gene_id="both_1"),
             _gene(gene_id="tx_1", provenance="transcriptome_only",
                   genome_locus="unplaced")]
    p1, p2 = tmp_path / "r1.md", tmp_path / "r2.md"
    rc.write_report(genes, str(p1))
    rc.write_report(genes, str(p2))
    assert p1.read_text() == p2.read_text()


# ---- main CLI -------------------------------------------------------------

def test_main_help_exits_zero(capsys):
    with pytest.raises(SystemExit) as excinfo:
        rc.main(["--help"])
    assert excinfo.value.code == 0


# ---- pure CLI helpers -----------------------------------------------------

@pytest.mark.parametrize("token", [
    "1", "true", "TRUE", "True", "t", "T", "yes", "Yes", "y", "Y",
    "complete", "Complete", "COMPLETE", " complete ", "1 ", "\ttrue"])
def test_parse_bool_truthy(token):
    assert rc._parse_bool(token) is True


@pytest.mark.parametrize("token", [
    "0", "false", "False", "partial", "PARTIAL", "", "   ", "no", "n",
    "2", "completed", "truthy"])
def test_parse_bool_falsy(token):
    # substrings/near-misses ("completed", "truthy") must NOT read as True.
    assert rc._parse_bool(token) is False


def test_read_fasta_joins_multiline_and_uses_first_token(tmp_path):
    p = tmp_path / "in.faa"
    p.write_text(">ref_1 some description here\nMAAA\nAAKK\n"
                 ">tx_2\nMBBB\n")
    seqs = rc._read_fasta(str(p))
    # id = first whitespace token of the header; sequence lines joined.
    assert seqs == {"ref_1": "MAAAAAKK", "tx_2": "MBBB"}


def test_read_fasta_empty_header(tmp_path):
    p = tmp_path / "in.faa"
    p.write_text(">\nMXYZ\n")
    seqs = rc._read_fasta(str(p))
    assert seqs == {"": "MXYZ"}


def _placement(query, method, chrom="chr1", start=100, end=500, strand="+"):
    return rc.Placement(query=query, chrom=chrom, start=start, end=end,
                        strand=strand, pct_identity=99.0, pct_coverage=100.0,
                        method=method, score=60.0)


def test_group_placements_buckets_by_query_and_method():
    a1, a2 = _placement("tx_A", "minimap2"), _placement("tx_A", "gmap")
    b1 = _placement("tx_B", "minimap2")
    grouped = rc._group_placements(("minimap2", [a1, b1]), ("gmap", [a2]),
                                   ("miniprot", []), ("rbh", []))
    assert set(grouped) == {"tx_A", "tx_B"}
    assert grouped["tx_A"]["minimap2"] == [a1]
    assert grouped["tx_A"]["gmap"] == [a2]
    assert grouped["tx_B"]["minimap2"] == [b1]
    assert "gmap" not in grouped["tx_B"]          # only present methods bucketed


def test_group_placements_routes_blastp_to_rbh_key():
    # The BLASTp parser tags method="blastp"; the caller buckets it under the
    # cascade's "rbh" key. The bucket key is the caller's pairing, NOT
    # Placement.method (reconcile reads bundle["rbh"]).
    hit = _placement("tx_C", "blastp")
    grouped = rc._group_placements(("minimap2", []), ("gmap", []),
                                   ("miniprot", []), ("rbh", [hit]))
    assert grouped["tx_C"]["rbh"] == [hit]
    assert "blastp" not in grouped["tx_C"]


# ---- metadata readers: fail loud on malformed rows ------------------------

def test_read_candidates_raises_on_missing_field(tmp_path):
    # tx_B is missing n_tm -> MUST raise, never silently vanish (a dropped
    # candidate would not even surface as unplaced -> "never drop" violated).
    p = tmp_path / "txome.tsv"
    p.write_text("tx_A\t1\t300\t7\ntx_B\t0\t200\n")
    with pytest.raises(ValueError):
        rc._read_candidates(str(p))


def test_read_refseq_models_raises_on_short_row(tmp_path):
    p = tmp_path / "refseq.tsv"
    p.write_text("ref_1\tchr1\t100\t500\t+\t1\t310\n")   # 7 fields, need >= 8
    with pytest.raises(ValueError):
        rc._read_refseq_models(str(p))


def test_read_loci_raises_on_short_row(tmp_path):
    # A short loci row must raise; silently dropping it empties refseq_loci
    # and the chimeric flag would never be raised.
    p = tmp_path / "loci.tsv"
    p.write_text("chr1\t100\t500\n")                     # missing strand
    with pytest.raises(ValueError):
        rc._read_loci(str(p))


def test_read_candidates_wellformed_tolerates_extra_fields_and_comments(tmp_path):
    p = tmp_path / "txome.tsv"
    p.write_text("#query\tcomplete\tlength\tn_tm\n"       # header skipped
                 "tx_A\t1\t300\t7\n"
                 "tx_B\tpartial\t200\t6\textra\n")        # extra col tolerated
    cands = rc._read_candidates(str(p))
    assert [c.query for c in cands] == ["tx_A", "tx_B"]
    assert cands[0].complete is True and cands[0].n_tm == 7
    assert cands[1].complete is False and cands[1].length == 200


def test_read_candidates_is_crlf_safe(tmp_path):
    # Without rstrip("\r\n") the trailing \r rides into int("7\r") -> crash.
    p = tmp_path / "txome.tsv"
    p.write_bytes(b"tx_A\t1\t300\t7\r\n")
    cands = rc._read_candidates(str(p))
    assert cands[0].query == "tx_A" and cands[0].n_tm == 7
