"""Tests for scripts/build_refseq_model_table.py.

This module was extracted from the inline RefSeq GFF-walk heredoc in
``02c_genome_reconcile.sh`` (Task-6 review follow-up). It walks the Berghia
RefSeq GFF3 to build the reconcile module's two inputs:

* ``--refseq-models`` — per-candidate ``query chrom start end strand complete
  length n_tm`` (``#``-commented header, values ``complete``/``partial``), and
* ``--refseq-loci`` — ALL gene intervals ``chrom start end strand`` (headerless).

These tests pin the coordinate walk (both the CDS ``gene=`` route and the
``Parent=rna- -> mRNA -> gene`` chain route), partial detection, residue-count
length, ``n_tm`` lookup, loci = every gene (not just candidates), and — the
Task-6 ``da7142e`` guard — the FAIL-LOUD-on-unresolvable / never-silently-drop
contract. The join key is the VERSIONED ``XP_<ver>`` accession throughout.

conftest.py puts ``scripts/`` on ``sys.path``.
"""
from __future__ import annotations

import build_refseq_model_table as brmt


# --- Synthetic RefSeq-format fixtures --------------------------------------

def _gff_line(chrom, ftype, start, end, strand, attrs):
    """One RefSeq-style 9-column GFF3 row (score/frame filled with '.')."""
    return "\t".join([chrom, "RefSeq", ftype, str(start), str(end),
                      ".", strand, ".", attrs])


# LOC1 resolves via the CDS ``gene=`` attribute (route A); LOC2 via the
# ``Parent=rna- -> mRNA -> gene`` chain (route B, CDS has NO gene=); LOC3 is
# partial; LOC9 is a NON-candidate gene (its interval must still reach loci).
# ``gene_biotype=`` is a decoy attribute (must not be mistaken for ``gene=``);
# the ``region`` feature must not land in loci.
GFF_ROWS = [
    "##gff-version 3",
    _gff_line("chr1", "region", 1, 100000, "+", "ID=chr1;Dbxref=taxon:1"),
    _gff_line("chr1", "gene", 1000, 2000, "+",
              "ID=gene-LOC1;Name=LOC1;gene_biotype=protein_coding;gene=LOC1"),
    _gff_line("chr1", "mRNA", 1000, 2000, "+",
              "ID=rna-XM1;Parent=gene-LOC1;gene=LOC1"),
    _gff_line("chr1", "CDS", 1000, 2000, "+",
              "ID=cds-XP_001.1;Parent=rna-XM1;gene=LOC1;protein_id=XP_001.1"),
    _gff_line("chr1", "gene", 5000, 6000, "-",
              "ID=gene-LOC2;Name=LOC2;gene_biotype=protein_coding;gene=LOC2"),
    _gff_line("chr1", "mRNA", 5000, 6000, "-", "ID=rna-XM2;Parent=gene-LOC2"),
    _gff_line("chr1", "CDS", 5000, 6000, "-",
              "ID=cds-XP_002.1;Parent=rna-XM2;protein_id=XP_002.1"),
    _gff_line("chr2", "gene", 3000, 4000, "+",
              "ID=gene-LOC3;Name=LOC3;gene=LOC3;partial=true;start_range=.,3000"),
    _gff_line("chr2", "mRNA", 3000, 4000, "+",
              "ID=rna-XM3;Parent=gene-LOC3;gene=LOC3;partial=true"),
    _gff_line("chr2", "CDS", 3000, 4000, "+",
              "ID=cds-XP_003.1;Parent=rna-XM3;gene=LOC3;protein_id=XP_003.1"),
    _gff_line("chr2", "gene", 8000, 9000, "-",
              "ID=gene-LOC9;Name=LOC9;gene_biotype=protein_coding;gene=LOC9"),
    _gff_line("chr2", "mRNA", 8000, 9000, "-", "ID=rna-XM9;Parent=gene-LOC9"),
    _gff_line("chr2", "CDS", 8000, 9000, "-",
              "ID=cds-XP_009.1;Parent=rna-XM9;gene=LOC9;protein_id=XP_009.1"),
]


def _make_fixture(tmp_path, cand_ids=("XP_001.1", "XP_002.1", "XP_003.1"),
                  gff_rows=None):
    """Write gff / ids / proteins / prediction fixtures; return their paths."""
    gff = tmp_path / "berghia.gff3"
    gff.write_text("\n".join(gff_rows if gff_rows is not None else GFF_ROWS) + "\n")

    ids = tmp_path / "cand_ids.txt"
    ids.write_text("\n".join(cand_ids) + "\n")

    # XP_001.1 spans two sequence lines (30 residues) to pin multi-line
    # accumulation; the header first token is the VERSIONED accession.
    proteins = tmp_path / "cands.faa"
    proteins.write_text("".join([
        ">XP_001.1 alpha desc\n", "A" * 15 + "\n", "A" * 15 + "\n",
        ">XP_002.1 beta\n", "A" * 20 + "\n",
        ">XP_003.1\n", "A" * 12 + "\n",
    ]))

    # DeepTMHMM prediction: col1 = versioned id, col5 = TM-region count.
    pred = tmp_path / "prediction"
    pred.write_text("XP_001.1\tTM\t0.98\tfoo\t7\n"
                    "XP_002.1\tTM\t0.95\tfoo\t6\n"
                    "XP_003.1\tTM\t0.80\tfoo\t7\n")

    return {
        "gff": str(gff), "ids": str(ids), "proteins": str(proteins),
        "pred": str(pred),
        "models_out": str(tmp_path / "refseq_models.tsv"),
        "loci_out": str(tmp_path / "refseq_loci.tsv"),
    }


# --- attr / is_partial (attribute parsing) ---------------------------------

def test_attr_extracts_gene_not_gene_biotype_decoy():
    a = "ID=gene-LOC1;Name=LOC1;gene_biotype=protein_coding;gene=LOC1"
    assert brmt.attr(a, "gene") == "LOC1"           # NOT "protein_coding"
    assert brmt.attr(a, "gene_biotype") == "protein_coding"
    assert brmt.attr(a, "protein_id") is None


def test_is_partial_detects_incomplete_tokens():
    assert brmt.is_partial("ID=x;partial=true") is True
    assert brmt.is_partial("ID=x;start_range=.,100") is True
    assert brmt.is_partial("ID=x;end_range=200,.") is True
    assert brmt.is_partial("ID=x;Name=y") is False


# --- read_protein_lengths / read_n_tm --------------------------------------

def test_read_protein_lengths_accumulates_multiline_and_keeps_version(tmp_path):
    fx = _make_fixture(tmp_path)
    plen = brmt.read_protein_lengths(fx["proteins"])
    assert plen["XP_001.1"] == 30   # 15 + 15 across two sequence lines
    assert plen["XP_002.1"] == 20
    assert plen["XP_003.1"] == 12


def test_read_n_tm_reads_col5_and_skips_short_rows(tmp_path):
    p = tmp_path / "pred"
    p.write_text("XP_1.1\ta\tb\tc\t6.0\n"   # int(float("6.0")) == 6
                 "SHORT\tonly\ttwo\n"        # < 5 fields -> skipped
                 "XP_2.1\ta\tb\tc\t7\n")
    assert brmt.read_n_tm(str(p)) == {"XP_1.1": 6, "XP_2.1": 7}


# --- resolve (the coordinate walk, both routes) ----------------------------

def _resolve(tmp_path, pid, cand_ids, gff_rows=None):
    fx = _make_fixture(tmp_path, cand_ids=cand_ids, gff_rows=gff_rows)
    gbi, gbl, mp, _loci, cds = brmt.parse_gff(fx["gff"], set(cand_ids))
    return brmt.resolve(pid, gbi, gbl, mp, cds)


def test_resolve_route_a_via_cds_gene_attribute(tmp_path):
    # CDS carries gene=LOC1 -> genes_by_loc route.
    assert _resolve(tmp_path, "XP_001.1", ("XP_001.1",)) == \
        ("chr1", 1000, 2000, "+", False)


def test_resolve_route_b_via_parent_mrna_gene_chain(tmp_path):
    # CDS has no gene= -> Parent=rna-XM2 -> mRNA -> gene-LOC2 route.
    assert _resolve(tmp_path, "XP_002.1", ("XP_002.1",)) == \
        ("chr1", 5000, 6000, "-", False)


def test_resolve_marks_partial_from_gene_flag(tmp_path):
    assert _resolve(tmp_path, "XP_003.1", ("XP_003.1",)) == \
        ("chr2", 3000, 4000, "+", True)


def test_resolve_ors_in_cds_level_partial(tmp_path):
    # gene + mRNA complete, but the CDS carries start_range= -> partial.
    rows = [
        _gff_line("c", "gene", 10, 20, "+", "ID=gene-LOCX;gene=LOCX"),
        _gff_line("c", "mRNA", 10, 20, "+", "ID=rna-X;Parent=gene-LOCX"),
        _gff_line("c", "CDS", 10, 20, "+",
                  "ID=cds;Parent=rna-X;gene=LOCX;protein_id=XP_x.1;start_range=.,10"),
    ]
    r = _resolve(tmp_path, "XP_x.1", ("XP_x.1",), gff_rows=rows)
    assert r == ("c", 10, 20, "+", True)


def test_resolve_ors_in_mrna_level_partial(tmp_path):
    # gene + CDS complete, but the mRNA is partial -> partial (even on route A).
    rows = [
        _gff_line("c", "gene", 10, 20, "+", "ID=gene-LOCY;gene=LOCY"),
        _gff_line("c", "mRNA", 10, 20, "+", "ID=rna-Y;Parent=gene-LOCY;partial=true"),
        _gff_line("c", "CDS", 10, 20, "+",
                  "ID=cds;Parent=rna-Y;gene=LOCY;protein_id=XP_y.1"),
    ]
    r = _resolve(tmp_path, "XP_y.1", ("XP_y.1",), gff_rows=rows)
    assert r == ("c", 10, 20, "+", True)


def test_resolve_returns_none_when_no_cds(tmp_path):
    assert _resolve(tmp_path, "XP_404.1", ("XP_404.1",)) is None


def test_resolve_returns_none_when_parent_chain_broken(tmp_path):
    # CDS with no gene= and a Parent that has no mRNA feature -> unresolvable.
    rows = [
        _gff_line("c", "gene", 10, 20, "+", "ID=gene-LOCZ;gene=LOCZ"),
        _gff_line("c", "CDS", 10, 20, "+",
                  "ID=cds;Parent=rna-MISSING;protein_id=XP_z.1"),
    ]
    assert _resolve(tmp_path, "XP_z.1", ("XP_z.1",), gff_rows=rows) is None


def test_resolve_route_a_carries_even_when_parent_chain_broken(tmp_path):
    # Fault-isolates route A: the CDS carries gene=LOCA (route A) BUT its Parent
    # points at a MISSING mRNA (route B is dead, exactly as the broken-chain case
    # above). The candidate must STILL resolve via route A. If route A were
    # disabled and everything fell through to route B, this would return None —
    # so this test is the route-A counterpart of route B's isolation (XP_002.1,
    # which has no gene=).
    rows = [
        _gff_line("c", "gene", 10, 20, "+", "ID=gene-LOCA;gene=LOCA"),
        _gff_line("c", "CDS", 10, 20, "+",
                  "ID=cds;Parent=rna-MISSING;gene=LOCA;protein_id=XP_a.1"),
    ]
    assert _resolve(tmp_path, "XP_a.1", ("XP_a.1",), gff_rows=rows) == \
        ("c", 10, 20, "+", False)


# --- build_tables (byte-level output contract) -----------------------------

def test_build_tables_writes_models_tsv_field_for_field(tmp_path):
    fx = _make_fixture(tmp_path)
    n_written, all_loci, unresolved = brmt.build_tables(
        fx["ids"], fx["gff"], fx["proteins"], fx["pred"],
        fx["models_out"], fx["loci_out"])
    assert unresolved == []
    assert n_written == 3
    with open(fx["models_out"]) as fh:
        assert fh.read() == (
            "#query\tchrom\tstart\tend\tstrand\tcomplete\tlength\tn_tm\n"
            "XP_001.1\tchr1\t1000\t2000\t+\tcomplete\t30\t7\n"
            "XP_002.1\tchr1\t5000\t6000\t-\tcomplete\t20\t6\n"
            "XP_003.1\tchr2\t3000\t4000\t+\tpartial\t12\t7\n")


def test_build_tables_loci_lists_all_genes_including_non_candidate(tmp_path):
    fx = _make_fixture(tmp_path)
    brmt.build_tables(fx["ids"], fx["gff"], fx["proteins"], fx["pred"],
                      fx["models_out"], fx["loci_out"])
    # LOC9 is NOT a candidate but its interval MUST be present (chimeric flag);
    # the region feature must NOT appear. Order follows the GFF gene order.
    with open(fx["loci_out"]) as fh:
        assert fh.read() == (
            "chr1\t1000\t2000\t+\n"
            "chr1\t5000\t6000\t-\n"
            "chr2\t3000\t4000\t+\n"
            "chr2\t8000\t9000\t-\n")


# --- main (CLI + the fail-loud never-drop guard) ---------------------------

def _argv(fx):
    return ["--ids", fx["ids"], "--gff", fx["gff"], "--proteins", fx["proteins"],
            "--prediction", fx["pred"], "--models-out", fx["models_out"],
            "--loci-out", fx["loci_out"]]


def test_main_returns_zero_and_reports_counts(tmp_path, capsys):
    fx = _make_fixture(tmp_path)
    assert brmt.main(_argv(fx)) == 0
    assert "refseq models: 3 written; refseq loci: 4" in capsys.readouterr().err


def test_main_fails_loud_and_nonzero_on_unresolvable_candidate(tmp_path, capsys):
    # XP_404.1 is a candidate with NO CDS: it must NOT vanish silently — the
    # stage must exit nonzero and name it (the da7142e never-drop guard).
    fx = _make_fixture(
        tmp_path, cand_ids=("XP_001.1", "XP_002.1", "XP_003.1", "XP_404.1"))
    rc = brmt.main(_argv(fx))
    assert rc != 0
    err = capsys.readouterr().err
    assert "XP_404.1" in err and "did not resolve" in err
    # The resolvable candidates are still written (loud failure, not a silent
    # loss of the whole table); the unresolvable one is absent-but-reported.
    with open(fx["models_out"]) as fh:
        models = fh.read()
    assert all(x in models for x in ("XP_001.1", "XP_002.1", "XP_003.1"))
    assert "XP_404.1" not in models
