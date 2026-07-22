"""Tests for scripts/reconcile_candidates.py — Task 1.

Genome-track + candidate reconciliation feature (bead: genome-track
stage 02c). Task 1 covers only the placement model, the four aligner
parsers, and the %match gate / best-vs-second-margin functions. Locus
grouping, representative selection, the placement cascade, and output
writers are built in later tasks and are out of scope here.

Fixtures use minimal-but-realistic aligner output built directly in
each test (no golden files) so the parsing contract is visible at the
call site.
"""
from __future__ import annotations

import reconcile_candidates as rc


# ---- Placement dataclass -------------------------------------------------

def test_placement_fields_positional():
    p = rc.Placement("q", "chr1", 1, 2, "+", 96.0, 95.0, "minimap2", 100.0)
    assert p.query == "q"
    assert p.chrom == "chr1"
    assert p.start == 1 and p.end == 2
    assert p.strand == "+"
    assert p.pct_identity == 96.0
    assert p.pct_coverage == 95.0
    assert p.method == "minimap2"
    assert p.score == 100.0


# ---- parse_minimap2_paf ---------------------------------------------------

def test_paf_identity_and_coverage(tmp_path):
    # minimap2 PAF: qname qlen qstart qend strand tname tlen tstart tend nmatch
    # alen mapq  tags...  Identity now comes from the intron-aware cs tag (NOT
    # nmatch/alen): cs :294*ac -> 294 matches + 1 substitution.
    paf = tmp_path / "m.paf"
    paf.write_text("BersteEVm1\t300\t0\t300\t+\tchr1\t10000\t500\t800\t294\t300\t60\t"
                   "tp:A:P\tcs:Z::294*ac\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert p.chrom == "chr1" and p.strand == "+"
    assert abs(p.pct_identity - 100.0 * 294 / 295) < 1e-6  # 294/(294 match + 1 sub)
    assert p.pct_coverage == 100.0         # 100 * (300-0)/300


def test_paf_coverage_uses_qlen_denominator(tmp_path):
    """Coverage uses qlen (f[1]) as its denominator, independent of the
    identity source. With qlen != (qend-qstart) the coverage is a proper
    fraction; identity comes from the de:f divergence tag (intron-aware), so a
    coverage regression cannot masquerade as an identity change and vice versa.
    Correct: identity = 100*(1-0.02) = 98.0 ; coverage = 100*300/305 = 98.36..
    """
    paf = tmp_path / "m.paf"
    # qname qlen=305 qstart=0 qend=300 strand tname tlen tstart tend nmatch alen mapq de:f
    paf.write_text("q1\t305\t0\t300\t+\tchr1\t10000\t500\t800\t294\t300\t60\t"
                   "tp:A:P\tde:f:0.02\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert abs(p.pct_identity - 98.0) < 1e-6                   # 100 * (1 - 0.02)
    assert abs(p.pct_coverage - 98.36065573770492) < 1e-9      # 100 * 300/305


def test_paf_query_and_score_and_method(tmp_path):
    paf = tmp_path / "m.paf"
    paf.write_text("BersteEVm1\t300\t0\t300\t+\tchr1\t10000\t500\t800\t294\t300\t60\t"
                   "tp:A:P\tde:f:0.01\n")
    p = rc.parse_minimap2_paf(str(paf))[0]
    assert p.query == "BersteEVm1"
    assert p.method == "minimap2"
    assert p.score == 60.0
    assert p.start == 500 and p.end == 800


def test_paf_multiple_lines_multiple_placements(tmp_path):
    paf = tmp_path / "m.paf"
    paf.write_text(
        "q1\t300\t0\t300\t+\tchr1\t10000\t500\t800\t294\t300\t60\tde:f:0.02\n"
        "q2\t200\t0\t200\t-\tchr2\t5000\t100\t300\t190\t200\t40\tde:f:0.05\n"
    )
    placements = rc.parse_minimap2_paf(str(paf))
    assert len(placements) == 2
    assert {p.query for p in placements} == {"q1", "q2"}
    q2 = [p for p in placements if p.query == "q2"][0]
    assert q2.chrom == "chr2" and q2.strand == "-"


def test_paf_skips_short_lines(tmp_path):
    paf = tmp_path / "m.paf"
    paf.write_text(
        "too\tshort\n"
        "\n"
        "q1\t300\t0\t300\t+\tchr1\t10000\t500\t800\t294\t300\t60\tde:f:0.02\n"
    )
    placements = rc.parse_minimap2_paf(str(paf))
    assert len(placements) == 1
    assert placements[0].query == "q1"


def test_paf_skips_non_numeric_fields(tmp_path):
    paf = tmp_path / "m.paf"
    paf.write_text(
        "bad\tNOTANUMBER\t0\t300\t+\tchr1\t10000\t500\t800\t294\t300\t60\tde:f:0.02\n"
        "q1\t300\t0\t300\t+\tchr1\t10000\t500\t800\t294\t300\t60\tde:f:0.02\n"
    )
    placements = rc.parse_minimap2_paf(str(paf))
    assert len(placements) == 1
    assert placements[0].query == "q1"


def test_paf_zero_qlen_coverage_zero_no_crash(tmp_path):
    """A zero qlen must not crash the coverage division (-> 0.0). Identity
    comes from the divergence tag, so alen (no longer used for identity) is
    irrelevant here."""
    paf = tmp_path / "m.paf"
    paf.write_text("q1\t0\t0\t0\t+\tchr1\t10000\t500\t800\t0\t0\t60\tde:f:0.0\n")
    placements = rc.parse_minimap2_paf(str(paf))
    assert len(placements) == 1
    assert placements[0].pct_identity == 100.0   # 100 * (1 - 0.0)
    assert placements[0].pct_coverage == 0.0     # qlen 0 -> 0.0, no ZeroDivisionError


def test_paf_missing_file_returns_empty():
    assert rc.parse_minimap2_paf("/no/such/path.paf") == []


def test_paf_empty_file_returns_empty(tmp_path):
    paf = tmp_path / "empty.paf"
    paf.write_text("")
    assert rc.parse_minimap2_paf(str(paf)) == []


# ---- parse_gmap_gff --------------------------------------------------------

def _gmap_gff_line(chrom="chr1", start=500, end=800, strand="+",
                    name="BersteEVm1", coverage=99.5, identity=97.2,
                    feature="mRNA"):
    attrs = f"ID={name}.path1.mrna1;Name={name};Parent={name}.path1"
    if coverage is not None:
        attrs += f";coverage={coverage}"
    if identity is not None:
        attrs += f";identity={identity}"
    return f"{chrom}\tgmap\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


def test_gmap_gff_identity_and_coverage(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        + _gmap_gff_line()
    )
    p = rc.parse_gmap_gff(str(gff))[0]
    assert p.query == "BersteEVm1"
    assert p.chrom == "chr1" and p.start == 500 and p.end == 800
    assert p.strand == "+"
    assert p.pct_identity == 97.2
    assert p.pct_coverage == 99.5
    assert p.method == "gmap"


def test_gmap_gff_gene_row_without_attrs_is_skipped(tmp_path):
    """A typical gmap gff3_gene file has a 'gene' summary row with no
    coverage=/identity= attributes and an 'mRNA' row that carries them.
    Only the row(s) that actually carry identity/coverage should
    produce a Placement — this naturally selects the mRNA row without
    needing to special-case feature-type ordering."""
    gene_line = "chr1\tgmap\tgene\t500\t800\t.\t+\t.\tID=BersteEVm1.path1;Name=BersteEVm1\n"
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n" + gene_line + _gmap_gff_line())
    placements = rc.parse_gmap_gff(str(gff))
    assert len(placements) == 1
    assert placements[0].pct_identity == 97.2


def test_gmap_gff_exon_rows_ignored(tmp_path):
    exon_line = ("chr1\tgmap\texon\t500\t600\t.\t+\t.\t"
                 "ID=BersteEVm1.path1.mrna1.exon1;Parent=BersteEVm1.path1.mrna1\n")
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n" + _gmap_gff_line() + exon_line)
    placements = rc.parse_gmap_gff(str(gff))
    assert len(placements) == 1


def test_gmap_gff_multiple_mrna_rows(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        + _gmap_gff_line(chrom="chr1", name="q1", identity=97.2, coverage=99.5)
        + _gmap_gff_line(chrom="chr2", name="q2", strand="-", identity=88.0, coverage=90.0)
    )
    placements = rc.parse_gmap_gff(str(gff))
    assert len(placements) == 2
    q2 = [p for p in placements if p.query == "q2"][0]
    assert q2.chrom == "chr2" and q2.strand == "-" and q2.pct_identity == 88.0


def test_gmap_gff_malformed_identity_value_skipped(tmp_path):
    bad = "chr1\tgmap\tmRNA\t500\t800\t.\t+\t.\tID=x;Name=q1;coverage=99.5;identity=NOTANUMBER\n"
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n" + bad + _gmap_gff_line(name="q2"))
    placements = rc.parse_gmap_gff(str(gff))
    assert len(placements) == 1
    assert placements[0].query == "q2"


def test_gmap_gff_short_and_blank_lines_skipped(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text("##comment\n\ntoo\tshort\n" + _gmap_gff_line())
    placements = rc.parse_gmap_gff(str(gff))
    assert len(placements) == 1


def test_gmap_gff_score_from_matches_attribute(tmp_path):
    line = ("chr1\tgmap\tmRNA\t500\t800\t.\t+\t.\t"
            "ID=x;Name=q1;coverage=99.5;identity=97.2;matches=291;mismatches=8\n")
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n" + line)
    p = rc.parse_gmap_gff(str(gff))[0]
    assert p.score == 291.0


def test_gmap_gff_score_defaults_zero_without_matches_attribute(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n" + _gmap_gff_line())  # no matches= attr
    p = rc.parse_gmap_gff(str(gff))[0]
    assert p.score == 0.0


def test_gmap_gff_query_from_target_when_name_absent(tmp_path):
    line = ("chr1\tgmap\tmRNA\t500\t800\t.\t+\t.\t"
            "ID=x;Parent=y;coverage=99.5;identity=97.2;Target=BersteEVm7 1 300 +\n")
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n" + line)
    p = rc.parse_gmap_gff(str(gff))[0]
    assert p.query == "BersteEVm7"


def test_gmap_gff_no_name_or_target_is_skipped(tmp_path):
    """Fail closed on the join key: an mRNA row that carries identity=/
    coverage= but has neither Name= nor Target= is SKIPPED rather than
    falling back to ID= (which is `<id>.pathN.mrnaN`, not the bare
    transcript id — a silently-wrong key that would break Task 4
    concordance). Losing a placement degrades gracefully (the cascade
    falls through to miniprot); a wrong key does not."""
    line = ("chr1\tgmap\tmRNA\t500\t800\t.\t+\t.\t"
            "ID=BersteEVm1.path1.mrna1;Parent=BersteEVm1.path1;"
            "coverage=99.5;identity=97.2\n")
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n" + line)
    assert rc.parse_gmap_gff(str(gff)) == []


def test_gmap_gff_missing_file_returns_empty():
    assert rc.parse_gmap_gff("/no/such/path.gff3") == []


def test_gmap_gff_empty_file_returns_empty(tmp_path):
    gff = tmp_path / "empty.gff3"
    gff.write_text("")
    assert rc.parse_gmap_gff(str(gff)) == []


# ---- parse_miniprot_gff -----------------------------------------------------

def _miniprot_gff_line(chrom="chr1", start=500, end=800, strand="+",
                        target="BersteEVm1 1 300", score="1450",
                        identity=0.972, extra_attrs=""):
    attrs = f"ID=MP00001;Rank=1;Identity={identity};Positive=0.985;Target={target}"
    if extra_attrs:
        attrs += ";" + extra_attrs
    return f"{chrom}\tminiprot\tmRNA\t{start}\t{end}\t{score}\t{strand}\t.\t{attrs}\n"


def test_miniprot_gff_identity_fraction_converted_to_percent(tmp_path):
    gff = tmp_path / "mp.gff3"
    gff.write_text("##gff-version 3\n" + _miniprot_gff_line(identity=0.972))
    p = rc.parse_miniprot_gff(str(gff))[0]
    assert p.query == "BersteEVm1"
    assert p.chrom == "chr1" and p.start == 500 and p.end == 800
    assert p.strand == "+"
    assert abs(p.pct_identity - 97.2) < 1e-6
    assert p.method == "miniprot"
    assert p.score == 1450.0


def test_miniprot_gff_missing_coverage_uses_documented_sentinel(tmp_path):
    """miniprot's native GFF3 mRNA row has no coverage-equivalent
    attribute (the Target range alone doesn't reveal full protein
    length), so pct_coverage falls back to a documented sentinel that
    always fails a positive-min_cov gate (conservative reject ->
    cascade falls through to the next placement method)."""
    gff = tmp_path / "mp.gff3"
    gff.write_text("##gff-version 3\n" + _miniprot_gff_line())
    p = rc.parse_miniprot_gff(str(gff))[0]
    assert p.pct_coverage == rc.MINIPROT_COVERAGE_SENTINEL
    assert rc.pass_gate(p, min_id=90.0, min_cov=1.0) is False


def test_miniprot_gff_explicit_coverage_attr_used_when_present(tmp_path):
    gff = tmp_path / "mp.gff3"
    gff.write_text("##gff-version 3\n" + _miniprot_gff_line(extra_attrs="Coverage=95.0"))
    p = rc.parse_miniprot_gff(str(gff))[0]
    assert p.pct_coverage == 95.0


def test_miniprot_gff_multiple_rows(tmp_path):
    gff = tmp_path / "mp.gff3"
    gff.write_text(
        "##gff-version 3\n"
        + _miniprot_gff_line(target="q1 1 300", identity=0.972)
        + _miniprot_gff_line(chrom="chr2", strand="-", target="q2 1 200", identity=0.80)
    )
    placements = rc.parse_miniprot_gff(str(gff))
    assert len(placements) == 2
    q2 = [p for p in placements if p.query == "q2"][0]
    assert q2.chrom == "chr2" and q2.strand == "-"
    assert abs(q2.pct_identity - 80.0) < 1e-6


def test_miniprot_gff_non_mrna_rows_ignored(tmp_path):
    cds_line = "chr1\tminiprot\tCDS\t500\t600\t.\t+\t0\tParent=MP00001\n"
    gff = tmp_path / "mp.gff3"
    gff.write_text("##gff-version 3\n" + _miniprot_gff_line() + cds_line)
    placements = rc.parse_miniprot_gff(str(gff))
    assert len(placements) == 1


def test_miniprot_gff_malformed_identity_skipped(tmp_path):
    bad = "chr1\tminiprot\tmRNA\t500\t800\t1450\t+\t.\tID=x;Target=q1 1 300;Identity=NOTANUMBER\n"
    gff = tmp_path / "mp.gff3"
    gff.write_text("##gff-version 3\n" + bad + _miniprot_gff_line(target="q2 1 300"))
    placements = rc.parse_miniprot_gff(str(gff))
    assert len(placements) == 1
    assert placements[0].query == "q2"


def test_miniprot_gff_missing_file_returns_empty():
    assert rc.parse_miniprot_gff("/no/such/path.gff3") == []


def test_miniprot_gff_empty_file_returns_empty(tmp_path):
    gff = tmp_path / "empty.gff3"
    gff.write_text("")
    assert rc.parse_miniprot_gff(str(gff)) == []


# ---- parse_blastp_tab -------------------------------------------------------

def _blast_row(qseqid="q1", sseqid="XP_000001.1", pident="96.500",
               length="300", mismatch="10", gapopen="0",
               qstart="1", qend="300", sstart="1", send="300",
               evalue="1e-150", bitscore="580", qcovs="98"):
    return "\t".join([qseqid, sseqid, pident, length, mismatch, gapopen,
                      qstart, qend, sstart, send, evalue, bitscore, qcovs]) + "\n"


def test_blastp_tab_basic_with_coords(tmp_path):
    tab = tmp_path / "b.tsv"
    tab.write_text(_blast_row())
    coords = {"XP_000001.1": ("chr3", 1000, 2000, "+")}
    p = rc.parse_blastp_tab(str(tab), coords)[0]
    assert p.query == "q1"
    assert p.chrom == "chr3" and p.start == 1000 and p.end == 2000
    assert p.strand == "+"
    assert p.pct_identity == 96.5
    assert p.pct_coverage == 98.0
    assert p.score == 580.0
    assert p.method == "blastp"


def test_blastp_tab_subject_absent_from_coords_skipped(tmp_path):
    tab = tmp_path / "b.tsv"
    tab.write_text(_blast_row(sseqid="XP_999999.1"))
    coords = {"XP_000001.1": ("chr3", 1000, 2000, "+")}
    placements = rc.parse_blastp_tab(str(tab), coords)
    assert placements == []


def test_blastp_tab_multiple_rows_mixed_coords(tmp_path):
    tab = tmp_path / "b.tsv"
    tab.write_text(
        _blast_row(qseqid="q1", sseqid="XP_000001.1")
        + _blast_row(qseqid="q2", sseqid="XP_999999.1")   # no coords -> skipped
        + _blast_row(qseqid="q3", sseqid="XP_000002.1", pident="88.0", qcovs="70")
    )
    coords = {
        "XP_000001.1": ("chr3", 1000, 2000, "+"),
        "XP_000002.1": ("chr4", 5000, 6000, "-"),
    }
    placements = rc.parse_blastp_tab(str(tab), coords)
    assert len(placements) == 2
    assert {p.query for p in placements} == {"q1", "q3"}
    q3 = [p for p in placements if p.query == "q3"][0]
    assert q3.chrom == "chr4" and q3.strand == "-"
    assert q3.pct_identity == 88.0 and q3.pct_coverage == 70.0


def test_blastp_tab_skips_short_or_garbage_lines(tmp_path):
    tab = tmp_path / "b.tsv"
    tab.write_text("too\tfew\tcolumns\n\n" + _blast_row())
    coords = {"XP_000001.1": ("chr3", 1000, 2000, "+")}
    placements = rc.parse_blastp_tab(str(tab), coords)
    assert len(placements) == 1


def test_blastp_tab_non_numeric_fields_skipped(tmp_path):
    tab = tmp_path / "b.tsv"
    tab.write_text(
        _blast_row(qseqid="bad", pident="NOTANUMBER")
        + _blast_row(qseqid="q1")
    )
    coords = {"XP_000001.1": ("chr3", 1000, 2000, "+")}
    placements = rc.parse_blastp_tab(str(tab), coords)
    assert len(placements) == 1
    assert placements[0].query == "q1"


def test_blastp_tab_empty_coords_dict_skips_everything(tmp_path):
    tab = tmp_path / "b.tsv"
    tab.write_text(_blast_row())
    assert rc.parse_blastp_tab(str(tab), {}) == []


def test_blastp_tab_missing_file_returns_empty():
    coords = {"XP_000001.1": ("chr3", 1000, 2000, "+")}
    assert rc.parse_blastp_tab("/no/such/path.tsv", coords) == []


def test_blastp_tab_empty_file_returns_empty(tmp_path):
    tab = tmp_path / "empty.tsv"
    tab.write_text("")
    assert rc.parse_blastp_tab(str(tab), {}) == []


# ---- pass_gate ---------------------------------------------------------------

def test_gate_absolute_and_margin():
    a = rc.Placement("q","chr1",1,2,"+",96.0,95.0,"minimap2",100)
    b = rc.Placement("q","chr2",1,2,"+",84.0,95.0,"minimap2",90)
    assert rc.pass_gate(a, 95.0, 90.0) is True
    assert rc.pass_gate(b, 95.0, 90.0) is False       # below id bar
    best, margin = rc.best_and_margin([a, b])
    assert best is a and abs(margin - 12.0) < 1e-6


def test_gate_below_coverage_bar_fails():
    a = rc.Placement("q", "chr1", 1, 2, "+", 99.0, 50.0, "minimap2", 100)
    assert rc.pass_gate(a, 95.0, 90.0) is False


def test_gate_exact_boundary_passes():
    a = rc.Placement("q", "chr1", 1, 2, "+", 95.0, 90.0, "minimap2", 100)
    assert rc.pass_gate(a, 95.0, 90.0) is True


def test_gate_none_placement_returns_false():
    assert rc.pass_gate(None, 95.0, 90.0) is False


# ---- best_and_margin -----------------------------------------------------------

def test_single_placement_margin_inf():
    a = rc.Placement("q","chr1",1,2,"+",99.0,99.0,"gmap",100)
    _, margin = rc.best_and_margin([a])
    assert margin == float("inf")


def test_best_and_margin_empty_list():
    best, margin = rc.best_and_margin([])
    assert best is None
    assert margin == 0.0


def test_best_and_margin_filters_none_entries():
    a = rc.Placement("q", "chr1", 1, 2, "+", 96.0, 95.0, "minimap2", 100)
    best, margin = rc.best_and_margin([a, None])
    assert best is a
    assert margin == float("inf")


def test_best_and_margin_three_way_picks_top_two():
    a = rc.Placement("q", "chr1", 1, 2, "+", 96.0, 95.0, "minimap2", 100)
    b = rc.Placement("q", "chr2", 1, 2, "+", 90.0, 95.0, "minimap2", 90)
    c = rc.Placement("q", "chr3", 1, 2, "+", 70.0, 95.0, "minimap2", 80)
    best, margin = rc.best_and_margin([c, a, b])
    assert best is a
    assert abs(margin - 6.0) < 1e-6
