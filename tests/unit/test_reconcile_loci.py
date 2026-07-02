"""Tests for scripts/reconcile_candidates.py — Task 2 (locus grouping).

Genome-track + candidate reconciliation feature (bead: genome-track
stage 02c). Task 2 covers only the ``PlacedModel`` / ``Locus`` value
records and ``group_into_loci`` — grouping placed candidate gene models
into loci (= genes) by same-strand genomic overlap and collapsing
isoforms. Representative selection, provenance labeling, the placement
cascade, and output writers belong to later tasks and are out of scope.

The overlap contract is half-open ``[start, end)``: two intervals share
a locus iff ``a.start < b.end and b.start < a.end`` AND they are on the
same chromosome and strand. Bookended intervals (touching at a boundary,
0 shared bases) never merge.

group_into_loci may shell out to ``bedtools merge`` when the binary is
present, but falls back to a pure-Python interval sweep otherwise. This
machine has no bedtools, so every test below (except the explicitly
monkeypatched parity test) exercises the pure-Python path.
"""
from __future__ import annotations

from types import SimpleNamespace

import reconcile_candidates as rc


# ---- PlacedModel dataclass -----------------------------------------------

def test_placed_model_positional_query_chrom_coords_strand_source_kwarg():
    pm = rc.PlacedModel("t1", "chr1", 100, 500, "+", source="txome")
    assert pm.query == "t1"
    assert pm.chrom == "chr1"
    assert pm.start == 100 and pm.end == 500
    assert pm.strand == "+"
    assert pm.source == "txome"


def test_placed_model_source_is_free_form_string():
    # downstream uses "genome" and "txome"; the field itself is free-form.
    assert rc.PlacedModel("g1", "chr1", 1, 2, "+", source="genome").source == "genome"
    assert rc.PlacedModel("x", "chr1", 1, 2, "+", source="anything").source == "anything"


# ---- Locus dataclass ------------------------------------------------------

def test_locus_fields_and_member_list():
    pm1 = rc.PlacedModel("a", "chr1", 100, 400, "+", source="genome")
    pm2 = rc.PlacedModel("b", "chr1", 300, 600, "+", source="txome")
    locus = rc.Locus("chr1", 100, 600, "+", [pm1, pm2])
    assert locus.chrom == "chr1"
    assert locus.start == 100 and locus.end == 600
    assert locus.strand == "+"
    assert locus.members == [pm1, pm2]


# ---- group_into_loci: degenerate inputs -----------------------------------

def test_group_empty_returns_empty_list():
    assert rc.group_into_loci([]) == []


def test_group_single_record_is_one_locus_one_member():
    pm = rc.PlacedModel("solo", "chr1", 100, 500, "+", source="genome")
    loci = rc.group_into_loci([pm])
    assert len(loci) == 1
    (locus,) = loci
    assert locus.members == [pm]
    assert locus.chrom == "chr1" and locus.strand == "+"
    assert locus.start == 100 and locus.end == 500


# ---- group_into_loci: isoform collapse ------------------------------------

def test_two_overlapping_same_strand_isoforms_collapse_to_one_locus():
    a = rc.PlacedModel("iso_a", "chr1", 100, 500, "+", source="genome")
    b = rc.PlacedModel("iso_b", "chr1", 400, 600, "+", source="txome")   # overlaps a
    loci = rc.group_into_loci([a, b])
    assert len(loci) == 1
    assert len(loci[0].members) == 2


def test_isoform_collapse_spans_min_start_max_end():
    a = rc.PlacedModel("iso_a", "chr1", 100, 500, "+", source="genome")
    b = rc.PlacedModel("iso_b", "chr1", 400, 600, "+", source="txome")
    (locus,) = rc.group_into_loci([a, b])
    assert locus.start == 100   # min(member starts)
    assert locus.end == 600     # max(member ends)


def test_isoform_collapse_preserves_both_source_labels():
    # a genome model + a transcriptome model at the same locus -> both
    # sources survive on the locus's members (Task 3 reads this to label
    # provenance "both"). Grouping must not drop or overwrite source.
    g = rc.PlacedModel("gmodel", "chr1", 100, 500, "+", source="genome")
    t = rc.PlacedModel("tmodel", "chr1", 300, 700, "+", source="txome")
    (locus,) = rc.group_into_loci([g, t])
    assert {m.source for m in locus.members} == {"genome", "txome"}


def test_transitive_overlap_chain_forms_single_locus():
    # A-B overlap, B-C overlap, A and C are disjoint: the bridging model B
    # ties all three into one connected-component locus.
    a = rc.PlacedModel("A", "chr1", 100, 300, "+", source="genome")
    b = rc.PlacedModel("B", "chr1", 250, 550, "+", source="genome")
    c = rc.PlacedModel("C", "chr1", 500, 700, "+", source="genome")
    loci = rc.group_into_loci([a, b, c])
    assert len(loci) == 1
    assert len(loci[0].members) == 3
    assert loci[0].start == 100 and loci[0].end == 700


# ---- group_into_loci: things that must NOT merge --------------------------

def test_disjoint_same_strand_paralogs_are_two_loci():
    a = rc.PlacedModel("a", "chr1", 100, 500, "+", source="genome")
    b = rc.PlacedModel("b", "chr1", 600, 900, "+", source="genome")
    loci = rc.group_into_loci([a, b])
    assert len(loci) == 2
    assert all(len(l.members) == 1 for l in loci)


def test_overlapping_opposite_strand_are_two_loci_ordered_plus_then_minus():
    # same chrom, same coordinates, opposite strand -> separate loci.
    # also pins the (chrom, start, strand) sort: "+" sorts before "-".
    m = rc.PlacedModel("m", "chr1", 100, 500, "-", source="txome")
    p = rc.PlacedModel("p", "chr1", 100, 500, "+", source="genome")
    loci = rc.group_into_loci([m, p])   # deliberately minus-first input
    assert len(loci) == 2
    assert loci[0].strand == "+" and loci[0].members[0].query == "p"
    assert loci[1].strand == "-" and loci[1].members[0].query == "m"


def test_same_coords_different_chrom_are_two_loci():
    a = rc.PlacedModel("a", "chr1", 100, 500, "+", source="genome")
    b = rc.PlacedModel("b", "chr2", 100, 500, "+", source="genome")
    assert len(rc.group_into_loci([a, b])) == 2


def test_bookended_same_strand_do_not_merge_gap_guard():
    # [100,500) and [500,900) touch at 500 with 0 shared bases (half-open):
    # they must stay separate. This is the gap guard.
    a = rc.PlacedModel("a", "chr1", 100, 500, "+", source="genome")
    b = rc.PlacedModel("b", "chr1", 500, 900, "+", source="txome")
    loci = rc.group_into_loci([a, b])
    assert len(loci) == 2


# ---- group_into_loci: determinism -----------------------------------------

def test_loci_sorted_by_chrom_start_strand():
    # unsorted input across two chromosomes / strands -> deterministic order.
    z = rc.PlacedModel("z", "chr2", 100, 200, "+", source="genome")
    a = rc.PlacedModel("a", "chr1", 500, 900, "-", source="txome")
    m = rc.PlacedModel("m", "chr1", 100, 400, "+", source="genome")
    loci = rc.group_into_loci([z, a, m])
    ordered = [(l.chrom, l.start, l.strand) for l in loci]
    assert ordered == [("chr1", 100, "+"), ("chr1", 500, "-"), ("chr2", 100, "+")]


def test_locus_members_sorted_by_query():
    # queries supplied out of lexical order collapse into one locus whose
    # members come back sorted by query.
    t9 = rc.PlacedModel("t9", "chr1", 100, 500, "+", source="txome")
    t1 = rc.PlacedModel("t1", "chr1", 200, 600, "+", source="genome")
    (locus,) = rc.group_into_loci([t9, t1])
    assert [m.query for m in locus.members] == ["t1", "t9"]


# ---- bedtools path parity (monkeypatched; no binary needed) ---------------

def test_bedtools_path_matches_python_path_and_uses_stranded_overlap_flags(monkeypatch):
    """When bedtools IS available the shell-out path must return loci
    identical to the pure-Python sweep, and must invoke `bedtools merge`
    with the same-strand (`-s`) + require-overlap (`-d -1`) flags that
    forbid book-ended merges. We fake both `shutil.which` and
    `subprocess.run` so no real binary is needed: the fake returns the
    collapse output bedtools would emit for these three records
    (indices 0,1 overlap; index 2 is a disjoint paralog).
    """
    records = [
        rc.PlacedModel("iso_a", "chr1", 100, 500, "+", source="genome"),
        rc.PlacedModel("iso_b", "chr1", 400, 600, "+", source="txome"),
        rc.PlacedModel("para", "chr1", 900, 1200, "+", source="txome"),
    ]

    # ground truth from the pure-Python path
    monkeypatch.setattr(rc.shutil, "which", lambda _binary: None)
    python_loci = rc.group_into_loci(records)

    # now force the bedtools branch with a canned merge -c 4 -o collapse
    # output (collapsed index list is the LAST tab field, per parser).
    recorded = {}

    def fake_run(cmd, *args, **kwargs):
        recorded["argv"] = cmd
        recorded["input"] = kwargs.get("input")
        out = "chr1\t100\t600\t+\t0,1\nchr1\t900\t1200\t+\t2\n"
        return SimpleNamespace(stdout=out, stderr="", returncode=0)

    monkeypatch.setattr(rc.shutil, "which", lambda _binary: "/usr/bin/bedtools")
    monkeypatch.setattr(rc.subprocess, "run", fake_run)
    bedtools_loci = rc.group_into_loci(records)

    assert bedtools_loci == python_loci
    argv = recorded["argv"]
    assert "merge" in argv and "-s" in argv
    assert "-d" in argv and "-1" in argv          # forbids book-ended merges


def test_clusters_from_merge_output_maps_collapsed_indices_back_to_records():
    # unit-level guard on the pure parse-back helper the bedtools branch
    # relies on: the last tab field is a comma-separated list of record
    # indices; each output line becomes one cluster.
    records = [
        rc.PlacedModel("r0", "chr1", 100, 500, "+", source="genome"),
        rc.PlacedModel("r1", "chr1", 400, 600, "+", source="txome"),
        rc.PlacedModel("r2", "chr2", 10, 90, "-", source="genome"),
    ]
    out = "chr1\t100\t600\t+\t0,1\nchr2\t10\t90\t-\t2\n"
    clusters = rc._clusters_from_merge_output(out, records)
    assert [sorted(m.query for m in c) for c in clusters] == [["r0", "r1"], ["r2"]]
