"""Unit tests for scripts/compute_tandem_clusters.py (bead -ar8)."""
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from compute_tandem_clusters import (
    find_tandem_clusters,
    load_candidate_ids,
)


class TestFindTandemClusters:
    def test_basic_three_gene_cluster(self):
        # Three candidates within 50kb on scaffold S1
        genes = [
            ("g1", "S1", 1000, 2000),
            ("g2", "S1", 30000, 31000),
            ("g3", "S1", 60000, 61000),
        ]
        out = find_tandem_clusters(genes, window_kb=100, min_size=3)
        assert all(out[g][0] == 3 for g in ["g1", "g2", "g3"])
        assert all(out[g][1] is not None for g in ["g1", "g2", "g3"])
        assert len({out[g][1] for g in ["g1", "g2", "g3"]}) == 1  # same cluster id

    def test_isolated_gene_size_1(self):
        genes = [("solo", "S1", 1000, 2000)]
        out = find_tandem_clusters(genes, window_kb=100, min_size=3)
        assert out["solo"] == (1, None)

    def test_two_genes_min_size_3_remain_isolated(self):
        genes = [("a", "S1", 1000, 2000), ("b", "S1", 50000, 51000)]
        out = find_tandem_clusters(genes, window_kb=100, min_size=3)
        assert out["a"] == (1, None)
        assert out["b"] == (1, None)

    def test_distance_too_far_breaks_cluster(self):
        # Gap > window between g2 and g3 -> g3 isolated
        genes = [
            ("g1", "S1", 1000, 2000),
            ("g2", "S1", 50000, 51000),
            ("g3", "S1", 200000, 201000),
        ]
        out = find_tandem_clusters(genes, window_kb=100, min_size=3)
        # Only 2 genes in window for the first attempt -> all should be size 1
        assert out["g1"] == (1, None)
        assert out["g2"] == (1, None)
        assert out["g3"] == (1, None)

    def test_two_clusters_on_same_scaffold(self):
        genes = [
            ("a1", "S1", 1000, 2000),
            ("a2", "S1", 50000, 51000),
            ("a3", "S1", 80000, 81000),
            # Gap of 500kb separates the second cluster
            ("b1", "S1", 600000, 601000),
            ("b2", "S1", 650000, 651000),
            ("b3", "S1", 680000, 681000),
        ]
        out = find_tandem_clusters(genes, window_kb=100, min_size=3)
        a_ids = {out[g][1] for g in ["a1", "a2", "a3"]}
        b_ids = {out[g][1] for g in ["b1", "b2", "b3"]}
        assert len(a_ids) == 1 and len(b_ids) == 1
        assert a_ids != b_ids
        assert all(out[g][0] == 3 for g in ["a1", "a2", "a3", "b1", "b2", "b3"])

    def test_cluster_size_six(self):
        genes = [(f"g{i}", "S1", i * 20000, i * 20000 + 1000) for i in range(6)]
        out = find_tandem_clusters(genes, window_kb=100, min_size=3)
        for g in [f"g{i}" for i in range(6)]:
            assert out[g][0] == 6

    def test_different_scaffolds_independent(self):
        genes = [
            ("a", "S1", 1000, 2000),
            ("b", "S1", 30000, 31000),
            ("c", "S2", 1000, 2000),
            ("d", "S2", 30000, 31000),
        ]
        # Min size 2 to test
        out = find_tandem_clusters(genes, window_kb=100, min_size=2)
        s1_ids = {out[g][1] for g in ["a", "b"]}
        s2_ids = {out[g][1] for g in ["c", "d"]}
        assert len(s1_ids) == 1 and len(s2_ids) == 1
        assert s1_ids != s2_ids

    def test_empty_input(self):
        assert find_tandem_clusters([]) == {}

    def test_min_size_param_respected(self):
        genes = [("g1", "S1", 1000, 2000), ("g2", "S1", 30000, 31000)]
        # min_size=2 should cluster these two
        out_2 = find_tandem_clusters(genes, window_kb=100, min_size=2)
        assert out_2["g1"][0] == 2
        # min_size=3 should leave them isolated
        out_3 = find_tandem_clusters(genes, window_kb=100, min_size=3)
        assert out_3["g1"][0] == 1


class TestLoadCandidateIds:
    def test_fasta_format(self, tmp_path):
        f = tmp_path / "cands.fa"
        f.write_text(">gene_1 descr\nACGT\n>gene_2\nACGT\n")
        ids = load_candidate_ids(str(f))
        assert ids == {"gene_1", "gene_2"}

    def test_plain_text(self, tmp_path):
        f = tmp_path / "cands.txt"
        f.write_text("gene_a\ngene_b\n\ngene_c\n")
        ids = load_candidate_ids(str(f))
        assert ids == {"gene_a", "gene_b", "gene_c"}

    def test_fasta_with_sequence_lines_skips_sequences(self, tmp_path):
        # A file starting with '>' is treated as FASTA. Sequence lines are
        # NOT IDs. This is the canonical pipeline-produced candidates file.
        f = tmp_path / "cands.fa"
        f.write_text(">gene_x\nMACGT\n>gene_y\nMNNNNNNN\n")
        ids = load_candidate_ids(str(f))
        assert ids == {"gene_x", "gene_y"}

    def test_plain_text_tolerates_stray_gt(self, tmp_path):
        # Plain text (no leading >) tolerates a > prefix on individual lines.
        f = tmp_path / "cands.txt"
        f.write_text("gene_a\n>gene_b\ngene_c\n")
        ids = load_candidate_ids(str(f))
        assert ids == {"gene_a", "gene_b", "gene_c"}


class TestEndToEndWithRealGFF:
    """Test against a small synthetic GFF using gffutils."""

    def test_gff_to_clusters(self, tmp_path):
        gff_text = textwrap.dedent("""\
            ##gff-version 3
            S1\ttest\tgene\t1000\t2000\t.\t+\t.\tID=g1
            S1\ttest\tgene\t30000\t31000\t.\t+\t.\tID=g2
            S1\ttest\tgene\t60000\t61000\t.\t+\t.\tID=g3
            S1\ttest\tgene\t500000\t501000\t.\t+\t.\tID=g4
            S2\ttest\tgene\t1000\t2000\t.\t+\t.\tID=g5
        """)
        gff = tmp_path / "test.gff3"
        gff.write_text(gff_text)
        cand = tmp_path / "cand.txt"
        cand.write_text("g1\ng2\ng3\ng4\ng5\n")

        from compute_tandem_clusters import build_or_load_gff_db, iter_genes

        db = build_or_load_gff_db(str(gff), str(tmp_path / "test.gff3.db"))
        cand_set = load_candidate_ids(str(cand))
        candidate_genes = [g for g in iter_genes(db) if g[0] in cand_set]
        out = find_tandem_clusters(candidate_genes, window_kb=100, min_size=3)
        # g1, g2, g3 should be in a cluster of size 3; g4 isolated; g5 isolated
        assert out["g1"][0] == 3
        assert out["g2"][0] == 3
        assert out["g3"][0] == 3
        assert out["g4"][0] == 1
        assert out["g5"][0] == 1
