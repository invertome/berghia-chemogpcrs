"""Unit tests for scripts/parse_jcvi_anchors.py (bead -e59)."""
import pytest

from parse_jcvi_anchors import parse_anchors, per_gene_summary


def write_anchors(tmp_path, blocks):
    """Helper: write a list-of-blocks to a synthetic .anchors file."""
    p = tmp_path / "x.anchors"
    lines = []
    for i, block in enumerate(blocks):
        lines.append(f"### block {i}")
        for a, b, e in block:
            lines.append(f"{a}\t{b}\t{e}")
    p.write_text("\n".join(lines) + "\n")
    return str(p)


class TestParseAnchors:
    def test_single_block(self, tmp_path):
        path = write_anchors(tmp_path, [[("g1", "h1", 1e-50), ("g2", "h2", 1e-40)]])
        blocks = parse_anchors(path)
        assert len(blocks) == 1
        assert blocks[0] == [("g1", "h1", 1e-50), ("g2", "h2", 1e-40)]

    def test_multiple_blocks(self, tmp_path):
        path = write_anchors(tmp_path, [
            [("g1", "h1", 1e-50)],
            [("g2", "h2", 1e-30), ("g3", "h3", 1e-20)],
        ])
        blocks = parse_anchors(path)
        assert len(blocks) == 2
        assert len(blocks[1]) == 2

    def test_empty_file(self, tmp_path):
        p = tmp_path / "empty.anchors"
        p.write_text("")
        assert parse_anchors(str(p)) == []


class TestPerGeneSummary:
    def test_basic_counts(self, tmp_path):
        path = write_anchors(tmp_path, [
            [("g1", "h1", 1e-50), ("g2", "h2", 1e-40)],   # block size 2
            [("g1", "h3", 1e-30), ("g3", "h4", 1e-20)],   # block size 2
            [("g4", "h5", 1e-10)],                          # block size 1
        ])
        blocks = parse_anchors(path)
        df = per_gene_summary(blocks, "test_target")
        # g1 appears in 2 blocks with sizes 2+2 = 4 total genes
        row_g1 = df[df["candidate_id"] == "g1"].iloc[0]
        assert row_g1["n_anchor_blocks"] == 2
        assert row_g1["total_anchor_genes"] == 4
        # g4 is only in one size-1 block
        row_g4 = df[df["candidate_id"] == "g4"].iloc[0]
        assert row_g4["n_anchor_blocks"] == 1
        assert row_g4["total_anchor_genes"] == 1
        # target_species column present everywhere
        assert (df["target_species"] == "test_target").all()

    def test_dedup_within_block(self, tmp_path):
        # Same gene-A appearing twice in a block should not double-count
        path = write_anchors(tmp_path, [
            [("g1", "h1", 1e-50), ("g1", "h2", 1e-40), ("g2", "h3", 1e-20)],
        ])
        blocks = parse_anchors(path)
        df = per_gene_summary(blocks, "tgt")
        row_g1 = df[df["candidate_id"] == "g1"].iloc[0]
        # n_blocks=1 (one block contains g1 once after dedup)
        assert row_g1["n_anchor_blocks"] == 1
        # total_anchor_genes still equals block size (3)
        assert row_g1["total_anchor_genes"] == 3

    def test_no_blocks(self, tmp_path):
        df = per_gene_summary([], "tgt")
        assert df.empty
        assert list(df.columns) == ["candidate_id", "n_anchor_blocks",
                                    "total_anchor_genes", "target_species"]
