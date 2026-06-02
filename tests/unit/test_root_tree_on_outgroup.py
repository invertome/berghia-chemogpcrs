"""Tests for scripts/root_tree_on_outgroup.py.

Downstream rooting of the per-class trees (bead vo8.3): IQ-TREE leaves them
unrooted; this roots on the swap-map sister-class outgroup taxa, falling back
to midpoint when the outgroup is absent or non-monophyletic.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

pytest.importorskip("ete3")
import root_tree_on_outgroup as rt


# ---------------------------------------------------------------------------
# outgroup id parsing
# ---------------------------------------------------------------------------

class TestParseOutgroupIds:
    def test_plain_id_file(self, tmp_path: Path) -> None:
        f = tmp_path / "ids.txt"
        f.write_text("O1\nO2\n\n")          # blank line ignored
        assert rt.parse_outgroup_ids(f) == ["O1", "O2"]

    def test_fasta_headers(self, tmp_path: Path) -> None:
        f = tmp_path / "og.fa"
        f.write_text(">O1 desc here\nACGT\n>O2\nACGT\n")
        # FASTA header → id is first whitespace-delimited token after '>'
        assert rt.parse_outgroup_ids(f) == ["O1", "O2"]


# ---------------------------------------------------------------------------
# rooting logic
# ---------------------------------------------------------------------------

class TestRootTree:
    def test_single_outgroup(self) -> None:
        from ete3 import Tree
        t = Tree("((A1:1,A2:1):1,O1:1);", format=1)
        mode = rt.root_tree(t, ["O1"])
        assert mode == "outgroup_single"
        # O1 is a direct child of the (new) root
        assert "O1" in [c.name for c in t.children if c.is_leaf()]

    def test_monophyletic_outgroup_clade(self) -> None:
        from ete3 import Tree
        t = Tree("((A1:1,A2:1):1,(O1:1,O2:1):1);", format=1)
        mode = rt.root_tree(t, ["O1", "O2"])
        assert mode == "outgroup_clade"
        # one side of the root is exactly {O1,O2}
        sides = [set(c.get_leaf_names()) for c in t.children]
        assert {"O1", "O2"} in sides

    def test_non_monophyletic_outgroup_falls_back_to_midpoint(self) -> None:
        from ete3 import Tree
        # O1 and O2 are on opposite sides → not monophyletic
        t = Tree("((A1:1,O1:1):1,(A2:1,O2:1):1);", format=1)
        mode = rt.root_tree(t, ["O1", "O2"], fallback="midpoint")
        assert mode == "midpoint"

    def test_absent_outgroup_falls_back_to_midpoint(self) -> None:
        from ete3 import Tree
        t = Tree("((A1:1,A2:1):1,(A3:1,A4:1):1);", format=1)
        mode = rt.root_tree(t, ["O1", "O2"], fallback="midpoint")
        assert mode == "midpoint"

    def test_absent_outgroup_fallback_none_leaves_unrooted(self) -> None:
        from ete3 import Tree
        t = Tree("((A1:1,A2:1):1,(A3:1,A4:1):1);", format=1)
        mode = rt.root_tree(t, ["O1"], fallback="none")
        assert mode == "unrooted"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

class TestCLI:
    def test_writes_rooted_treefile(self, tmp_path: Path) -> None:
        from ete3 import Tree
        tree_f = tmp_path / "class_A.treefile"
        tree_f.write_text("((A1:1,A2:1):1,(O1:1,O2:1):1);")
        og_f = tmp_path / "outgroup.fa"
        og_f.write_text(">O1\nACGT\n>O2\nACGT\n")
        out = tmp_path / "class_A.rooted.treefile"
        rc = rt.main(["--tree", str(tree_f), "--outgroup-ids", str(og_f),
                      "--out", str(out)])
        assert rc == 0
        assert out.exists()
        # output is valid newick and rooted on the outgroup clade
        rooted = Tree(str(out), format=1)
        sides = [set(c.get_leaf_names()) for c in rooted.children]
        assert {"O1", "O2"} in sides

    def test_missing_tree_returns_nonzero(self, tmp_path: Path) -> None:
        out = tmp_path / "out.treefile"
        rc = rt.main(["--tree", str(tmp_path / "nope.treefile"),
                      "--outgroup-ids", str(tmp_path / "og.fa"),
                      "--out", str(out)])
        assert rc != 0
        assert not out.exists()
