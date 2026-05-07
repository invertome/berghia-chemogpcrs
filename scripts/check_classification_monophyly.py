"""Family/subfamily monophyly check on the classification reference trees.

Tree leaves carry metadata as pipe-delimited fields:
    accession|family|subfamily|species
We parse the leaf name directly — no TSV lookup needed.
"""
from collections import defaultdict
from ete3 import Tree


def parse_leaf(name: str) -> tuple[str, str, str, str]:
    parts = name.split("|")
    while len(parts) < 4:
        parts.append("")
    return parts[0], parts[1], parts[2], parts[3]


def check_tree_monophyly(tree_path: str, label_idx: int) -> dict:
    """label_idx: 1 for family, 2 for subfamily."""
    t = Tree(tree_path, format=1)  # IQ-TREE 3 internal node labels
    by_label = defaultdict(list)
    for leaf in t.get_leaves():
        fields = parse_leaf(leaf.name)
        lbl = fields[label_idx]
        if not lbl:
            continue
        by_label[lbl].append(leaf.name)
    results = {}
    for label, leaves in by_label.items():
        if len(leaves) < 2:
            results[label] = (None, len(leaves))
            continue
        is_mono, _, _ = t.check_monophyly(values=leaves, target_attr="name")
        results[label] = (is_mono, len(leaves))
    return results


def report(tree_name, results, kind):
    print(f"\n=== {tree_name} ({kind} monophyly) ===")
    n_mono = sum(1 for v, _ in results.values() if v is True)
    n_poly = sum(1 for v, _ in results.values() if v is False)
    n_single = sum(1 for v, _ in results.values() if v is None)
    print(f"  monophyletic: {n_mono}, polyphyletic: {n_poly}, "
          f"singletons: {n_single}")
    for label in sorted(results.keys()):
        status, n = results[label]
        if status is None:
            mark = "  ·"
        elif status:
            mark = "  ✓"
        else:
            mark = "  ✗"
        print(f"{mark} {label:<35} (n={n})")


if __name__ == "__main__":
    base = ("/scratch3/workspace/jperezmoreno_umass_edu-jorge/"
            "chemogpcrs_2026-05/results/classification/trees")
    bb = check_tree_monophyly(f"{base}/backbone.treefile", label_idx=1)
    report("backbone", bb, "family")
    am = check_tree_monophyly(f"{base}/aminergic_subtree.treefile", label_idx=2)
    report("aminergic_subtree", am, "subfamily")
    pe = check_tree_monophyly(f"{base}/peptide_subtree.treefile", label_idx=2)
    report("peptide_subtree", pe, "subfamily")
