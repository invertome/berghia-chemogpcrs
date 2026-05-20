#!/bin/bash
# fetch_audino2024_chemoreceptors.sh — Pull the Audino et al. 2024 Sci Rep
# 14:31539 lophotrochozoan chemoGPCR alignment + treefile from Figshare.
# Identify the verified close-cluster of the 11 Anomia simplex tentacle
# chemoreceptor candidates by tree analysis (NOT by alignment membership).
# Strip alignment gaps. Reformat headers. Write to
# references/curated_chemoreceptors_lophotrochozoa/apgr/apgr_seqs_audino2024.fa
#
# Why tree-based filtering (verified vs assumed):
# Audino's alignment file (Asim_chemoGPCR_ali.fasta) contains 81 sequences
# spanning the full chemoGPCR clade definition (incl. distant references
# and an Nvec outgroup). "Closely clustered to Anomia" requires actually
# parsing the published tree (Asim_chemoGPCR.treefile) — alignment
# membership alone is insufficient.
#
# Per-Asim sister-clade analysis (LCA walk with clade-size cap ≤15)
# yields 23 close-cluster tips:
#   11 TRINITY_DN (Anomia simplex)        — chemoreceptor candidates
#    5 Lven (Lottia ventricosa)           — verified clade members
#    3 Smed (Schmidtea mediterranea)      — verified clade members
#    3 Pcan (Pomacea canaliculata)        — verified clade members
#    1 Cele (Caenorhabditis elegans)      — Cele_NP_502887.1_Rho, sister of
#                                           TRINITY_DN29133 per Audino's tree
#                                           (included despite ecdysozoan
#                                           scope because it IS the tightest
#                                           phylogenetic neighbor)
#
# Source: Figshare 10.6084/m9.figshare.26058316
#   - file 47120083 (Asim_chemoGPCR_ali.fasta, 73 KB)
#   - file 47120071 (Asim_chemoGPCR.treefile, 13 KB, NEXUS w/ FigTree annotations)
#
# Bead -k0g v2, 2026-05-19.

set -eo pipefail

OUT_FASTA="${OUT_FASTA:-references/curated_chemoreceptors_lophotrochozoa/apgr/apgr_seqs_audino2024.fa}"
TMP="$(mktemp -d /tmp/audino2024.XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

ALI_URL="https://ndownloader.figshare.com/files/47120083"
TREE_URL="https://ndownloader.figshare.com/files/47120071"

echo "Downloading Audino 2024 chemoGPCR alignment + treefile from Figshare..."
curl -sL --max-time 60 -o "$TMP/audino_ali.fasta" "$ALI_URL"
curl -sL --max-time 60 -o "$TMP/audino_tree.nex" "$TREE_URL"
if [ ! -s "$TMP/audino_ali.fasta" ] || [ ! -s "$TMP/audino_tree.nex" ]; then
    echo "ERROR: Audino download failed or returned empty file" >&2
    exit 1
fi

ali_count=$(grep -c '^>' "$TMP/audino_ali.fasta")
echo "Alignment: $ali_count sequences."

mkdir -p "$(dirname "$OUT_FASTA")"

# Tree analysis: identify per-Asim close-cluster (LCA walk, clade size ≤15)
# then filter the alignment to only those tips.
python3 - "$TMP/audino_tree.nex" "$TMP/close_cluster.txt" <<'PYEOF'
import sys, re
from ete3 import Tree

tree_file, out_file = sys.argv[1], sys.argv[2]
with open(tree_file) as f:
    nexus = f.read()
m = re.search(r"tree\s+\w+\s*=\s*(?:\[[^\]]*\])?\s*(.+?;)\s*$",
              nexus, re.IGNORECASE | re.DOTALL | re.MULTILINE)
if not m:
    m = re.search(r"(\(.+\));", nexus, re.DOTALL)
newick = m.group(1)
# strip FigTree comments [&...] and single quotes
newick = re.sub(r"\[&[^\]]*\]", "", newick).replace("'", "") + ";"
t = Tree(newick, format=1)
asim_tips = [n for n in t.get_leaf_names() if n.startswith("TRINITY_DN")]

close = set()
for asim in asim_tips:
    node = t.search_nodes(name=asim)[0]
    cur = node
    while cur.up:
        cur = cur.up
        leaves = cur.get_leaf_names()
        if 2 <= len(leaves) <= 10:
            close.update(leaves)
            break
        if len(leaves) > 10:
            for child in cur.children:
                child_leaves = child.get_leaf_names()
                if asim in child_leaves and len(child_leaves) <= 15:
                    close.update(child_leaves)
            break

with open(out_file, "w") as f:
    for n in sorted(close):
        f.write(n + "\n")
print(f"close-cluster size: {len(close)}", file=sys.stderr)
PYEOF

# Filter the alignment FASTA to only close-cluster tips + standardize headers.
python3 - "$TMP/audino_ali.fasta" "$TMP/close_cluster.txt" "$OUT_FASTA" <<'PYEOF'
import sys, re

src, cluster_file, dst = sys.argv[1], sys.argv[2], sys.argv[3]

# Alignment headers use Asim_DN36691_c1_g1_i1.p1 form;
# treefile uses TRINITY_DN36691_c1_g1_i1.p1 form (different prefix).
# Normalize both to a comparable token.
def norm(name: str) -> str:
    return re.sub(r"^(Asim_|TRINITY_)", "", name)

close = set()
with open(cluster_file) as f:
    for line in f:
        line = line.strip()
        if line:
            close.add(norm(line))

# species abbreviation -> (full organism, taxid)
SPECIES_MAP = {
    "Asim":    ("Anomia_simplex",          138454),
    "Acal":    ("Aplysia_californica",     6500),
    "Acal5A":  ("Aplysia_californica",     6500),
    "Acal4B":  ("Aplysia_californica",     6500),
    "Acal10C": ("Aplysia_californica",     6500),
    "Lven":    ("Lottia_ventricosa",       1396875),
    "Pcan":    ("Pomacea_canaliculata",    400727),
    "Lgig":    ("Lottia_gigantea",         225164),
    "Smed":    ("Schmidtea_mediterranea",  79327),
    "Cele":    ("Caenorhabditis_elegans",  6239),  # included via tree analysis
}

# Idempotent: skip already-present sequences in OUT_FASTA
existing = set()
try:
    with open(dst) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].split("__")
                if parts:
                    existing.add(parts[-1].strip())
except FileNotFoundError:
    pass

def process(src_path: str, dst_path: str) -> None:
    counts = {"kept": 0, "skipped_not_in_cluster": 0, "skipped_cached": 0}
    with open(src_path) as fin, open(dst_path, "a") as fout:
        header = None
        seq_lines: list[str] = []
        def flush(header, seq_lines):
            if header is None:
                return
            if norm(header) not in close:
                counts["skipped_not_in_cluster"] += 1
                return
            sp_prefix = header.split("_", 1)[0]
            organism_taxid = SPECIES_MAP.get(sp_prefix)
            if organism_taxid is None:
                return
            organism, taxid = organism_taxid
            token = "audino2024_" + re.sub(r"[^A-Za-z0-9._-]", "_", header)
            if token in existing:
                counts["skipped_cached"] += 1
                return
            seq = "".join(seq_lines).replace("-", "").replace(".", "").upper()
            if len(seq) < 50:
                return
            gene_name = re.sub(r"[^A-Za-z0-9._-]", "_", header)
            new_header = f"apgr__{gene_name}__{taxid}__{token}"
            fout.write(f">{new_header}\n{seq}\n")
            counts["kept"] += 1

        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                flush(header, seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        flush(header, seq_lines)
    print(f"kept={counts['kept']} "
          f"skipped_not_in_cluster={counts['skipped_not_in_cluster']} "
          f"skipped_cached={counts['skipped_cached']}", file=sys.stderr)

process(src, dst)
PYEOF

n_total=$(grep -c '^>' "$OUT_FASTA")
echo "$OUT_FASTA now has $n_total sequences total."
