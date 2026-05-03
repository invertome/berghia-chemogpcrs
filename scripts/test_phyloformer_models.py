#!/usr/bin/env python3
# test_phyloformer_models.py
# Purpose: Select the best model for Phyloformer using IQ-TREE and infer a tree.
# Inputs: FASTA file ($1), output prefix ($2), threads ($3)
# Outputs: Phyloformer tree (${output_prefix}.tre)
# Author: Jorge L. Perez-Moreno, Ph.D.

import subprocess
import sys

fasta_file = sys.argv[1]
output_prefix = sys.argv[2]
threads = sys.argv[3]

models = ['LG', 'WAG', 'JTT']
# Bead -ryr / user direction: use iqtree3 (config.sh's IQTREE var) — both
# 'iqtree2' and 'iqtree3' alias to IQ-TREE 3 in the active env, but invoke
# the canonical name explicitly. Allow override via $IQTREE.
import os
iqtree_bin = os.environ.get("IQTREE", "iqtree3")
cmd = [iqtree_bin, "-s", fasta_file, "-m", "TESTONLY", "-mset", ",".join(models), "-nt", threads, "-pre", f"{output_prefix}_modeltest"]
subprocess.run(cmd, check=True)

with open(f"{output_prefix}_modeltest.log", 'r') as f:
    for line in f:
        if "Best-fit model according to" in line:
            best_model = line.split()[3]
            break

cmd = ["phyloformer", "-i", fasta_file, "-o", f"{output_prefix}.tre", "--model", best_model, "--threads", threads]
subprocess.run(cmd, check=True)
