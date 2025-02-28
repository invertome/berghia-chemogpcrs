#!/usr/bin/env python3
# prune_alignment.py
# Purpose: Prune sequences from an alignment based on length and gap thresholds.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from Bio import SeqIO

input_file = sys.argv[1]
output_file = sys.argv[2]
min_length = int(sys.argv[3])
max_gap_percent = float(sys.argv[4])

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for record in SeqIO.parse(infile, 'fasta'):
        seq_len = len(record.seq)
        gap_count = record.seq.count('-')
        gap_percent = (gap_count / seq_len) * 100 if seq_len > 0 else 0
        if seq_len >= min_length and gap_percent <= max_gap_percent:
            SeqIO.write(record, outfile, 'fasta')
