#!/usr/bin/env python3
# prune_alignment.py
# Purpose: Prune sequences from an alignment based on length and gap thresholds after trimming.
# Inputs: Input alignment FASTA ($1), output FASTA ($2), min length ($3), max gap percent ($4)
# Outputs: Pruned alignment FASTA
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

import sys
from Bio import SeqIO

input_file = sys.argv[1]       # Input alignment FASTA file
output_file = sys.argv[2]      # Output pruned FASTA file
min_length = int(sys.argv[3])  # Minimum sequence length threshold
max_gap_percent = float(sys.argv[4])  # Maximum gap percentage threshold

# Open input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Iterate through sequences in the alignment
    for record in SeqIO.parse(infile, 'fasta'):
        seq_len = len(record.seq)
        gap_count = record.seq.count('-')
        gap_percent = (gap_count / seq_len) * 100 if seq_len > 0 else 0
        # Write sequence if it meets criteria
        if seq_len >= min_length and gap_percent <= max_gap_percent:
            SeqIO.write(record, outfile, 'fasta')
