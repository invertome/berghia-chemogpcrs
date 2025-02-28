#!/usr/bin/env python3
# update_headers.py
# Purpose: Update FASTA headers with short IDs and generate an ID mapping file.
# Inputs: FASTA file ($1), output ID map CSV ($2)
# Outputs: Updated FASTA file (${fasta_file}_updated.fa), ID map CSV
# Author: Jorge L. Perez-Moreno, Ph.D.

import sys

fasta_file = sys.argv[1]
id_map_file = sys.argv[2]
output_file = f"{fasta_file}_updated.fa"

with open(fasta_file, 'r') as f, open(output_file, 'w') as out, open(id_map_file, 'w') as map_out:
    counter = 1
    map_out.write("original_id,short_id\n")
    for line in f:
        if line.startswith('>'):
            original_id = line.strip()[1:]
            short_id = f"ref_{counter}"
            map_out.write(f"{original_id},{short_id}\n")
            out.write(f">{short_id}\n")
            counter += 1
        else:
            out.write(line)
