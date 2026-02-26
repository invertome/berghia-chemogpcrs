#!/bin/bash
# 02b_interproscan_classification.sh
# Purpose: Classify GPCR candidates by Pfam domain (7tm_1 = Class A rhodopsin)
# Optional step — set RUN_INTERPROSCAN=true in config.sh to enable
# Inputs: GPCR FASTA from step 02
# Outputs: gpcr_classes.csv with gene_id, gpcr_class, evidence columns
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=interproscan_classification
#SBATCH --output=${LOGS_DIR}/02b_interproscan_%j.out
#SBATCH --error=${LOGS_DIR}/02b_interproscan_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

source config.sh
source functions.sh

if [ "${RUN_INTERPROSCAN}" != "true" ]; then
    log "InterProScan classification skipped (RUN_INTERPROSCAN=${RUN_INTERPROSCAN})"
    # Create a default classification file (all Class A)
    mkdir -p "${RESULTS_DIR}/interproscan"
    CANDIDATES_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
    if [ -f "$CANDIDATES_FA" ]; then
        python3 -c "
from Bio import SeqIO
with open('${RESULTS_DIR}/interproscan/gpcr_classes.csv', 'w') as f:
    f.write('gene_id,gpcr_class,evidence\n')
    for rec in SeqIO.parse('${CANDIDATES_FA}', 'fasta'):
        f.write(f'{rec.id},Class_A,no_interproscan\n')
"
    fi
    touch "${RESULTS_DIR}/step_completed_02b.txt"
    exit 0
fi

mkdir -p "${RESULTS_DIR}/interproscan" "${LOGS_DIR}"
check_file "${RESULTS_DIR}/step_completed_02.txt"

CANDIDATES_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
check_file "$CANDIDATES_FA"

log "Running InterProScan classification..."

run_command "interproscan_berghia" ${INTERPROSCAN} \
    -i "$CANDIDATES_FA" \
    -o "${RESULTS_DIR}/interproscan/interproscan_results.tsv" \
    -f tsv \
    -appl Pfam \
    --goterms \
    -cpu "${CPUS}"

# Parse: extract Pfam domain hits and classify
python3 -c "
import sys

# Pfam domains for GPCR classification
CLASS_DOMAINS = {
    '7tm_1': 'Class_A',        # Rhodopsin-like
    '7tm_2': 'Class_B1',       # Secretin-like
    '7tm_3': 'Class_C',        # Glutamate
    '7TM_GPCR_Srv': 'Class_A',
    '7TM_GPCR_Srsx': 'Class_A',
    '7TM_GPCR_Srw': 'Class_A',
    'Frizzled': 'Class_F',
    'GPS': 'Adhesion',
}

# Parse InterProScan TSV (standard 11-column format)
gene_classes = {}
try:
    with open('${RESULTS_DIR}/interproscan/interproscan_results.tsv') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                gene_id = parts[0]
                pfam_id = parts[4] if len(parts) > 4 else ''
                pfam_desc = parts[5] if len(parts) > 5 else ''
                for domain, gpcr_class in CLASS_DOMAINS.items():
                    if domain in pfam_id or domain in pfam_desc:
                        gene_classes[gene_id] = gpcr_class
                        break
except FileNotFoundError:
    print('Warning: InterProScan results file not found', file=sys.stderr)

# Write classification
from Bio import SeqIO
with open('${RESULTS_DIR}/interproscan/gpcr_classes.csv', 'w') as f:
    f.write('gene_id,gpcr_class,evidence\n')
    for rec in SeqIO.parse('${CANDIDATES_FA}', 'fasta'):
        cls = gene_classes.get(rec.id, 'unknown')
        evidence = 'interproscan' if rec.id in gene_classes else 'unclassified'
        f.write(f'{rec.id},{cls},{evidence}\n')

classified = len(gene_classes)
total = sum(1 for _ in SeqIO.parse('${CANDIDATES_FA}', 'fasta'))
print(f'Classified {classified}/{total} candidates by Pfam domain', file=sys.stderr)
"

touch "${RESULTS_DIR}/step_completed_02b.txt"
log "InterProScan classification completed."
