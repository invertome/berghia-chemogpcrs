#!/bin/bash
# 02a_cluster_sequences.sh
# Purpose: Cluster identified GPCR candidates to remove redundant sequences (splice variants, alleles).
# Inputs: GPCR candidates from step 02
# Outputs: Non-redundant candidate set with cluster membership mapping
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=cluster_sequences
#SBATCH --output=${LOGS_DIR}/02a_cluster_sequences_%j.out
#SBATCH --error=${LOGS_DIR}/02a_cluster_sequences_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Initialize pipeline (provenance tracking)
init_pipeline

# Create output directories
CLUSTER_DIR="${RESULTS_DIR}/clustering"
mkdir -p "${CLUSTER_DIR}" "${LOGS_DIR}" || { log --level=ERROR "Cannot create directories"; exit 1; }

# Check for --force flag
FORCE_RERUN=false
if [ "$1" = "--force" ]; then
    FORCE_RERUN=true
    shift
fi

# Skip if already completed (unless --force)
if [ "$FORCE_RERUN" = false ]; then
    skip_if_completed "cluster_sequences"
fi

log "Starting sequence clustering (CD-HIT at ${CDHIT_IDENTITY} identity)"

# --- Check Dependencies ---
# Step 02 creates step_completed_02.txt
check_file "${RESULTS_DIR}/step_completed_02.txt"

# Find input FASTA file
INPUT_FASTA=""
CANDIDATES_DIR="${RESULTS_DIR}/candidates"

# Try different possible input locations
for candidate_file in \
    "${CANDIDATES_DIR}/chemogpcr_candidates.fa" \
    "${CANDIDATES_DIR}/chemogpcr_candidates.fasta" \
    "${CANDIDATES_DIR}/all_candidates.fa" \
    "${RESULTS_DIR}/chemogpcrs/candidates.fa"; do
    if [ -f "$candidate_file" ]; then
        INPUT_FASTA="$candidate_file"
        break
    fi
done

if [ -z "$INPUT_FASTA" ]; then
    log --level=ERROR "No candidate FASTA file found in ${CANDIDATES_DIR}"
    exit 1
fi

check_file "$INPUT_FASTA"

# Count input sequences
INPUT_COUNT=$(grep -c "^>" "$INPUT_FASTA")
log "Input sequences: ${INPUT_COUNT}"

if [ "$INPUT_COUNT" -eq 0 ]; then
    log --level=ERROR "No sequences in input file"
    exit 1
fi

# --- Detect Resources ---
detect_resources
THREADS=${DETECTED_CPUS:-8}

# --- Run CD-HIT ---
OUTPUT_FASTA="${CLUSTER_DIR}/candidates_nr${CDHIT_IDENTITY/./}.fa"
CLUSTER_FILE="${OUTPUT_FASTA}.clstr"

log "Running CD-HIT with ${THREADS} threads"
log "  Identity threshold: ${CDHIT_IDENTITY}"
log "  Word size: ${CDHIT_WORDSIZE}"
log "  Memory limit: ${CDHIT_MEMORY}MB"

run_command "cdhit_cluster" \
    --inputs="${INPUT_FASTA}" \
    --outputs="${OUTPUT_FASTA},${CLUSTER_FILE}" \
    ${CDHIT} \
    -i "$INPUT_FASTA" \
    -o "$OUTPUT_FASTA" \
    -c ${CDHIT_IDENTITY} \
    -n ${CDHIT_WORDSIZE} \
    -M ${CDHIT_MEMORY} \
    -d 0 \
    -T ${THREADS}

# --- Parse Cluster File and Create Mapping ---
log "Parsing cluster results"

python3 << 'PYTHON_SCRIPT'
import sys
import os
import re
from collections import defaultdict

# Get paths from environment
cluster_dir = os.environ.get('CLUSTER_DIR', 'results/clustering')
cluster_file = f"{cluster_dir}/candidates_nr{os.environ.get('CDHIT_IDENTITY', '0.98').replace('.', '')}.fa.clstr"
output_mapping = f"{cluster_dir}/cluster_mapping.tsv"
output_summary = f"{cluster_dir}/cluster_summary.tsv"

# Parse CD-HIT cluster file
clusters = defaultdict(list)
current_cluster = None
representative = None

with open(cluster_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>Cluster'):
            current_cluster = int(line.split()[1])
            representative = None
        elif line:
            # Parse sequence entry
            # Format: 0	123aa, >seq_id... *  or  0	123aa, >seq_id... at 98.00%
            match = re.search(r'>([^\s.]+)', line)
            if match:
                seq_id = match.group(1)
                is_rep = line.endswith('*')
                identity = 100.0 if is_rep else float(re.search(r'at ([\d.]+)%', line).group(1)) if 'at' in line else 0

                clusters[current_cluster].append({
                    'id': seq_id,
                    'is_representative': is_rep,
                    'identity': identity
                })

                if is_rep:
                    representative = seq_id

# Write detailed mapping
with open(output_mapping, 'w') as f:
    f.write("cluster_id\trepresentative\tmember_id\tidentity\n")
    for cluster_id, members in sorted(clusters.items()):
        rep = next((m['id'] for m in members if m['is_representative']), members[0]['id'])
        for member in members:
            f.write(f"{cluster_id}\t{rep}\t{member['id']}\t{member['identity']:.2f}\n")

# Write summary
with open(output_summary, 'w') as f:
    f.write("cluster_id\trepresentative\tcluster_size\n")
    for cluster_id, members in sorted(clusters.items()):
        rep = next((m['id'] for m in members if m['is_representative']), members[0]['id'])
        f.write(f"{cluster_id}\t{rep}\t{len(members)}\n")

# Print statistics
total_seqs = sum(len(m) for m in clusters.values())
num_clusters = len(clusters)
singletons = sum(1 for m in clusters.values() if len(m) == 1)
max_cluster = max(len(m) for m in clusters.values()) if clusters else 0
avg_cluster = total_seqs / num_clusters if num_clusters > 0 else 0

print(f"\nCluster Statistics:")
print(f"  Total input sequences: {total_seqs}")
print(f"  Number of clusters: {num_clusters}")
print(f"  Reduction: {total_seqs - num_clusters} sequences ({(total_seqs - num_clusters) / total_seqs * 100:.1f}%)")
print(f"  Singletons: {singletons}")
print(f"  Largest cluster: {max_cluster}")
print(f"  Average cluster size: {avg_cluster:.2f}")
PYTHON_SCRIPT

# --- Count Output Sequences ---
OUTPUT_COUNT=$(grep -c "^>" "$OUTPUT_FASTA")
REDUCTION=$((INPUT_COUNT - OUTPUT_COUNT))
REDUCTION_PCT=$(echo "scale=1; $REDUCTION * 100 / $INPUT_COUNT" | bc)

log "Clustering complete:"
log "  Input sequences: ${INPUT_COUNT}"
log "  Output clusters: ${OUTPUT_COUNT}"
log "  Sequences removed: ${REDUCTION} (${REDUCTION_PCT}%)"

# --- Create Symlink for Downstream Steps ---
# This allows downstream scripts to use a consistent filename
CANONICAL_OUTPUT="${CANDIDATES_DIR}/candidates_clustered.fa"
ln -sf "${OUTPUT_FASTA}" "${CANONICAL_OUTPUT}"
log "Created symlink: ${CANONICAL_OUTPUT}"

# --- Save Cluster Statistics ---
cat > "${CLUSTER_DIR}/cluster_stats.json" <<EOF
{
    "input_file": "${INPUT_FASTA}",
    "output_file": "${OUTPUT_FASTA}",
    "identity_threshold": ${CDHIT_IDENTITY},
    "word_size": ${CDHIT_WORDSIZE},
    "input_sequences": ${INPUT_COUNT},
    "output_clusters": ${OUTPUT_COUNT},
    "sequences_removed": ${REDUCTION},
    "reduction_percent": ${REDUCTION_PCT}
}
EOF

# Mark step as completed
create_checkpoint "cluster_sequences"

# Also create legacy completion flag for downstream steps
touch "${RESULTS_DIR}/step_completed_02a.txt"

log "Sequence clustering completed successfully"
log "Output files:"
log "  Clustered FASTA: ${OUTPUT_FASTA}"
log "  Cluster mapping: ${CLUSTER_DIR}/cluster_mapping.tsv"
log "  Cluster summary: ${CLUSTER_DIR}/cluster_summary.tsv"
