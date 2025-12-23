#!/bin/bash
# 02b_classify_gpcrs.sh
# Purpose: Classify identified GPCRs using InterProScan and sub-family HMMs.
# Inputs: GPCR candidates from clustering step
# Outputs: Classified GPCRs with family/sub-family annotations
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=classify_gpcrs
#SBATCH --output=${LOGS_DIR}/02b_classify_gpcrs_%j.out
#SBATCH --error=${LOGS_DIR}/02b_classify_gpcrs_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Initialize pipeline
init_pipeline

# Create output directories
CLASSIFY_DIR="${RESULTS_DIR}/classification"
mkdir -p "${CLASSIFY_DIR}" "${LOGS_DIR}" || { log --level=ERROR "Cannot create directories"; exit 1; }

# Check for --force flag
FORCE_RERUN=false
if [ "$1" = "--force" ]; then
    FORCE_RERUN=true
    shift
fi

if [ "$FORCE_RERUN" = false ]; then
    skip_if_completed "classify_gpcrs"
fi

log "Starting GPCR classification"

# --- Check Dependencies ---
# Step 02a creates step_completed_02a.txt (or checkpoint)
check_file "${RESULTS_DIR}/step_completed_02a.txt"

# Find input FASTA file (clustered candidates)
INPUT_FASTA=""
for candidate_file in \
    "${RESULTS_DIR}/candidates/candidates_clustered.fa" \
    "${RESULTS_DIR}/clustering/candidates_nr098.fa" \
    "${RESULTS_DIR}/candidates/chemogpcr_candidates.fa"; do
    if [ -f "$candidate_file" ]; then
        INPUT_FASTA="$candidate_file"
        break
    fi
done

if [ -z "$INPUT_FASTA" ]; then
    log --level=ERROR "No candidate FASTA file found"
    exit 1
fi

check_file "$INPUT_FASTA"
INPUT_COUNT=$(grep -c "^>" "$INPUT_FASTA")
log "Input sequences: ${INPUT_COUNT}"

# --- Detect Resources ---
detect_resources
THREADS=${DETECTED_CPUS:-16}

# --- Step 1: Run InterProScan (if available) ---
log "Step 1: Running InterProScan domain annotation"

INTERPRO_OUTPUT="${CLASSIFY_DIR}/interproscan_results.tsv"

if command -v interproscan.sh &> /dev/null; then
    log "Running InterProScan with ${THREADS} threads..."

    run_command "interproscan" \
        --inputs="$INPUT_FASTA" \
        --outputs="$INTERPRO_OUTPUT" \
        --allow-fail \
        interproscan.sh \
        -i "$INPUT_FASTA" \
        -f TSV \
        -o "$INTERPRO_OUTPUT" \
        -appl Pfam,PANTHER,CDD,SUPERFAMILY \
        -goterms \
        -pa \
        -cpu ${THREADS}

    if [ -f "$INTERPRO_OUTPUT" ]; then
        log "InterProScan completed successfully"
    else
        log --level=WARN "InterProScan did not produce output"
    fi
else
    log --level=WARN "InterProScan not available, skipping domain annotation"
    log "  Install with: conda install -c bioconda interproscan"
    touch "$INTERPRO_OUTPUT"  # Create empty file
fi

# --- Step 2: Run GPCR Sub-family Classification ---
log "Step 2: Running GPCR sub-family classification"

python3 << 'PYTHON_SCRIPT'
import os
import sys
import json
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Get paths
classify_dir = os.environ.get('CLASSIFY_DIR', 'results/classification')
input_fasta = os.environ.get('INPUT_FASTA', '')
interpro_output = f"{classify_dir}/interproscan_results.tsv"
local_db_dir = os.environ.get('LOCAL_DB_DIR', '')

# GPCR domain patterns for classification
GPCR_DOMAINS = {
    'Class_A': {
        'pfam': ['PF00001', 'PF10320', 'PF10321', 'PF10322', 'PF10323', 'PF10324'],
        'interpro': ['IPR000276', 'IPR017452', 'IPR019422'],
        'keywords': ['7tm_1', 'rhodopsin', 'gpcr_a']
    },
    'Class_B1': {
        'pfam': ['PF00002'],
        'interpro': ['IPR000832', 'IPR001879'],
        'keywords': ['7tm_2', 'secretin']
    },
    'Class_B2': {
        'pfam': ['PF12430'],
        'interpro': ['IPR001879'],
        'keywords': ['adhesion', 'gpcr_b2']
    },
    'Class_C': {
        'pfam': ['PF00003', 'PF01094'],
        'interpro': ['IPR000337', 'IPR001828'],
        'keywords': ['7tm_3', 'metabotropic', 'glutamate']
    },
    'Frizzled': {
        'pfam': ['PF01534', 'PF01392'],
        'interpro': ['IPR000539', 'IPR020067'],
        'keywords': ['frizzled', 'smoothened', 'fz']
    },
    'Taste2': {
        'pfam': ['PF05296'],
        'interpro': ['IPR007960'],
        'keywords': ['taste', 'tas2r', 't2r']
    }
}

# Chemoreceptor-specific sub-families for Class A
CLASS_A_SUBFAMILIES = {
    'Olfactory': ['olfr', 'or_', 'odorant', 'olfactory'],
    'Gustatory': ['taste', 'gust', 'tas1r', 'sweet', 'umami'],
    'Vomeronasal': ['vomero', 'v1r', 'v2r', 'vmn'],
    'Trace_amine': ['taar', 'trace_amine', 'tar_'],
    'Formyl_peptide': ['fpr', 'formyl'],
    'Chemokine': ['cxcr', 'ccr', 'chemokine'],
    'Neuropeptide': ['npr', 'npy', 'neuropeptide'],
    'Opsin': ['opsin', 'rhodopsin', 'rh1', 'rh2'],
    'Aminergic': ['drd', 'htr', '5ht', 'adra', 'adrb', 'chrm', 'hrh'],
    'Peptide': ['sstr', 'tacr', 'ednr', 'avpr', 'oxtr'],
    'Lipid': ['s1pr', 'lpar', 'ptger', 'ptgfr'],
    'Other_chemosensory': ['chemoreceptor', 'chem', 'srx', 'srw', 'str']
}


def parse_fasta_headers(fasta_file):
    """Extract sequence IDs and descriptions from FASTA file."""
    sequences = {}
    current_id = None

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split(None, 1)
                seq_id = parts[0]
                description = parts[1] if len(parts) > 1 else ''
                sequences[seq_id] = {'description': description, 'domains': [], 'class': None, 'subfamily': None}
                current_id = seq_id

    return sequences


def parse_interpro_results(interpro_file, sequences):
    """Parse InterProScan TSV output and annotate sequences."""
    if not os.path.exists(interpro_file) or os.path.getsize(interpro_file) == 0:
        return sequences

    try:
        # InterProScan TSV format: seq_id, md5, length, analysis, accession, description, start, end, score, status, date, interpro_acc, interpro_desc, go_terms, pathways
        with open(interpro_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                seq_id = parts[0]
                analysis = parts[3] if len(parts) > 3 else ''
                accession = parts[4] if len(parts) > 4 else ''
                description = parts[5] if len(parts) > 5 else ''
                interpro_acc = parts[11] if len(parts) > 11 else ''

                if seq_id in sequences:
                    sequences[seq_id]['domains'].append({
                        'analysis': analysis,
                        'accession': accession,
                        'description': description,
                        'interpro': interpro_acc
                    })
    except Exception as e:
        print(f"Warning: Error parsing InterPro results: {e}", file=sys.stderr)

    return sequences


def classify_gpcr_class(sequence_info):
    """Classify GPCR into major class based on domains."""
    domains = sequence_info['domains']
    description = sequence_info['description'].lower()

    # Collect all domain accessions
    pfam_ids = set()
    interpro_ids = set()
    keywords = set()

    for domain in domains:
        if domain['accession']:
            if domain['accession'].startswith('PF'):
                pfam_ids.add(domain['accession'])
            keywords.add(domain['accession'].lower())
        if domain['interpro']:
            interpro_ids.add(domain['interpro'])
        if domain['description']:
            keywords.update(domain['description'].lower().split())

    # Add description keywords
    keywords.update(description.split())

    # Check each class
    scores = {}
    for gpcr_class, patterns in GPCR_DOMAINS.items():
        score = 0
        # Check Pfam
        score += len(pfam_ids & set(patterns['pfam'])) * 3
        # Check InterPro
        score += len(interpro_ids & set(patterns['interpro'])) * 2
        # Check keywords
        for kw in patterns['keywords']:
            if any(kw in k for k in keywords):
                score += 1
        scores[gpcr_class] = score

    # Return highest scoring class (default to Class_A for GPCRs)
    if max(scores.values()) > 0:
        return max(scores, key=scores.get)
    return 'Class_A'  # Default for unclassified GPCRs


def classify_subfamily(sequence_info, gpcr_class):
    """Classify Class A GPCRs into sub-families."""
    if gpcr_class != 'Class_A':
        return None

    description = sequence_info['description'].lower()
    seq_id = sequence_info.get('id', '').lower()

    # Collect keywords from domains
    keywords = set(description.split())
    keywords.update(seq_id.split('_'))

    for domain in sequence_info['domains']:
        if domain['description']:
            keywords.update(domain['description'].lower().split())

    # Check each subfamily
    for subfamily, patterns in CLASS_A_SUBFAMILIES.items():
        for pattern in patterns:
            if any(pattern in k for k in keywords):
                return subfamily

    return 'Unclassified_Class_A'


def load_gpcrdb_classifications(local_db_dir):
    """Load reference classifications from local GPCRdb cache."""
    classifications = {}

    if not local_db_dir:
        return classifications

    gpcrdb_file = Path(local_db_dir) / "gpcrdb" / "receptor_list.json"
    if gpcrdb_file.exists():
        try:
            with open(gpcrdb_file) as f:
                data = json.load(f)
                for receptor in data:
                    entry_name = receptor.get('entry_name', '')
                    family = receptor.get('family', '')
                    classifications[entry_name] = family
            print(f"Loaded {len(classifications)} GPCRdb classifications", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Could not load GPCRdb classifications: {e}", file=sys.stderr)

    return classifications


# Main execution
print("=== GPCR Classification ===", file=sys.stderr)

# Parse sequences
sequences = parse_fasta_headers(input_fasta)
print(f"Loaded {len(sequences)} sequences", file=sys.stderr)

# Parse InterPro results
sequences = parse_interpro_results(interpro_output, sequences)
domain_count = sum(1 for s in sequences.values() if s['domains'])
print(f"Sequences with domain annotations: {domain_count}", file=sys.stderr)

# Load reference classifications
ref_classifications = load_gpcrdb_classifications(local_db_dir)

# Classify each sequence
results = []
class_counts = defaultdict(int)
subfamily_counts = defaultdict(int)

for seq_id, info in sequences.items():
    info['id'] = seq_id

    # Classify GPCR class
    gpcr_class = classify_gpcr_class(info)
    info['class'] = gpcr_class
    class_counts[gpcr_class] += 1

    # Classify subfamily (for Class A)
    subfamily = classify_subfamily(info, gpcr_class)
    info['subfamily'] = subfamily
    if subfamily:
        subfamily_counts[subfamily] += 1

    # Get domain summary
    domain_list = [d['accession'] for d in info['domains'] if d['accession']]
    interpro_list = [d['interpro'] for d in info['domains'] if d['interpro']]

    results.append({
        'id': seq_id,
        'gpcr_class': gpcr_class,
        'subfamily': subfamily or '',
        'n_domains': len(info['domains']),
        'pfam_domains': ','.join(set(d for d in domain_list if d.startswith('PF'))),
        'interpro_domains': ','.join(set(interpro_list)),
        'description': info['description'][:100]
    })

# Write results
output_file = f"{classify_dir}/gpcr_classifications.tsv"
df = pd.DataFrame(results)
df.to_csv(output_file, sep='\t', index=False)
print(f"\nClassifications written to: {output_file}", file=sys.stderr)

# Write summary
summary = {
    'total_sequences': len(sequences),
    'with_domains': domain_count,
    'class_distribution': dict(class_counts),
    'subfamily_distribution': dict(subfamily_counts)
}

summary_file = f"{classify_dir}/classification_summary.json"
with open(summary_file, 'w') as f:
    json.dump(summary, f, indent=2)

# Print summary
print("\n=== Classification Summary ===", file=sys.stderr)
print(f"Total sequences: {len(sequences)}", file=sys.stderr)
print("\nGPCR Class Distribution:", file=sys.stderr)
for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
    print(f"  {cls}: {count}", file=sys.stderr)

if subfamily_counts:
    print("\nClass A Sub-family Distribution:", file=sys.stderr)
    for sf, count in sorted(subfamily_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {sf}: {count}", file=sys.stderr)
PYTHON_SCRIPT

# Export variables for Python script
export CLASSIFY_DIR INPUT_FASTA

check_file "${CLASSIFY_DIR}/gpcr_classifications.tsv"

# --- Step 3: Create Classification Summary for Downstream ---
log "Step 3: Creating classification files for downstream analysis"

# Create simple mapping file for ranking
python3 << 'PYTHON_SCRIPT'
import pandas as pd
import os

classify_dir = os.environ.get('CLASSIFY_DIR', 'results/classification')

# Read classifications
df = pd.read_csv(f"{classify_dir}/gpcr_classifications.tsv", sep='\t')

# Create mapping for ranking (id -> class + subfamily)
mapping = df[['id', 'gpcr_class', 'subfamily']].copy()
mapping['classification'] = mapping.apply(
    lambda r: f"{r['gpcr_class']}:{r['subfamily']}" if r['subfamily'] else r['gpcr_class'],
    axis=1
)

# Write mapping
mapping[['id', 'classification']].to_csv(
    f"{classify_dir}/classification_mapping.tsv",
    sep='\t', index=False, header=False
)

# Write chemoreceptor candidates (subfamilies likely to be chemosensory)
chemosensory_subfamilies = ['Olfactory', 'Gustatory', 'Vomeronasal', 'Trace_amine',
                           'Formyl_peptide', 'Other_chemosensory']
chemoreceptors = df[df['subfamily'].isin(chemosensory_subfamilies)]
chemoreceptors['id'].to_csv(
    f"{classify_dir}/putative_chemoreceptors.txt",
    index=False, header=False
)

print(f"Putative chemoreceptors identified: {len(chemoreceptors)}")
PYTHON_SCRIPT

# Mark step as completed
create_checkpoint "classify_gpcrs"

# Also create legacy completion flag for downstream steps
touch "${RESULTS_DIR}/step_completed_02b.txt"

log "GPCR classification completed"
log "Output files:"
log "  Classifications: ${CLASSIFY_DIR}/gpcr_classifications.tsv"
log "  Summary: ${CLASSIFY_DIR}/classification_summary.json"
log "  Mapping: ${CLASSIFY_DIR}/classification_mapping.tsv"
log "  Putative chemoreceptors: ${CLASSIFY_DIR}/putative_chemoreceptors.txt"
