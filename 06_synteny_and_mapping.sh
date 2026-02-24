#!/bin/bash
# 06_synteny_and_mapping.sh
# Purpose: Analyze synteny across available genomes using MCScanX; mapping is optional if transcriptomes exist.
# Inputs: Genomes in ${GENOME_DIR}/*.fasta, transcriptomes for mapping if available
# Outputs: Synteny results in ${RESULTS_DIR}/synteny/, optional mapping BAMs in ${RESULTS_DIR}/mapping/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=synteny_mapping
#SBATCH --output=${LOGS_DIR}/06_synteny_and_mapping_%j.out
#SBATCH --error=${LOGS_DIR}/06_synteny_and_mapping_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/mapping" "${RESULTS_DIR}/synteny" "${RESULTS_DIR}/synteny/blast" "${RESULTS_DIR}/synteny/gff" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_extract_berghia.txt"

log "Starting synteny and mapping analysis."

# --- Validate genome count for synteny analysis ---
genome_count=$(find "${GENOME_DIR}" -name "*.fasta" -type f 2>/dev/null | wc -l)
if [ "$genome_count" -lt 2 ]; then
    log --level=WARN "Synteny analysis requires at least 2 genomes, found ${genome_count}"
    log --level=WARN "Skipping synteny analysis. Set SYNTENY_WEIGHT=0 in config.sh to disable synteny scoring."
    # Create empty output files to prevent downstream failures
    touch "${RESULTS_DIR}/synteny/synteny_ids.txt"
    touch "${RESULTS_DIR}/step_completed_synteny.txt"
    log "Synteny step completed (skipped due to insufficient genomes)"
    exit 0
fi

log "Found ${genome_count} genomes for synteny analysis"

# --- Optional mapping of transcriptomes to genomes ---
# Note: minimap2 -ax splice requires NUCLEOTIDE sequences (mRNA/cDNA), not protein
for genome in "${GENOME_DIR}"/*.fasta; do
    [ -f "$genome" ] || continue
    taxid_sample=$(basename "$genome" .fasta)

    # Look for nucleotide transcriptome (.mrna, .cds, .fa, .fna)
    nuc_trans=""
    for ext in mrna cds fa fna; do
        if [ -f "${TRANSCRIPTOME_DIR}/${taxid_sample}.${ext}" ]; then
            nuc_trans="${TRANSCRIPTOME_DIR}/${taxid_sample}.${ext}"
            break
        fi
    done

    if [ -n "$nuc_trans" ]; then
        log "Mapping ${taxid_sample} transcriptome to genome"
        # Validate minimap2 is available before running
        if ! command -v "${MINIMAP2}" &>/dev/null; then
            log --level=WARN "minimap2 not found, skipping mapping for ${taxid_sample}"
            continue
        fi
        run_command "minimap2_${taxid_sample}" --stdout="${RESULTS_DIR}/mapping/${taxid_sample}.sam" ${MINIMAP2} -ax splice -uf -k14 "$genome" "$nuc_trans"
        run_command "samtools_${taxid_sample}" ${SAMTOOLS} view -bS "${RESULTS_DIR}/mapping/${taxid_sample}.sam" | ${SAMTOOLS} sort -o "${RESULTS_DIR}/mapping/${taxid_sample}.bam"
        run_command "samtools_index_${taxid_sample}" ${SAMTOOLS} index "${RESULTS_DIR}/mapping/${taxid_sample}.bam"
    else
        log "Note: No nucleotide transcriptome for $taxid_sample, skipping mapping"
    fi
done

# --- Prepare MCScanX input files ---
# MCScanX requires:
# 1. All-vs-all BLAST results (xyz.blast)
# 2. GFF/BED files with gene positions (xyz.gff)

log "Preparing MCScanX input files"

# Create protein BLAST database for each genome's predicted proteins
for genome in "${GENOME_DIR}"/*.fasta; do
    [ -f "$genome" ] || continue
    taxid=$(basename "$genome" .fasta)

    # Look for protein predictions (from gene annotation)
    prot_file="${GENOME_DIR}/${taxid}.proteins.fa"
    gff_file="${GENOME_DIR}/${taxid}.gff"

    if [ -f "$prot_file" ]; then
        # Create BLAST database
        run_command "makeblastdb_${taxid}" makeblastdb -in "$prot_file" -dbtype prot -out "${RESULTS_DIR}/synteny/blast/${taxid}_db"

        # Convert GFF to MCScanX format if available
        if [ -f "$gff_file" ]; then
            awk -F'\t' '$3=="gene" || $3=="CDS" {
                split($9, attrs, ";");
                for (i in attrs) {
                    if (attrs[i] ~ /^ID=/ || attrs[i] ~ /^Name=/) {
                        gsub(/ID=|Name=/, "", attrs[i]);
                        gene_id = attrs[i];
                        break;
                    }
                }
                print $1, gene_id, $4, $5
            }' "$gff_file" > "${RESULTS_DIR}/synteny/gff/${taxid}.gff"
        fi
    else
        log "Warning: No protein file for ${taxid}, using GPCR candidates"
        # Use GPCR candidates if no genome-wide proteins
        if [ -f "${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid}.fa" ]; then
            run_command "makeblastdb_${taxid}" makeblastdb -in "${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid}.fa" -dbtype prot -out "${RESULTS_DIR}/synteny/blast/${taxid}_db"
        fi
    fi
done

# --- Run all-vs-all BLASTP for synteny detection ---
log "Running pairwise BLASTP for synteny analysis"

genome_list=("${GENOME_DIR}"/*.fasta)
for ((i=0; i<${#genome_list[@]}; i++)); do
    [ -f "${genome_list[$i]}" ] || continue
    taxid1=$(basename "${genome_list[$i]}" .fasta)

    for ((j=i+1; j<${#genome_list[@]}; j++)); do
        [ -f "${genome_list[$j]}" ] || continue
        taxid2=$(basename "${genome_list[$j]}" .fasta)

        # Check if BLAST databases exist
        if [ -f "${RESULTS_DIR}/synteny/blast/${taxid1}_db.phr" ] && [ -f "${RESULTS_DIR}/synteny/blast/${taxid2}_db.phr" ]; then
            # Get query proteins
            query_prot="${GENOME_DIR}/${taxid1}.proteins.fa"
            [ -f "$query_prot" ] || query_prot="${RESULTS_DIR}/chemogpcrs/chemogpcrs_${taxid1}.fa"

            if [ -f "$query_prot" ]; then
                log "Running: blastp_${taxid1}_vs_${taxid2}"
                blastp -query "$query_prot" -db "${RESULTS_DIR}/synteny/blast/${taxid2}_db" -outfmt 6 -evalue 1e-5 -num_threads "${CPUS}" >> "${RESULTS_DIR}/synteny/${taxid1}_${taxid2}.blast" 2> "${LOGS_DIR}/blastp_${taxid1}_vs_${taxid2}.err"
            fi
        fi
    done
done

# --- Run MCScanX for synteny detection ---
log "Running MCScanX synteny analysis"

for blast_file in "${RESULTS_DIR}/synteny/"*.blast; do
    [ -f "$blast_file" ] || continue
    base=$(basename "$blast_file" .blast)

    # Extract taxids from filename (pipeline-controlled format: taxid1_taxid2.blast)
    # Note: This is parsing pipeline-generated filenames, not external data headers
    taxid1=$(echo "$base" | cut -d'_' -f1)
    taxid2=$(echo "$base" | cut -d'_' -f2-)

    # Check for GFF files
    gff1="${RESULTS_DIR}/synteny/gff/${taxid1}.gff"
    gff2="${RESULTS_DIR}/synteny/gff/${taxid2}.gff"

    if [ -f "$gff1" ] && [ -f "$gff2" ]; then
        # Combine GFF files for MCScanX
        cat "$gff1" "$gff2" > "${RESULTS_DIR}/synteny/${base}.gff"

        # Run MCScanX
        run_command "mcscanx_${base}" ${MCSCANX} -a -b 2 "${RESULTS_DIR}/synteny/${base}"
    else
        log "Warning: Missing GFF files for ${base}, skipping MCScanX"
    fi
done

# --- Extract synteny-supported gene IDs ---
log "Extracting synteny-supported genes"
> "${RESULTS_DIR}/synteny/synteny_ids.txt"

for collinearity in "${RESULTS_DIR}/synteny/"*.collinearity; do
    [ -f "$collinearity" ] || continue
    # Extract gene IDs from collinearity blocks
    grep -v "^#" "$collinearity" | awk '{print $2; print $3}' >> "${RESULTS_DIR}/synteny/synteny_ids.txt"
done

# Deduplicate
sort -u "${RESULTS_DIR}/synteny/synteny_ids.txt" > "${RESULTS_DIR}/synteny/synteny_ids_unique.txt"
mv "${RESULTS_DIR}/synteny/synteny_ids_unique.txt" "${RESULTS_DIR}/synteny/synteny_ids.txt"

# --- Plot synteny at different taxonomic levels ---
python3 "${SCRIPTS_DIR}/plot_synteny.py" "${RESULTS_DIR}/synteny" "${RESULTS_DIR}/synteny/synteny_plot" || log "Warning: Synteny plotting failed"

touch "${RESULTS_DIR}/step_completed_synteny.txt"
log "Synteny and mapping analysis completed."
