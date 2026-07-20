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

# Bead -e59: SYNTENY_BACKEND=jcvi (default) uses JCVI MCscan against close
# mollusc genomes (Aplysia first-tier); SYNTENY_BACKEND=mcscanx falls back
# to the legacy MCScanX path. Both paths produce a per-Berghia-gene CSV
# at ${RESULTS_DIR}/synteny/jcvi_anchors_per_gene.csv consumed by
# rank_candidates.py via TANDEM_CLUSTERS_FILE-style env var SYNTENY_ANCHORS_CSV.
SYNTENY_BACKEND="${SYNTENY_BACKEND:-jcvi}"
SYNTENY_ANCHORS_CSV="${RESULTS_DIR}/synteny/jcvi_anchors_per_gene.csv"

# --- JCVI MCscan branch (preferred when annotated genomes available) ---
if [ "$SYNTENY_BACKEND" = "jcvi" ]; then
    if [ ! -f "${GENOME_GFF:-}" ] || [ ! -f "${GENOME_PROTEIN:-}" ]; then
        log --level=WARN "JCVI synteny requires GENOME_GFF and GENOME_PROTEIN; falling back to MCScanX backend"
        SYNTENY_BACKEND="mcscanx"
    fi
fi

if [ "$SYNTENY_BACKEND" = "jcvi" ]; then
    # Discover target species: any *.proteins.fa in GENOME_DIR that is NOT
    # the Berghia file. Pair each with its corresponding .gff3.
    JCVI_OUT="${RESULTS_DIR}/synteny/jcvi"
    mkdir -p "$JCVI_OUT"
    > "$SYNTENY_ANCHORS_CSV"
    JCVI_HEADER_WRITTEN=0
    for tgt_pep in "${GENOME_DIR}"/*.proteins.fa; do
        [ -f "$tgt_pep" ] || continue
        # Skip Berghia itself (the symlink and the canonical name both point at it)
        if [ "$(realpath "$tgt_pep")" = "$(realpath "$GENOME_PROTEIN" 2>/dev/null || echo $tgt_pep)" ]; then
            continue
        fi
        tgt_base=$(basename "$tgt_pep" .proteins.fa)
        tgt_gff="${GENOME_DIR}/${tgt_base}.gff3"
        if [ ! -f "$tgt_gff" ]; then
            log --level=WARN "Skipping JCVI for $tgt_base — no matching .gff3"
            continue
        fi
        log "JCVI MCscan: ${BERGHIA_FILE_PREFIX} vs ${tgt_base}"
        bash "${SCRIPTS_DIR}/run_jcvi_synteny.sh" \
            --berghia-prefix="berghia" \
            --berghia-gff="$GENOME_GFF" --berghia-pep="$GENOME_PROTEIN" \
            --target="$tgt_pep" --target-gff="$tgt_gff" \
            --target-prefix="$tgt_base" \
            --output-dir="$JCVI_OUT/berghia_vs_${tgt_base}" \
            2> "${LOGS_DIR}/jcvi_berghia_vs_${tgt_base}.err" || \
            log --level=WARN "JCVI failed for ${tgt_base} (see logs)"
        anc="$JCVI_OUT/berghia_vs_${tgt_base}/berghia.${tgt_base}.lifted.anchors"
        [ -s "$anc" ] || anc="$JCVI_OUT/berghia_vs_${tgt_base}/berghia.${tgt_base}.anchors"
        if [ -s "$anc" ]; then
            tmp_csv="${JCVI_OUT}/berghia_vs_${tgt_base}/per_gene.csv"
            python3 "${SCRIPTS_DIR}/parse_jcvi_anchors.py" \
                --anchors "$anc" --target-species "$tgt_base" --out "$tmp_csv" \
                && {
                    if [ "$JCVI_HEADER_WRITTEN" -eq 0 ]; then
                        cp "$tmp_csv" "$SYNTENY_ANCHORS_CSV"
                        JCVI_HEADER_WRITTEN=1
                    else
                        tail -n +2 "$tmp_csv" >> "$SYNTENY_ANCHORS_CSV"
                    fi
                }
        fi
    done
    if [ ! -s "$SYNTENY_ANCHORS_CSV" ]; then
        log --level=WARN "No JCVI anchors produced (no target genomes? all failed?). " \
            "Drop SYNTENY_WEIGHT=0 in config.sh if you want to disable synteny scoring entirely."
        # Touch empty placeholder so downstream stages don't fail on file-not-found
        echo "candidate_id,n_anchor_blocks,total_anchor_genes,target_species" > "$SYNTENY_ANCHORS_CSV"
    fi
    touch "${RESULTS_DIR}/step_completed_synteny.txt"
    log "JCVI synteny step completed → $SYNTENY_ANCHORS_CSV"
    exit 0
fi

# --- Validate genome count for legacy MCScanX backend ---
genome_count=$(find "${GENOME_DIR}" -name "*.fasta" -type f 2>/dev/null | wc -l)
if [ "$genome_count" -lt 2 ]; then
    log --level=WARN "MCScanX synteny analysis requires at least 2 genomes, found ${genome_count}"
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
        run_command "samtools_view_${taxid_sample}" ${SAMTOOLS} view -bS "${RESULTS_DIR}/mapping/${taxid_sample}.sam" -o "${RESULTS_DIR}/mapping/${taxid_sample}_unsorted.bam"
        run_command "samtools_sort_${taxid_sample}" ${SAMTOOLS} sort -o "${RESULTS_DIR}/mapping/${taxid_sample}.bam" "${RESULTS_DIR}/mapping/${taxid_sample}_unsorted.bam"
        run_command "samtools_index_${taxid_sample}" ${SAMTOOLS} index "${RESULTS_DIR}/mapping/${taxid_sample}.bam"
        rm -f "${RESULTS_DIR}/mapping/${taxid_sample}_unsorted.bam"
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
    # The pipeline writes .gff3 (scripts/fetch_berghia_genome.sh:32, and the
    # JCVI branch above reads .gff3); probing only the bare .gff meant this
    # branch never found an annotation. .gff is kept as a second choice so a
    # hand-staged file still works.
    gff_file=""
    for gff_candidate in "${GENOME_DIR}/${taxid}.gff3" "${GENOME_DIR}/${taxid}.gff"; do
        if [ -f "$gff_candidate" ]; then
            gff_file="$gff_candidate"
            break
        fi
    done

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

# Bead -ar8: intra-genome tandem-cluster detection on Berghia GPCR candidates.
# Reads $GENOME_GFF + the candidate FASTA, writes per-candidate cluster size +
# ID to $TANDEM_CLUSTERS_FILE (consumed by rank_candidates.py in stage 07).
# Gated on RUN_TANDEM_DETECTION (default 1). Cheap (~seconds) given gffutils
# SQLite caching; safe to re-run.
if [ "${RUN_TANDEM_DETECTION:-1}" = "1" ] && [ -f "${GENOME_GFF:-}" ]; then
    CANDIDATES_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
    if [ -f "$CANDIDATES_FA" ]; then
        log "Running intra-genome tandem-cluster detection on Berghia candidates..."
        bash "${SCRIPTS_DIR}/run_tandem_detection.sh" "$CANDIDATES_FA" \
            2>> "${LOGS_DIR}/tandem_detection.err" \
            || log --level=WARN "Tandem-cluster detection failed (rank_candidates will skip the tandem axis)"
    else
        log --level=WARN "No candidate FASTA at $CANDIDATES_FA — skipping tandem-cluster detection"
    fi
fi

touch "${RESULTS_DIR}/step_completed_synteny.txt"
log "Synteny and mapping analysis completed."
