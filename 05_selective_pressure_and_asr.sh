#!/bin/bash
# 05_selective_pressure_and_asr.sh
# Purpose: Analyze selective pressure with HyPhy's aBSREL and reconstruct ancestral sequences with FastML.
# Inputs: Orthogroups from step 03, alignments from step 04
# Outputs: Selective pressure results in ${RESULTS_DIR}/selective_pressure/, ASR sequences in ${RESULTS_DIR}/asr/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=selective_pressure_asr
#SBATCH --output=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.out
#SBATCH --error=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --array=0-999%50  # Adjusted for orthogroup processing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/selective_pressure" "${RESULTS_DIR}/selective_pressure/nucleotide" "${RESULTS_DIR}/asr" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_lse_classification.txt"

log "Starting selective pressure and ASR analysis."

# Get orthogroups
ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)

# Handle case where no orthogroups exist
if [ ${#ORTHOGROUPS[@]} -eq 0 ] || [ ! -f "${ORTHOGROUPS[0]}" ]; then
    log "Error: No orthogroups found"
    exit 1
fi

# Handle SLURM array indexing - skip if index exceeds available orthogroups
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    if [ "$SLURM_ARRAY_TASK_ID" -ge ${#ORTHOGROUPS[@]} ]; then
        log "Skipping: Array index ${SLURM_ARRAY_TASK_ID} exceeds available orthogroups (${#ORTHOGROUPS[@]})"
        exit 0
    fi
    og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
else
    # Non-array mode: process first orthogroup (for testing)
    og="${ORTHOGROUPS[0]}"
fi

[ -z "$og" ] || [ ! -f "$og" ] && { log "Skipping missing orthogroup"; exit 0; }

base=$(basename "$og" .fa)
taxids=$(grep "^>" "$og" | sed 's/>//' | cut -d'_' -f1 | sort -u)
taxa_count=$(echo "$taxids" | wc -l)

# --- Helper function to find nucleotide sequences for an orthogroup ---
find_nucleotide_sequences() {
    local protein_file="$1"
    local output_file="$2"

    > "$output_file"

    # Extract sequence IDs from protein alignment
    while IFS= read -r header; do
        seq_id=$(echo "$header" | sed 's/>//')
        taxid_sample=$(echo "$seq_id" | cut -d'_' -f1,2)
        protein_id=$(echo "$seq_id" | cut -d'_' -f3-)

        # Try to find nucleotide sequence in multiple locations
        found=false
        for ext in mrna cds fna fa; do
            nuc_file="${TRANSCRIPTOME_DIR}/${taxid_sample}.${ext}"
            if [ -f "$nuc_file" ]; then
                # Extract the corresponding nucleotide sequence
                # Try exact match first, then partial match
                if grep -q "^>${seq_id}" "$nuc_file"; then
                    grep -A1 "^>${seq_id}" "$nuc_file" | head -2 >> "$output_file"
                    found=true
                    break
                elif grep -q "^>${protein_id}" "$nuc_file"; then
                    # Rename header to match protein
                    grep -A1 "^>${protein_id}" "$nuc_file" | head -2 | sed "1s/.*/>$seq_id/" >> "$output_file"
                    found=true
                    break
                fi
            fi
        done

        if [ "$found" = false ]; then
            log "Warning: No nucleotide sequence found for ${seq_id}"
        fi
    done < <(grep "^>" "$protein_file")

    # Return success if we found at least some sequences
    [ -s "$output_file" ]
}

# --- Selective Pressure with aBSREL ---
if [ "$taxa_count" -gt 1 ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    tree="${RESULTS_DIR}/phylogenies/protein/${base}.treefile"

    if [ -f "$protein_align" ] && [ -f "$tree" ]; then
        # Find nucleotide sequences for this orthogroup
        nuc_align="${RESULTS_DIR}/selective_pressure/nucleotide/${base}_nuc.fa"

        if find_nucleotide_sequences "$protein_align" "$nuc_align"; then
            # Create codon alignment
            run_command "${base}_codon" pal2nal.pl "$protein_align" "$nuc_align" -output paml > "${RESULTS_DIR}/selective_pressure/${base}_codon.phy" 2>/dev/null

            if [ -s "${RESULTS_DIR}/selective_pressure/${base}_codon.phy" ]; then
                # Run aBSREL
                run_command "${base}_absrel" hyphy aBSREL --alignment "${RESULTS_DIR}/selective_pressure/${base}_codon.phy" --tree "$tree" --output "${RESULTS_DIR}/selective_pressure/${base}_absrel.json"

                # Parse results
                if [ -f "${RESULTS_DIR}/selective_pressure/${base}_absrel.json" ]; then
                    python3 "${SCRIPTS_DIR}/parse_absrel.py" "${RESULTS_DIR}/selective_pressure/${base}_absrel.json" "${RESULTS_DIR}/selective_pressure/absrel_results.csv" || log "Warning: Failed to parse aBSREL for $base"
                fi
            else
                log "Warning: pal2nal failed for ${base}, skipping aBSREL"
            fi
        else
            log "Warning: Could not find nucleotide sequences for ${base}, skipping dN/dS analysis"
        fi
    else
        log "Warning: Missing alignment or tree for ${base}"
    fi
fi

# --- Berghia-specific LSEs with ASR ---
# Check if this orthogroup contains Berghia sequences AND has multiple sequences (for meaningful ASR)
berghia_count=$(echo "$taxids" | grep -c "${BERGHIA_TAXID}" || echo 0)
seq_count=$(grep -c "^>" "$og")

if [ "$berghia_count" -gt 0 ] && [ "$seq_count" -gt 2 ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    tree="${RESULTS_DIR}/phylogenies/protein/${base}.treefile"

    if [ -f "$protein_align" ] && [ -f "$tree" ]; then
        # Find nucleotide sequences
        nuc_align="${RESULTS_DIR}/selective_pressure/nucleotide/${base}_nuc.fa"

        if [ ! -f "$nuc_align" ]; then
            find_nucleotide_sequences "$protein_align" "$nuc_align"
        fi

        if [ -f "$nuc_align" ] && [ -s "$nuc_align" ]; then
            # Run LSE-specific aBSREL if not already done
            if [ ! -f "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json" ]; then
                codon_file="${RESULTS_DIR}/selective_pressure/${base}_codon.phy"
                if [ ! -f "$codon_file" ]; then
                    run_command "${base}_codon_lse" pal2nal.pl "$protein_align" "$nuc_align" -output paml > "$codon_file" 2>/dev/null
                fi

                if [ -s "$codon_file" ]; then
                    run_command "${base}_absrel_lse" hyphy aBSREL --alignment "$codon_file" --tree "$tree" --output "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json"
                    [ -f "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json" ] && \
                        python3 "${SCRIPTS_DIR}/parse_absrel.py" "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json" "${RESULTS_DIR}/selective_pressure/absrel_results_lse.csv"
                fi
            fi

            # ASR for deep nodes
            deep_nodes=$(python3 "${SCRIPTS_DIR}/select_deep_nodes.py" "$tree" "${BERGHIA_TAXID}" "${MIN_ASR_DISTANCE}" 2>/dev/null)

            if [ -n "$deep_nodes" ]; then
                for node in $deep_nodes; do
                    run_command "${base}_${node}_asr" ${FASTML} --seq "$nuc_align" --tree "$tree" \
                        --out_seq "${RESULTS_DIR}/asr/${base}_${node}_asr.fa" \
                        --out_tree "${RESULTS_DIR}/asr/${base}_${node}_asr.tree" \
                        --node "$node" -t "${CPUS}" --verbose 2>/dev/null || log "Warning: FastML failed for ${base} node ${node}"
                done

                # Plot ASR for first deep node
                first_node="${deep_nodes%% *}"
                if [ -f "${RESULTS_DIR}/asr/${base}_${first_node}_asr.fa" ]; then
                    python3 "${SCRIPTS_DIR}/plot_asr.py" "$tree" "${RESULTS_DIR}/asr/${base}_${first_node}_asr.fa" "${RESULTS_DIR}/asr/${base}_asr_plot" || log "Warning: ASR plotting failed for ${base}"
                fi
            fi
        fi
    fi
fi

# Create completion flag for this orthogroup
touch "${RESULTS_DIR}/selective_pressure/step_completed_${base}.txt"

log "Selective pressure and ASR completed for ${base}."

# If this is the last array job, create overall completion flag
if [ -n "$SLURM_ARRAY_TASK_ID" ] && [ "$SLURM_ARRAY_TASK_ID" -eq $((${#ORTHOGROUPS[@]} - 1)) ]; then
    touch "${RESULTS_DIR}/step_completed_selective_pressure.txt"
fi
