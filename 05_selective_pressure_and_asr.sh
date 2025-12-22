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

# Get orthogroup count using manifest (preferred) or fallback to globbing
MANIFEST_FILE="${RESULTS_DIR}/orthogroup_manifest.tsv"
OG_COUNT=$(get_orthogroup_count "$MANIFEST_FILE")

if [ "$OG_COUNT" -eq 0 ]; then
    log "Error: No orthogroups found"
    exit 1
fi

# Handle SLURM array indexing using helper function
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    validate_array_index "$SLURM_ARRAY_TASK_ID" "$OG_COUNT"

    # Get orthogroup name from manifest or fallback to globbing
    if [ -f "$MANIFEST_FILE" ]; then
        base=$(get_orthogroup_by_index "$SLURM_ARRAY_TASK_ID" "$MANIFEST_FILE")
    else
        ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)
        og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
        base=$(basename "$og" .fa)
    fi
else
    # Non-array mode: process first orthogroup (for testing)
    if [ -f "$MANIFEST_FILE" ]; then
        base=$(get_orthogroup_by_index 0 "$MANIFEST_FILE")
    else
        ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)
        base=$(basename "${ORTHOGROUPS[0]}" .fa)
    fi
fi

# Find the orthogroup FASTA file
og=$(find "${RESULTS_DIR}/orthogroups" -name "${base}.fa" -type f 2>/dev/null | head -1)
[ -z "$og" ] || [ ! -f "$og" ] && { log "Skipping missing orthogroup: ${base}"; exit 0; }
# Use centralized metadata lookup to get taxids (excludes references)
taxids=$(get_taxids_from_fasta "$og" --exclude-refs)
taxa_count=$(echo "$taxids" | wc -w)

# --- Define reference CDS file location ---
REFERENCE_CDS_FILE="${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna"

# --- Helper function to find nucleotide sequences for an orthogroup ---
# Searches both transcriptome directory (for target species) and reference CDS
find_nucleotide_sequences() {
    local protein_file="$1"
    local output_file="$2"

    > "$output_file"
    local found_count=0
    local missing_count=0
    local ref_count=0

    # Extract sequence IDs from protein alignment
    while IFS= read -r header; do
        seq_id=$(echo "$header" | sed 's/>//')

        # Check if this is a reference sequence
        is_ref=$(is_reference_seq "$seq_id" && echo "yes" || echo "no")

        # Use metadata lookup to get taxid (with fallback to header parsing)
        taxid=$(get_taxid_for_seq "$seq_id")

        # Extract protein ID portion (after taxid prefix)
        # For ref_TAXID_N format: protein_id is N
        # For TAXID_N format: protein_id is N
        if [[ "$seq_id" == ref_* ]]; then
            protein_id=$(echo "$seq_id" | cut -d'_' -f3-)
        else
            protein_id=$(echo "$seq_id" | cut -d'_' -f2-)
        fi

        found=false

        # --- Check reference CDS file first for reference sequences ---
        if [ "$is_ref" = "yes" ] && [ -f "$REFERENCE_CDS_FILE" ]; then
            # Try exact match with the seq_id
            if grep -q "^>${seq_id}" "$REFERENCE_CDS_FILE"; then
                # Extract sequence (handle multi-line FASTA)
                awk -v id="$seq_id" '
                    /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print; found=1; next } else { found=0 } }
                    found { print }
                ' "$REFERENCE_CDS_FILE" >> "$output_file"
                found=true
                ((ref_count++))
            fi
        fi

        # --- Search transcriptome directory for non-reference sequences ---
        if [ "$found" = false ]; then
            for ext in mrna cds fna fa; do
                # Try different naming conventions for nucleotide files
                for nuc_file in "${TRANSCRIPTOME_DIR}/${taxid}.${ext}" \
                               "${TRANSCRIPTOME_DIR}/taxid_${taxid}.${ext}" \
                               "${TRANSCRIPTOME_DIR}/${taxid}_${taxid}.${ext}"; do
                    if [ -f "$nuc_file" ]; then
                        # Extract the corresponding nucleotide sequence
                        # Try exact match first, then partial match
                        if grep -q "^>${seq_id}" "$nuc_file"; then
                            # Handle multi-line FASTA
                            awk -v id="$seq_id" '
                                /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print; found=1; next } else { found=0 } }
                                found { print }
                            ' "$nuc_file" >> "$output_file"
                            found=true
                            break 2
                        elif grep -q "^>${protein_id}" "$nuc_file"; then
                            # Extract and rename header to match protein
                            awk -v id="$protein_id" -v new_id="$seq_id" '
                                /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print ">" new_id; found=1; next } else { found=0 } }
                                found { print }
                            ' "$nuc_file" >> "$output_file"
                            found=true
                            break 2
                        fi
                    fi
                done
            done
        fi

        # --- Fallback: check all reference CDS for non-reference sequences too ---
        if [ "$found" = false ] && [ -f "$REFERENCE_CDS_FILE" ]; then
            if grep -q "^>${seq_id}" "$REFERENCE_CDS_FILE"; then
                awk -v id="$seq_id" '
                    /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print; found=1; next } else { found=0 } }
                    found { print }
                ' "$REFERENCE_CDS_FILE" >> "$output_file"
                found=true
            fi
        fi

        if [ "$found" = true ]; then
            ((found_count++))
        else
            ((missing_count++))
            log "Warning: No nucleotide sequence found for ${seq_id}"
        fi
    done < <(grep "^>" "$protein_file")

    # Log summary
    log "Nucleotide lookup: found=${found_count} (${ref_count} from refs), missing=${missing_count}"

    # Return success if we found at least some sequences (majority required for meaningful analysis)
    local total=$((found_count + missing_count))
    if [ "$found_count" -gt 0 ] && [ "$found_count" -ge $((total / 2)) ]; then
        return 0
    else
        return 1
    fi
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

# Create array checkpoint for resume capability
[ -n "$SLURM_ARRAY_TASK_ID" ] && create_array_checkpoint "05_selective" "$SLURM_ARRAY_TASK_ID"

log "Selective pressure and ASR completed for ${base}."

# Create overall step completion flag
# Note: We create this flag after each task completes, since downstream steps
# can start processing as results become available
touch "${RESULTS_DIR}/step_completed_05.txt"
