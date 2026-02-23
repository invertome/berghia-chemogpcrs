#!/bin/bash
# 01_reference_processing.sh
# Purpose: Process reference FASTA files and build HMMs if custom HMMs are not provided.
#          Generates an ID mapping file for reference sequences.
#          Supports both legacy format and new nath_et_al directory structure.
# Inputs: Reference FASTA files in ${REFERENCE_DIR}
#         Legacy: taxid_conserved_refs.aa, taxid_lse_refs.aa
#         New: nath_et_al/{lse,one_to_one_ortholog}/{taxonomic_group}/Species_name.faa
# Outputs: Combined reference FASTA, HMMs, ID map CSV, optional CDS files
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=reference_processing
#SBATCH --output=${LOGS_DIR}/01_reference_processing_%j.out
#SBATCH --error=${LOGS_DIR}/01_reference_processing_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/reference_sequences" "${RESULTS_DIR}/hmms" "${LOGS_DIR}" "${SCRIPTS_DIR}" || { log "Error: Cannot create directories"; exit 1; }
mkdir -p "${RESULTS_DIR}/reference_sequences/cds" || { log "Error: Cannot create CDS directory"; exit 1; }

log "Starting reference processing."

# --- Detect reference file structure ---
NATH_ET_AL_DIR="${REFERENCE_DIR}/nath_et_al"
USE_NATH_ET_AL=false

if [ -d "${NATH_ET_AL_DIR}" ]; then
    USE_NATH_ET_AL=true
    log "Detected nath_et_al reference structure"
fi

# --- Step 1: Rename reference files with taxid if needed (nath_et_al structure) ---
if [ "${USE_NATH_ET_AL}" = true ]; then
    # Check if files already have taxid prefix
    has_taxid_prefix=$(find "${NATH_ET_AL_DIR}" -name "*.faa" | head -1 | xargs basename 2>/dev/null | grep -E "^[0-9]+_" || true)

    if [ -z "${has_taxid_prefix}" ]; then
        log "Renaming reference files to include taxid prefix..."

        if [ -f "${SCRIPTS_DIR}/rename_references_with_taxid.py" ]; then
            run_command "ref_rename_taxid" python3 "${SCRIPTS_DIR}/rename_references_with_taxid.py" \
                "${NATH_ET_AL_DIR}" \
                --output-map "${RESULTS_DIR}/reference_sequences/taxid_rename_map.csv" \
                ${LOCAL_DB_DIR:+--db-dir "${LOCAL_DB_DIR}"}
        else
            log "Warning: rename_references_with_taxid.py not found, skipping taxid rename"
        fi
    else
        log "Reference files already have taxid prefix"
    fi
fi

# --- Step 2: Fetch CDS for reference proteins if needed (nath_et_al structure) ---
if [ "${USE_NATH_ET_AL}" = true ] && [ "${FETCH_REFERENCE_CDS:-true}" = true ]; then
    # Check if CDS files already exist
    cds_count=$(find "${RESULTS_DIR}/reference_sequences/cds" -name "*_cds.fna" 2>/dev/null | wc -l)

    if [ "${cds_count}" -eq 0 ]; then
        log "Fetching CDS sequences for reference proteins..."

        if [ -f "${SCRIPTS_DIR}/fetch_reference_cds.py" ]; then
            run_command "ref_fetch_cds" python3 "${SCRIPTS_DIR}/fetch_reference_cds.py" \
                "${NATH_ET_AL_DIR}" \
                -o "${RESULTS_DIR}/reference_sequences/cds" \
                --log-failed "${RESULTS_DIR}/reference_sequences/failed_cds_accessions.csv"

            log "CDS fetch complete. Check failed_cds_accessions.csv for any missing sequences."
        else
            log "Warning: fetch_reference_cds.py not found, skipping CDS fetch"
        fi
    else
        log "CDS files already exist (${cds_count} files found)"
    fi
fi

# --- Step 3: Combine reference sequences by category ---
if [ "${USE_NATH_ET_AL}" = true ]; then
    # New structure: combine by LSE vs conserved (one_to_one_ortholog)
    log "Combining reference sequences from nath_et_al structure..."

    # Combine LSE references
    find "${NATH_ET_AL_DIR}/lse" -name "*.faa" -exec cat {} + > "${RESULTS_DIR}/reference_sequences/lse_references.fa" 2>/dev/null || true
    lse_count=$(grep -c "^>" "${RESULTS_DIR}/reference_sequences/lse_references.fa" 2>/dev/null | tail -1)
    lse_count=${lse_count:-0}
    log "Combined ${lse_count} LSE reference sequences"

    # Combine conserved (one_to_one_ortholog) references
    find "${NATH_ET_AL_DIR}/one_to_one_ortholog" -name "*.faa" -exec cat {} + > "${RESULTS_DIR}/reference_sequences/conserved_references.fa" 2>/dev/null || true
    conserved_count=$(grep -c "^>" "${RESULTS_DIR}/reference_sequences/conserved_references.fa" 2>/dev/null | tail -1)
    conserved_count=${conserved_count:-0}
    log "Combined ${conserved_count} conserved reference sequences"

    # Combine all references
    cat "${RESULTS_DIR}/reference_sequences/lse_references.fa" "${RESULTS_DIR}/reference_sequences/conserved_references.fa" \
        > "${RESULTS_DIR}/reference_sequences/all_references.fa" 2>/dev/null

    # Combine CDS files similarly
    if [ -d "${RESULTS_DIR}/reference_sequences/cds/lse" ]; then
        find "${RESULTS_DIR}/reference_sequences/cds/lse" -name "*.fna" -exec cat {} + \
            > "${RESULTS_DIR}/reference_sequences/cds/lse_references_cds.fna" 2>/dev/null || true
    fi
    if [ -d "${RESULTS_DIR}/reference_sequences/cds/one_to_one_ortholog" ]; then
        find "${RESULTS_DIR}/reference_sequences/cds/one_to_one_ortholog" -name "*.fna" -exec cat {} + \
            > "${RESULTS_DIR}/reference_sequences/cds/conserved_references_cds.fna" 2>/dev/null || true
    fi
    cat "${RESULTS_DIR}/reference_sequences/cds/"*.fna > "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" 2>/dev/null || true

else
    # Legacy structure: use taxid_*_refs.aa files
    log "Using legacy reference structure..."

    # --- Process Conserved HMMs ---
    if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
        cp "${CONSERVED_HMM}" "${RESULTS_DIR}/hmms/conserved.hmm" || { log "Error: Failed to copy ${CONSERVED_HMM}"; exit 1; }
        log "Using custom conserved HMM: ${CONSERVED_HMM}"
    else
        for taxid in "${TAXA[@]}"; do
            conserved_aa="${REFERENCE_DIR}/${taxid}_conserved_refs.aa"
            if [ -f "$conserved_aa" ]; then
                run_command "ref_hmm_conserved_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_conserved.hmm" "$conserved_aa"
            fi
        done
        cat "${RESULTS_DIR}/hmms/"*_conserved.hmm > "${RESULTS_DIR}/hmms/conserved.hmm" 2>/dev/null || log "Warning: No conserved HMMs generated"
    fi

    # --- Process LSE HMMs ---
    if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
        cp "${LSE_HMM}" "${RESULTS_DIR}/hmms/lse.hmm" || { log "Error: Failed to copy ${LSE_HMM}"; exit 1; }
        log "Using custom LSE HMM: ${LSE_HMM}"
    else
        for taxid in "${TAXA[@]}"; do
            lse_aa="${REFERENCE_DIR}/${taxid}_lse_refs.aa"
            if [ -f "$lse_aa" ]; then
                run_command "ref_hmm_lse_${taxid}" ${HMMBUILD} "${RESULTS_DIR}/hmms/${taxid}_lse.hmm" "$lse_aa"
            fi
        done
        cat "${RESULTS_DIR}/hmms/"*_lse.hmm > "${RESULTS_DIR}/hmms/lse.hmm" 2>/dev/null || log "Warning: No LSE HMMs generated"
    fi

    # Combine all reference sequences (legacy)
    cat "${REFERENCE_DIR}"/*_conserved_refs.aa "${REFERENCE_DIR}"/*_lse_refs.aa \
        > "${RESULTS_DIR}/reference_sequences/all_references.fa" 2>/dev/null || { log "Error: No reference sequences found"; exit 1; }
fi

# --- Step 4: Build HMMs from references (nath_et_al structure) ---
#
# HMM Strategy: ${HMM_BUILD_STRATEGY:-per_species}
#
#   per_species  (default, preliminary) — Build one HMM per species FASTA file.
#       Each species has ~200-300 GPCR sequences, small enough for fast MAFFT
#       alignment. Multiple species-level HMMs are more sensitive than one giant
#       alignment because a Berghia GPCR similar to Aplysia's repertoire hits
#       that species' HMM strongly rather than being diluted by distant phyla.
#       Limitation: mixes functionally distinct receptor subtypes within a species.
#
#   per_orthogroup  (full analysis) — Cluster ALL 109K reference sequences into
#       ortholog groups via DIAMOND + MCL (or OrthoFinder), then build one HMM
#       per ortholog group. This produces functionally coherent profiles (e.g.,
#       "all lophotrochozoan rhodopsin-type receptors") and is the gold standard
#       for remote homolog detection. Requires all-vs-all DIAMOND (~109K seqs),
#       which needs a machine with >=32GB RAM and several hours of compute.
#       Set HMM_BUILD_STRATEGY=per_orthogroup in config.sh when resources allow.
#
if [ "${USE_NATH_ET_AL}" = true ]; then

    HMM_STRATEGY="${HMM_BUILD_STRATEGY:-per_species}"
    log "HMM build strategy: ${HMM_STRATEGY}"

    MAX_ALIGN_SEQS=500  # CD-HIT + subsample if a group exceeds this

    # Helper: align sequences and build one HMM from a FASTA file
    build_hmm_from_fasta() {
        local faa_file="$1"
        local hmm_out="$2"
        local label="$3"

        local seq_count
        seq_count=$(grep -c "^>" "$faa_file" 2>/dev/null | tail -1)
        seq_count=${seq_count:-0}
        [ "$seq_count" -lt 2 ] && return 1  # Need at least 2 sequences for alignment

        local align_input="$faa_file"

        if [ "$seq_count" -gt "$MAX_ALIGN_SEQS" ]; then
            local rep_fa="${hmm_out%.hmm}_reps.fa"
            ${CDHIT} -i "$faa_file" -o "$rep_fa" -c 0.7 -n 5 \
                -M "${CDHIT_MEMORY:-8000}" -T "${CPUS}" -d 0 > /dev/null 2>&1
            align_input="$rep_fa"
        fi

        local aligned="${hmm_out%.hmm}_aligned.fa"
        ${MAFFT} --auto --thread "${CPUS}" "$align_input" > "$aligned" 2>/dev/null
        [ -s "$aligned" ] || return 1

        ${HMMBUILD} --cpu "${CPUS}" "$hmm_out" "$aligned" > /dev/null 2>&1
    }

    # ---- per_orthogroup strategy (full analysis, gold standard) ----
    if [ "${HMM_STRATEGY}" = "per_orthogroup" ]; then
        log "Building per-orthogroup HMMs (DIAMOND + MCL clustering)..."
        DIAMOND="${DIAMOND:-diamond}"
        MCL="${MCL:-mcl}"

        OG_WORKDIR="${RESULTS_DIR}/hmms/orthogroup_build"
        mkdir -p "${OG_WORKDIR}"

        # Build combined FASTA per category for clustering
        for category in conserved lse; do
            if [ "$category" = "conserved" ]; then
                src_dir="${NATH_ET_AL_DIR}/one_to_one_ortholog"
            else
                src_dir="${NATH_ET_AL_DIR}/lse"
            fi

            combined="${OG_WORKDIR}/${category}_all.fa"
            find "$src_dir" -name "*.faa" -exec cat {} + > "$combined" 2>/dev/null

            total_seqs=$(grep -c "^>" "$combined" 2>/dev/null | tail -1)
            total_seqs=${total_seqs:-0}
            log "${category}: ${total_seqs} sequences for orthogroup clustering"

            # All-vs-all DIAMOND blastp
            db="${OG_WORKDIR}/${category}_db"
            hits="${OG_WORKDIR}/${category}_hits.tsv"
            run_command "diamond_makedb_${category}" ${DIAMOND} makedb --in "$combined" -d "$db"
            run_command "diamond_blastp_${category}" ${DIAMOND} blastp \
                -d "$db" -q "$combined" -o "$hits" \
                --very-sensitive -e 1e-5 --max-target-seqs 500 \
                --outfmt 6 qseqid sseqid pident length evalue bitscore \
                --threads "${CPUS}"

            # MCL clustering on DIAMOND hits (use -log10(evalue) as edge weight)
            abc="${OG_WORKDIR}/${category}_hits.abc"
            awk '{print $1, $2, $5}' "$hits" > "$abc"
            mci="${OG_WORKDIR}/${category}.mci"
            tab="${OG_WORKDIR}/${category}.tab"
            run_command "mcxload_${category}" mcxload -abc "$abc" --stream-mirror -o "$mci" -write-tab "$tab"
            run_command "mcl_${category}" ${MCL} "$mci" -I 2.0 -o "${OG_WORKDIR}/${category}_clusters.mcl" -use-tab "$tab"

            # Parse MCL clusters → one FASTA per orthogroup → align → hmmbuild
            hmm_dir="${RESULTS_DIR}/hmms/${category}_orthogroups"
            mkdir -p "$hmm_dir"

            python3 "${SCRIPTS_DIR}/split_mcl_to_fasta.py" \
                "${OG_WORKDIR}/${category}_clusters.mcl" "$combined" "$hmm_dir" \
                --min-seqs 3 --prefix "OG_${category}"

            og_hmm_count=0
            for og_fasta in "$hmm_dir"/OG_*.fa; do
                [ -f "$og_fasta" ] || continue
                og_name=$(basename "$og_fasta" .fa)
                hmm_file="${hmm_dir}/${og_name}.hmm"
                if build_hmm_from_fasta "$og_fasta" "$hmm_file" "$og_name"; then
                    og_hmm_count=$((og_hmm_count + 1))
                fi
            done

            cat "$hmm_dir"/*.hmm > "${RESULTS_DIR}/hmms/${category}.hmm" 2>/dev/null || true
            log "Built ${og_hmm_count} ${category} orthogroup HMMs"
        done

    # ---- per_species strategy (preliminary, fast) ----
    else
        log "Building per-species HMMs from reference sequences..."

        # Build conserved HMMs (one per species)
        if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
            cp "${CONSERVED_HMM}" "${RESULTS_DIR}/hmms/conserved.hmm"
            log "Using custom conserved HMM: ${CONSERVED_HMM}"
        else
            mkdir -p "${RESULTS_DIR}/hmms/conserved_species"
            conserved_hmm_count=0
            for faa in "${NATH_ET_AL_DIR}/one_to_one_ortholog"/*/*.faa; do
                [ -f "$faa" ] || continue
                species=$(basename "$faa" .faa)
                hmm_file="${RESULTS_DIR}/hmms/conserved_species/${species}.hmm"
                if build_hmm_from_fasta "$faa" "$hmm_file" "$species"; then
                    conserved_hmm_count=$((conserved_hmm_count + 1))
                fi
            done
            # Concatenate all species HMMs into one file for hmmsearch
            cat "${RESULTS_DIR}/hmms/conserved_species"/*.hmm > "${RESULTS_DIR}/hmms/conserved.hmm" 2>/dev/null || true
            log "Built ${conserved_hmm_count} conserved species HMMs"
        fi

        # Build LSE HMMs (one per species)
        if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
            cp "${LSE_HMM}" "${RESULTS_DIR}/hmms/lse.hmm"
            log "Using custom LSE HMM: ${LSE_HMM}"
        else
            mkdir -p "${RESULTS_DIR}/hmms/lse_species"
            lse_hmm_count=0
            for faa in "${NATH_ET_AL_DIR}/lse"/*/*.faa; do
                [ -f "$faa" ] || continue
                species=$(basename "$faa" .faa)
                hmm_file="${RESULTS_DIR}/hmms/lse_species/${species}.hmm"
                if build_hmm_from_fasta "$faa" "$hmm_file" "$species"; then
                    lse_hmm_count=$((lse_hmm_count + 1))
                fi
            done
            # Concatenate all species HMMs
            cat "${RESULTS_DIR}/hmms/lse_species"/*.hmm > "${RESULTS_DIR}/hmms/lse.hmm" 2>/dev/null || true
            log "Built ${lse_hmm_count} LSE species HMMs"
        fi
    fi
fi

# --- Step 5: Update headers with short IDs and generate ID map ---
# Check if we have any references
ref_count=$(grep -c "^>" "${RESULTS_DIR}/reference_sequences/all_references.fa" 2>/dev/null | tail -1 || echo 0)
ref_count=${ref_count:-0}
if [ "${ref_count}" -eq 0 ]; then
    log "Error: No reference sequences found in all_references.fa"
    exit 1
fi

log "Processing ${ref_count} total reference sequences..."

run_command "ref_id_map" python3 "${SCRIPTS_DIR}/update_headers.py" \
    "${RESULTS_DIR}/reference_sequences/all_references.fa" "${ID_MAP}" \
    --source-type reference
mv "${RESULTS_DIR}/reference_sequences/all_references.fa_updated.fa" "${RESULTS_DIR}/reference_sequences/all_references.fa" || { log "Error: Failed to move updated reference file"; exit 1; }

# --- Step 6: Update CDS headers to match protein headers ---
if [ -f "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" ]; then
    log "Updating CDS headers to match protein ID mapping..."
    python3 "${SCRIPTS_DIR}/update_headers.py" \
        "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" "${ID_MAP}" \
        --source-type reference --append
    mv "${RESULTS_DIR}/reference_sequences/cds/all_references_cds_updated.fna" \
       "${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna" 2>/dev/null || true
fi

# Create completion flag
touch "${RESULTS_DIR}/step_completed_01.txt"
log "Reference processing completed. Processed ${ref_count} sequences."
