#!/bin/bash
# 03a_busco_species_tree.sh
# Purpose: Generate a species tree using BUSCO single-copy orthologs with IQ-TREE for gene trees.
# Inputs: Transcriptome files in ${TRANSCRIPTOME_DIR}/*.aa, GPCR output from step 02
# Outputs: Species tree in ${SPECIES_TREE}, gene trees in ${RESULTS_DIR}/busco/gene_trees/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=busco_species_tree
#SBATCH --output=${LOGS_DIR}/03a_busco_species_tree_%j.out
#SBATCH --error=${LOGS_DIR}/03a_busco_species_tree_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# Create output directories
mkdir -p "${RESULTS_DIR}/busco" "${RESULTS_DIR}/busco/single_copy" "${RESULTS_DIR}/busco/alignments" "${RESULTS_DIR}/busco/gene_trees" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency
check_file "${RESULTS_DIR}/step_completed_extract_berghia.txt"

log "Starting BUSCO species tree generation."

# --- Build dedup'd input-sample list ---
# Same realpath = same data (the transcriptomes/ dir has Berghia .aa as
# both the real file and a lowercase-symlink alias; running BUSCO on
# both wastes ~11 min and produces a fake "second species" in the
# combine loop). Resolve realpaths and skip duplicates. 2026-05-18.
declare -A _seen_realpath
unique_inputs=()
for trans in "${TRANSCRIPTOME_DIR}"/*.aa "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"; do
    [ -f "$trans" ] || continue
    real=$(realpath "$trans" 2>/dev/null) || continue
    if [ -z "${_seen_realpath[$real]:-}" ]; then
        _seen_realpath[$real]=1
        unique_inputs+=("$trans")
    fi
done
log "BUSCO inputs (after realpath dedup): ${#unique_inputs[@]} unique sample(s)"

# --- Run BUSCO on all unique transcriptomes ---
for trans in "${unique_inputs[@]}"; do
    sample=$(basename "$trans" .fa | sed 's/chemogpcrs_//')
    taxid_sample="${sample}"
    # `-m proteins`: input is amino-acid sequences (transcriptomes/*.aa and
    # chemogpcrs_*.fa are both protein). The prior `-m transcriptome` mode
    # required nucleotide input; passing .aa files made BUSCO abort with
    # "The input file does not contain nucleotide sequences" (stage 03a
    # job 57824881, 2026-05-18). Reference taxa only have .aa files
    # (miniprot-recovered protein translations), so proteins mode is the
    # correct choice for both Berghia and references.
    #
    # BUSCO's `-o` argument is a run NAME (no slashes allowed); the prior
    # full-path form `-o /scratch3/.../busco_<sample>` produced misleading
    # "A run with the name scratch3/... already exists" errors (job
    # 57824889). Use `--out_path` for the directory, `-o` for just the
    # run name.
    run_command "busco_${taxid_sample}" ${BUSCO} -i "$trans" \
        --out_path "${RESULTS_DIR}/busco" \
        -o "busco_${taxid_sample}" \
        -m proteins -l mollusca_odb10 -c "${CPUS}" -f
done

# --- Extract single-copy BUSCOs ---
for busco_dir in "${RESULTS_DIR}/busco/busco_"*; do
    [ -d "$busco_dir" ] || continue
    taxid_sample=$(basename "$busco_dir" | sed 's/busco_//')

    # Create per-sample directory for BUSCO sequences
    mkdir -p "${RESULTS_DIR}/busco/single_copy/${taxid_sample}"

    # Copy single-copy BUSCO sequences, prefixing with sample name.
    # BUSCO 5.x layout: run_<lineage>/busco_sequences/single_copy_busco_sequences/
    # (the `busco_sequences/` parent dir was added in BUSCO 5; older
    # versions had single_copy_busco_sequences/ directly under run_*).
    # Job 57824897 silently extracted 0 single-copies because the
    # `busco_sequences/` level was missing from the path. 2026-05-18.
    for busco_faa in "${busco_dir}/run_mollusca_odb10/busco_sequences/single_copy_busco_sequences/"*.faa; do
        [ -f "$busco_faa" ] || continue
        busco_name=$(basename "$busco_faa")
        cp "$busco_faa" "${RESULTS_DIR}/busco/single_copy/${taxid_sample}/${busco_name}"
    done

    [ -z "$(ls -A "${RESULTS_DIR}/busco/single_copy/${taxid_sample}" 2>/dev/null)" ] && log "Warning: No single-copy BUSCOs for ${taxid_sample}"
done

# --- Combine orthologous BUSCOs across species and align ---
# First, identify all unique BUSCO IDs present in multiple species
declare -A busco_files
for sample_dir in "${RESULTS_DIR}/busco/single_copy/"*/; do
    [ -d "$sample_dir" ] || continue
    sample=$(basename "$sample_dir")
    for busco_faa in "${sample_dir}"*.faa; do
        [ -f "$busco_faa" ] || continue
        busco_id=$(basename "$busco_faa" .faa)
        busco_files["$busco_id"]+="${busco_faa} "
    done
done

# Concatenate same-ID BUSCOs from different species and align
for busco_id in "${!busco_files[@]}"; do
    files=(${busco_files[$busco_id]})
    # Only process if present in at least 2 species
    if [ ${#files[@]} -ge 2 ]; then
        # Combine sequences, adding sample name to headers
        combined_file="${RESULTS_DIR}/busco/alignments/${busco_id}_combined.fa"
        > "$combined_file"
        for faa in "${files[@]}"; do
            sample=$(basename "$(dirname "$faa")")
            # Add sample prefix to sequence headers
            sed "s/^>/>${sample}_/" "$faa" >> "$combined_file"
        done

        # Bead -align: regime-based aligner. BUSCO per-gene alignments have
        # one sequence per species (~hundreds), so the wrapper picks
        # MAFFT L-INS-i (<200) or MAFFT --auto (200-999) automatically.
        run_command "align_busco_${busco_id}" bash "${SCRIPTS_DIR}/run_aligner.sh" \
            --input="$combined_file" \
            --output="${RESULTS_DIR}/busco/alignments/${busco_id}.afa" \
            --threads="${CPUS}"
        run_command "trim_busco_${busco_id}" ${TRIMAL} -in "${RESULTS_DIR}/busco/alignments/${busco_id}.afa" -out "${RESULTS_DIR}/busco/alignments/${busco_id}_trimmed.afa" -automated1
    fi
done

# --- Infer gene trees with IQ-TREE 3 ---
# Three bugs from IQ-TREE 1.x carryover (caught in stage 03a audit, 2026-05-08):
#   - `-m "${IQTREE_MODEL}"` collapsed the composite "MFP -mset ..." into a
#     single quoted argument that IQ-TREE 3 rejects as an unknown model.
#     Pass the discriminating flags separately like stage 04 does.
#   - `-nt` (IQ-TREE 1.x threads) → `-T` (IQ-TREE 3.x).
#   - `-pre` (IQ-TREE 1.x prefix) → `--prefix` (IQ-TREE 3.x).
# Also force `-st AA` and run the defensive near-all-gap filter to match
# stages 04 and the classification tree builder (avoids "Unknown sequence
# type" on TrimAl outputs where one species is entirely missing).
shopt -s nullglob
trimmed_alns=("${RESULTS_DIR}/busco/alignments/"*_trimmed.afa)
shopt -u nullglob
# Graceful degradation conditions:
#   (a) Zero multi-species BUSCO alignments produced (combine loop saw
#       no busco_id in >=2 samples).
#   (b) <3 distinct sample inputs in total. Even if (a) didn't fire (the
#       combine loop saw 2 samples that happen to be the same species —
#       e.g., chemogpcrs_berghia.fa shares a BUSCO with transcriptomes/
#       *.aa), there's no meaningful tree to build from 2 species.
# The species tree is only consumed by 03c (CAFE) and 03d (Notung),
# which aren't on the chemoreceptor-ranking critical path (stages 04-09
# don't use SPECIES_TREE). Mark 03a complete so downstream gates can
# proceed; flag the limitation explicitly. Bead -? (P3): rebuild with
# full reference proteomes from RefSeq.
if [ "${#trimmed_alns[@]}" -eq 0 ] || [ "${#unique_inputs[@]}" -lt 3 ]; then
    log --level=WARN "03a: degenerate species tree case (${#unique_inputs[@]} unique sample(s), ${#trimmed_alns[@]} multi-sample alignments). Skipping IQ-TREE + ASTRAL; species tree empty. OK for the chemoreceptor pipeline (gates 4-7); only blocks CAFE/Notung (03c/03d)."
    touch "${RESULTS_DIR}/step_completed_busco_species_tree.txt"
    log "BUSCO species tree generation completed (degenerate mode — insufficient species)."
    exit 0
fi
for aln in "${trimmed_alns[@]}"; do
    base=$(basename "$aln" _trimmed.afa)
    python3 "${SCRIPTS_DIR}/drop_near_all_gap_rows.py" \
        --input "$aln" --output "$aln" \
        2>>"${LOGS_DIR}/busco_${base}_drop_gappy.log" || true
    run_command "tree_busco_${base}" ${IQTREE} \
        -s "$aln" \
        -st AA \
        -m "${IQTREE_MODEL_FIND}" -msub "${IQTREE_MSUB}" -mset "${IQTREE_MODEL_SET}" \
        -B "${IQTREE_BOOTSTRAP}" -alrt 1000 \
        -seed "${IQTREE_SEED}" \
        -T "${CPUS}" \
        --prefix "${RESULTS_DIR}/busco/gene_trees/${base}"
done

# --- Generate species tree with ASTRAL ---
find "${RESULTS_DIR}/busco/gene_trees/" -name "*.treefile" > "${RESULTS_DIR}/busco/gene_trees_list.txt" || { log "Error: No gene trees found"; exit 1; }
run_command "astral_species_tree" ${ASTRAL} -i "${RESULTS_DIR}/busco/gene_trees_list.txt" -o "${SPECIES_TREE}"

touch "${RESULTS_DIR}/step_completed_busco_species_tree.txt"
log "BUSCO species tree generation completed."
