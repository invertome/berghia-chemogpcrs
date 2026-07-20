#!/bin/bash
# 08_structural_analysis.sh
# Purpose: Predict structures with AlphaFold, build structural phylogenies including references, and compare with sequence trees.
# Inputs: Ranked candidates from step 07, reference structures from GPCRdb
# Outputs: Predicted structures in ${RESULTS_DIR}/structural_analysis/alphafold/, structural trees
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=structural_analysis
#SBATCH --output=${LOGS_DIR}/08_structural_analysis_%j.out
#SBATCH --error=${LOGS_DIR}/08_structural_analysis_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
# Bead -05o (2026-05-20): switched to GPU partition because AF3 inference
# requires CUDA (Ampere+ for Triton flash-attention; XLA fallback is
# portable but ~3x slower). CPU-only AF3 is impractical (hours per
# protein) — the wrapper aborts if no GPU is allocated.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

# --- Collect stage 05's ancestral-sequence reconstructions ---
# Arguments: $1 - ASR directory (${RESULTS_DIR}/asr), $2 - output FASTA
# Always creates $2 (possibly empty) and returns 0: ASR is optional, so a
# missing/empty results/asr/ must never fail the stage.
#
# Why not the bare `cat "${dir}"/*_asr.fa > out 2>/dev/null` idiom (bead 444):
# with the default (nullglob off) an unmatched glob is passed through
# literally, so `cat` fails on a nonexistent '*_asr.fa' path while the
# redirection has already truncated the output. find keeps the no-match case
# honest and empty.
collect_asr_sequences() {
    local asr_dir="$1"
    local out="$2"
    local -a asr_files=()

    if [ -d "$asr_dir" ]; then
        while IFS= read -r f; do
            [ -n "$f" ] && asr_files+=("$f")
        done < <(find "$asr_dir" -maxdepth 1 -type f -name '*_asr.fa' 2>/dev/null | sort)
    fi

    # Truncate unconditionally so a re-run cannot inherit a previous run's ASR.
    : > "$out"
    if [ "${#asr_files[@]}" -eq 0 ]; then
        return 0
    fi
    cat "${asr_files[@]}" >> "$out"
    return 0
}

# --- Build the sequence source that seqtk subseq extracts from ---
# Arguments: $1 - extant candidate FASTA, $2 - ASR FASTA, $3 - output FASTA
# Returns 0 on success; logs ERROR and returns 1 on empty input or duplicate ids.
#
# Why (bead 444): top_ids.txt mixes extant candidate ids with ancestral node
# ids, but the extraction previously read only the extant candidate FASTA.
# Ancestral reconstructions are by definition absent from it, so seqtk subseq
# silently dropped every ASR id, the per-candidate `grep -A1` produced a
# zero-byte input and hit `continue`, and nothing ancestral was ever folded.
# Concatenating both sources makes ids from both resolve. A duplicate id here
# would be silent corruption (seqtk emits both records; the per-candidate grep
# folds whichever comes first), hence the fail-loud guard.
build_structural_seq_source() {
    local extant="$1"
    local asr="$2"
    local out="$3"

    # Inline rather than via assert_fasta_non_empty: that helper is in-flight
    # work in functions.sh, and this stage should not depend on it landing.
    if [ ! -f "$extant" ]; then
        log --level=ERROR "extant candidate FASTA not found: ${extant}"
        return 1
    fi
    local n_extant
    n_extant=$(grep -c '^>' "$extant") || true
    if [ "${n_extant:-0}" -eq 0 ]; then
        log --level=ERROR "extant candidate FASTA ${extant} contains 0 FASTA records - refusing to continue"
        return 1
    fi

    cat "$extant" > "$out"
    # A candidate FASTA with no trailing newline would otherwise fuse its last
    # sequence line onto the first ASR header, corrupting both records.
    if [ -n "$(tail -c 1 "$out")" ]; then
        printf '\n' >> "$out"
    fi

    # ASR is optional. Require an actual record, not just a non-empty file.
    if [ -s "$asr" ] && grep -q '^>' "$asr"; then
        local n_asr
        n_asr=$(grep -c '^>' "$asr") || true
        cat "$asr" >> "$out"
        log "Added ${n_asr} ancestral (ASR) sequence(s) to the structural extraction source"
    else
        log "No ASR sequences available - folding extant candidates only"
    fi

    assert_no_duplicate_fasta_ids "$out" "structural extraction source" || return 1
    return 0
}

# --- Count what will ACTUALLY be folded ---
# Arguments: $1 - extracted FASTA (top_seqs.fa), $2 - requested id list
# Sets FOLDABLE_COUNT / REQUESTED_COUNT; warns when ids failed to resolve.
# Reporting the id-list length as the candidate count was how the dropped ASR
# ids stayed invisible: they inflated the total while never being folded.
count_foldable_sequences() {
    local top_seqs="$1"
    local top_ids="$2"
    local requested folded

    requested=$(grep -c '[^[:space:]]' "$top_ids" 2>/dev/null) || true
    folded=$(grep -c '^>' "$top_seqs" 2>/dev/null) || true
    REQUESTED_COUNT=${requested:-0}
    FOLDABLE_COUNT=${folded:-0}

    if [ "${FOLDABLE_COUNT}" -lt "${REQUESTED_COUNT}" ]; then
        log --level=WARN "Only ${FOLDABLE_COUNT} of ${REQUESTED_COUNT} requested id(s) resolved in the extraction source - the rest cannot be folded"
    fi
    return 0
}

# Create output directories
mkdir -p "${RESULTS_DIR}/structural_analysis/alphafold" "${RESULTS_DIR}/structural_analysis/references" "${RESULTS_DIR}/structural_analysis/all_pdb" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# Check dependency (step 07 creates step_completed_07.txt)
check_file "${RESULTS_DIR}/step_completed_07.txt"

log "Starting structural analysis."

# --- Phase 6: Select top N candidates for AlphaFold based on ranking ---
log "Selecting top ${NUM_STRUCTURAL_CANDIDATES} candidates for AlphaFold..."

# Method 1: Use ranked candidates directly (top N by rank score)
if [ -f "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" ]; then
    # Extract top N candidate IDs (skip header, get first column)
    head -n $((NUM_STRUCTURAL_CANDIDATES + 1)) "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" | \
        tail -n ${NUM_STRUCTURAL_CANDIDATES} | \
        cut -d',' -f1 > "${RESULTS_DIR}/structural_analysis/candidates_for_alphafold.txt"
    log "Selected ${NUM_STRUCTURAL_CANDIDATES} top-ranked candidates"
else
    log "Warning: Ranked candidates not found, using diverse selection"
    # Fallback: use diverse selection method
    run_command "select_diverse" python3 "${SCRIPTS_DIR}/select_diverse_candidates.py" \
        "${RESULTS_DIR}/phylogenies/protein/class_A/class_A.treefile" \
        "${RESULTS_DIR}/ranking/ranked_candidates_sorted.csv" \
        "${NUM_STRUCTURAL_CANDIDATES}" \
        "${RESULTS_DIR}/structural_analysis/candidates_for_alphafold.txt"
fi

# Create top_ids.txt from selected candidates (for compatibility)
cp "${RESULTS_DIR}/structural_analysis/candidates_for_alphafold.txt" \
   "${RESULTS_DIR}/structural_analysis/top_ids.txt"

# --- Include ASR sequences ---
# Ancestral reconstructions are folded alongside the extant candidates so the
# structural phylogeny can be compared ancestor-vs-extant. Stage 05 may not
# have run, so this is best-effort and never fatal.
ASR_SEQS="${RESULTS_DIR}/structural_analysis/asr_seqs.fa"
collect_asr_sequences "${RESULTS_DIR}/asr" "${ASR_SEQS}"
if [ -s "${ASR_SEQS}" ]; then
    # Bare id only (header up to the first whitespace), matching how
    # assert_no_duplicate_fasta_ids and the per-candidate `grep -A1 "^>id$"`
    # below identify records.
    grep "^>" "${ASR_SEQS}" | sed 's/^>//; s/[[:space:]].*//' \
        >> "${RESULTS_DIR}/structural_analysis/top_ids.txt"
fi

# --- Deduplicate IDs ---
sort -u "${RESULTS_DIR}/structural_analysis/top_ids.txt" > "${RESULTS_DIR}/structural_analysis/top_ids_unique.txt"
mv "${RESULTS_DIR}/structural_analysis/top_ids_unique.txt" "${RESULTS_DIR}/structural_analysis/top_ids.txt"

# --- Build the extraction source: extant candidates + ASR reconstructions ---
# top_ids.txt now mixes both id namespaces, so the source seqtk reads must
# cover both or the ASR ids resolve to nothing (bead 444).
SEQ_SOURCE="${RESULTS_DIR}/structural_analysis/extraction_source.fa"
build_structural_seq_source \
    "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" \
    "${ASR_SEQS}" \
    "${SEQ_SOURCE}" \
    || { log --level=ERROR "Could not build the structural extraction source"; exit 1; }

# --- Extract sequences for AlphaFold ---
run_command "seqtk_top" --stdout="${RESULTS_DIR}/structural_analysis/top_seqs.fa" ${SEQTK} subseq "${SEQ_SOURCE}" "${RESULTS_DIR}/structural_analysis/top_ids.txt"

# Count what will actually be folded (not merely what was requested)
count_foldable_sequences "${RESULTS_DIR}/structural_analysis/top_seqs.fa" \
                         "${RESULTS_DIR}/structural_analysis/top_ids.txt"
final_count=${FOLDABLE_COUNT}
log "Total candidates for structural prediction: ${final_count}"

# --- Predict structures with AlphaFold 3 (only if not already done) ---
# AF3 emits mmCIF (*_model.cif). foldseek (and therefore FoldTree) ingests
# mmCIF natively, so NO CIF->PDB conversion is performed; the CIF is fed
# straight into the structural-phylogeny step alongside reference structures.
log "Running AlphaFold 3 predictions..."

# AF3 needs the Google-gated model weights AND the sequence/structure
# databases. Abort loudly if either is missing rather than silently
# producing no structures.
if [ -z "${ALPHAFOLD3_MODEL_DIR}" ] || [ ! -e "${ALPHAFOLD3_MODEL_DIR}/af3.bin" ]; then
    log --level=ERROR "ALPHAFOLD3_MODEL_DIR='${ALPHAFOLD3_MODEL_DIR}' has no af3.bin — cannot run AlphaFold 3. Set it in config.sh (see scripts/unity/download_af3_weights.sh)."
    exit 1
fi
if [ ! -d "${ALPHAFOLD3_DB_DIR}" ]; then
    log --level=ERROR "ALPHAFOLD3_DB_DIR='${ALPHAFOLD3_DB_DIR}' not found — cannot run AlphaFold 3."
    exit 1
fi

# AF3 ships as an apptainer module on the cluster (provides the `af3` helper
# + run_alphafold.py; GPU via --nv). Keep the image cache off the home quota.
module load alphafold3/latest 2>/dev/null || log --level=WARN "could not 'module load alphafold3/latest'; assuming 'af3' is already on PATH"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${RESULTS_DIR}/.apptainer_cache}"
mkdir -p "${APPTAINER_CACHEDIR}"

while read -r candidate_id; do
    cand_dir="${RESULTS_DIR}/structural_analysis/alphafold/${candidate_id}"

    # Skip if a model CIF already exists for this candidate
    if find "${cand_dir}" -name '*_model.cif' 2>/dev/null | grep -q .; then
        log "  Skipping ${candidate_id} (structure exists)"
        continue
    fi
    mkdir -p "${cand_dir}"

    # Extract the single candidate sequence
    grep -A1 "^>${candidate_id}$" "${RESULTS_DIR}/structural_analysis/top_seqs.fa" > \
        "${cand_dir}/${candidate_id}_input.fasta" 2>/dev/null || continue
    [ -s "${cand_dir}/${candidate_id}_input.fasta" ] || continue

    log "  Running AlphaFold 3 for ${candidate_id}..."
    # 1) FASTA -> AF3 input JSON (via the cluster `af3` helper)
    ${AF3} convert --fasta "${cand_dir}/${candidate_id}_input.fasta" \
                --output_dir "${cand_dir}/json" \
        || { log --level=WARN "  af3 convert failed for ${candidate_id}"; continue; }
    af3_json=$(find "${cand_dir}/json" -maxdepth 1 -name '*.json' 2>/dev/null | head -1)
    [ -n "${af3_json}" ] || { log --level=WARN "  no AF3 JSON produced for ${candidate_id}"; continue; }

    # 2) Structure inference (data pipeline + model; GPU via the module's --nv)
    ${AF3} run --json_path "${af3_json}" \
            --output_dir "${cand_dir}" \
            --model_dir "${ALPHAFOLD3_MODEL_DIR}" \
            --db_dir "${ALPHAFOLD3_DB_DIR}" \
        || log --level=WARN "  AlphaFold 3 failed for ${candidate_id}"
done < "${RESULTS_DIR}/structural_analysis/top_ids.txt"

# --- Fetch reference structures and metadata ---
run_command "fetch_references" python3 "${SCRIPTS_DIR}/fetch_ligands.py" \
    "${RESULTS_DIR}/structural_analysis/references" \
    "${RESULTS_DIR}/structural_analysis/reference_ligands.csv" \
    "${GPCRDB_SEARCH_TERMS}" \
    "${GPCRDB_SPECIES}"

# --- Combine predicted (mmCIF) and reference structures for phylogeny ---
# foldseek/FoldTree read mmCIF and PDB from the same directory, so AF3 CIFs
# and GPCRdb reference structures are pooled without any conversion.
# AlphaFold 3 predictions (mmCIF)
alphafold_struct_count=$(find "${RESULTS_DIR}/structural_analysis/alphafold/" -name "*_model.cif" 2>/dev/null | wc -l)
if [ "$alphafold_struct_count" -gt 0 ]; then
    find "${RESULTS_DIR}/structural_analysis/alphafold/" -name "*_model.cif" -exec cp {} "${RESULTS_DIR}/structural_analysis/all_pdb/" \;
    log "Copied ${alphafold_struct_count} AlphaFold 3 structures (mmCIF)"
else
    log --level=WARN "No AlphaFold 3 structures found - structural phylogeny may be limited"
fi

# Reference structures (PDB and/or mmCIF)
ref_pdb_count=$(find "${RESULTS_DIR}/structural_analysis/references/" \( -name "*.pdb" -o -name "*.cif" \) 2>/dev/null | wc -l)
if [ "$ref_pdb_count" -gt 0 ]; then
    find "${RESULTS_DIR}/structural_analysis/references/" \( -name "*.pdb" -o -name "*.cif" \) -exec cp {} "${RESULTS_DIR}/structural_analysis/all_pdb/" \;
    log "Copied ${ref_pdb_count} reference structures"
else
    log --level=WARN "No reference structures found - structural comparison may be limited"
fi

# Verify we have enough structures for phylogeny
total_struct_count=$(find "${RESULTS_DIR}/structural_analysis/all_pdb/" \( -name "*.pdb" -o -name "*.cif" \) 2>/dev/null | wc -l)
if [ "$total_struct_count" -lt 3 ]; then
    log --level=WARN "Only ${total_struct_count} structures available - need at least 3 for meaningful structural phylogeny"
fi

# --- Render publication figures of the predicted structures ---
# PyMOL (headless OSMesa software render via the `pymol` skill) draws each AF3
# prediction as a cartoon colored by pLDDT confidence (mmCIF B-factor), then a
# montage grid is assembled for slides/manuscripts. STRUCT_RENDERER defaults to
# `uv run scripts/render_structure_figure.py` (uv pulls pymol-open-source-whl).
fig_dir="${RESULTS_DIR}/structural_analysis/figures"
mkdir -p "${fig_dir}"
export PATH="${HOME}/.local/bin:${PATH}"          # uv default install location
export UV_CACHE_DIR="${UV_CACHE_DIR:-${RESULTS_DIR}/.uv_cache}"
rendered_pngs=()
while IFS= read -r cif; do
    [ -n "$cif" ] || continue
    cand=$(basename "$(dirname "$(dirname "$cif")")")   # alphafold/<cand>/<name>/<name>_model.cif
    png="${fig_dir}/${cand}.png"
    if ${STRUCT_RENDERER} "$cif" "$png" --session "${fig_dir}/${cand}.pse" --color-mode plddt; then
        [ -s "$png" ] && rendered_pngs+=("$png")
    else
        log --level=WARN "  structure figure render failed for ${cand}"
    fi
done < <(find "${RESULTS_DIR}/structural_analysis/alphafold/" -name "*_model.cif" 2>/dev/null)

if [ "${#rendered_pngs[@]}" -gt 0 ]; then
    python3 "${SCRIPTS_DIR}/montage_structure_figures.py" \
        --out "${fig_dir}/all_structures_montage.png" "${rendered_pngs[@]}" \
        || log --level=WARN "  structure-figure montage assembly failed"
    log "Rendered ${#rendered_pngs[@]} structure figures (color-by-pLDDT) -> ${fig_dir}"
else
    log --level=WARN "No structure figures rendered (no AF3 predictions or renderer unavailable)"
fi

# --- Build structural phylogeny with FoldTree ---
run_command "foldtree" ${FOLDTREE} --input_dir "${RESULTS_DIR}/structural_analysis/all_pdb" --output "${RESULTS_DIR}/structural_analysis/foldtree.tre" --method "${FOLDTREE_METHOD}"

# --- Compare structural and sequence trees ---
python3 "${SCRIPTS_DIR}/plot_struct_vs_seq.py" "${RESULTS_DIR}/structural_analysis/foldtree.tre" "${RESULTS_DIR}/phylogenies/protein/class_A/class_A.treefile" "${RESULTS_DIR}/structural_analysis/struct_vs_seq_plot"

touch "${RESULTS_DIR}/step_completed_foldtree.txt"
log "Structural analysis completed."
