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
cat "${RESULTS_DIR}/asr/"*_asr.fa > "${RESULTS_DIR}/structural_analysis/asr_seqs.fa" 2>/dev/null
grep "^>" "${RESULTS_DIR}/structural_analysis/asr_seqs.fa" | sed 's/>//' >> "${RESULTS_DIR}/structural_analysis/top_ids.txt" 2>/dev/null

# --- Deduplicate IDs ---
sort -u "${RESULTS_DIR}/structural_analysis/top_ids.txt" > "${RESULTS_DIR}/structural_analysis/top_ids_unique.txt"
mv "${RESULTS_DIR}/structural_analysis/top_ids_unique.txt" "${RESULTS_DIR}/structural_analysis/top_ids.txt"

# Count final candidates
final_count=$(wc -l < "${RESULTS_DIR}/structural_analysis/top_ids.txt")
log "Total candidates for structural prediction: ${final_count}"

# --- Extract sequences for AlphaFold ---
run_command "seqtk_top" --stdout="${RESULTS_DIR}/structural_analysis/top_seqs.fa" ${SEQTK} subseq "${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa" "${RESULTS_DIR}/structural_analysis/top_ids.txt"

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
