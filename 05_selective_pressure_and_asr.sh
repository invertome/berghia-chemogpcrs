#!/bin/bash
# 05_selective_pressure_and_asr.sh
# Purpose: Analyze selective pressure and reconstruct ancestral sequences.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

#SBATCH --job-name=selective_pressure_asr
#SBATCH --output=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.out
#SBATCH --error=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --cpus-per-task=8  # Cap for FastML
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --array=0-999%100  # Increased concurrency
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh

mkdir -p "${RESULTS_DIR}/selective_pressure" "${RESULTS_DIR}/asr" "${LOGS_DIR}"

if [ ! -f "${RESULTS_DIR}/step_completed_lse_classification.txt" ]; then
    log "Error: LSE classification step not completed."
    exit 1
fi

log "Starting selective pressure and ASR analysis."

ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results*/Orthogroups/OG"*.fa)
if [ ${#ORTHOGROUPS[@]} -eq 0 ]; then
    log "Error: No orthogroups found."
    exit 1
fi

og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$og" ] || [ ! -f "$og" ]; then
    log "Skipping empty or missing orthogroup at index ${SLURM_ARRAY_TASK_ID}."
    exit 0
fi

base=$(basename "$og" .fa")
taxids=$(awk '{print substr($1, 1, index($1, "_")-1)}' "$og" | sort -u)
taxa_count=$(echo "$taxids" | wc -l)

if [ "$taxa_count" -gt 1 ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    nuc_align="${RESULTS_DIR}/phylogenies/nucleotide/${base}_trimmed.fa"
    tree="${RESULTS_DIR}/phylogenies/nucleotide/${base}.treefile"
    run_command "${base}_codon" pal2nal.pl "$protein_align" "$nuc_align" -output paml > "${RESULTS_DIR}/selective_pressure/${base}_codon.phy"
    cat > "${RESULTS_DIR}/selective_pressure/${base}_codeml_null.ctl" <<EOF
seqfile = ${RESULTS_DIR}/selective_pressure/${base}_codon.phy
treefile = $tree
outfile = ${RESULTS_DIR}/selective_pressure/${base}_codeml_null.out
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
EOF
    run_command "${base}_codeml_null" ${CODEML} "${RESULTS_DIR}/selective_pressure/${base}_codeml_null.ctl"
    cat > "${RESULTS_DIR}/selective_pressure/${base}_codeml_alt.ctl" <<EOF
seqfile = ${RESULTS_DIR}/selective_pressure/${base}_codon.phy
treefile = $tree
outfile = ${RESULTS_DIR}/selective_pressure/${base}_codeml_alt.out
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 2
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
EOF
    run_command "${base}_codeml_alt" ${CODEML} "${RESULTS_DIR}/selective_pressure/${base}_codeml_alt.ctl"
    python3 "${SCRIPTS_DIR}/compute_lrt.py" "$base" "${RESULTS_DIR}/selective_pressure/${base}_codeml_null.out" "${RESULTS_DIR}/selective_pressure/${base}_codeml_alt.out" "${RESULTS_DIR}/selective_pressure/lrt_results.csv"
    python3 "${SCRIPTS_DIR}/plot_selective_pressure.py" "${RESULTS_DIR}/selective_pressure/lrt_results.csv" "${RESULTS_DIR}/selective_pressure/lrt_plot_${base}.png"
fi

if [ "$(echo "$taxids" | grep -c "${BERGHIA_TAXID}")" -eq "$taxa_count" ] && [ "$taxa_count" -eq 1 ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    nuc_align="${RESULTS_DIR}/phylogenies/nucleotide/${base}_trimmed.fa"
    tree="${RESULTS_DIR}/phylogenies/nucleotide/${base}.treefile"
    run_command "${base}_codon_lse" pal2nal.pl "$protein_align" "$nuc_align" -output paml > "${RESULTS_DIR}/selective_pressure/${base}_codon_lse.phy"
    cat > "${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_null.ctl" <<EOF
seqfile = ${RESULTS_DIR}/selective_pressure/${base}_codon_lse.phy
treefile = $tree
outfile = ${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_null.out
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
EOF
    run_command "${base}_codeml_lse_null" ${CODEML} "${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_null.ctl"
    cat > "${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_alt.ctl" <<EOF
seqfile = ${RESULTS_DIR}/selective_pressure/${base}_codon_lse.phy
treefile = $tree
outfile = ${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_alt.out
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 2
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
EOF
    run_command "${base}_codeml_lse_alt" ${CODEML} "${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_alt.ctl"
    python3 "${SCRIPTS_DIR}/compute_lrt.py" "$base" "${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_null.out" "${RESULTS_DIR}/selective_pressure/${base}_codeml_lse_alt.out" "${RESULTS_DIR}/selective_pressure/lrt_results_lse.csv"
    python3 "${SCRIPTS_DIR}/plot_selective_pressure.py" "${RESULTS_DIR}/selective_pressure/lrt_results_lse.csv" "${RESULTS_DIR}/selective_pressure/lrt_plot_lse_${base}.png"

    deep_nodes=$(python3 "${SCRIPTS_DIR}/select_deep_nodes.py" "$tree" "${BERGHIA_TAXID}")
    for node in $deep_nodes; do
        run_command "${base}_${node}_asr" ${FASTML} --seq "$nuc_align" --tree "$tree" --out_seq "${RESULTS_DIR}/asr/${base}_${node}_asr.fa" --out_tree "${RESULTS_DIR}/asr/${base}_${node}_asr.tree" --node "$node" -t 8 --verbose
    done
    if [ -f "${RESULTS_DIR}/asr/${base}_${deep_nodes%% *}_asr.fa" ]; then
        python3 "${SCRIPTS_DIR}/plot_asr.py" "$tree" "${RESULTS_DIR}/asr/${base}_${deep_nodes%% *}_asr.fa" "${RESULTS_DIR}/asr/${base}_asr_plot.png"
    fi
fi

log "Selective pressure and ASR analysis completed for ${base}."
