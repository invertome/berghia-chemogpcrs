#!/bin/bash
# run_pipeline.sh
# Purpose: Submit pipeline steps to SLURM with proper configuration
# Usage: ./run_pipeline.sh [options] <step|all|range>
#
# Examples:
#   ./run_pipeline.sh 02                    # Submit step 02
#   ./run_pipeline.sh 02 03 04              # Submit steps 02, 03, 04
#   ./run_pipeline.sh 02-05                 # Submit steps 02 through 05
#   ./run_pipeline.sh all                   # Submit all steps sequentially
#   ./run_pipeline.sh --local 07            # Run step 07 locally (no SLURM)
#   ./run_pipeline.sh --dry-run all         # Show what would be submitted
#   ./run_pipeline.sh --wait 02 03          # Wait for 02 to complete before submitting 03
#
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

set -euo pipefail

# --- Source Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/functions.sh"

# --- Default Options ---
DRY_RUN=false
LOCAL_MODE=false
WAIT_MODE=false
VALIDATE_FIRST=true
VERBOSE=false

# --- Pipeline Steps Definition ---
# Format: step_number:script_name:description:is_array_job
declare -A PIPELINE_STEPS=(
    ["01"]="01_reference_processing.sh:Reference Processing:no"
    ["02"]="02_chemogpcrs_identification.sh:ChemoGPCR Identification:no"
    ["02a"]="02a_cluster_sequences.sh:Sequence Clustering:no"
    ["02b"]="02b_classify_gpcrs.sh:GPCR Classification:no"
    ["03"]="03_orthology_clustering.sh:Orthology Clustering:no"
    ["03a"]="03a_busco_species_tree.sh:BUSCO Species Tree:no"
    ["03b"]="03b_lse_classification.sh:LSE Classification:no"
    ["03c"]="03c_cafe_analysis.sh:CAFE5 Analysis:no"
    ["03d"]="03d_notung_reconciliation.sh:NOTUNG Reconciliation:no"
    ["04"]="04_phylogenetic_analysis.sh:Phylogenetic Analysis:yes"
    ["05"]="05_selective_pressure_and_asr.sh:Selective Pressure & ASR:yes"
    ["06"]="06_synteny_and_mapping.sh:Synteny & Mapping:no"
    ["07"]="07_candidate_ranking.sh:Candidate Ranking:no"
    ["08"]="08_structural_analysis.sh:Structural Analysis:no"
    ["09"]="09_report_generation.sh:Report Generation:no"
)

# Ordered list for sequential execution
STEP_ORDER=("01" "02" "02a" "02b" "03" "03a" "03b" "03c" "03d" "04" "05" "06" "07" "08" "09")

# --- Usage ---
usage() {
    cat << 'EOF'
Usage: ./run_pipeline.sh [options] <step|all|range>

Submit pipeline steps to SLURM with proper configuration from config.sh.

STEPS:
  01        Reference Processing
  02        ChemoGPCR Identification
  02a       Sequence Clustering (CD-HIT)
  02b       GPCR Classification (InterProScan)
  03        Orthology Clustering (OrthoFinder)
  03a       BUSCO Species Tree
  03b       LSE Classification
  03c       CAFE5 Analysis
  03d       NOTUNG Reconciliation
  04        Phylogenetic Analysis (array job)
  05        Selective Pressure & ASR (array job)
  06        Synteny & Mapping
  07        Candidate Ranking
  08        Structural Analysis
  09        Report Generation

OPTIONS:
  -l, --local         Run locally instead of submitting to SLURM
  -d, --dry-run       Show commands without executing
  -w, --wait          Wait for each job to complete before submitting next
  -n, --no-validate   Skip configuration validation
  -v, --verbose       Show detailed output
  -h, --help          Show this help message

EXAMPLES:
  ./run_pipeline.sh all                 # Run all steps sequentially
  ./run_pipeline.sh 02 03 04            # Submit specific steps
  ./run_pipeline.sh 02-05               # Submit range of steps
  ./run_pipeline.sh --wait 04 05        # Wait for 04 before starting 05
  ./run_pipeline.sh --local 07          # Run step 07 locally

SLURM CONFIGURATION (set in config.sh):
  SLURM_PARTITION     Partition/queue name
  SLURM_ACCOUNT       Account/allocation
  SLURM_QOS           Quality of service
  SLURM_CONSTRAINT    Node constraints
  SLURM_RESERVATION   Reservation name
  SLURM_EXTRA_ARGS    Additional sbatch flags

EOF
    exit 0
}

# --- Build SLURM Options ---
build_slurm_opts() {
    local opts=""

    # Add conditional options only if set
    [ -n "${SLURM_PARTITION:-}" ] && opts+=" --partition=${SLURM_PARTITION}"
    [ -n "${SLURM_ACCOUNT:-}" ] && opts+=" --account=${SLURM_ACCOUNT}"
    [ -n "${SLURM_QOS:-}" ] && opts+=" --qos=${SLURM_QOS}"
    [ -n "${SLURM_CONSTRAINT:-}" ] && opts+=" --constraint=${SLURM_CONSTRAINT}"
    [ -n "${SLURM_RESERVATION:-}" ] && opts+=" --reservation=${SLURM_RESERVATION}"
    [ -n "${SLURM_EXTRA_ARGS:-}" ] && opts+=" ${SLURM_EXTRA_ARGS}"

    echo "$opts"
}

# --- Submit a Step ---
submit_step() {
    local step="$1"
    local prev_job_id="${2:-}"

    # Get step info
    local step_info="${PIPELINE_STEPS[$step]:-}"
    if [ -z "$step_info" ]; then
        echo "Error: Unknown step '$step'" >&2
        echo "Valid steps: ${!PIPELINE_STEPS[*]}" >&2
        return 1
    fi

    IFS=':' read -r script_name description is_array <<< "$step_info"
    local script_path="${SCRIPT_DIR}/${script_name}"

    # Check script exists
    if [ ! -f "$script_path" ]; then
        echo "Error: Script not found: $script_path" >&2
        return 1
    fi

    # Build command
    local slurm_opts
    slurm_opts=$(build_slurm_opts)

    # Add dependency if waiting and previous job exists
    if [ "$WAIT_MODE" = true ] && [ -n "$prev_job_id" ]; then
        slurm_opts+=" --dependency=afterok:${prev_job_id}"
    fi

    if [ "$LOCAL_MODE" = true ]; then
        # Run locally
        echo "[$step] Running locally: $description" >&2
        if [ "$DRY_RUN" = true ]; then
            echo "  Would run: bash $script_path" >&2
        else
            bash "$script_path"
        fi
        echo "" >&2
    else
        # Submit to SLURM
        echo "[$step] Submitting: $description" >&2
        [ "$VERBOSE" = true ] && echo "  Script: $script_name" >&2
        [ "$VERBOSE" = true ] && echo "  Options: $slurm_opts" >&2

        if [ "$DRY_RUN" = true ]; then
            echo "  Would run: sbatch${slurm_opts} $script_path" >&2
            echo "12345"  # Fake job ID for dry run (stdout for capture)
        else
            local output
            output=$(sbatch${slurm_opts} "$script_path" 2>&1)
            local job_id
            job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+' || echo "")

            if [ -n "$job_id" ]; then
                echo "  Job ID: $job_id" >&2
                echo "$job_id"  # stdout for capture
            else
                echo "  Warning: Could not parse job ID from: $output" >&2
                echo ""
            fi
        fi
    fi
}

# --- Expand Step Range ---
expand_range() {
    local range="$1"
    local start end

    if [[ "$range" =~ ^([0-9]+[a-z]?)-([0-9]+[a-z]?)$ ]]; then
        start="${BASH_REMATCH[1]}"
        end="${BASH_REMATCH[2]}"

        local in_range=false
        local result=()

        for step in "${STEP_ORDER[@]}"; do
            if [ "$step" = "$start" ]; then
                in_range=true
            fi
            if [ "$in_range" = true ]; then
                result+=("$step")
            fi
            if [ "$step" = "$end" ]; then
                break
            fi
        done

        echo "${result[@]}"
    else
        echo "$range"
    fi
}

# --- Main ---
main() {
    local steps=()

    # Parse options
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                usage
                ;;
            -l|--local)
                LOCAL_MODE=true
                shift
                ;;
            -d|--dry-run)
                DRY_RUN=true
                shift
                ;;
            -w|--wait)
                WAIT_MODE=true
                shift
                ;;
            -n|--no-validate)
                VALIDATE_FIRST=false
                shift
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            all)
                steps=("${STEP_ORDER[@]}")
                shift
                ;;
            *-*)
                # Range like "02-05"
                expanded=$(expand_range "$1")
                read -ra expanded_arr <<< "$expanded"
                steps+=("${expanded_arr[@]}")
                shift
                ;;
            *)
                steps+=("$1")
                shift
                ;;
        esac
    done

    # Check if any steps specified
    if [ ${#steps[@]} -eq 0 ]; then
        echo "Error: No steps specified" >&2
        echo "Usage: ./run_pipeline.sh [options] <step|all|range>" >&2
        echo "Try './run_pipeline.sh --help' for more information." >&2
        exit 1
    fi

    # Validate configuration first
    if [ "$VALIDATE_FIRST" = true ] && [ -f "${SCRIPT_DIR}/validate_config.sh" ]; then
        echo "Validating configuration..."
        if ! "${SCRIPT_DIR}/validate_config.sh" --quiet; then
            echo "Configuration validation failed. Run './validate_config.sh' for details." >&2
            echo "Use '--no-validate' to skip this check." >&2
            exit 1
        fi
        echo "Configuration OK"
        echo ""
    fi

    # Show SLURM settings if verbose
    if [ "$VERBOSE" = true ]; then
        echo "SLURM Configuration:"
        echo "  Partition:   ${SLURM_PARTITION:-<not set>}"
        echo "  Account:     ${SLURM_ACCOUNT:-<not set>}"
        echo "  QOS:         ${SLURM_QOS:-<not set>}"
        echo "  Constraint:  ${SLURM_CONSTRAINT:-<not set>}"
        echo "  Time:        ${DEFAULT_TIME}"
        echo "  Memory:      ${DEFAULT_MEM}"
        echo "  CPUs:        ${CPUS}"
        echo ""
    fi

    # Submit steps
    local prev_job_id=""
    for step in "${steps[@]}"; do
        prev_job_id=$(submit_step "$step" "$prev_job_id")
    done

    echo ""
    if [ "$DRY_RUN" = true ]; then
        echo "Dry run complete. No jobs were submitted."
    elif [ "$LOCAL_MODE" = true ]; then
        echo "Local execution complete."
    else
        echo "All jobs submitted. Monitor with: squeue -u $USER"
    fi
}

main "$@"
