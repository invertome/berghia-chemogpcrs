#!/bin/bash
# resume_array_job.sh
# Purpose: Resume incomplete SLURM array jobs by identifying failed tasks
# Usage: ./resume_array_job.sh <step_name> <script> [--dry-run]
# Example: ./resume_array_job.sh 04_phylo 04_phylogenetic_analysis.sh --dry-run
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

source config.sh
source functions.sh

# Parse arguments
STEP_NAME="${1:-}"
SCRIPT="${2:-}"
DRY_RUN=false

if [ "$3" = "--dry-run" ]; then
    DRY_RUN=true
fi

# Validate arguments
if [ -z "$STEP_NAME" ] || [ -z "$SCRIPT" ]; then
    echo "Usage: $0 <step_name> <script> [--dry-run]"
    echo ""
    echo "Examples:"
    echo "  $0 04_phylo 04_phylogenetic_analysis.sh          # Resume step 04 phylogenetic analysis"
    echo "  $0 05_selective 05_selective_pressure_and_asr.sh # Resume step 05 selective pressure"
    echo "  $0 04_phylo 04_phylogenetic_analysis.sh --dry-run  # Show what would be resubmitted"
    echo ""
    echo "Available checkpoints:"
    if [ -d "${RESULTS_DIR}/checkpoints" ]; then
        ls -d "${RESULTS_DIR}/checkpoints"/*/ 2>/dev/null | xargs -n1 basename
    else
        echo "  (none found)"
    fi
    exit 1
fi

# Get orthogroup count from manifest
MANIFEST_FILE="${RESULTS_DIR}/orthogroup_manifest.tsv"
OG_COUNT=$(get_orthogroup_count "$MANIFEST_FILE")

if [ "$OG_COUNT" -eq 0 ]; then
    echo "Error: No orthogroups found. Have you run step 03?"
    exit 1
fi

# Get incomplete task IDs
INCOMPLETE=$(get_incomplete_tasks "$STEP_NAME" "$OG_COUNT")

if [ -z "$INCOMPLETE" ]; then
    echo "All $OG_COUNT tasks for ${STEP_NAME} have completed."
    exit 0
fi

# Count incomplete tasks
INCOMPLETE_COUNT=$(echo "$INCOMPLETE" | tr ',' '\n' | wc -l)

echo "=== Array Job Resume Utility ==="
echo "Step: ${STEP_NAME}"
echo "Script: ${SCRIPT}"
echo "Total tasks: ${OG_COUNT}"
echo "Completed: $((OG_COUNT - INCOMPLETE_COUNT))"
echo "Incomplete: ${INCOMPLETE_COUNT}"
echo ""

# Show first few incomplete task IDs
PREVIEW=$(echo "$INCOMPLETE" | tr ',' '\n' | head -10 | tr '\n' ',' | sed 's/,$//')
if [ "$INCOMPLETE_COUNT" -gt 10 ]; then
    echo "Incomplete tasks (first 10): ${PREVIEW}..."
else
    echo "Incomplete tasks: ${INCOMPLETE}"
fi
echo ""

# Build sbatch command
SBATCH_CMD="sbatch --array=${INCOMPLETE} ${SCRIPT}"

if [ "$DRY_RUN" = true ]; then
    echo "[DRY RUN] Would execute:"
    echo "  $SBATCH_CMD"
    echo ""
    echo "To actually submit, run without --dry-run flag"
else
    echo "Submitting array job..."
    echo "  $SBATCH_CMD"
    $SBATCH_CMD
    echo ""
    echo "Job submitted. Monitor with: squeue -u \$USER"
fi
