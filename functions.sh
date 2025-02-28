#!/bin/bash
# functions.sh
# Purpose: Define reusable functions for logging, file checking, and command execution across the pipeline.
# Notes: Sourced by all scripts to standardize common operations.

# Log messages with timestamps to both file and stdout
log() {
    local message="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" >> "${LOGS_DIR}/pipeline.log"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message"
}

# Check if a file exists and is non-empty, exit with error if not
check_file() {
    local file="$1"
    if [ ! -s "$file" ]; then
        log "Error: File $file is missing or empty."
        exit 1
    fi
}

# Run a command with error handling and logging
# Arguments: job_name, command, command_args...
run_command() {
    local job_name="$1"
    shift
    log "Running: $*"
    "$@" 2>> "${LOGS_DIR}/${job_name}.err" >> "${LOGS_DIR}/${job_name}.out"
    if [ $? -ne 0 ]; then
        log "Error: Command failed - $*"
        exit 1
    fi
    touch "${RESULTS_DIR}/step_completed_${job_name}.txt"  # Mark step completion
}

# Scale SLURM resources dynamically based on config settings
scale_resources() {
    echo "--cpus-per-task=${DEFAULT_CPUS} --mem=${DEFAULT_MEM}"
}
