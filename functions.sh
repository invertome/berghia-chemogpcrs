#!/bin/bash
# functions.sh
# Purpose: Define reusable functions for the pipeline (logging, command execution, file checks).
# Usage: Source this file in all pipeline scripts (source functions.sh).

# --- Logging Function ---
# Arguments: $1 - Message to log
log() {
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    echo "[$timestamp] $1" >> "${LOGS_DIR}/pipeline.log"
    echo "[$timestamp] $1"
}

# --- Command Execution Function ---
# Arguments: $1 - Command name, $2+ - Command and arguments
run_command() {
    local name="$1"
    shift
    log "Running: $*"
    "$@" > "${LOGS_DIR}/${name}.log" 2> "${LOGS_DIR}/${name}.err"
    if [ $? -ne 0 ]; then
        log "Error: Command failed: $* (see ${LOGS_DIR}/${name}.err)"
        exit 1
    fi
    touch "${RESULTS_DIR}/step_completed_${name}.txt"
}

# --- File Check Function ---
# Arguments: $1+ - File paths
check_file() {
    for file in "$@"; do
        if [ ! -f "$file" ]; then
            log "Error: File not found: $file"
            exit 1
        fi
    done
}

# --- Resource Scaling Function ---
# Returns: SLURM resource directives based on CPUS and DEFAULT_MEM
scale_resources() {
    echo "--cpus-per-task=${CPUS} --mem=${DEFAULT_MEM}"
}
