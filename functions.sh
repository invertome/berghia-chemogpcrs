#!/bin/bash
# functions.sh
# Purpose: Define reusable functions for the pipeline, including file checking, logging, and command execution.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst.

# --- Check if a file exists ---
# Arguments: $1 - File path to check
# Returns: Exits with error if file is missing
check_file() {
    local file="$1"
    if [ ! -f "$file" ]; then
        log "Error: File $file not found."
        echo "Error: File $file not found." >&2
        exit 1
    fi
}

# --- Log messages to file ---
# Arguments: $1 - Message to log
# Returns: Appends message with timestamp to pipeline.log
log() {
    local message="$1"
    echo "$(date): $message" >> "${LOGS_DIR}/pipeline.log"
}

# --- Run a command with error handling and checkpointing ---
# Arguments: $1 - Step name for checkpointing, $2+ - Command and its arguments
# Returns: Executes command, logs output/errors, and creates checkpoint on success
run_command() {
    local step="$1"
    shift
    if [ ! -f "${RESULTS_DIR}/step_completed_${step}.txt" ]; then
        log "Running command for $step: $*"
        "$@" > "${LOGS_DIR}/${step}.log" 2>> "${LOGS_DIR}/${step}.err"
        local exit_code=$?
        if [ $exit_code -ne 0 ]; then
            log "Error: Command failed for $step with exit code $exit_code. Check ${LOGS_DIR}/${step}.err for details."
            echo "Error: Command failed - $*" >&2
            cat "${LOGS_DIR}/${step}.err" >&2
            exit 1
        fi
        log "Completed $step successfully."
        touch "${RESULTS_DIR}/step_completed_${step}.txt"
    else
        log "Step $step already completed, skipping."
    fi
}

# --- Run commands in parallel ---
# Arguments: $1 - Command to run, $2 - Input list file, $3 - Number of CPUs (optional)
# Returns: Executes command in parallel across inputs
run_parallel() {
    local cmd="$1"
    local input_list="$2"
    local cpus="${3:-$DEFAULT_CPUS}"
    log "Running parallel command: $cmd with $cpus jobs"
    parallel --jobs "$cpus" "$cmd" :::: "$input_list" 2>> "${LOGS_DIR}/parallel.err"
    if [ $? -ne 0 ]; then
        log "Error: Parallel execution failed. Check ${LOGS_DIR}/parallel.err"
        exit 1
    fi
}

# --- Scale resources based on SLURM node limits ---
# Returns: SLURM resource string with CPU and memory allocation
scale_resources() {
    local max_cpus=$(sinfo -N -o "%C" | grep -m 1 "[0-9]" | awk '{print $1}' | cut -d'/' -f2)
    local max_mem=$(sinfo -N -o "%m" | grep -m 1 "[0-9]" | awk '{print $1}')
    CPUS=$((max_cpus < DEFAULT_CPUS ? max_cpus : DEFAULT_CPUS))
    MEM=$((max_mem / 1024 < DEFAULT_MEM ? max_mem / 1024 : DEFAULT_MEM))
    echo "--cpus-per-task=$CPUS --mem=${MEM}G"
}
