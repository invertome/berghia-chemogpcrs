#!/bin/bash
# functions.sh
# Purpose: Define reusable functions for the pipeline (logging, command execution,
#          checkpointing, provenance tracking, and resource management).
# Usage: Source this file in all pipeline scripts (source functions.sh).

# --- Initialize Pipeline Environment ---
init_pipeline() {
    # Create required directories
    mkdir -p "${LOGS_DIR}" "${RESULTS_DIR}/checkpoints" "${RESULTS_DIR}/provenance"

    # Initialize provenance file for this run
    export PIPELINE_RUN_ID=$(date "+%Y%m%d_%H%M%S")_$$
    export PROVENANCE_FILE="${RESULTS_DIR}/provenance/run_${PIPELINE_RUN_ID}.json"

    # Start provenance tracking
    cat > "$PROVENANCE_FILE" <<EOF
{
  "run_id": "${PIPELINE_RUN_ID}",
  "start_time": "$(date -Iseconds)",
  "hostname": "$(hostname)",
  "user": "${USER}",
  "working_directory": "${BASE_DIR}",
  "pipeline_version": "$(get_pipeline_version)",
  "environment": {
    "CPUS": ${CPUS:-1},
    "DEFAULT_MEM": "${DEFAULT_MEM:-4G}",
    "SLURM_JOB_ID": "${SLURM_JOB_ID:-local}"
  },
  "tool_versions": $(get_tool_versions_json),
  "steps": []
}
EOF

    log "Pipeline initialized: run_id=${PIPELINE_RUN_ID}"
}

# --- Get Pipeline Version ---
get_pipeline_version() {
    if [ -d "${BASE_DIR}/.git" ]; then
        git -C "${BASE_DIR}" describe --tags --always 2>/dev/null || echo "unknown"
    else
        echo "unknown"
    fi
}

# --- Get Tool Versions as JSON ---
get_tool_versions_json() {
    local versions="{"
    local first=true

    # Check each tool and get version
    declare -A tool_cmds=(
        ["hmmsearch"]="hmmsearch -h | head -2 | tail -1"
        ["hhblits"]="hhblits -h 2>&1 | head -1"
        ["mafft"]="mafft --version 2>&1"
        ["iqtree"]="iqtree2 --version 2>&1 | head -1"
        ["fasttree"]="FastTree 2>&1 | head -1"
        ["orthofinder"]="orthofinder -h 2>&1 | grep -i version | head -1"
        ["python"]="python3 --version 2>&1"
        ["R"]="R --version 2>&1 | head -1"
    )

    for tool in "${!tool_cmds[@]}"; do
        if command -v "${tool%% *}" &>/dev/null; then
            local ver=$(eval "${tool_cmds[$tool]}" 2>/dev/null | tr -d '\n' | tr '"' "'" | head -c 100)
            if [ "$first" = true ]; then
                first=false
            else
                versions+=","
            fi
            versions+="\"$tool\":\"$ver\""
        fi
    done

    versions+="}"
    echo "$versions"
}

# --- Logging Function ---
# Arguments: $1 - Message to log
# Options: --level=INFO|WARN|ERROR
log() {
    local level="INFO"
    local message="$1"

    # Parse options
    if [[ "$1" == --level=* ]]; then
        level="${1#--level=}"
        shift
        message="$1"
    fi

    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    local log_entry="[$timestamp] [$level] $message"

    # Write to log file
    echo "$log_entry" >> "${LOGS_DIR}/pipeline.log"

    # Write to console with color
    case "$level" in
        ERROR) echo -e "\033[0;31m$log_entry\033[0m" >&2 ;;
        WARN)  echo -e "\033[0;33m$log_entry\033[0m" ;;
        *)     echo "$log_entry" ;;
    esac
}

# --- Compute File Checksum ---
# Arguments: $1 - File path
# Returns: MD5 checksum or "missing" if file doesn't exist
compute_checksum() {
    local file="$1"
    if [ -f "$file" ]; then
        md5sum "$file" 2>/dev/null | cut -d' ' -f1
    else
        echo "missing"
    fi
}

# --- Atomic Write Function ---
# Arguments: $1 - Target file, $2 - Content
# Writes content atomically (write to temp, then rename)
# Uses SLURM job info + hostname + PID for unique temp files in array jobs
atomic_write() {
    local target="$1"
    local content="$2"
    # Create unique temp file name to avoid race conditions in SLURM array jobs
    # Different nodes could have same PID, so include hostname and SLURM task ID
    local unique_id="${SLURM_ARRAY_TASK_ID:-0}_${SLURM_JOB_ID:-local}_$(hostname -s 2>/dev/null || echo local)_$$"
    local temp_file="${target}.tmp.${unique_id}"

    echo "$content" > "$temp_file"
    mv "$temp_file" "$target"
    # Clean up any stale temp files from previous failed runs (older than 1 hour)
    find "$(dirname "$target")" -name "$(basename "$target").tmp.*" -mmin +60 -delete 2>/dev/null || true
}

# --- Create Checkpoint ---
# Arguments: $1 - Step name
# Creates a checkpoint with metadata for resume capability
create_checkpoint() {
    local step_name="$1"
    local checkpoint_dir="${RESULTS_DIR}/checkpoints"
    local checkpoint_file="${checkpoint_dir}/${step_name}.checkpoint"

    # Create checkpoint with metadata
    local checkpoint_data=$(cat <<EOF
{
    "step": "${step_name}",
    "completed_at": "$(date -Iseconds)",
    "run_id": "${PIPELINE_RUN_ID:-unknown}",
    "duration_seconds": ${STEP_DURATION:-0},
    "exit_code": 0
}
EOF
)

    atomic_write "$checkpoint_file" "$checkpoint_data"

    # Also create legacy marker for backwards compatibility
    touch "${RESULTS_DIR}/step_completed_${step_name}.txt"

    log "Checkpoint created: ${step_name}"
}

# --- Check if Step Completed ---
# Arguments: $1 - Step name
# Returns: 0 if completed, 1 if not
is_step_completed() {
    local step_name="$1"
    local checkpoint_file="${RESULTS_DIR}/checkpoints/${step_name}.checkpoint"
    local legacy_file="${RESULTS_DIR}/step_completed_${step_name}.txt"

    [ -f "$checkpoint_file" ] || [ -f "$legacy_file" ]
}

# --- Skip if Completed ---
# Arguments: $1 - Step name
# Exits with 0 if step already completed (for resume capability)
skip_if_completed() {
    local step_name="$1"

    if is_step_completed "$step_name"; then
        log "Skipping ${step_name}: already completed (use --force to rerun)"
        exit 0
    fi
}

# --- Record Step in Provenance ---
# Arguments: $1 - Step name, $2 - Command, $3 - Input files (comma-separated), $4 - Output files (comma-separated)
record_provenance() {
    local step_name="$1"
    local command="$2"
    local inputs="$3"
    local outputs="$4"

    [ -z "$PROVENANCE_FILE" ] && return

    # Save IFS to restore after parsing comma-separated values
    local OLD_IFS="$IFS"

    # Compute input checksums
    local input_checksums="{"
    local first=true
    local INPUT_FILES
    IFS=',' read -ra INPUT_FILES <<< "$inputs"
    IFS="$OLD_IFS"  # Restore IFS immediately after use
    for f in "${INPUT_FILES[@]}"; do
        f=$(echo "$f" | xargs)  # Trim whitespace
        if [ -n "$f" ]; then
            [ "$first" = false ] && input_checksums+=","
            input_checksums+="\"$f\":\"$(compute_checksum "$f")\""
            first=false
        fi
    done
    input_checksums+="}"

    # Compute output checksums
    local output_checksums="{"
    first=true
    local OUTPUT_FILES
    IFS=',' read -ra OUTPUT_FILES <<< "$outputs"
    IFS="$OLD_IFS"  # Restore IFS immediately after use
    for f in "${OUTPUT_FILES[@]}"; do
        f=$(echo "$f" | xargs)
        if [ -n "$f" ]; then
            [ "$first" = false ] && output_checksums+=","
            output_checksums+="\"$f\":\"$(compute_checksum "$f")\""
            first=false
        fi
    done
    output_checksums+="}"

    # Create step record
    local step_record=$(cat <<EOF
{
    "name": "${step_name}",
    "command": "$(echo "$command" | tr '"' "'" | head -c 500)",
    "started_at": "$(date -Iseconds)",
    "input_checksums": ${input_checksums},
    "output_checksums": ${output_checksums}
}
EOF
)

    # Append to provenance file (this is simplified; production would use jq)
    log "Provenance recorded for: ${step_name}"
}

# --- Command Execution Function (Enhanced) ---
# Arguments: $1 - Command name, $2+ - Command and arguments
# Options: --inputs=file1,file2 --outputs=file1,file2 --allow-fail
#          --stdout=FILE  Redirect command stdout to FILE instead of log
#                         (use for commands like seqtk/mafft/fasttree that
#                         write results to stdout)
run_command() {
    local name="$1"
    shift

    local inputs=""
    local outputs=""
    local allow_fail=false
    local stdout_file=""

    # Parse options
    while [[ "$1" == --* ]]; do
        case "$1" in
            --inputs=*) inputs="${1#--inputs=}"; shift ;;
            --outputs=*) outputs="${1#--outputs=}"; shift ;;
            --allow-fail) allow_fail=true; shift ;;
            --stdout=*) stdout_file="${1#--stdout=}"; shift ;;
            *) break ;;
        esac
    done

    local start_time=$(date +%s)
    log "Running: $name"
    log "  Command: $*"

    # Record in provenance
    record_provenance "$name" "$*" "$inputs" "$outputs"

    # Execute command — redirect stdout to file if --stdout given, else to log
    if [ -n "$stdout_file" ]; then
        "$@" > "$stdout_file" 2> "${LOGS_DIR}/${name}.err"
    else
        "$@" > "${LOGS_DIR}/${name}.log" 2> "${LOGS_DIR}/${name}.err"
    fi
    local exit_code=$?

    local end_time=$(date +%s)
    export STEP_DURATION=$((end_time - start_time))

    if [ $exit_code -ne 0 ]; then
        log --level=ERROR "Command failed: $name (exit code: $exit_code)"
        log --level=ERROR "  See: ${LOGS_DIR}/${name}.err"

        if [ "$allow_fail" = false ]; then
            exit 1
        fi
        return $exit_code
    fi

    # Create checkpoint
    create_checkpoint "$name"
    log "Completed: $name (${STEP_DURATION}s)"
}

# --- File Check Function (Enhanced) ---
# Arguments: $1+ - File paths
# Options: --warn-only (don't exit on missing file)
check_file() {
    local warn_only=false

    if [ "$1" = "--warn-only" ]; then
        warn_only=true
        shift
    fi

    for file in "$@"; do
        if [ ! -f "$file" ]; then
            if [ "$warn_only" = true ]; then
                log --level=WARN "File not found (non-critical): $file"
            else
                log --level=ERROR "File not found: $file"
                exit 1
            fi
        fi
    done
}

# --- Check Directory ---
# Arguments: $1+ - Directory paths
check_dir() {
    for dir in "$@"; do
        if [ ! -d "$dir" ]; then
            log --level=ERROR "Directory not found: $dir"
            exit 1
        fi
    done
}

# --- Detect Available Resources ---
# Sets DETECTED_CPUS and DETECTED_MEM based on system
detect_resources() {
    # Detect CPUs
    if [ -n "$SLURM_CPUS_PER_TASK" ]; then
        export DETECTED_CPUS=$SLURM_CPUS_PER_TASK
    elif [ -f /proc/cpuinfo ]; then
        export DETECTED_CPUS=$(grep -c ^processor /proc/cpuinfo)
    else
        export DETECTED_CPUS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
    fi

    # Detect Memory (in GB)
    if [ -n "$SLURM_MEM_PER_NODE" ]; then
        # SLURM memory is in MB
        export DETECTED_MEM_GB=$((SLURM_MEM_PER_NODE / 1024))
    elif [ -f /proc/meminfo ]; then
        local mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        export DETECTED_MEM_GB=$((mem_kb / 1024 / 1024))
    else
        export DETECTED_MEM_GB=$(sysctl -n hw.memsize 2>/dev/null | awk '{print int($1/1024/1024/1024)}' || echo 8)
    fi

    log "Detected resources: ${DETECTED_CPUS} CPUs, ${DETECTED_MEM_GB}GB RAM"
}

# --- Resource Scaling Function (Enhanced) ---
# Arguments: $1 - Scale factor (small/medium/large/auto)
# Returns: SLURM resource directives based on input size or explicit scale
scale_resources() {
    local scale="${1:-auto}"
    local cpus=$CPUS
    local mem=$DEFAULT_MEM

    case "$scale" in
        small)
            cpus=$((CPUS / 4))
            mem="16G"
            ;;
        medium)
            cpus=$((CPUS / 2))
            mem="32G"
            ;;
        large)
            cpus=$CPUS
            mem=$DEFAULT_MEM
            ;;
        auto)
            # Use detected resources if available
            if [ -n "$DETECTED_CPUS" ]; then
                cpus=$DETECTED_CPUS
            fi
            if [ -n "$DETECTED_MEM_GB" ]; then
                mem="${DETECTED_MEM_GB}G"
            fi
            ;;
    esac

    # Ensure minimum values
    [ $cpus -lt 1 ] && cpus=1

    echo "--cpus-per-task=${cpus} --mem=${mem}"
}

# --- Generate Conditional SLURM Options ---
# Returns SLURM directives for partition, account, qos, etc. only if they are set
# Usage: In scripts, call this function and eval the output, or use in sbatch command
get_slurm_opts() {
    local opts=""

    # Add partition if set
    [ -n "$SLURM_PARTITION" ] && opts+=" --partition=${SLURM_PARTITION}"

    # Add account if set
    [ -n "$SLURM_ACCOUNT" ] && opts+=" --account=${SLURM_ACCOUNT}"

    # Add QOS if set
    [ -n "$SLURM_QOS" ] && opts+=" --qos=${SLURM_QOS}"

    # Add constraint if set
    [ -n "$SLURM_CONSTRAINT" ] && opts+=" --constraint=${SLURM_CONSTRAINT}"

    # Add reservation if set
    [ -n "$SLURM_RESERVATION" ] && opts+=" --reservation=${SLURM_RESERVATION}"

    # Add any extra args
    [ -n "$SLURM_EXTRA_ARGS" ] && opts+=" ${SLURM_EXTRA_ARGS}"

    echo "$opts"
}

# --- Generate SBATCH Header Block ---
# Arguments: $1 - job name, $2 - log prefix (optional, defaults to job name)
#            $3 - array spec (optional, e.g., "0-999")
# Prints a complete SBATCH header block for sourcing or embedding
generate_sbatch_header() {
    local job_name="$1"
    local log_prefix="${2:-$1}"
    local array_spec="$3"

    echo "#!/bin/bash"
    echo "#SBATCH --job-name=${job_name}"

    # Handle array job logs differently
    if [ -n "$array_spec" ]; then
        echo "#SBATCH --output=${LOGS_DIR}/${log_prefix}_%j_%a.out"
        echo "#SBATCH --error=${LOGS_DIR}/${log_prefix}_%j_%a.err"
        echo "#SBATCH --array=${array_spec}%${SLURM_ARRAY_LIMIT:-50}"
    else
        echo "#SBATCH --output=${LOGS_DIR}/${log_prefix}_%j.out"
        echo "#SBATCH --error=${LOGS_DIR}/${log_prefix}_%j.err"
    fi

    echo "#SBATCH --time=${DEFAULT_TIME}"
    echo "#SBATCH --cpus-per-task=${CPUS}"
    echo "#SBATCH --mem=${DEFAULT_MEM}"

    # Conditional options (only if set)
    [ -n "$SLURM_PARTITION" ] && echo "#SBATCH --partition=${SLURM_PARTITION}"
    [ -n "$SLURM_ACCOUNT" ] && echo "#SBATCH --account=${SLURM_ACCOUNT}"
    [ -n "$SLURM_QOS" ] && echo "#SBATCH --qos=${SLURM_QOS}"
    [ -n "$SLURM_CONSTRAINT" ] && echo "#SBATCH --constraint=${SLURM_CONSTRAINT}"
    [ -n "$SLURM_RESERVATION" ] && echo "#SBATCH --reservation=${SLURM_RESERVATION}"

    # Mail notifications
    [ -n "$SLURM_EMAIL" ] && echo "#SBATCH --mail-type=ALL" && echo "#SBATCH --mail-user=${SLURM_EMAIL}"

    # Extra args (parse into separate lines if multiple)
    if [ -n "$SLURM_EXTRA_ARGS" ]; then
        for arg in $SLURM_EXTRA_ARGS; do
            echo "#SBATCH ${arg}"
        done
    fi
}

# --- Check if Running in SLURM ---
is_slurm() {
    [ -n "$SLURM_JOB_ID" ]
}

# --- Run Locally or via SLURM ---
# Arguments: $1 - Script path, $2+ - Script arguments
# Uses --local flag to force local execution
run_step() {
    local script="$1"
    shift

    local force_local=false
    if [ "$1" = "--local" ]; then
        force_local=true
        shift
    fi

    if [ "$force_local" = true ] || [ "${RUN_MODE:-slurm}" = "local" ]; then
        log "Running locally: $script"
        bash "$script" "$@"
    else
        log "Submitting to SLURM: $script"
        sbatch "$script" "$@"
    fi
}

# --- Wait for File with Timeout ---
# Arguments: $1 - File path, $2 - Timeout in seconds (default 3600)
wait_for_file() {
    local file="$1"
    local timeout="${2:-3600}"
    local elapsed=0
    local interval=10

    while [ ! -f "$file" ] && [ $elapsed -lt $timeout ]; do
        sleep $interval
        elapsed=$((elapsed + interval))
    done

    if [ ! -f "$file" ]; then
        log --level=ERROR "Timeout waiting for file: $file"
        return 1
    fi
    return 0
}

# --- Cleanup Temporary Files ---
# Arguments: $1 - Pattern to match (e.g., "*.tmp")
cleanup_temp() {
    local pattern="${1:-*.tmp}"
    find "${RESULTS_DIR}" -name "$pattern" -type f -mmin +60 -delete 2>/dev/null
    log "Cleaned up temporary files matching: $pattern"
}

# --- Finalize Pipeline Run ---
# Guard variable to prevent re-entrance
_PIPELINE_FINALIZED=false

finalize_pipeline() {
    # Re-entrance guard: prevent multiple calls during cleanup
    if [ "$_PIPELINE_FINALIZED" = true ]; then
        return
    fi
    _PIPELINE_FINALIZED=true

    local exit_code="${1:-0}"

    if [ -n "${PROVENANCE_FILE:-}" ] && [ -f "${PROVENANCE_FILE:-}" ]; then
        # Update provenance with end time
        local end_time=$(date -Iseconds)
        # Note: In production, use jq to properly update JSON
        log "Pipeline finalized: exit_code=${exit_code}"
    fi

    # Cleanup old temp files
    cleanup_temp
}

# --- Trap for Clean Exit ---
trap 'finalize_pipeline $?' EXIT

# =============================================================================
# METADATA LOOKUP FUNCTIONS
# =============================================================================
# These functions provide centralized sequence metadata lookups, eliminating
# the need for fragile FASTA header parsing throughout the pipeline.

# --- Get Taxid for Sequence ID ---
# Arguments: $1 - Sequence ID
# Returns: Taxid string
get_taxid_for_seq() {
    local seq_id="$1"
    local metadata_file="${ID_MAP:-${RESULTS_DIR}/reference_sequences/id_map.csv}"

    if [ -f "$metadata_file" ]; then
        python3 "${SCRIPTS_DIR}/get_metadata.py" --metadata "$metadata_file" --seq-id "$seq_id" --field taxid
    else
        # Fallback to header parsing if metadata file doesn't exist
        if [[ "$seq_id" == ref_* ]]; then
            echo "$seq_id" | cut -d'_' -f2
        else
            echo "$seq_id" | cut -d'_' -f1
        fi
    fi
}

# --- Get Source Type for Sequence ID ---
# Arguments: $1 - Sequence ID
# Returns: source_type (reference/target/outgroup/unknown)
get_source_type_for_seq() {
    local seq_id="$1"
    local metadata_file="${ID_MAP:-${RESULTS_DIR}/reference_sequences/id_map.csv}"

    if [ -f "$metadata_file" ]; then
        python3 "${SCRIPTS_DIR}/get_metadata.py" --metadata "$metadata_file" --seq-id "$seq_id" --field source_type
    else
        # Fallback
        if [[ "$seq_id" == ref_* ]]; then
            echo "reference"
        else
            echo "unknown"
        fi
    fi
}

# --- Check if Sequence is Reference ---
# Arguments: $1 - Sequence ID
# Returns: 0 if reference, 1 otherwise
is_reference_seq() {
    local seq_id="$1"
    local source_type=$(get_source_type_for_seq "$seq_id")
    [ "$source_type" = "reference" ]
}

# --- Get Unique Taxids from FASTA ---
# Arguments: $1 - FASTA file, $2 - (optional) --exclude-refs
# Returns: Space-separated list of unique taxids
get_taxids_from_fasta() {
    local fasta_file="$1"
    local exclude_refs="${2:-}"
    local metadata_file="${ID_MAP:-${RESULTS_DIR}/reference_sequences/id_map.csv}"

    if [ -f "$metadata_file" ]; then
        if [ "$exclude_refs" = "--exclude-refs" ]; then
            python3 "${SCRIPTS_DIR}/get_metadata.py" --metadata "$metadata_file" --fasta "$fasta_file" --unique-taxids --exclude-refs
        else
            python3 "${SCRIPTS_DIR}/get_metadata.py" --metadata "$metadata_file" --fasta "$fasta_file" --unique-taxids
        fi
    else
        # Fallback: parse headers directly
        grep "^>" "$fasta_file" | sed 's/>//' | while read -r header; do
            if [[ "$header" == ref_* ]]; then
                echo "$header" | cut -d'_' -f2
            else
                echo "$header" | cut -d'_' -f1
            fi
        done | sort -u | tr '\n' ' ' | sed 's/ $//'
    fi
}

# --- Get Non-Reference Taxids from FASTA ---
# Arguments: $1 - FASTA file
# Returns: Space-separated list of unique taxids (excluding references)
get_non_ref_taxids_from_fasta() {
    get_taxids_from_fasta "$1" --exclude-refs
}

# =============================================================================
# RESOURCE ESTIMATION FUNCTIONS
# =============================================================================
# These functions estimate memory requirements based on data size to prevent
# OOM failures on large datasets.

# Tool-specific memory multipliers (bytes per unit)
# These can be overridden in config.sh
: ${MEM_MULT_IQTREE:=100}        # bytes per site per taxon
: ${MEM_MULT_ORTHOFINDER:=50}    # bytes per sequence pair
: ${MEM_MULT_MAFFT:=20}          # bytes per residue pair
: ${MEM_MULT_HYPHY:=200}         # bytes per site per branch
: ${MEM_MULT_FASTTREE:=50}       # bytes per site per taxon
: ${MEM_BASE_GB:=4}              # base memory overhead in GB
: ${MEM_MAX_GB:=128}             # maximum memory to request

# --- Count Sequences in FASTA ---
# Arguments: $1 - FASTA file
# Returns: Number of sequences
count_sequences() {
    local fasta="$1"
    grep -c "^>" "$fasta" 2>/dev/null || echo 0
}

# --- Get Average Sequence Length ---
# Arguments: $1 - FASTA file
# Returns: Average sequence length (integer)
get_avg_seq_length() {
    local fasta="$1"
    awk '/^>/{if(seq)print length(seq);seq=""} !/^>/{seq=seq$0} END{if(seq)print length(seq)}' "$fasta" 2>/dev/null | \
        awk '{s+=$1;n++}END{if(n>0)print int(s/n);else print 0}'
}

# --- Get Alignment Length ---
# Arguments: $1 - Alignment file (FASTA format)
# Returns: Alignment length (number of columns), always returns "0" for empty/missing files
get_alignment_length() {
    local aln="$1"

    # Return 0 if file doesn't exist or is empty
    if [ ! -f "$aln" ] || [ ! -s "$aln" ]; then
        echo 0
        return
    fi

    # Get length of first sequence (all should be same length in alignment)
    local len
    len=$(awk '/^>/{if(seq){print length(seq);exit}seq=""} !/^>/{seq=seq$0} END{if(seq)print length(seq)}' "$aln" 2>/dev/null)

    # Ensure we always return a valid number
    if [ -z "$len" ] || ! [[ "$len" =~ ^[0-9]+$ ]]; then
        echo 0
    else
        echo "$len"
    fi
}

# --- Estimate Memory for Alignment ---
# Arguments: $1 - FASTA file
# Returns: Estimated memory in GB (e.g., "8G")
estimate_memory_for_alignment() {
    local fasta="$1"
    local num_seqs=$(count_sequences "$fasta")
    local avg_len=$(get_avg_seq_length "$fasta")

    # Use awk for safe arithmetic to prevent integer overflow on large datasets
    # Distance matrix: O(n^2) with ~10 bytes per cell
    local mem_gb=$(awk -v n="$num_seqs" -v l="$avg_len" -v mult="$MEM_MULT_MAFFT" -v base="$MEM_BASE_GB" -v max="$MEM_MAX_GB" '
        BEGIN {
            mem_bytes = (n * n * 10) + (n * l * mult)
            mem_gb = int(mem_bytes / 1073741824) + base
            if (mem_gb > max) mem_gb = max
            if (mem_gb < base) mem_gb = base
            print mem_gb
        }
    ')

    echo "${mem_gb}G"
}

# --- Estimate Memory for Tree Building ---
# Arguments: $1 - Alignment file
# Returns: Estimated memory in GB (e.g., "16G")
estimate_memory_for_tree() {
    local aln="$1"
    local num_seqs=$(count_sequences "$aln")
    local aln_len=$(get_alignment_length "$aln")

    # Use awk for safe arithmetic to prevent integer overflow on large datasets
    # IQ-TREE memory: sites * taxa * multiplier
    local mem_gb=$(awk -v n="$num_seqs" -v l="$aln_len" -v mult="$MEM_MULT_IQTREE" -v base="$MEM_BASE_GB" -v max="$MEM_MAX_GB" '
        BEGIN {
            mem_bytes = n * l * mult
            mem_gb = int(mem_bytes / 1073741824) + base
            if (mem_gb > max) mem_gb = max
            if (mem_gb < base) mem_gb = base
            print mem_gb
        }
    ')

    echo "${mem_gb}G"
}

# --- Estimate Memory for OrthoFinder ---
# Arguments: $1 - Directory with proteomes or total sequence count
# Returns: Estimated memory in GB
estimate_memory_for_orthofinder() {
    local input="$1"
    local total_seqs=0

    if [ -d "$input" ]; then
        # Count sequences across all FASTA files in directory
        # Use nullglob to handle case where no files match
        local old_nullglob=$(shopt -p nullglob)
        shopt -s nullglob
        for f in "$input"/*.fa "$input"/*.faa "$input"/*.fasta; do
            [ -f "$f" ] && total_seqs=$((total_seqs + $(count_sequences "$f")))
        done
        eval "$old_nullglob"  # Restore previous nullglob setting
    else
        # Assume input is sequence count
        total_seqs="$input"
    fi

    # Use awk for safe arithmetic to prevent integer overflow on large datasets
    # OrthoFinder: O(n^2) for all-vs-all DIAMOND
    local mem_gb=$(awk -v n="$total_seqs" -v mult="$MEM_MULT_ORTHOFINDER" -v base="$MEM_BASE_GB" -v max="$MEM_MAX_GB" '
        BEGIN {
            mem_bytes = n * n * mult
            mem_gb = int(mem_bytes / 1073741824) + base
            if (mem_gb > max) mem_gb = max
            if (mem_gb < 8) mem_gb = 8  # OrthoFinder needs at least 8GB
            print mem_gb
        }
    ')

    echo "${mem_gb}G"
}

# --- Estimate Memory for HyPhy ---
# Arguments: $1 - Alignment file, $2 - Tree file (optional, for branch count)
# Returns: Estimated memory in GB
estimate_memory_for_hyphy() {
    local aln="$1"
    local tree="${2:-}"
    local num_seqs=$(count_sequences "$aln")
    local aln_len=$(get_alignment_length "$aln")

    # Estimate branches: 2n-2 for binary tree
    local num_branches=$((2 * num_seqs - 2))

    # Use awk for safe arithmetic to prevent integer overflow on large datasets
    # HyPhy aBSREL: sites * branches * multiplier
    local mem_gb=$(awk -v l="$aln_len" -v b="$num_branches" -v mult="$MEM_MULT_HYPHY" -v base="$MEM_BASE_GB" -v max="$MEM_MAX_GB" '
        BEGIN {
            mem_bytes = l * b * mult
            mem_gb = int(mem_bytes / 1073741824) + base
            if (mem_gb > max) mem_gb = max
            if (mem_gb < base) mem_gb = base
            print mem_gb
        }
    ')

    echo "${mem_gb}G"
}

# --- Check Resource Requirements (Pre-flight) ---
# Arguments: $1 - Input file, $2 - Tool name (alignment/tree/orthofinder/hyphy)
# Returns: 0 if OK, 1 if insufficient resources (with warning)
check_resource_requirements() {
    local input="$1"
    local tool="$2"
    local estimated=""

    case "$tool" in
        alignment|mafft)
            estimated=$(estimate_memory_for_alignment "$input")
            ;;
        tree|iqtree|fasttree)
            estimated=$(estimate_memory_for_tree "$input")
            ;;
        orthofinder)
            estimated=$(estimate_memory_for_orthofinder "$input")
            ;;
        hyphy|absrel)
            estimated=$(estimate_memory_for_hyphy "$input")
            ;;
        *)
            log --level=WARN "Unknown tool for resource estimation: $tool"
            return 0
            ;;
    esac

    local estimated_gb=${estimated%G}
    local available_gb=${DETECTED_MEM_GB:-64}

    if [ "$estimated_gb" -gt "$available_gb" ]; then
        log --level=WARN "Resource warning for $tool:"
        log --level=WARN "  Estimated memory: ${estimated}"
        log --level=WARN "  Available memory: ${available_gb}G"
        log --level=WARN "  Consider: (1) reducing input size with CD-HIT"
        log --level=WARN "            (2) requesting more memory via SLURM"
        log --level=WARN "            (3) using a high-memory node"
        return 1
    fi

    if [ "$VERBOSE" = true ] || [ "${DEBUG:-}" = true ]; then
        log "Resource check for $tool: ${estimated} estimated, ${available_gb}G available"
    fi

    return 0
}

# --- Data-Aware Resource Scaling ---
# Arguments: $1 - Input file, $2 - Tool name
# Returns: SLURM resource directives
scale_resources_for_data() {
    local input="$1"
    local tool="$2"
    local estimated=""

    case "$tool" in
        alignment|mafft)
            estimated=$(estimate_memory_for_alignment "$input")
            ;;
        tree|iqtree|fasttree)
            estimated=$(estimate_memory_for_tree "$input")
            ;;
        orthofinder)
            estimated=$(estimate_memory_for_orthofinder "$input")
            ;;
        hyphy|absrel)
            estimated=$(estimate_memory_for_hyphy "$input")
            ;;
        *)
            estimated="${DEFAULT_MEM:-32G}"
            ;;
    esac

    echo "--cpus-per-task=${CPUS:-8} --mem=${estimated}"
}

# --- Get Dataset Statistics ---
# Arguments: $1 - FASTA file or directory
# Returns: Prints statistics to stdout
get_dataset_stats() {
    local input="$1"

    if [ -f "$input" ]; then
        local num_seqs=$(count_sequences "$input")
        local avg_len=$(get_avg_seq_length "$input")
        local total_residues=$((num_seqs * avg_len))
        echo "Sequences: $num_seqs"
        echo "Avg length: $avg_len"
        echo "Total residues: $total_residues"
    elif [ -d "$input" ]; then
        local total_seqs=0
        local num_files=0
        for f in "$input"/*.fa "$input"/*.faa "$input"/*.fasta; do
            if [ -f "$f" ]; then
                total_seqs=$((total_seqs + $(count_sequences "$f")))
                num_files=$((num_files + 1))
            fi
        done
        echo "Files: $num_files"
        echo "Total sequences: $total_seqs"
    fi
}

# =============================================================================
# STEP INITIALIZATION AND VALIDATION
# =============================================================================
# These functions reduce boilerplate in pipeline step scripts.

# --- Initialize Pipeline Step ---
# Arguments: $1 - Step number (e.g., "01", "03b")
#            $2 - Step name (e.g., "reference_processing")
#            $3+ - Prerequisite step numbers (e.g., "01" "02")
# Creates output directories, checks prerequisites, and sets up logging.
# Usage: init_step "04" "phylogenetic_analysis" "03"
init_step() {
    local step_num="$1"
    local step_name="$2"
    shift 2
    local prereqs=("$@")

    # Export step info for use by other functions
    export CURRENT_STEP_NUM="$step_num"
    export CURRENT_STEP_NAME="$step_name"

    # Create required directories
    mkdir -p "${LOGS_DIR}" "${RESULTS_DIR}" || {
        echo "Error: Cannot create directories" >&2
        exit 1
    }

    # Initialize pipeline if not already done
    [ -z "$PIPELINE_RUN_ID" ] && init_pipeline

    log "=== Step ${step_num}: ${step_name} ==="

    # Check prerequisites
    for prereq in "${prereqs[@]}"; do
        local prereq_file="${RESULTS_DIR}/step_completed_${prereq}.txt"
        if [ ! -f "$prereq_file" ]; then
            log --level=ERROR "Prerequisite not met: step ${prereq} has not completed"
            log --level=ERROR "Expected file: ${prereq_file}"
            exit 1
        fi
    done

    # Detect available resources
    detect_resources

    log "Starting ${step_name} (step ${step_num})"
}

# --- Complete Pipeline Step ---
# Arguments: $1 - Step number (optional, uses CURRENT_STEP_NUM if not provided)
# Creates completion flag and logs success.
complete_step() {
    local step_num="${1:-$CURRENT_STEP_NUM}"
    local step_name="${CURRENT_STEP_NAME:-unknown}"

    touch "${RESULTS_DIR}/step_completed_${step_num}.txt"
    log "Step ${step_num} (${step_name}) completed successfully."
}

# --- Validate Output Files ---
# Arguments: $1+ - Expected output files
# Options: --warn-only (don't exit on missing files)
# Returns: 0 if all files exist, 1 otherwise (exits unless --warn-only)
validate_outputs() {
    local warn_only=false
    local files=()

    # Parse arguments
    while [ $# -gt 0 ]; do
        case "$1" in
            --warn-only) warn_only=true; shift ;;
            *) files+=("$1"); shift ;;
        esac
    done

    local missing=()
    local empty=()

    for file in "${files[@]}"; do
        if [ ! -e "$file" ]; then
            missing+=("$file")
        elif [ -f "$file" ] && [ ! -s "$file" ]; then
            empty+=("$file")
        fi
    done

    if [ ${#missing[@]} -gt 0 ] || [ ${#empty[@]} -gt 0 ]; then
        if [ ${#missing[@]} -gt 0 ]; then
            log --level=ERROR "Missing output files:"
            for f in "${missing[@]}"; do
                log --level=ERROR "  - $f"
            done
        fi
        if [ ${#empty[@]} -gt 0 ]; then
            log --level=WARN "Empty output files:"
            for f in "${empty[@]}"; do
                log --level=WARN "  - $f"
            done
        fi

        if [ "$warn_only" = true ]; then
            return 1
        else
            exit 1
        fi
    fi

    return 0
}

# --- Validate Directory with Expected Files ---
# Arguments: $1 - Directory path
#            $2 - Minimum file count (default: 1)
#            $3 - File pattern (default: "*")
# Returns: 0 if valid, 1 otherwise
validate_output_dir() {
    local dir="$1"
    local min_count="${2:-1}"
    local pattern="${3:-*}"

    if [ ! -d "$dir" ]; then
        log --level=ERROR "Output directory does not exist: $dir"
        return 1
    fi

    local count=$(find "$dir" -maxdepth 1 -name "$pattern" -type f | wc -l)
    if [ "$count" -lt "$min_count" ]; then
        log --level=ERROR "Output directory has insufficient files: $dir"
        log --level=ERROR "  Expected at least $min_count files matching '$pattern', found $count"
        return 1
    fi

    return 0
}

# =============================================================================
# ORTHOGROUP MANIFEST FUNCTIONS
# =============================================================================
# These functions create and read manifest files for orthogroup tracking.

# --- Create Orthogroup Manifest ---
# Arguments: $1 - OrthoFinder results directory
#            $2 - Output manifest file (default: orthogroup_manifest.tsv)
# Creates a TSV with orthogroup metadata for array job coordination.
create_orthogroup_manifest() {
    local orthofinder_dir="$1"
    local manifest_file="${2:-${RESULTS_DIR}/orthogroup_manifest.tsv}"

    # Find the latest OrthoFinder results
    local results_dir=$(find "$orthofinder_dir" -maxdepth 2 -type d -name "Results_*" | sort -r | head -1)
    if [ -z "$results_dir" ]; then
        log --level=ERROR "No OrthoFinder Results directory found in $orthofinder_dir"
        return 1
    fi

    local orthogroups_file="${results_dir}/Orthogroups/Orthogroups.tsv"
    if [ ! -f "$orthogroups_file" ]; then
        log --level=ERROR "Orthogroups.tsv not found: $orthogroups_file"
        return 1
    fi

    # Create manifest header
    echo -e "index\torthogroup\tnum_seqs\thas_berghia\thas_reference" > "$manifest_file"

    # Parse orthogroups and create manifest
    # Note: Use process substitution to avoid subshell - variables persist correctly
    local index=0
    while IFS=$'\t' read -r og rest; do
        # Count sequences (tab-separated columns after OG name)
        local num_seqs=$(echo "$rest" | tr '\t' '\n' | grep -v "^$" | wc -l)

        # Check for Berghia and reference sequences
        local has_berghia="false"
        local has_reference="false"

        if echo "$rest" | grep -qi "berghia\|${BERGHIA_TAXID:-taxid_berghia}"; then
            has_berghia="true"
        fi

        if echo "$rest" | grep -q "ref_"; then
            has_reference="true"
        fi

        echo -e "${index}\t${og}\t${num_seqs}\t${has_berghia}\t${has_reference}" >> "$manifest_file"
        index=$((index + 1))
    done < <(tail -n +2 "$orthogroups_file")

    local total_ogs=$(tail -n +2 "$manifest_file" | wc -l)
    log "Created orthogroup manifest: $manifest_file ($total_ogs orthogroups)"
    echo "$manifest_file"
}

# --- Get Orthogroup Count from Manifest ---
# Arguments: $1 - Manifest file (optional, uses default location)
# Returns: Number of orthogroups
get_orthogroup_count() {
    local manifest_file="${1:-${RESULTS_DIR}/orthogroup_manifest.tsv}"

    if [ ! -f "$manifest_file" ]; then
        # Fall back to counting directories
        local og_dirs=$(find "${RESULTS_DIR}/phylogenies" -maxdepth 1 -type d -name "OG*" 2>/dev/null | wc -l)
        echo "$og_dirs"
        return
    fi

    tail -n +2 "$manifest_file" | wc -l
}

# --- Get Orthogroup by Index ---
# Arguments: $1 - Index (0-based)
#            $2 - Manifest file (optional)
# Returns: Orthogroup name (e.g., "OG0000001")
get_orthogroup_by_index() {
    local index="$1"
    local manifest_file="${2:-${RESULTS_DIR}/orthogroup_manifest.tsv}"

    if [ ! -f "$manifest_file" ]; then
        log --level=ERROR "Manifest file not found: $manifest_file"
        return 1
    fi

    # Validate index is a non-negative integer
    if ! [[ "$index" =~ ^[0-9]+$ ]]; then
        log --level=ERROR "Invalid index: $index (must be non-negative integer)"
        return 1
    fi

    # Check bounds: count data lines (excluding header)
    local total_count=$(tail -n +2 "$manifest_file" | wc -l)
    if [ "$index" -ge "$total_count" ]; then
        log --level=ERROR "Index $index out of bounds (max: $((total_count - 1)))"
        return 1
    fi

    # Add 2 to skip header (awk is 1-based, header is line 1)
    local line_num=$((index + 2))
    local og_name=$(awk -F'\t' "NR==$line_num {print \$2}" "$manifest_file")

    if [ -z "$og_name" ]; then
        log --level=ERROR "No orthogroup found at index $index"
        return 1
    fi

    echo "$og_name"
}

# --- Filter Orthogroups for Processing ---
# Arguments: $1 - Manifest file
#            $2+ - Filter options
# Options: --with-berghia (only OGs containing Berghia)
#          --with-reference (only OGs containing references)
#          --min-seqs=N (minimum sequence count)
# Returns: Space-separated list of orthogroup names
filter_orthogroups() {
    local manifest_file="$1"
    shift

    local with_berghia=false
    local with_reference=false
    local min_seqs=0

    # Parse options
    while [ $# -gt 0 ]; do
        case "$1" in
            --with-berghia) with_berghia=true; shift ;;
            --with-reference) with_reference=true; shift ;;
            --min-seqs=*) min_seqs="${1#--min-seqs=}"; shift ;;
            *) shift ;;
        esac
    done

    awk -F'\t' -v berghia="$with_berghia" -v ref="$with_reference" -v min="$min_seqs" '
        NR > 1 {
            if ($3 >= min) {
                if (berghia == "true" && $4 != "true") next
                if (ref == "true" && $5 != "true") next
                print $2
            }
        }
    ' "$manifest_file" | tr '\n' ' ' | sed 's/ $//'
}

# =============================================================================
# ARRAY JOB HELPERS
# =============================================================================
# Functions to help with SLURM array job coordination.

# --- Validate Array Index ---
# Arguments: $1 - Array task ID (SLURM_ARRAY_TASK_ID)
#            $2 - Maximum valid index
# Returns: 0 if valid, exits with error if invalid
validate_array_index() {
    local task_id="$1"
    local max_index="$2"

    if [ -z "$task_id" ]; then
        log --level=ERROR "SLURM_ARRAY_TASK_ID is not set"
        exit 1
    fi

    if [ "$task_id" -ge "$max_index" ]; then
        log "Array task $task_id exceeds orthogroup count ($max_index). Exiting gracefully."
        exit 0
    fi
}

# --- Create Array Job Checkpoint ---
# Arguments: $1 - Step name
#            $2 - Array task ID
# Creates a per-task checkpoint for resume capability
create_array_checkpoint() {
    local step_name="$1"
    local task_id="$2"
    local checkpoint_dir="${RESULTS_DIR}/checkpoints/${step_name}"

    # Validate checkpoint directory can be created
    if ! mkdir -p "$checkpoint_dir" 2>/dev/null; then
        log --level=WARN "Could not create checkpoint directory: $checkpoint_dir"
        return 1
    fi

    if ! touch "${checkpoint_dir}/task_${task_id}.done" 2>/dev/null; then
        log --level=WARN "Could not create checkpoint file for task $task_id"
        return 1
    fi

    return 0
}

# --- Check Array Job Completion ---
# Arguments: $1 - Step name
#            $2 - Expected task count
# Returns: 0 if all tasks completed, 1 otherwise
check_array_completion() {
    local step_name="$1"
    local expected_count="$2"
    local checkpoint_dir="${RESULTS_DIR}/checkpoints/${step_name}"

    if [ ! -d "$checkpoint_dir" ]; then
        return 1
    fi

    local completed=$(find "$checkpoint_dir" -name "task_*.done" | wc -l)
    if [ "$completed" -ge "$expected_count" ]; then
        return 0
    fi

    log "Array job ${step_name}: $completed / $expected_count tasks completed"
    return 1
}

# --- Get Incomplete Array Tasks ---
# Arguments: $1 - Step name
#            $2 - Total task count
# Returns: Comma-separated list of incomplete task IDs (for --array resubmission)
get_incomplete_tasks() {
    local step_name="$1"
    local total_count="$2"
    local checkpoint_dir="${RESULTS_DIR}/checkpoints/${step_name}"

    local incomplete=()
    for i in $(seq 0 $((total_count - 1))); do
        if [ ! -f "${checkpoint_dir}/task_${i}.done" ]; then
            incomplete+=("$i")
        fi
    done

    # Output as comma-separated for SLURM --array format
    # Use subshell to avoid polluting IFS in calling context
    ( IFS=','; echo "${incomplete[*]}" )
}
