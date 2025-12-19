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
atomic_write() {
    local target="$1"
    local content="$2"
    local temp_file="${target}.tmp.$$"

    echo "$content" > "$temp_file"
    mv "$temp_file" "$target"
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

    # Compute input checksums
    local input_checksums="{"
    local first=true
    IFS=',' read -ra INPUT_FILES <<< "$inputs"
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
    IFS=',' read -ra OUTPUT_FILES <<< "$outputs"
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
run_command() {
    local name="$1"
    shift

    local inputs=""
    local outputs=""
    local allow_fail=false

    # Parse options
    while [[ "$1" == --* ]]; do
        case "$1" in
            --inputs=*) inputs="${1#--inputs=}"; shift ;;
            --outputs=*) outputs="${1#--outputs=}"; shift ;;
            --allow-fail) allow_fail=true; shift ;;
            *) break ;;
        esac
    done

    local start_time=$(date +%s)
    log "Running: $name"
    log "  Command: $*"

    # Record in provenance
    record_provenance "$name" "$*" "$inputs" "$outputs"

    # Execute command
    "$@" > "${LOGS_DIR}/${name}.log" 2> "${LOGS_DIR}/${name}.err"
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
finalize_pipeline() {
    local exit_code="${1:-0}"

    if [ -n "$PROVENANCE_FILE" ] && [ -f "$PROVENANCE_FILE" ]; then
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
# Returns: Alignment length (number of columns)
get_alignment_length() {
    local aln="$1"
    # Get length of first sequence (all should be same length in alignment)
    awk '/^>/{if(seq){print length(seq);exit}seq=""} !/^>/{seq=seq$0} END{if(seq)print length(seq)}' "$aln" 2>/dev/null || echo 0
}

# --- Estimate Memory for Alignment ---
# Arguments: $1 - FASTA file
# Returns: Estimated memory in GB (e.g., "8G")
estimate_memory_for_alignment() {
    local fasta="$1"
    local num_seqs=$(count_sequences "$fasta")
    local avg_len=$(get_avg_seq_length "$fasta")

    # Distance matrix: O(n^2) with ~10 bytes per cell
    local mem_bytes=$(( num_seqs * num_seqs * 10 + num_seqs * avg_len * MEM_MULT_MAFFT ))
    local mem_gb=$(( mem_bytes / 1073741824 + MEM_BASE_GB ))

    # Cap at maximum
    [ $mem_gb -gt $MEM_MAX_GB ] && mem_gb=$MEM_MAX_GB
    [ $mem_gb -lt $MEM_BASE_GB ] && mem_gb=$MEM_BASE_GB

    echo "${mem_gb}G"
}

# --- Estimate Memory for Tree Building ---
# Arguments: $1 - Alignment file
# Returns: Estimated memory in GB (e.g., "16G")
estimate_memory_for_tree() {
    local aln="$1"
    local num_seqs=$(count_sequences "$aln")
    local aln_len=$(get_alignment_length "$aln")

    # IQ-TREE memory: sites * taxa * multiplier
    local mem_bytes=$(( num_seqs * aln_len * MEM_MULT_IQTREE ))
    local mem_gb=$(( mem_bytes / 1073741824 + MEM_BASE_GB ))

    # Cap at maximum
    [ $mem_gb -gt $MEM_MAX_GB ] && mem_gb=$MEM_MAX_GB
    [ $mem_gb -lt $MEM_BASE_GB ] && mem_gb=$MEM_BASE_GB

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
        for f in "$input"/*.fa "$input"/*.faa "$input"/*.fasta; do
            [ -f "$f" ] && total_seqs=$((total_seqs + $(count_sequences "$f")))
        done
    else
        # Assume input is sequence count
        total_seqs="$input"
    fi

    # OrthoFinder: O(n^2) for all-vs-all DIAMOND
    local mem_bytes=$(( total_seqs * total_seqs * MEM_MULT_ORTHOFINDER ))
    local mem_gb=$(( mem_bytes / 1073741824 + MEM_BASE_GB ))

    # Cap at maximum
    [ $mem_gb -gt $MEM_MAX_GB ] && mem_gb=$MEM_MAX_GB
    [ $mem_gb -lt 8 ] && mem_gb=8  # OrthoFinder needs at least 8GB

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

    # HyPhy aBSREL: sites * branches * multiplier
    local mem_bytes=$(( aln_len * num_branches * MEM_MULT_HYPHY ))
    local mem_gb=$(( mem_bytes / 1073741824 + MEM_BASE_GB ))

    # Cap at maximum
    [ $mem_gb -gt $MEM_MAX_GB ] && mem_gb=$MEM_MAX_GB
    [ $mem_gb -lt $MEM_BASE_GB ] && mem_gb=$MEM_BASE_GB

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
