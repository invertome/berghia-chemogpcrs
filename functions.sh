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
        ["iqtree"]="${IQTREE:-iqtree3} --version 2>&1 | head -1"
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
    sync "$temp_file" 2>/dev/null || true
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

    mkdir -p "$checkpoint_dir"
    atomic_write "$checkpoint_file" "$checkpoint_data"

    # Also create legacy marker for backwards compatibility
    touch "${RESULTS_DIR}/step_completed_${step_name}.txt"

    # Flush to disk so checkpoints survive crashes
    sync "${RESULTS_DIR}/step_completed_${step_name}.txt" 2>/dev/null || true

    log "Checkpoint created: ${step_name}"
}

# filter_fasta_by_length INPUT OUTPUT [MAX_LEN]
#
# Write to OUTPUT a copy of INPUT (FASTA) containing only sequences whose
# length is <= MAX_LEN. MAX_LEN defaults to the MAX_AA_LENGTH env var, or
# 1500 if neither is set. Inclusive on the upper bound.
#
# Why this exists (bead -m1f follow-up, stage 02 job 57653730): TMbed sorts
# inputs by length and processes longest last. ProtT5 has quadratic GPU
# memory in sequence length, so a few transcript-assembly outliers (chimeras,
# unspliced run-throughs) at thousands of aa take 130+ s/seq on a 2080 Ti
# and burn the wallclock. Real chemoreceptor GPCRs are 300-500 aa, so
# dropping >1500 aa is biologically safe and computationally essential.
filter_fasta_by_length() {
    local input="$1"
    local output="$2"
    local max_len="${3:-${MAX_AA_LENGTH:-1500}}"

    if [[ ! -f "$input" ]]; then
        echo "filter_fasta_by_length: input not found: $input" >&2
        return 1
    fi

    python3 - "$input" "$output" "$max_len" <<'PYEOF'
import sys
from Bio import SeqIO
src, dst, max_len = sys.argv[1], sys.argv[2], int(sys.argv[3])
kept = total = 0
with open(dst, 'w') as out:
    for rec in SeqIO.parse(src, 'fasta'):
        total += 1
        if len(rec.seq) <= max_len:
            SeqIO.write(rec, out, 'fasta')
            kept += 1
print(f'filter_fasta_by_length: kept {kept}/{total} sequences (max_len={max_len})',
      file=sys.stderr)
PYEOF
}

# identify_gpcr_candidates INPUT_FA OUTPUT_FA [CENSUS_TSV]
#
# HMM-first GPCR identification, detection-only (bead -m1f v3).
#
# Earlier versions added a Stage B annotation pass via lse.hmm (14k OG
# HMMs from Nath et al. lineage-specific expansion data) + conserved.hmm
# (34k OG HMMs from Nath et al. 1:1 orthologs). On Berghia (893 GPCR-
# positive candidates) Stage B took >16h to scan both libraries because
# the OG-level granularity meant each chemoreceptor-like protein matched
# hundreds of similar OG HMMs. Dropped: orthogroup placement is what
# stage 03's OrthoFinder does properly, so we don't need to re-derive it
# via massive hmmsearch in stage 02.
#
# Current detection layer (small + fast):
#   - classify_via_hmm.py against curated Swiss-Prot family HMMs +
#     Pfam fallback. LOO-thresholded, high precision; catches the
#     non-chemoreceptor families (bioamine, peptide, opsin, lipid, ...).
#   - Direct hmmsearch against TIAMMAT_GPCR_HMM at E <
#     GPCR_HMM_EVALUE (default 1e-5). 17 TIAMMAT-revised HMMs:
#     5x 7tm_N (Pfam GPCR class A-D retrained by TIAMMAT)
#     + 11x 7TM_GPCR_Sr* (invertebrate Sr chemoreceptor families,
#     retrained) + ABA_GPCR. Replaces the plain 4-HMM Pfam fallback —
#     better invertebrate-GPCR coverage at no extra runtime cost.
#
# Output:
#   - OUTPUT_FA: GPCR-positive subset of INPUT_FA (passed to TMbed
#     downstream; the >=6 TM filter then selects chemoreceptors).
#   - CENSUS_TSV (Berghia only): per-sequence family/subfamily/evalue/
#     source for the follow-up Berghia GPCR/brain-expression paper.
#     OG-level annotation lands later via stage 03 OrthoFinder.
identify_gpcr_candidates() {
    local input_fa="$1"
    local output_fa="$2"
    local census_tsv="${3:-}"
    local threads="${SLURM_CPUS_PER_TASK:-${CPUS:-4}}"
    local evalue="${GPCR_HMM_EVALUE:-1e-5}"

    if [[ ! -f "$input_fa" ]]; then
        echo "identify_gpcr_candidates: input not found: $input_fa" >&2
        return 1
    fi

    local out_dir; out_dir=$(dirname "$output_fa")
    local stem; stem=$(basename "$output_fa")
    stem="${stem%.fa}"; stem="${stem%.fasta}"
    local work="$out_dir/_hmm_filter_${stem}"
    mkdir -p "$work"

    local hmm_dir="${RESULTS_DIR}/classification/hmms"
    local pfam_dir="${hmm_dir}/pfam_fallback"
    local loo="${RESULTS_DIR}/classification/loo/loo_metrics.tsv"
    local tiammat_hmm="${TIAMMAT_GPCR_HMM:-${BASE_DIR}/references/tiammat_mollusca_gpcr.hmm}"
    local curated_hmm="${CURATED_CHEMO_HMM:-${RESULTS_DIR}/classification/hmms_curated_chemo/curated_chemo.hmm}"
    local class_tsv="$work/classification.tsv"
    local tiammat_tbl="$work/tiammat_gpcr_hits.tbl"
    local curated_tbl="$work/curated_chemo_hits.tbl"
    local ids_out="$work/gpcr_ids.txt"

    # --- Detection layer 1: classify_via_hmm.py (curated Swiss-Prot
    #     families + Pfam fallback, LOO-thresholded)
    if [[ -d "$hmm_dir" ]]; then
        local pfam_args=()
        [[ -d "$pfam_dir" ]] && pfam_args=(--pfam-fallback-dir "$pfam_dir")
        local loo_args=()
        [[ -f "$loo" ]] && loo_args=(--loo-metrics "$loo")
        python3 "${SCRIPTS_DIR:-./scripts}/classify_via_hmm.py" \
            --candidate-fasta "$input_fa" \
            --hmm-dir "$hmm_dir" \
            "${pfam_args[@]}" \
            "${loo_args[@]}" \
            --output-tsv "$class_tsv" \
            --threads "$threads" \
            >&2 || {
                echo "identify_gpcr_candidates: classify_via_hmm.py failed" >&2
                return 2
            }
    else
        echo "identify_gpcr_candidates: classification HMM dir not found ($hmm_dir); skipping curated scan" >&2
        : > "$class_tsv"
    fi

    # --- Detection layer 2: TIAMMAT-revised GPCR HMMs.
    #     17 HMMs (Pfam 7tm_1-7 retrained + invertebrate Sr-family
    #     chemoreceptors + ABA_GPCR). Single small library -> fast scan.
    if [[ -f "$tiammat_hmm" && -s "$tiammat_hmm" ]]; then
        hmmsearch --cpu "$threads" -E "$evalue" \
            --tblout "$tiammat_tbl" \
            "$tiammat_hmm" "$input_fa" > /dev/null 2>&1 || {
                echo "identify_gpcr_candidates: hmmsearch on TIAMMAT GPCR HMMs failed" >&2
                return 6
            }
    else
        echo "identify_gpcr_candidates: TIAMMAT GPCR HMM not found ($tiammat_hmm); skipping TIAMMAT scan" >&2
        : > "$tiammat_tbl"
    fi

    # --- Detection layer 3: curated invertebrate chemoreceptor HMMs.
    #     HIGHEST-confidence detection: HMMs built from functionally
    #     validated lophotrochozoan chemoreceptors (Cummins 2009 Aplysia
    #     + deorphanized Schistosoma + ...). Each HMM is one family
    #     (apgr, schisto_chemo, annelid_chemo, other_gastropod_chemo).
    #     Bead -k0g, 2026-05-19.
    if [[ -f "$curated_hmm" && -s "$curated_hmm" ]]; then
        hmmsearch --cpu "$threads" -E "$evalue" \
            --tblout "$curated_tbl" \
            "$curated_hmm" "$input_fa" > /dev/null 2>&1 || {
                echo "identify_gpcr_candidates: hmmsearch on curated chemoreceptor HMMs failed" >&2
                return 7
            }
    else
        echo "identify_gpcr_candidates: curated chemoreceptor HMM not found or empty ($curated_hmm); skipping curated_chemo scan (run scripts/build_curated_chemo_hmms.sh if you want this layer)" >&2
        : > "$curated_tbl"
    fi

    # --- Merge detection evidence -> GPCR-positive IDs + (optional) census
    local census_args=()
    [[ -n "$census_tsv" ]] && census_args=(--census-out "$census_tsv")
    python3 "${SCRIPTS_DIR:-./scripts}/identify_gpcrs.py" \
        --classification-tsv "$class_tsv" \
        --tiammat-tblout "$tiammat_tbl" \
        --curated-chemo-tblout "$curated_tbl" \
        --ids-out "$ids_out" \
        "${census_args[@]}" || {
            echo "identify_gpcr_candidates: identify_gpcrs.py failed" >&2
            return 4
        }

    # Extract GPCR-positive sequences from input
    "${SEQTK:-seqtk}" subseq "$input_fa" "$ids_out" > "$output_fa"

    local n_in n_gpcr
    n_in=$(grep -c '^>' "$input_fa" 2>/dev/null || echo 0)
    n_gpcr=$(grep -c '^>' "$output_fa" 2>/dev/null || echo 0)
    echo "identify_gpcr_candidates: $n_gpcr / $n_in proteins flagged as GPCR (input=$(basename "$input_fa"))" >&2
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

    # Bead -ryr: actually append step_record into the provenance JSON.
    # Use Python (always available) since jq may not be installed on the HPC node.
    if [ -f "$PROVENANCE_FILE" ]; then
        python3 - "$PROVENANCE_FILE" <<PYEOF
import json, sys
prov_path = sys.argv[1]
new_step_json = '''${step_record}'''
try:
    with open(prov_path) as f:
        prov = json.load(f)
except (FileNotFoundError, json.JSONDecodeError):
    prov = {"steps": []}
if "steps" not in prov or not isinstance(prov["steps"], list):
    prov["steps"] = []
try:
    prov["steps"].append(json.loads(new_step_json))
except json.JSONDecodeError as e:
    sys.stderr.write(f"WARNING: provenance step record could not be parsed: {e}\n")
    sys.exit(0)
with open(prov_path, "w") as f:
    json.dump(prov, f, indent=2)
PYEOF
    fi
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

    # Flush output files to disk so results survive crashes
    if [ $exit_code -eq 0 ]; then
        if [ -n "$stdout_file" ]; then
            sync "$stdout_file" 2>/dev/null || true
        fi
        # Sync declared output files
        if [ -n "$outputs" ]; then
            local IFS=','
            for outf in $outputs; do
                outf=$(echo "$outf" | xargs)
                [ -f "$outf" ] && sync "$outf" 2>/dev/null || true
            done
        fi
    fi

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

# --- Assert no duplicate FASTA sequence IDs (fail-loud) ---
# Arguments: $1 - FASTA file; $2 - optional context label for the message.
# Returns 1 (and logs an ERROR) if any sequence ID — the first whitespace-
# delimited token of a header — occurs more than once. IQ-TREE aborts on
# duplicate taxon names, so this guards the stage-04 per-class combined pool
# against silent double-adds (bead berghia-chemogpcrs-3sd).
assert_no_duplicate_fasta_ids() {
    local fasta="$1"
    local label="${2:-$fasta}"
    if [ ! -f "$fasta" ]; then
        log --level=ERROR "assert_no_duplicate_fasta_ids: file not found: $fasta"
        return 1
    fi
    local dups
    # `|| true`: a header-less file makes grep exit 1, which would trip the
    # caller's `set -o pipefail`; we only care about the (possibly empty) output.
    dups=$(grep '^>' "$fasta" | sed 's/^>//; s/[[:space:]].*//' | sort | uniq -d) || true
    if [ -n "$dups" ]; then
        local n examples
        n=$(printf '%s\n' "$dups" | sed '/^$/d' | wc -l | tr -d ' ')
        examples=$(printf '%s\n' "$dups" | sed '/^$/d' | head -3 | tr '\n' ' ')
        log --level=ERROR "Duplicate FASTA IDs in ${label}: ${n} id(s) occur more than once (e.g. ${examples})"
        return 1
    fi
    return 0
}

# --- Assert a FASTA holds at least one record ---
# Arguments: $1 - FASTA path
#            $2 - human-readable label for the error message (default: the path)
# Returns 0 when the file has >=1 '>' record; logs ERROR and returns 1 otherwise.
#
# Why (bead 444): `[ -f ]` alone is satisfied by a present-but-empty file. An
# empty per-class pool then flows through cat/length-filter/MAFFT/ClipKit/
# FastTree/IQ-TREE — all of which "succeed" on empty input — and the stage
# emits a degenerate tree with no error. Fail at the boundary instead.
assert_fasta_non_empty() {
    local fasta="$1"
    local label="${2:-$fasta}"
    if [ ! -f "$fasta" ]; then
        log --level=ERROR "${label}: FASTA not found: ${fasta}"
        return 1
    fi
    local n
    # `|| true`: grep -c exits 1 on a header-less file; we want the count, not
    # the status (and must not trip a caller's `set -o pipefail`).
    n=$(grep -c '^>' "$fasta") || true
    if [ "${n:-0}" -eq 0 ]; then
        log --level=ERROR "${label}: ${fasta} contains 0 FASTA records — refusing to continue"
        return 1
    fi
    return 0
}

# --- Draw a per-class outgroup, excluding Berghia's own paralogs ---
# Arguments: $1 - source FASTA (the SISTER class's P2 pool, per the outgroup
#                 swap-map A<-C, B<-A, C<-A, F<-A; locked decision 2026-05-28)
#            $2 - that same sister class's pool_members_class_<C>.tsv manifest
#                 (columns: seq_id, taxid, source in {anchor, berghia, ref})
#            $3 - number of records to draw
#            $4 - output FASTA
# Writes $4 (always created, possibly empty) and returns 0; returns 1 only when
# the source FASTA is missing.
#
# Why (bead 444): stage 04 used to draw the outgroup as the first N records of
# the sister pool. But build_per_class_reference_pools.py writes
#   selected = berghia_pairs + anchor_pairs + selected_refs
# so Berghia records come FIRST in every pool file, and the head-N draw took
# exactly Berghia's own class-B/C/F candidates. Every per-class tree was
# therefore rooted on unannotated Berghia paralogs whose class assignment came
# from this pipeline's own classifier — circular (the classifier's output roots
# the tree that tests it), single-species, and not the characterized reference
# GPCRs the design intended. Anchors and refs are both valid outgroup material;
# only `berghia` is excluded.
#
# Records absent from the manifest are treated as non-Berghia: the manifest and
# the pool FASTA are co-written from the same `selected` list, so an absent id
# means the two files are out of sync, not that provenance is Berghia.
draw_outgroup_fasta() {
    local source_fa="$1"
    local members_tsv="$2"
    local n="${3:-10}"
    local out_fa="$4"

    if [ ! -f "$source_fa" ]; then
        log --level=ERROR "draw_outgroup_fasta: outgroup source FASTA not found: ${source_fa}"
        return 1
    fi

    # `-s`, not `-f`: a zero-byte manifest carries no provenance either, and it
    # would additionally break the NR==FNR two-pass idiom below (with an empty
    # first file, NR==FNR stays true and the FASTA gets parsed as the manifest).
    if [ ! -s "$members_tsv" ]; then
        # Provenance is undeterminable without the manifest. Fall back to the
        # previous unfiltered behaviour rather than hard-failing the stage, but
        # say so plainly — the outgroup may be Berghia paralogs again.
        log --level=WARN "draw_outgroup_fasta: pool-members manifest missing or empty: ${members_tsv} — cannot determine provenance; falling back to the first ${n} records UNFILTERED, so the outgroup MAY CONTAIN Berghia paralogs"
        awk -v n="$n" '/^>/{c++} c>n{exit} {print}' "$source_fa" > "$out_fa"
    else
        # Strict consistency gate (user decision, bead 444). The pool FASTA and
        # its manifest are written back-to-back from the same `selected` list in
        # one loop of build_per_class_reference_pools.py, so a record with no
        # manifest row means the two files are OUT OF SYNC (a run interrupted
        # between the writes, or a stale manifest beside a fresh FASTA).
        #
        # Do not guess. Assuming such a record is non-Berghia readmits the exact
        # circular-rooting bug this helper exists to prevent; dropping it quietly
        # shrinks the outgroup instead. Both are SILENT. Refuse and let the
        # caller skip the class, which is the only outcome that is visible.
        local unmapped n_unmapped
        unmapped=$(awk -F'\t' '
            NR == FNR {
                if (FNR == 1 && $1 == "seq_id") next
                known[$1] = 1
                next
            }
            /^>/ {
                id = substr($0, 2)
                sub(/[ \t].*/, "", id)
                if (!(id in known)) print id
            }
        ' "$members_tsv" "$source_fa")
        if [ -n "$unmapped" ]; then
            n_unmapped=$(printf '%s\n' "$unmapped" | wc -l | tr -d ' ')
            log --level=ERROR "draw_outgroup_fasta: ${source_fa} and ${members_tsv} are OUT OF SYNC — ${n_unmapped} pool record(s) have no provenance row, e.g. $(printf '%s' "$unmapped" | head -5 | tr '\n' ' '). Refusing to draw an outgroup rather than guess whether they are Berghia."
            : > "$out_fa"
            return 2
        fi

        # Pass 1 (NR==FNR): collect the Berghia ids from the manifest.
        # Pass 2: stream the FASTA, copying complete records (header + the whole
        # multi-line sequence block, via the sticky `emit` flag) for the first
        # $n non-Berghia ids. The manifest id is the bare first header token, so
        # any description after whitespace is stripped before the lookup.
        awk -F'\t' -v n="$n" '
            NR == FNR {
                if (FNR == 1 && $1 == "seq_id") next
                if ($3 == "berghia") berghia[$1] = 1
                next
            }
            /^>/ {
                id = substr($0, 2)
                sub(/[ \t].*/, "", id)
                if (kept >= n) exit
                emit = (id in berghia) ? 0 : 1
                if (emit) kept++
            }
            emit
        ' "$members_tsv" "$source_fa" > "$out_fa"
    fi

    local drawn
    # `|| true`: grep -c exits 1 on an empty draw; we want the count, not the
    # status (and must not trip a caller's `set -o pipefail`).
    drawn=$(grep -c '^>' "$out_fa") || true
    if [ "${drawn:-0}" -eq 0 ]; then
        log --level=WARN "draw_outgroup_fasta: 0 non-Berghia records drawn from ${source_fa} — the tree rooted with this outgroup will be effectively unrooted"
    fi
    return 0
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

# --- Resolve the Berghia Candidate FASTA Source (genome-track aware) ---
# Echoes the Berghia candidate FASTA that downstream stages (03 orthology
# clustering) should consume. When the genome track is enabled (RUN_GENOME_TRACK
# != 0, mirroring stage 02c's toggle) AND stage 02c has written the reconciled
# set, that reconciled FASTA is the source; otherwise the legacy stage-02
# transcriptome candidate set — byte-identical to pre-genome-track behavior.
# Pure: echoes a path only, no side effects.
berghia_candidate_fasta() {
    local reconciled="${RESULTS_DIR}/reconciliation/reconciled_candidates.faa"
    local legacy="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
    if [ "${RUN_GENOME_TRACK:-1}" != "0" ] && [ -f "$reconciled" ]; then
        printf '%s\n' "$reconciled"
    else
        printf '%s\n' "$legacy"
    fi
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

# --- Count exact occurrences of a taxid in a taxid list ---
# Arguments: $1 - whitespace-separated taxid list (as emitted by
#                 get_taxids_from_fasta)
#            $2 - taxid to count
# Returns: echoes exactly one integer line (0 when absent)
#
# Replaces the `echo "$list" | grep -c "$taxid" || echo 0` idiom, which had two
# defects (bead 444):
#   (a) `grep -c` prints "0" AND exits 1 on no match, so `|| echo 0` appended a
#       SECOND "0". The captured value was the two-line string "0\n0", and
#       `[ "$n" -gt 0 ]` then died with "integer expected" instead of comparing.
#   (b) grep matched substrings, so BERGHIA_TAXID=1287507 was "found" inside an
#       unrelated taxid such as 12875070 — a silent false positive.
# -x (whole line) + -F (literal) after splitting on whitespace fixes both.
count_taxid_occurrences() {
    local taxid_list="$1"
    local taxid="$2"
    local n
    n=$(printf '%s' "$taxid_list" | tr -s '[:space:]' '\n' | grep -c -x -F -- "$taxid") || true
    echo "${n:-0}"
}

# --- Find nucleotide sequences for an orthogroup's protein alignment ---
# Arguments: $1 - protein alignment (headers give seq_ids); $2 - output FASTA
# Reads globals TRANSCRIPTOME_DIR and REFERENCE_CDS_FILE. Returns 0 if a
# majority of sequences were found, 1 otherwise.
find_nucleotide_sequences() {
    local protein_file="$1"
    local output_file="$2"

    > "$output_file"
    local found_count=0
    local missing_count=0
    local ref_count=0

    # Extract sequence IDs from protein alignment
    while IFS= read -r header; do
        seq_id=$(echo "$header" | sed 's/>//')

        # Check if this is a reference sequence
        is_ref=$(is_reference_seq "$seq_id" && echo "yes" || echo "no")

        # Use metadata lookup to get taxid (with fallback to header parsing)
        taxid=$(get_taxid_for_seq "$seq_id")

        # Extract protein ID portion (after taxid prefix)
        # For ref_TAXID_N format: protein_id is N
        # For TAXID_N format: protein_id is N
        if [[ "$seq_id" == ref_* ]]; then
            protein_id=$(echo "$seq_id" | cut -d'_' -f3-)
        else
            protein_id=$(echo "$seq_id" | cut -d'_' -f2-)
        fi

        found=false

        # --- Check reference CDS file first for reference sequences ---
        if [ "$is_ref" = "yes" ] && [ -f "$REFERENCE_CDS_FILE" ]; then
            # Try exact match with the seq_id
            if grep -q "^>${seq_id}" "$REFERENCE_CDS_FILE"; then
                # Extract sequence (handle multi-line FASTA)
                awk -v id="$seq_id" '
                    /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print; found=1; next } else { found=0 } }
                    found { print }
                ' "$REFERENCE_CDS_FILE" >> "$output_file"
                found=true
                ((ref_count++))
            fi
        fi

        # --- Search transcriptome directory for non-reference sequences ---
        if [ "$found" = false ]; then
            # nullglob so the per-species glob below drops out cleanly when it
            # has no match, leaving the literal patterns to run.
            shopt -s nullglob
            for ext in mrna cds fna fa; do
                # Naming conventions: bare taxid, taxid_taxid, and the canonical
                # per-species <taxid>_<genus>_<species> layout (config.sh
                # BERGHIA_FILE_PREFIX) — the last matched with a glob.
                for nuc_file in "${TRANSCRIPTOME_DIR}/${taxid}.${ext}" \
                               "${TRANSCRIPTOME_DIR}/taxid_${taxid}.${ext}" \
                               "${TRANSCRIPTOME_DIR}/${taxid}_${taxid}.${ext}" \
                               "${TRANSCRIPTOME_DIR}/${taxid}_"*".${ext}"; do
                    if [ -f "$nuc_file" ]; then
                        # Extract the corresponding nucleotide sequence
                        # Try exact match first, then partial match
                        if grep -q "^>${seq_id}" "$nuc_file"; then
                            # Handle multi-line FASTA
                            awk -v id="$seq_id" '
                                /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print; found=1; next } else { found=0 } }
                                found { print }
                            ' "$nuc_file" >> "$output_file"
                            found=true
                            break 2
                        elif grep -q "^>${protein_id}" "$nuc_file"; then
                            # Extract and rename header to match protein
                            awk -v id="$protein_id" -v new_id="$seq_id" '
                                /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print ">" new_id; found=1; next } else { found=0 } }
                                found { print }
                            ' "$nuc_file" >> "$output_file"
                            found=true
                            break 2
                        fi
                    fi
                done
            done
            shopt -u nullglob
        fi

        # --- Fallback: check all reference CDS for non-reference sequences too ---
        if [ "$found" = false ] && [ -f "$REFERENCE_CDS_FILE" ]; then
            if grep -q "^>${seq_id}" "$REFERENCE_CDS_FILE"; then
                awk -v id="$seq_id" '
                    /^>/ { if (match($0, "^>" id "($|[[:space:]])")) { print; found=1; next } else { found=0 } }
                    found { print }
                ' "$REFERENCE_CDS_FILE" >> "$output_file"
                found=true
            fi
        fi

        if [ "$found" = true ]; then
            ((found_count++))
        else
            ((missing_count++))
            log "Warning: No nucleotide sequence found for ${seq_id}"
        fi
    done < <(grep "^>" "$protein_file")

    # Log summary
    log "Nucleotide lookup: found=${found_count} (${ref_count} from refs), missing=${missing_count}"

    # Return success if we found at least some sequences (majority required for meaningful analysis)
    local total=$((found_count + missing_count))
    if [ "$found_count" -gt 0 ] && [ "$found_count" -ge $((total / 2)) ]; then
        return 0
    else
        return 1
    fi
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

    # Find Orthogroups.tsv — may be in original or resumed Results_* directory
    local orthogroups_file=$(find "$orthofinder_dir" -maxdepth 5 -path "*/Orthogroups/Orthogroups.tsv" -type f 2>/dev/null | head -1)
    if [ -z "$orthogroups_file" ]; then
        log --level=ERROR "Orthogroups.tsv not found under $orthofinder_dir"
        return 1
    fi
    local results_dir=$(dirname "$(dirname "$orthogroups_file")")

    # Create manifest using Python for reliable TSV parsing
    export ORTHOGROUPS_FILE="$orthogroups_file"
    export MANIFEST_FILE="$manifest_file"
    export BERGHIA_TAXID="${BERGHIA_TAXID:-1287507}"

    python3 << 'PYEOF'
import os, csv

og_file = os.environ["ORTHOGROUPS_FILE"]
out_file = os.environ["MANIFEST_FILE"]
berghia_taxid = os.environ["BERGHIA_TAXID"]

with open(og_file) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    # Find Berghia column (contains taxid in column name)
    berghia_col = None
    for i, col in enumerate(header):
        if berghia_taxid in col or 'berghia' in col.lower():
            berghia_col = i
            break

    with open(out_file, 'w') as out:
        out.write("index\torthogroup\tnum_seqs\thas_berghia\thas_reference\n")
        for idx, row in enumerate(reader):
            og = row[0]
            # Count total sequences across all species (comma-separated in each cell)
            num_seqs = 0
            has_reference = False
            for i, cell in enumerate(row[1:], 1):
                cell = cell.strip()
                if cell:
                    seqs = [s.strip() for s in cell.split(',') if s.strip()]
                    num_seqs += len(seqs)
                    if 'ref_' in header[i].lower():
                        has_reference = True

            has_berghia = False
            if berghia_col is not None and berghia_col < len(row):
                has_berghia = bool(row[berghia_col].strip())

            out.write(f"{idx}\t{og}\t{num_seqs}\t{str(has_berghia).lower()}\t{str(has_reference).lower()}\n")

    print(f"Wrote {idx + 1} orthogroups to {out_file}")
PYEOF

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

# --- Assert the SLURM array range covers the whole manifest (fail-loud) ---
# Arguments: $1 - orthogroup count (get_orthogroup_count / OG_COUNT)
#            $2 - highest scheduled array index (default: SLURM_ARRAY_TASK_MAX)
# Returns 1 (and logs an ERROR) when the array's top index is below OG_COUNT-1,
# i.e. orthogroups at index > that top would be silently never scheduled. The
# hardcoded `#SBATCH --array=0-999` truncates any larger manifest without
# warning; the submit wrapper must size --array from OG_COUNT (chunk if it
# exceeds MaxArraySize). Bead berghia-chemogpcrs-fxx.
assert_array_covers_manifest() {
    local og_count="$1"
    local array_max="${2:-${SLURM_ARRAY_TASK_MAX}}"

    # Not an array job (local / non-array run) → nothing to check.
    if [ -z "$array_max" ]; then
        return 0
    fi
    # Count unknown → warn but do not block (fail-open on missing info).
    if [ -z "$og_count" ] || ! [ "$og_count" -gt 0 ] 2>/dev/null; then
        log --level=WARN "assert_array_covers_manifest: orthogroup count unknown ($og_count); skipping array-range check"
        return 0
    fi

    local needed=$((og_count - 1))
    if [ "$array_max" -lt "$needed" ]; then
        local missing=$((og_count - array_max - 1))
        log --level=ERROR "SLURM array range too small: top index ${array_max} < ${needed} (OG_COUNT=${og_count}); ${missing} orthogroup(s) at index >${array_max} would be silently skipped. Resize the submit array (--array=0-${needed}%K, chunk if > MaxArraySize) and resubmit."
        return 1
    fi
    return 0
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

# -----------------------------------------------------------------------------
# Alignment-cleanup stack orchestration (May 2026 v2 — replaces HmmCleaner)
#
# Runs the full filter chain on a protein FASTA:
#   PREQUAL (residue mask)
#     -> alignment ensemble (canonical + MAFFT variants + FAMSA)
#     -> CLOAK (consensus mask)
#     -> TAPER (residue-outlier mask)
# Each step is gated by RUN_PREQUAL / RUN_CLOAK / RUN_TAPER (defaults ON).
# Failures are logged but non-fatal: the chain falls through and the caller
# receives the best-effort output (e.g. canonical alignment if CLOAK fails).
#
# Arguments:
#   $1 - input FASTA (unaligned sequences)
#   $2 - output FASTA (cleaned alignment)
#   $3 - work dir prefix (used for ensemble + intermediates; deletable on success)
#   $4 - tag for logs (e.g. orthogroup base name)
#   $5 - threads (optional; defaults to $CPUS)
#
# Returns: 0 on success (output exists, non-empty), 1 if even the canonical
# alignment couldn't be produced.
# -----------------------------------------------------------------------------
run_alignment_filter_stack() {
    local input="$1"
    local output="$2"
    local workdir="$3"
    local tag="$4"
    local threads="${5:-${CPUS:-4}}"

    if [ ! -f "$input" ]; then
        log --level=ERROR "filter_stack: input not found: $input"
        return 1
    fi

    mkdir -p "$workdir" "$(dirname "$output")"
    local cur="$input"

    # --- Stage 1: PREQUAL ---
    # Reuse pre-computed PREQUAL output when present (idempotent +
    # avoids re-running PREQUAL's O(n^2) pairwise step on big inputs
    # like all_berghia_refs_filtered.fa which has 2829 seqs and takes
    # 6-24h). For per-OG calls with <100 seqs, PREQUAL runs in seconds
    # and the precompute step is no-op for them. The big-tree precompute
    # is shipped via scripts/unity/precompute_prequal_big.sh.
    if [ "${RUN_PREQUAL:-1}" = "1" ]; then
        local prequal_out="${workdir}/${tag}_prequal.fa"
        if [ -s "$prequal_out" ]; then
            log "filter_stack[$tag]: reusing existing PREQUAL output ($(grep -c '^>' "$prequal_out") seqs)"
            cur="$prequal_out"
        elif bash "${SCRIPTS_DIR}/run_prequal.sh" \
                --input="$cur" --output="$prequal_out" \
                2>"${LOGS_DIR:-/tmp}/prequal_${tag}.err"; then
            cur="$prequal_out"
        else
            log --level=WARN "filter_stack[$tag]: PREQUAL failed, continuing with un-masked input"
        fi
    fi

    # --- Stage 2: alignment ensemble (canonical + variants) ---
    local ensemble_dir="${workdir}/${tag}_ensemble"
    if ! bash "${SCRIPTS_DIR}/run_alignment_ensemble.sh" \
            --input="$cur" \
            --output-dir="$ensemble_dir" \
            --threads="$threads" \
            --canonical-aligner="${ENSEMBLE_CANONICAL_ALIGNER:-auto}" \
            2>"${LOGS_DIR:-/tmp}/ensemble_${tag}.err"; then
        log --level=ERROR "filter_stack[$tag]: ensemble failed; aborting"
        return 1
    fi
    cur="${ensemble_dir}/canonical.fa"
    [ -s "$cur" ] || { log --level=ERROR "filter_stack[$tag]: empty canonical alignment"; return 1; }

    # --- Stage 3: CLOAK consensus ---
    if [ "${RUN_CLOAK:-1}" = "1" ]; then
        local cloak_out="${workdir}/${tag}_cloaked.fa"
        local n_variants
        n_variants=$(find "$ensemble_dir" -maxdepth 1 -name '*.fa' -type f | wc -l)
        if [ "$n_variants" -lt 2 ]; then
            log --level=WARN "filter_stack[$tag]: only $n_variants alignment(s); skipping CLOAK"
        elif bash "${SCRIPTS_DIR}/run_cloak.sh" \
                --ensemble-dir="$ensemble_dir" \
                --output="$cloak_out" \
                2>"${LOGS_DIR:-/tmp}/cloak_${tag}.err"; then
            cur="$cloak_out"
        else
            log --level=WARN "filter_stack[$tag]: CLOAK failed, continuing with canonical alignment"
        fi
    fi

    # --- Stage 4: TAPER residue-outlier mask ---
    if [ "${RUN_TAPER:-1}" = "1" ]; then
        local taper_out="${workdir}/${tag}_tapered.fa"
        if bash "${SCRIPTS_DIR}/run_taper.sh" \
                --input="$cur" --output="$taper_out" \
                2>"${LOGS_DIR:-/tmp}/taper_${tag}.err"; then
            cur="$taper_out"
        else
            log --level=WARN "filter_stack[$tag]: TAPER failed, continuing with prior alignment"
        fi
    fi

    # --- Optional 5th pass: HmmCleaner (deprecated, off by default) ---
    if [ "${RUN_HMMCLEANER:-0}" = "1" ] && command -v HmmCleaner.pl >/dev/null 2>&1; then
        local hmmc_out="${workdir}/${tag}_hmmcleaned.fa"
        if bash "${SCRIPTS_DIR}/run_hmmcleaner.sh" \
                --input="$cur" --output="$hmmc_out" \
                2>"${LOGS_DIR:-/tmp}/hmmcleaner_${tag}.err"; then
            cur="$hmmc_out"
        else
            log --level=WARN "filter_stack[$tag]: HmmCleaner failed, continuing"
        fi
    fi

    # --- Finalize: copy final cleaned alignment to caller's output path ---
    cp "$cur" "$output"

    # Sibling: preserve the canonical (pre-CLOAK / pre-TAPER) alignment for
    # downstream consumers that need the variable columns intact. 04b ECL
    # divergence analysis specifically *wants* the variable ECL columns —
    # CLOAK would have masked exactly those (because they're alignment-
    # uncertain). Save the un-cleaned canonical MAFFT/FAMSA output so 04b
    # can read it instead of the filter-stack output.
    local canonical_sibling="${output%.fa}_canonical.fa"
    if [ -s "${ensemble_dir}/canonical.fa" ]; then
        cp "${ensemble_dir}/canonical.fa" "$canonical_sibling"
    fi

    # Provenance summary alongside the output
    {
        echo "tool=run_alignment_filter_stack"
        echo "input=$input"
        echo "output=$output"
        echo "canonical_sibling=$canonical_sibling"
        echo "stages_run=prequal=${RUN_PREQUAL:-1},cloak=${RUN_CLOAK:-1},taper=${RUN_TAPER:-1},hmmcleaner=${RUN_HMMCLEANER:-0}"
        echo "ensemble_dir=$ensemble_dir"
        echo "final_alignment=$cur"
        echo "n_seqs=$(grep -c '^>' "$output" 2>/dev/null || echo 0)"
        echo "run_at=$(date -Iseconds)"
    } > "${output}.filter_stack.txt"

    log "filter_stack[$tag]: $output ready (n=$(grep -c '^>' "$output"); canonical sibling: $(basename "$canonical_sibling"))"
    return 0
}
