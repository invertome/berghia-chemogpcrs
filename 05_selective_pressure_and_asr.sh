#!/bin/bash
# 05_selective_pressure_and_asr.sh
# Purpose: Analyze selective pressure with HyPhy's aBSREL and reconstruct ancestral sequences with FastML.
# Inputs: Orthogroups from step 03, alignments from step 04
# Outputs: Selective pressure results in ${RESULTS_DIR}/selective_pressure/, ASR sequences in ${RESULTS_DIR}/asr/
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=selective_pressure_asr
#SBATCH --output=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.out
#SBATCH --error=${LOGS_DIR}/05_selective_pressure_asr_%j_%a.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH --cpus-per-task=8
#SBATCH --mem=${DEFAULT_MEM}
#SBATCH --array=0-999%50  # Adjusted for orthogroup processing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

source config.sh
source functions.sh
activate_conda_env   # self-activate the project conda env before running tools
# ONE deterministic, chronologically-correct rule for which OrthoFinder
# run is authoritative (mtime of Orthogroups.tsv). Shared by stages
# 03/03b/04/05/06c/07 so they can no longer resolve different runs.
# shellcheck source=scripts/orthofinder_paths.sh
source "${SCRIPTS_DIR:-scripts}/orthofinder_paths.sh"

# Create output directories
mkdir -p "${RESULTS_DIR}/selective_pressure" "${RESULTS_DIR}/selective_pressure/nucleotide" "${RESULTS_DIR}/asr" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

# This stage runs as `#SBATCH --array=0-999%50`, so up to 50 tasks run on
# DIFFERENT nodes against the same shared results directory. Every path that
# is NOT namespaced by the orthogroup needs a per-task-unique staging name.
# $$ alone is insufficient: PIDs are unique per node, the staging directory
# is shared storage, so two tasks on two nodes can hold the same PID. Compose
# the SLURM array identity with the PID so the tag is unique cluster-wide.
# Identical to TASK_TAG in scripts/hpc/run_selection_stack.sh (bead -ih5u).
TASK_TAG="${SLURM_JOB_ID:-nojob}.${SLURM_ARRAY_TASK_ID:-0}.$$"

# Check dependency
check_file "${RESULTS_DIR}/step_completed_lse_classification.txt"

log "Starting selective pressure and ASR analysis."

# --- Validate required tools ---
HYPHY_AVAILABLE=true
if ! command -v hyphy &>/dev/null; then
    log --level=WARN "HyPhy not found - dN/dS analysis will be skipped"
    HYPHY_AVAILABLE=false
fi

FASTML_AVAILABLE=true
if ! command -v "${FASTML:-fastml}" &>/dev/null; then
    log --level=WARN "FastML not found - ASR analysis will be skipped"
    FASTML_AVAILABLE=false
fi

if [ "$HYPHY_AVAILABLE" = false ] && [ "$FASTML_AVAILABLE" = false ]; then
    log --level=WARN "Neither HyPhy nor FastML available - step 05 will produce minimal output"
fi

# Get orthogroup count using manifest (preferred) or fallback to globbing
MANIFEST_FILE="${RESULTS_DIR}/orthogroup_manifest.tsv"
OG_COUNT=$(get_orthogroup_count "$MANIFEST_FILE")

if [ "$OG_COUNT" -eq 0 ]; then
    log "Error: No orthogroups found"
    exit 1
fi

# Handle SLURM array indexing using helper function
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    validate_array_index "$SLURM_ARRAY_TASK_ID" "$OG_COUNT"
    # Fail loud if the submitted array is too small to cover the manifest —
    # otherwise OGs at index > the array top are silently skipped (bead fxx).
    assert_array_covers_manifest "$OG_COUNT" || exit 1

    # Get orthogroup name from manifest or fallback to globbing
    if [ -f "$MANIFEST_FILE" ]; then
        base=$(get_orthogroup_by_index "$SLURM_ARRAY_TASK_ID" "$MANIFEST_FILE")
    else
        # Fallback: enumerate the authoritative run's Orthogroup_Sequences/
        # (scripts/orthofinder_paths.sh). The old glob pointed at
        # orthogroups/OrthoFinder/Results*/Orthogroups/ — wrong on BOTH
        # components (the real root is orthogroups/input/OrthoFinder, and
        # per-OG FASTAs live in Orthogroup_Sequences/, not Orthogroups/), so it
        # never matched and basename turned the literal pattern into "$base".
        _og_seq_dir=$(resolve_orthogroup_sequences_dir "${RESULTS_DIR}/orthogroups" || true)
        [ -n "$_og_seq_dir" ] || { log "Error: no OrthoFinder Orthogroup_Sequences directory"; exit 1; }
        # Deterministic, locale-independent order. The array index selects
        # POSITIONALLY and the 50 concurrent tasks run on different nodes, so
        # index N must denote the SAME orthogroup everywhere; sort under LC_ALL=C
        # rather than trusting filesystem traversal order or the ambient locale.
        mapfile -t ORTHOGROUPS < <(LC_ALL=C find "${_og_seq_dir}" -maxdepth 1 -type f -name 'OG*.fa' | LC_ALL=C sort)
        # A zero-match enumeration must ABORT, never iterate zero times. Bash
        # leaves an unmatched glob as its literal pattern (nullglob is off), so
        # the previous array-glob produced the orthogroup name "OG*" and the
        # stage ran a full analysis against a nonexistent orthogroup, exiting 0.
        [ "${#ORTHOGROUPS[@]}" -gt 0 ] || { log "Error: ${_og_seq_dir} contains no OG*.fa"; exit 1; }
        [ "$SLURM_ARRAY_TASK_ID" -lt "${#ORTHOGROUPS[@]}" ] || { log "Error: array index ${SLURM_ARRAY_TASK_ID} is past the ${#ORTHOGROUPS[@]} orthogroups in ${_og_seq_dir}"; exit 1; }
        og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
        base=$(basename "$og" .fa)
    fi
else
    # Non-array mode: process first orthogroup (for testing)
    if [ -f "$MANIFEST_FILE" ]; then
        base=$(get_orthogroup_by_index 0 "$MANIFEST_FILE")
    else
        _og_seq_dir=$(resolve_orthogroup_sequences_dir "${RESULTS_DIR}/orthogroups" || true)
        [ -n "$_og_seq_dir" ] || { log "Error: no OrthoFinder Orthogroup_Sequences directory"; exit 1; }
        # Same deterministic order and loud zero-match guard as the array
        # branch above, so non-array test runs pick the same OG index 0 does.
        mapfile -t ORTHOGROUPS < <(LC_ALL=C find "${_og_seq_dir}" -maxdepth 1 -type f -name 'OG*.fa' | LC_ALL=C sort)
        [ "${#ORTHOGROUPS[@]}" -gt 0 ] || { log "Error: ${_og_seq_dir} contains no OG*.fa"; exit 1; }
        base=$(basename "${ORTHOGROUPS[0]}" .fa)
    fi
fi

# Find the orthogroup FASTA file.
# Anchored to the authoritative run's Orthogroup_Sequences/ (see
# scripts/orthofinder_paths.sh). The previous unsorted `find ... | head -1`
# could return the WorkingDirectory/ copy, whose headers are OrthoFinder's
# internal integer ids — get_taxids_from_fasta then yields nothing, taxa_count
# is 0, and the `-gt 1` gate below SILENTLY skips selection analysis for that
# orthogroup.
og=$(resolve_orthogroup_fasta "${base}" "${RESULTS_DIR}/orthogroups" || true)
[ -z "$og" ] || [ ! -f "$og" ] && { log "Skipping missing orthogroup: ${base}"; exit 0; }
# Use centralized metadata lookup to get taxids (excludes references)
taxids=$(get_taxids_from_fasta "$og" --exclude-refs)
taxa_count=$(echo "$taxids" | wc -w)

# --- Define reference CDS file location ---
REFERENCE_CDS_FILE="${RESULTS_DIR}/reference_sequences/cds/all_references_cds.fna"

# find_nucleotide_sequences() now lives in functions.sh (moved for unit-test
# coverage; bead sz6). It reads the REFERENCE_CDS_FILE global set just above.

# --- Codon alignment validation ---
# Validates pal2nal output for proper PAML format and codon structure
# Arguments: $1 - codon alignment file (PAML format)
#            $2 - expected sequence count (optional)
# Returns: 0 if valid, 1 if invalid
validate_codon_alignment() {
    local codon_file="$1"
    local expected_count="${2:-0}"

    if [ ! -f "$codon_file" ] || [ ! -s "$codon_file" ]; then
        return 1
    fi

    # Validate PAML format using Python for robust checking
    python3 << PYTHON_SCRIPT
import sys

codon_file = "${codon_file}"
expected_count = ${expected_count}

try:
    with open(codon_file, 'r') as f:
        first_line = f.readline().strip()

        # Bead 771: MACSE (default RUN_MACSE=1) emits a FASTA codon alignment;
        # pal2nal emits PAML. HyPhy reads FASTA natively, so accept either —
        # the old PAML-only check failed on MACSE output and silently skipped
        # the ENTIRE selection stack (GARD/BUSTED/aBSREL/MEME) for every OG.
        if first_line.startswith('>'):
            f.seek(0)
            seqs = {}
            name = None
            buf = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if name is not None:
                        seqs[name] = ''.join(buf)
                    name = line[1:].split()[0] if len(line) > 1 else ''
                    buf = []
                elif name is not None and line:
                    buf.append(line.replace(' ', ''))
            if name is not None:
                seqs[name] = ''.join(buf)
            if len(seqs) < 2:
                print(f"FASTA codon alignment has <2 sequences ({len(seqs)})", file=sys.stderr)
                sys.exit(1)
            lengths = {len(s) for s in seqs.values()}
            if len(lengths) != 1:
                print(f"FASTA codon alignment not aligned (lengths {sorted(lengths)})", file=sys.stderr)
                sys.exit(1)
            nalign = lengths.pop()
            if nalign == 0 or nalign % 3 != 0:
                print(f"FASTA codon alignment length {nalign} not a positive multiple of 3", file=sys.stderr)
                sys.exit(1)
            if expected_count > 0 and len(seqs) != expected_count:
                print(f"Sequence count mismatch: expected {expected_count}, got {len(seqs)}", file=sys.stderr)
                sys.exit(1)
            print(f"Codon alignment valid (FASTA): {len(seqs)} sequences, {nalign} bp ({nalign//3} codons)", file=sys.stderr)
            sys.exit(0)

        parts = first_line.split()

        # PAML format: first line is "nseqs nalign"
        if len(parts) != 2:
            print(f"Invalid PAML header: expected 'nseqs nalign' or FASTA, got '{first_line}'", file=sys.stderr)
            sys.exit(1)

        try:
            nseqs, nalign = int(parts[0]), int(parts[1])
        except ValueError:
            print(f"Invalid PAML header values: {parts}", file=sys.stderr)
            sys.exit(1)

        # Check alignment length is divisible by 3 (codon-based)
        if nalign % 3 != 0:
            print(f"Alignment length {nalign} not divisible by 3 (not codon-aligned)", file=sys.stderr)
            sys.exit(1)

        # Check expected count if provided
        if expected_count > 0 and nseqs != expected_count:
            print(f"Sequence count mismatch: expected {expected_count}, got {nseqs}", file=sys.stderr)
            sys.exit(1)

        # Count actual sequences and validate lengths
        content = f.read()
        sequences = {}
        current_name = None
        current_seq = []

        for line in content.split('\n'):
            line = line.strip()
            if not line:
                continue
            # In PAML, sequence name is on its own line, followed by sequence
            if line and not current_name:
                current_name = line
            elif current_name:
                current_seq.append(line.replace(' ', ''))
                # Check if we have accumulated the full sequence
                full_seq = ''.join(current_seq)
                if len(full_seq) >= nalign:
                    sequences[current_name] = full_seq[:nalign]
                    current_name = None
                    current_seq = []

        if len(sequences) != nseqs:
            print(f"Parsed {len(sequences)} sequences but header says {nseqs}", file=sys.stderr)
            sys.exit(1)

        # Check each sequence length
        stop_codons = {'TAA', 'TAG', 'TGA', 'taa', 'tag', 'tga'}
        for name, seq in sequences.items():
            if len(seq) != nalign:
                print(f"Sequence {name} length {len(seq)} != expected {nalign}", file=sys.stderr)
                sys.exit(1)

            # Check for internal stop codons (warning only, not failure)
            seq_upper = seq.upper().replace('-', 'N')
            for i in range(0, len(seq_upper) - 3, 3):  # Exclude last codon
                codon = seq_upper[i:i+3]
                if codon in {'TAA', 'TAG', 'TGA'}:
                    print(f"Warning: Internal stop codon {codon} at position {i} in {name}", file=sys.stderr)

        # All checks passed
        print(f"Codon alignment valid: {nseqs} sequences, {nalign} bp ({nalign//3} codons)", file=sys.stderr)
        sys.exit(0)

except Exception as e:
    print(f"Validation error: {e}", file=sys.stderr)
    sys.exit(1)
PYTHON_SCRIPT

    return $?
}

# Bead 1z7: rebuild a cumulative CSV from all per-OG CSVs in <dir> matching
# *<suffix>. Idempotent and safe under SLURM arrays — each task writes only its
# own per-OG csv then re-runs this concat, so the last task yields the complete
# cumulative. Header is taken from the first file only.
#
# Staging is namespaced by ${TASK_TAG}, NOT by $$. The original `.tmp.$$` was
# cited as this project's canonical concurrency fix but is insufficient: PIDs
# are unique per node while <dir> is shared storage, so two of the 50 concurrent
# array tasks can collide on the same staging path — task A's `: >` truncates
# the file task B is mid-append to, and B publishes the truncated result. The
# final `mv` is atomic within the directory, so readers of the cumulative never
# observe a partial file.
#
# The per-OG CSVs this globs need their OWN atomicity guarantee. Staging
# protects the OUTPUT of this loop; a non-atomically-written per-OG CSV
# corrupts its INPUT, and no amount of output staging can recover from that.
# A plain open(path,"w") truncates the published file at open and refills it
# incrementally, so a concurrent task globbing this directory can read an
# empty or half-written per-OG CSV and publish it into the cumulative.
# Verified status of the three writers (do not assume — this comment
# previously asserted an atomicity nobody had checked, and it was false):
#   parse_busted.py  ATOMIC (write_csv_atomic: unique temp + os.replace)
#   parse_meme.py    ATOMIC (write_csv_atomic: unique temp + os.replace)
#   parse_absrel.py  ATOMIC (write_csv_atomic: unique temp + os.replace).
#                    Repaired after this note first recorded it as unsafe; its
#                    docstring now distinguishes per-OG isolation from write
#                    atomicity, which is what the two words used to conflate.
concat_per_og_csv() {
    local dir="$1" suffix="$2" cumulative="$3"
    # The fallback recomposes the SLURM identity rather than degrading to a
    # bare $$, so the guarantee holds even if this function is sourced
    # somewhere TASK_TAG was never set.
    local tmp="${cumulative}.tmp.${TASK_TAG:-${SLURM_JOB_ID:-nojob}.${SLURM_ARRAY_TASK_ID:-0}.$$}"
    : > "$tmp"
    local first=1 csv
    for csv in "$dir"/*"$suffix"; do
        [ -f "$csv" ] || continue
        if [ "$first" = 1 ]; then cat "$csv" >> "$tmp"; first=0; else tail -n +2 "$csv" >> "$tmp"; fi
    done
    if [ -s "$tmp" ]; then mv "$tmp" "$cumulative"; else rm -f "$tmp"; fi
}

# --- Per-OG path resolution (mirrors stage 04's per-class routing) ---
# Bead F1. Stage 04's per-class refactor writes every per-orthogroup product
# TWO LEVELS deeper than stage 05 used to look:
#     ${RESULTS_DIR}/phylogenies/protein/class_<CLASS>/<base>/<base>_trimmed.fa
#     ${RESULTS_DIR}/phylogenies/protein/class_<CLASS>/<base>/<base>.treefile
# (see 04:535 OG_OUT_DIR, 04:555 ClipKit output, 04:583 IQ-TREE --prefix).
# Stage 05 read the pre-refactor FLAT pair, which has no writer anywhere in the
# repository, so its `[ -f ... ] && [ -f ... ]` guard was NEVER true: no aBSREL
# / GARD / BUSTED-S / BUSTED-MH / MEME / ASR output was ever produced, while
# step_completed_05.txt was still touched and stage 07 ran on regardless.
#
# The class is resolved exactly the way stage 04 does it (04:524-534): majority
# class from results/classification/og_class_majority.tsv, with stage 04's two
# fallbacks preserved verbatim -- OG absent from the TSV -> "unclassified",
# TSV absent entirely -> "A" (back-compat).
#
# Sets globals OG_CLASS, protein_align and tree. Returns 0 only when BOTH
# products exist; on failure protein_align/tree are left empty so no caller can
# accidentally run HyPhy against a stale path.
resolve_perog_paths() {
    local og_base="$1"
    local og_class_tsv="${RESULTS_DIR}/classification/og_class_majority.tsv"
    local og_dir cand found_dir
    local -a matches=()

    protein_align=""
    tree=""

    # --- class resolution: identical to stage 04:524-534 ---
    if [ -f "${og_class_tsv}" ]; then
        OG_CLASS=$(awk -F'\t' -v og="${og_base}" 'NR>1 && $1==og {print $2; exit}' "${og_class_tsv}")
        if [ -z "${OG_CLASS}" ]; then
            log --level=WARN "OG '${og_base}' not found in ${og_class_tsv}; reading from class_unclassified"
            OG_CLASS="unclassified"
        fi
    else
        log "WARN: ${og_class_tsv} not found; reading ${og_base} from class_A (back-compat)"
        OG_CLASS="A"
    fi

    og_dir="${RESULTS_DIR}/phylogenies/protein/class_${OG_CLASS}/${og_base}"
    if [ -f "${og_dir}/${og_base}_trimmed.fa" ] && [ -f "${og_dir}/${og_base}.treefile" ]; then
        protein_align="${og_dir}/${og_base}_trimmed.fa"
        tree="${og_dir}/${og_base}.treefile"
        return 0
    fi

    # --- the products are not where the class map says they should be ---
    # The old code logged one undifferentiated "Missing alignment or tree"
    # warning here. That warning fired for EVERY orthogroup, which is exactly
    # why a whole broken stage went unnoticed: a message that is always printed
    # carries no information. Split it by severity instead, so the log
    # distinguishes "stage 04 legitimately built no tree for this OG" (normal,
    # quiet) from "the path is wrong again" (loud).
    for cand in "${RESULTS_DIR}/phylogenies/protein/class_"*"/${og_base}/${og_base}.treefile"; do
        [ -f "$cand" ] && matches+=("$cand")
    done

    if [ "${#matches[@]}" -eq 0 ]; then
        # Nothing anywhere. Stage 04 skips orthogroups it cannot build a tree
        # for, so this is an expected, uninteresting outcome -- INFO, not WARN.
        log "No per-OG tree for ${og_base} (class ${OG_CLASS}); stage 04 built none. Skipping."
        return 1
    fi

    if [ "${#matches[@]}" -gt 1 ]; then
        # The same orthogroup under two classes. Choosing one would be a guess
        # about which run is authoritative, and a wrong guess silently attaches
        # selection results to the wrong tree. Refuse.
        log --level=ERROR "Orthogroup ${og_base} has trees under MULTIPLE classes (${matches[*]}); refusing to guess. Rebuild stage 04 outputs or prune the stale class directory."
        return 1
    fi

    found_dir=$(dirname "${matches[0]}")
    if [ "${found_dir}" != "${og_dir}" ]; then
        # Stage 04 and stage 05 disagree about this OG's class. Happens when
        # og_class_majority.tsv is (re)generated after stage 04 already ran
        # under the class_A back-compat fallback. The single match is
        # unambiguous, so use it rather than silently producing nothing -- but
        # shout, because the two stages disagreeing is a real inconsistency.
        log --level=ERROR "Class-routing mismatch for ${og_base}: class map says '${OG_CLASS}' (${og_dir}) but stage 04 wrote ${found_dir}. Using the found path; re-run stage 04 to make the layout consistent."
    fi

    if [ ! -f "${found_dir}/${og_base}_trimmed.fa" ]; then
        # Tree present, trimmed alignment gone: stage 04 got as far as IQ-TREE
        # but its ClipKit product is missing. A genuine per-OG defect.
        log --level=WARN "Orthogroup ${og_base}: tree found at ${matches[0]} but no ${og_base}_trimmed.fa beside it; skipping selection/ASR."
        return 1
    fi

    protein_align="${found_dir}/${og_base}_trimmed.fa"
    tree="${matches[0]}"
    return 0
}

# --- Selective Pressure with aBSREL ---
if [ "$taxa_count" -gt 1 ] && [ "$HYPHY_AVAILABLE" = true ]; then
    if resolve_perog_paths "${base}"; then
        # Find nucleotide sequences for this orthogroup
        nuc_align="${RESULTS_DIR}/selective_pressure/nucleotide/${base}_nuc.fa"

        if find_nucleotide_sequences "$protein_align" "$nuc_align"; then
            # Create codon alignment
            codon_file="${RESULTS_DIR}/selective_pressure/${base}_codon.phy"
            # Bead -i61: codon-aware MSA via MACSE v2 (preferred for
            # frame-fragile miniprot-recovered CDS), falling back to pal2nal
            # naive back-translation if MACSE unavailable / disabled.
            if [ "${RUN_MACSE:-1}" = "1" ] \
               && [ -n "${MACSE_JAR:-}" ] && [ -f "${MACSE_JAR:-}" ]; then
                bash "${SCRIPTS_DIR}/run_macse.sh" \
                    --input="$nuc_align" --output="$codon_file" \
                    2>> "${LOGS_DIR}/macse_${base}.err" \
                    || run_command "${base}_codon" --stdout="$codon_file" \
                        pal2nal.pl "$protein_align" "$nuc_align" -output paml
            else
                run_command "${base}_codon" --stdout="$codon_file" pal2nal.pl "$protein_align" "$nuc_align" -output paml
            fi

            # Validate codon alignment (checks PAML format, codon structure, internal stops)
            if validate_codon_alignment "$codon_file"; then
                # Bead -urk: SELECTION_BACKEND=stack (default) runs the modern
                # HyPhy stack: GARD → BUSTED-S → BUSTED-MH → aBSREL → MEME.
                # SELECTION_BACKEND=absrel falls back to legacy aBSREL-only
                # (~5x faster, but no recombination screen / synonymous-rate
                # variation / multi-hit correction / site-level inference).
                if [ "${SELECTION_BACKEND:-stack}" = "stack" ]; then
                    bash "${SCRIPTS_DIR}/hpc/run_selection_stack.sh" \
                        "$codon_file" "$tree" "${base}" \
                        2>> "${LOGS_DIR}/selection_stack_${base}.err" \
                        || log --level=WARN "Selection stack failed for ${base} (continuing — partial outputs may be present)"
                else
                    run_command "${base}_absrel" hyphy aBSREL --alignment "$codon_file" --tree "$tree" --output "${RESULTS_DIR}/selective_pressure/${base}_absrel.json"
                fi

                # Parse aBSREL JSON (always — both backends produce it)
                if [ -f "${RESULTS_DIR}/selective_pressure/${base}_absrel.json" ]; then
                    # Bead 1z7: write a PER-OG csv (not the shared cumulative) so
                    # parallel array tasks don't clobber each other (parse_absrel
                    # default 'atomic' opens the target with mode 'w'); rebuild the
                    # cumulative absrel_results.csv from all per-OG csvs afterward.
                    python3 "${SCRIPTS_DIR}/parse_absrel.py" "${RESULTS_DIR}/selective_pressure/${base}_absrel.json" "${RESULTS_DIR}/selective_pressure/${base}_absrel.csv" || log "Warning: Failed to parse aBSREL for $base"
                    concat_per_og_csv "${RESULTS_DIR}/selective_pressure" "_absrel.csv" "${RESULTS_DIR}/selective_pressure/absrel_results.csv"
                fi
            else
                log "Warning: pal2nal produced invalid codon alignment for ${base}, skipping selection analysis"
            fi
        else
            log "Warning: Could not find nucleotide sequences for ${base}, skipping dN/dS analysis"
        fi
    fi
    # No `else` here on purpose: resolve_perog_paths already logged a message
    # specific to WHY the products are unavailable (none built / class-routing
    # mismatch / ambiguous / incomplete). A second, undifferentiated "Missing
    # alignment or tree" line would just restore the noise that hid bead F1.
fi

# --- Berghia-specific LSEs with ASR ---
# Check if this orthogroup contains Berghia sequences AND has multiple sequences (for meaningful ASR)
berghia_count=$(count_taxid_occurrences "$taxids" "${BERGHIA_TAXID}")
seq_count=$(grep -c "^>" "$og")

if [ "$berghia_count" -gt 0 ] && [ "$seq_count" -gt 2 ]; then
    # Bead F1: same per-class resolution as the selection branch above. This
    # branch carried an identical copy of the broken flat pair, so ASR was a
    # no-op for every orthogroup too. Re-resolving (rather than reusing the
    # values from the selection branch) keeps the two independent: the
    # selection branch is gated on HYPHY_AVAILABLE and taxa_count, so it may
    # never have run.
    if resolve_perog_paths "${base}"; then
        # Find nucleotide sequences
        nuc_align="${RESULTS_DIR}/selective_pressure/nucleotide/${base}_nuc.fa"

        if [ ! -f "$nuc_align" ]; then
            find_nucleotide_sequences "$protein_align" "$nuc_align"
        fi

        if [ -f "$nuc_align" ] && [ -s "$nuc_align" ]; then
            # Run LSE-specific aBSREL if not already done
            if [ ! -f "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json" ]; then
                codon_file="${RESULTS_DIR}/selective_pressure/${base}_codon.phy"
                if [ ! -f "$codon_file" ]; then
                    run_command "${base}_codon_lse" --stdout="$codon_file" pal2nal.pl "$protein_align" "$nuc_align" -output paml
                fi

                # Validate codon alignment before running aBSREL
                if validate_codon_alignment "$codon_file"; then
                    run_command "${base}_absrel_lse" hyphy aBSREL --alignment "$codon_file" --tree "$tree" --output "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json"
                    if [ -f "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json" ]; then
                        # Bead 1z7: per-OG csv + cumulative rebuild (see main path above).
                        python3 "${SCRIPTS_DIR}/parse_absrel.py" "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.json" "${RESULTS_DIR}/selective_pressure/${base}_absrel_lse.csv"
                        concat_per_og_csv "${RESULTS_DIR}/selective_pressure" "_absrel_lse.csv" "${RESULTS_DIR}/selective_pressure/absrel_results_lse.csv"
                    fi
                else
                    log "Warning: Invalid codon alignment for ${base} LSE analysis"
                fi
            fi

            # ASR for deep nodes — bead -mqt: pass full ${BERGHIA_FILE_PREFIX}
            # (e.g. "1287507_berghia_stephanieae") so Berghia leaves match by
            # exact prefix, not just taxid alone (which could spuriously match
            # other species' references that happen to share an underscore-
            # delimited substring).
            # Bead -m6k: prefer null-calibrated cutoff if available.
            DEPTH_CALIB_FILE="${RESULTS_DIR}/calibration/depth_thresholds.json"
            DEPTH_CALIB_ARG=""
            if [ -f "$DEPTH_CALIB_FILE" ]; then
                DEPTH_CALIB_ARG="--null-threshold-file $DEPTH_CALIB_FILE"
            fi
            deep_nodes=$(python3 "${SCRIPTS_DIR}/select_deep_nodes.py" \
                "$tree" "${BERGHIA_FILE_PREFIX}" "${MIN_ASR_DISTANCE}" \
                $DEPTH_CALIB_ARG 2>/dev/null | head -1)

            if [ -n "$deep_nodes" ]; then
                # Bead 4z1: ASR requires a multiple-sequence ALIGNMENT, but
                # nuc_align is the UNALIGNED recovered CDS. Prefer the codon
                # alignment built earlier this task (same taxa as the tree),
                # falling back to nuc_align only if it is missing/invalid.
                asr_codon="${RESULTS_DIR}/selective_pressure/${base}_codon.phy"
                if [ -f "$asr_codon" ] && validate_codon_alignment "$asr_codon" 2>/dev/null; then
                    asr_input="$asr_codon"
                else
                    asr_input="$nuc_align"
                    log --level=WARN "ASR for ${base}: codon alignment unavailable; using unaligned nucleotides (ASR may be unreliable)"
                fi
                # Bead -j44: switch from FastML to IQ-TREE --ancestral (model-
                # consistent with the inference tree, scales to thousands of
                # taxa, avoids re-running ML under a different model). FastML
                # is kept as a fallback when ASR_BACKEND=fastml.
                ASR_BACKEND="${ASR_BACKEND:-iqtree}"
                if [ "$ASR_BACKEND" = "iqtree" ]; then
                    # IQ-TREE empirical-Bayes marginal ASR using the tree's
                    # checkpoint. The .state file contains per-site per-node
                    # marginal posterior probabilities for every state.
                    asr_prefix="${RESULTS_DIR}/asr/${base}_asr"
                    mkdir -p "${RESULTS_DIR}/asr"
                    run_command "${base}_asr_iqtree" ${IQTREE} \
                        -s "$asr_input" -te "$tree" \
                        --ancestral -seed "${IQTREE_SEED}" -T "${SLURM_CPUS_PER_TASK:-${CPUS}}" \
                        --prefix "$asr_prefix" 2>/dev/null \
                        || log "Warning: IQ-TREE ASR failed for ${base}"
                    # Extract per-deep-node ancestral sequences for downstream use.
                    if [ -f "${asr_prefix}.state" ]; then
                        for node in $deep_nodes; do
                            # "$node" is the .state lookup key (IQ-TREE writes
                            # the bare label). --record-id is the emitted FASTA
                            # header: IQ-TREE numbers internal nodes per tree,
                            # so a bare "NodeN" is only unique within THIS
                            # orthogroup, and stage 08 pools every *_asr.fa into
                            # one extraction source. Namespacing by ${base}
                            # matches the filename convention and keeps the
                            # pooled ids unambiguous.
                            python3 "${SCRIPTS_DIR}/extract_iqtree_asr.py" \
                                "${asr_prefix}.state" "$node" \
                                "${RESULTS_DIR}/asr/${base}_${node}_asr.fa" \
                                --record-id "${base}_${node}" \
                                2>/dev/null \
                                || log "Warning: ASR extract failed for ${base} ${node}"
                        done
                    fi
                else
                    # KNOWN LIMITATION (ASR_BACKEND=fastml only; iqtree is the
                    # default and is unaffected). Two separate problems here,
                    # neither fixed by the --record-id namespacing above --
                    # FastML writes its own headers, so the extractor that
                    # honours --record-id is not in this code path at all:
                    #
                    #  1. The invocation does not match the real FastML CLI.
                    #     Per the fastml(1) manual the binary takes single-dash
                    #     options only (-s sequence, -t tree, -x/-y tree out,
                    #     -j/-k sequence out); --seq/--tree/--out_seq/--out_tree/
                    #     --node/--verbose do not exist, and `-t` here is handed
                    #     a CPU count where FastML expects the tree file. The
                    #     call therefore fails, no *_asr.fa is written, and
                    #     stage 08 silently folds extant candidates only.
                    #  2. FastML has no single-node selection flag: it writes
                    #     EVERY internal node into one file (seq.marginal.txt /
                    #     seq.joint.txt) under its own labels (N1, N2, ...).
                    #     So if (1) were fixed naively, each iteration of this
                    #     per-node loop would write the same all-nodes file, and
                    #     those labels are per-tree -- ids would collide both
                    #     within an orthogroup and across orthogroups, which is
                    #     exactly what stage 08's duplicate-id guard rejects.
                    #
                    # Fixing this means choosing a FastML output contract (run
                    # once per orthogroup, then split seq.marginal.txt and
                    # namespace each record as ${base}_<fastml-label>), which is
                    # a pipeline-behaviour decision, not a mechanical fix.
                    for node in $deep_nodes; do
                        run_command "${base}_${node}_asr" ${FASTML} --seq "$asr_input" --tree "$tree" \
                            --out_seq "${RESULTS_DIR}/asr/${base}_${node}_asr.fa" \
                            --out_tree "${RESULTS_DIR}/asr/${base}_${node}_asr.tree" \
                            --node "$node" -t "${SLURM_CPUS_PER_TASK:-${CPUS}}" --verbose 2>/dev/null \
                            || log "Warning: FastML failed for ${base} node ${node}"
                    done
                fi

                # Plot ASR for first deep node
                first_node="${deep_nodes%% *}"
                if [ -f "${RESULTS_DIR}/asr/${base}_${first_node}_asr.fa" ]; then
                    python3 "${SCRIPTS_DIR}/plot_asr.py" "$tree" "${RESULTS_DIR}/asr/${base}_${first_node}_asr.fa" "${RESULTS_DIR}/asr/${base}_asr_plot" || log "Warning: ASR plotting failed for ${base}"
                fi
            fi
        fi
    fi
fi

# --- Coverage accounting (bead higb) ----------------------------------------
# Count PRODUCTS ON DISK rather than trusting control flow. Bead F1 (stage 05
# reading a path nothing wrote) stayed invisible for months precisely because
# the stage reported what it BELIEVED it had done: both completion markers were
# touched unconditionally, so a run in which every orthogroup resolved to
# nothing still reported success and let stage 07 proceed on an empty table.
# "Complete" without coverage numbers is an unverified claim.
selection_ran=0
[ -s "${RESULTS_DIR}/selective_pressure/${base}_absrel.csv" ] && selection_ran=1
asr_ran=0
ls "${RESULTS_DIR}/asr/${base}"_*_asr.fa >/dev/null 2>&1 && asr_ran=1
og_products=$((selection_ran + asr_ran))

# One record per orthogroup, never a shared append: stage 05 is a SLURM array
# and concurrent tasks appending to one file is a write race.
coverage_dir="${RESULTS_DIR}/selective_pressure/coverage"
mkdir -p "$coverage_dir"
{
    printf 'orthogroup\tselection_ran\tasr_ran\tproducts\ttaxa_count\n'
    printf '%s\t%s\t%s\t%s\t%s\n' "$base" "$selection_ran" "$asr_ran" \
        "$og_products" "$taxa_count"
} > "${coverage_dir}/${base}.tsv"

if [ "$og_products" -eq 0 ]; then
    # Explicit, and deliberately NOT a completion. Downstream must be able to
    # tell "this orthogroup produced nothing" from "this orthogroup was fine".
    log "ERROR: ${base} produced nothing (selection=0 asr=0); refusing to mark it complete"
else
    log "Coverage for ${base}: selection=${selection_ran} asr=${asr_ran}"
    # Create completion flag for this orthogroup
    touch "${RESULTS_DIR}/selective_pressure/step_completed_${base}.txt"
fi

# Create array checkpoint for resume capability. This records that the TASK
# ran, which is independent of whether the orthogroup yielded products, so it
# stays outside the gate: resume must not re-run a task that already completed.
[ -n "$SLURM_ARRAY_TASK_ID" ] && create_array_checkpoint "05_selective" "$SLURM_ARRAY_TASK_ID"

log "Selective pressure and ASR finished for ${base} (products=${og_products})."

# Create overall step completion flag
# Note: We create this flag after each task completes, since downstream steps
# can start processing as results become available. Gated on this task having
# produced something, so a run where every orthogroup fails leaves no marker
# at all and stage 07 cannot mistake emptiness for completion.
if [ "$og_products" -gt 0 ]; then
    touch "${RESULTS_DIR}/step_completed_05.txt"
fi
