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

# Create output directories
mkdir -p "${RESULTS_DIR}/selective_pressure" "${RESULTS_DIR}/selective_pressure/nucleotide" "${RESULTS_DIR}/asr" "${LOGS_DIR}" || { log "Error: Cannot create directories"; exit 1; }

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
        ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)
        og="${ORTHOGROUPS[$SLURM_ARRAY_TASK_ID]}"
        base=$(basename "$og" .fa)
    fi
else
    # Non-array mode: process first orthogroup (for testing)
    if [ -f "$MANIFEST_FILE" ]; then
        base=$(get_orthogroup_by_index 0 "$MANIFEST_FILE")
    else
        ORTHOGROUPS=("${RESULTS_DIR}/orthogroups/OrthoFinder/Results"*/Orthogroups/OG*.fa)
        base=$(basename "${ORTHOGROUPS[0]}" .fa)
    fi
fi

# Find the orthogroup FASTA file
og=$(find "${RESULTS_DIR}/orthogroups" -name "${base}.fa" -type f 2>/dev/null | head -1)
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
# cumulative. Header is taken from the first file only; per-task temp avoids
# clobbering between concurrent tasks.
concat_per_og_csv() {
    local dir="$1" suffix="$2" cumulative="$3"
    local tmp="${cumulative}.tmp.$$"
    : > "$tmp"
    local first=1 csv
    for csv in "$dir"/*"$suffix"; do
        [ -f "$csv" ] || continue
        if [ "$first" = 1 ]; then cat "$csv" >> "$tmp"; first=0; else tail -n +2 "$csv" >> "$tmp"; fi
    done
    if [ -s "$tmp" ]; then mv "$tmp" "$cumulative"; else rm -f "$tmp"; fi
}

# --- Selective Pressure with aBSREL ---
if [ "$taxa_count" -gt 1 ] && [ "$HYPHY_AVAILABLE" = true ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    tree="${RESULTS_DIR}/phylogenies/protein/${base}.treefile"

    if [ -f "$protein_align" ] && [ -f "$tree" ]; then
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
    else
        log "Warning: Missing alignment or tree for ${base}"
    fi
fi

# --- Berghia-specific LSEs with ASR ---
# Check if this orthogroup contains Berghia sequences AND has multiple sequences (for meaningful ASR)
berghia_count=$(count_taxid_occurrences "$taxids" "${BERGHIA_TAXID}")
seq_count=$(grep -c "^>" "$og")

if [ "$berghia_count" -gt 0 ] && [ "$seq_count" -gt 2 ]; then
    protein_align="${RESULTS_DIR}/phylogenies/protein/${base}_trimmed.fa"
    tree="${RESULTS_DIR}/phylogenies/protein/${base}.treefile"

    if [ -f "$protein_align" ] && [ -f "$tree" ]; then
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
                            python3 "${SCRIPTS_DIR}/extract_iqtree_asr.py" \
                                "${asr_prefix}.state" "$node" \
                                "${RESULTS_DIR}/asr/${base}_${node}_asr.fa" \
                                2>/dev/null \
                                || log "Warning: ASR extract failed for ${base} ${node}"
                        done
                    fi
                else
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

# Create completion flag for this orthogroup
touch "${RESULTS_DIR}/selective_pressure/step_completed_${base}.txt"

# Create array checkpoint for resume capability
[ -n "$SLURM_ARRAY_TASK_ID" ] && create_array_checkpoint "05_selective" "$SLURM_ARRAY_TASK_ID"

log "Selective pressure and ASR completed for ${base}."

# Create overall step completion flag
# Note: We create this flag after each task completes, since downstream steps
# can start processing as results become available
touch "${RESULTS_DIR}/step_completed_05.txt"
