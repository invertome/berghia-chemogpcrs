# orthofinder_paths.sh — ONE deterministic rule for "which OrthoFinder run is
# authoritative", and for the paths derived from it.
#
# Sourced by the stages that read OrthoFinder output (03, 03b, 04, 05, 06c, 07).
# Not executable on its own; defines functions only, with no side effects, so it
# is safe to source anywhere and testable in isolation
# (tests/unit/test_ogpath_orthofinder_run.py).
#
# ---------------------------------------------------------------------------
# WHY THIS EXISTS
#
# Five call sites each rolled their own resolution and, on a real results tree,
# resolved DIFFERENTLY — so the same orthogroup id denoted different gene sets
# in different stages while every join still succeeded. Three failure modes:
#
#   1. `find ... -name Orthogroups.tsv | head -1` — unsorted. `find` order is
#      filesystem-dependent, so with more than one Results_* run present the
#      winner was arbitrary and not even stable between invocations.
#
#   2. `find ... -name "Results_*" | sort -r | head -1` — lexicographic on
#      OrthoFinder's `Results_<Mon><DD>` stamp, which is neither chronological
#      nor numeric: "Results_Sep05" sorts ABOVE "Results_Dec01" (S > D) even
#      though December is later, and same-day rerun suffixes sort "_10" below
#      "_2".
#
#   3. Neither checked that the run it picked actually produced results, so a
#      crashed/half-written run — typically the NEWEST directory on disk —
#      could win.
#
# THE RULE: the run whose Orthogroups/Orthogroups.tsv has the newest mtime,
# ties broken by reverse path order.
#
# Why mtime is the chronologically correct key: `Results_<Mon><DD>` carries NO
# YEAR, so no name-based rule can order a December run against the following
# January's — the information simply is not in the name. mtime is the only
# signal on disk that is actually chronological, and it is monotonic across
# reruns, resumes (`orthofinder -fg`, which rewrites Orthogroups.tsv) and
# same-day repeats alike. Keying on Orthogroups.tsv rather than the directory
# means a run must have PRODUCED the table to be eligible, which rules out
# failure mode 3; the reverse-path tie-break keeps two runs written in the same
# second from resolving differently on different invocations.
#
# Requires GNU find (`-printf`), which is what the pipeline runs on.
# ---------------------------------------------------------------------------

# --- Resolve the authoritative OrthoFinder run directory ---
# Arguments: $1 - search root (default: ${RESULTS_DIR}/orthogroups)
# Prints the Results_<stamp> directory path; returns non-zero and prints
# nothing when no completed run exists.
resolve_orthofinder_run() {
    local root="${1:-${RESULTS_DIR:-results}/orthogroups}"
    [ -d "$root" ] || return 1

    local newest
    newest=$(find "$root" -maxdepth 6 -type f -path "*/Orthogroups/Orthogroups.tsv" \
                  -printf '%T@\t%p\n' 2>/dev/null \
             | sort -k1,1nr -k2,2r \
             | head -1 \
             | cut -f2-)
    [ -n "$newest" ] || return 1

    # <run>/Orthogroups/Orthogroups.tsv -> <run>
    printf '%s\n' "$(dirname "$(dirname "$newest")")"
}

# --- Resolve Orthogroups.tsv within the authoritative run ---
# Arguments: $1 - search root (default: ${RESULTS_DIR}/orthogroups)
resolve_orthogroups_tsv() {
    local run
    run=$(resolve_orthofinder_run "${1:-}") || return 1
    local tsv="${run}/Orthogroups/Orthogroups.tsv"
    [ -f "$tsv" ] || return 1
    printf '%s\n' "$tsv"
}

# --- Resolve Orthogroup_Sequences/ within the authoritative run ---
# Arguments: $1 - search root (default: ${RESULTS_DIR}/orthogroups)
# This is where OrthoFinder writes the per-OG UNALIGNED protein FASTAs with
# their original headers — the file every stage means when it says "the
# orthogroup FASTA".
resolve_orthogroup_sequences_dir() {
    local run
    run=$(resolve_orthofinder_run "${1:-}") || return 1
    local dir="${run}/Orthogroup_Sequences"
    [ -d "$dir" ] || return 1
    printf '%s\n' "$dir"
}

# --- Resolve one orthogroup's FASTA ---
# Arguments: $1 - orthogroup base name (e.g. OG0000000)
#            $2 - search root (default: ${RESULTS_DIR}/orthogroups)
# Anchored to the authoritative run's Orthogroup_Sequences/ so it can never
# return a same-basename decoy: the MultipleSequenceAlignments/ copy (an
# already-aligned MSA — re-aligning it yields a garbage tree and exits 0), a
# WorkingDirectory/ copy keyed on OrthoFinder's internal integer ids (taxid
# extraction yields nothing, so stage 05's `-gt 1` gate silently skips
# selection analysis), or a zero-byte leftover. Requires a NON-EMPTY file.
resolve_orthogroup_fasta() {
    local base="$1"
    [ -n "$base" ] || return 1
    local dir
    dir=$(resolve_orthogroup_sequences_dir "${2:-}") || return 1
    local fasta="${dir}/${base}.fa"
    [ -s "$fasta" ] || return 1
    printf '%s\n' "$fasta"
}
