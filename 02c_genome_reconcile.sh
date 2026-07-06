#!/bin/bash
# 02c_genome_reconcile.sh
# Purpose: Genome-track chemoreceptor-candidate reconciliation for Berghia.
#   Runs the SAME stage-02 detector (HMM + TMbed) on the Berghia RefSeq proteins
#   to build a genome-derived candidate track, places the transcriptome candidates
#   on the genome via a placement cascade (minimap2+GMAP concordance -> miniprot
#   -> BLASTp reciprocal-best-hit), and calls scripts/reconcile_candidates.py to
#   merge both tracks into one genome-anchored, provenance-labeled candidate set.
# Inputs:
#   $BERGHIA_GENOME_PROTEINS  RefSeq XP_* proteins (genome track)
#   $GENOME, $GENOME_GFF      RefSeq genome FASTA + GFF3 (coords / loci)
#   $BERGHIA_TRANSCRIPTOME_MRNA  EvidentialGene nucleotide models (.mrna, NEW input)
#   $TRANSCRIPTOME            EvidentialGene .aa (authoritative aalen completeness)
#   stage-02 outputs: results/chemogpcrs/chemogpcrs_berghia.fa + deeptmhmm_berghia/prediction
# Outputs:
#   ${RESULTS_DIR}/reconciliation/reconciled_candidates.{tsv,faa} + reconciliation_report.md
# Toggle: RUN_GENOME_TRACK (default 1). 0 = transcriptome-only legacy: the stage-02
#   candidate set is passed through UNCHANGED to reconciled_candidates.faa.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

#SBATCH --job-name=genome_reconcile
#SBATCH --output=${LOGS_DIR}/02c_genome_reconcile_%j.out
#SBATCH --error=${LOGS_DIR}/02c_genome_reconcile_%j.err
#SBATCH --time=${DEFAULT_TIME}
#SBATCH $(scale_resources)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${SLURM_EMAIL}

# The env is activated by the caller (sbatch wrapper / interactive shell), as with
# every other numbered stage; we do not self-activate. NOT using `set -u` here
# because config.sh/functions.sh reference many optionally-set vars.
set -eo pipefail
source config.sh
source functions.sh

RECON_DIR="${RESULTS_DIR}/reconciliation"
WORK="${RECON_DIR}/work"
OUT_FAA="${RECON_DIR}/reconciled_candidates.faa"
TXOME_CANDS_FA="${RESULTS_DIR}/chemogpcrs/chemogpcrs_berghia.fa"
TXOME_PRED="${RESULTS_DIR}/chemogpcrs/deeptmhmm_berghia/prediction"
# Honor the SLURM allocation for tool threads (falls back to CPUS for local
# runs) so a smaller grant isn't over-subscribed.
THREADS="${SLURM_CPUS_PER_TASK:-${CPUS}}"

mkdir -p "${RECON_DIR}" "${LOGS_DIR}" || { log "Error: cannot create ${RECON_DIR}"; exit 1; }

# ============================================================================
# Toggle OFF: transcriptome-only legacy behavior (byte-identical pass-through).
# ============================================================================
if [ "${RUN_GENOME_TRACK:-1}" = "0" ]; then
    log "RUN_GENOME_TRACK=0 — genome track disabled; passing the transcriptome candidate set through unchanged."
    check_file "$TXOME_CANDS_FA"
    cp "$TXOME_CANDS_FA" "$OUT_FAA"
    touch "${RESULTS_DIR}/step_completed_02c.txt"
    log "Genome track skipped; ${OUT_FAA} == the transcriptome-only candidate set."
    exit 0
fi

# ============================================================================
# Toggle ON: genome-track reconciliation.
# ============================================================================
log "Starting genome-track candidate reconciliation (RUN_GENOME_TRACK=1)."

# --- (a) Staging fail-fast: required inputs must exist BEFORE any alignment. ---
check_file "$BERGHIA_GENOME_PROTEINS"

# The EvidentialGene nucleotide models (.mrna, type=mRNA, 1:1 with the .aa
# proteins) are a NEW input this feature requires; no prior stage needs them and
# nothing auto-stages them. Fail loud with an explicit staging instruction.
if [ ! -f "$BERGHIA_TRANSCRIPTOME_MRNA" ]; then
    log --level=ERROR "Required nucleotide models not found: ${BERGHIA_TRANSCRIPTOME_MRNA}"
    log --level=ERROR "The genome track needs the EvidentialGene nucleotide models (type=mRNA, paired 1:1"
    log --level=ERROR "with the .aa proteins). Stage or symlink them (lowercase, mirroring how genomes/"
    log --level=ERROR "symlinks the RefSeq files) at:"
    log --level=ERROR "    ${TRANSCRIPTOME_DIR}/${BERGHIA_FILE_PREFIX}.mrna"
    log --level=ERROR "then re-run. (Set RUN_GENOME_TRACK=0 to reproduce the transcriptome-only behavior.)"
    exit 1
fi

# Genome coords, the authoritative aalen source, and the stage-02 dependency.
check_file "$GENOME" "$GENOME_GFF" "$TRANSCRIPTOME" "$TXOME_CANDS_FA" "$TXOME_PRED"

mkdir -p "$WORK"

# --- (b) Detector on the RefSeq proteins (same HMM+TMbed flow as stage 02). ---
# Identical detection on both sources = comparability (design §6). This mirrors
# stage 02's 5 steps: HMM-first filter -> length pre-filter -> TMbed -> >=TM
# filter -> extract.
log "Detecting chemoreceptor candidates on the Berghia RefSeq proteins (HMM + TMbed)."
REFSEQ_GPCRS_FA="${WORK}/gpcrs_refseq.fa"
REFSEQ_TMBED_IN="${WORK}/_tmbed_input_refseq.aa"
REFSEQ_TMBED_OUT="${WORK}/deeptmhmm_refseq"
REFSEQ_CAND_IDS="${WORK}/refseq_cand_ids.txt"
REFSEQ_CAND_FA="${RECON_DIR}/refseq_candidates.faa"
REFSEQ_PRED="${REFSEQ_TMBED_OUT}/prediction"

identify_gpcr_candidates "$BERGHIA_GENOME_PROTEINS" "$REFSEQ_GPCRS_FA"
if [ ! -s "$REFSEQ_GPCRS_FA" ]; then
    log --level=WARN "No GPCR-positive RefSeq proteins detected; the genome track will contribute no candidates."
    : > "$REFSEQ_CAND_FA"
    : > "$REFSEQ_CAND_IDS"
    mkdir -p "$REFSEQ_TMBED_OUT"; : > "$REFSEQ_PRED"
else
    filter_fasta_by_length "$REFSEQ_GPCRS_FA" "$REFSEQ_TMBED_IN"
    run_command "deeptmhmm_refseq" ${DEEPTMHMM} -f "$REFSEQ_TMBED_IN" -o "$REFSEQ_TMBED_OUT"
    awk -F"\\t" -v min_tm="${MIN_TM_REGIONS}" -v min_conf="${DEEPTMHMM_MIN_CONFIDENCE:-0.5}" \
        'NF >= 5 && $5+0 >= min_tm && $3+0 >= min_conf {print $1}' \
        "$REFSEQ_PRED" > "$REFSEQ_CAND_IDS" || { log --level=ERROR "Failed to parse RefSeq TMbed prediction"; exit 1; }
    if [ -s "$REFSEQ_CAND_IDS" ]; then
        run_command "extract_refseq_cands" --stdout="$REFSEQ_CAND_FA" \
            ${SEQTK} subseq "$REFSEQ_GPCRS_FA" "$REFSEQ_CAND_IDS"
    else
        log --level=WARN "No RefSeq chemoreceptor candidates after the >=${MIN_TM_REGIONS} TM filter."
        : > "$REFSEQ_CAND_FA"
    fi
fi
log "RefSeq candidates: $(grep -c '^>' "$REFSEQ_CAND_FA" 2>/dev/null || true)"

# --- Build the transcriptome-candidate module input (--txome-candidates). ---
# query, complete, length, n_tm. completeness+length come from the ORIGINAL
# EvidentialGene .aa `aalen=<len>,<pct>,<class>` attribute (authoritative,
# regardless of any header rewrite); n_tm from the stage-02 DeepTMHMM prediction.
log "Assembling transcriptome-candidate table..."
TXOME_CAND_IDS="${WORK}/txome_cand_ids.txt"
TXOME_CAND_TSV="${WORK}/txome_candidates.tsv"
export TXOME_CANDS_FA TRANSCRIPTOME TXOME_PRED TXOME_CAND_IDS TXOME_CAND_TSV
python3 <<'PYEOF'
import os, re, sys
cand_fa = os.environ["TXOME_CANDS_FA"]
aa_src  = os.environ["TRANSCRIPTOME"]
pred    = os.environ["TXOME_PRED"]
ids_out = os.environ["TXOME_CAND_IDS"]
tsv_out = os.environ["TXOME_CAND_TSV"]

# Candidate ids = first token of each header in the stage-02 candidate FASTA.
cand_ids, seen = [], set()
with open(cand_fa) as fh:
    for line in fh:
        if line.startswith(">"):
            cid = line[1:].split()[0]
            if cid not in seen:
                seen.add(cid); cand_ids.append(cid)

# aalen (completeness + aa length) from the ORIGINAL EvidentialGene .aa headers,
# e.g. `>BersteEVm036785t1 type=protein; aalen=104,97%,partial3; clen=321;`.
aalen = {}
with open(aa_src) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        hid = line[1:].split()[0]
        if hid not in seen:
            continue
        m = re.search(r'aalen=([0-9]+),[^,]*,([A-Za-z0-9]+)', line)
        if m:
            aalen[hid] = (m.group(2).lower() == "complete", int(m.group(1)))

# n_TM from the stage-02 DeepTMHMM prediction (col1 = id, col5 = TM-region count).
ntm = {}
if os.path.exists(pred):
    with open(pred) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                ntm[f[0]] = int(float(f[4]))
            except ValueError:
                continue

missing = 0
with open(ids_out, "w") as idf, open(tsv_out, "w") as tf:
    tf.write("#query\tcomplete\tlength\tn_tm\n")
    for cid in cand_ids:
        idf.write(cid + "\n")
        complete, length = aalen.get(cid, (False, 0))
        if cid not in aalen:
            missing += 1
        tf.write(f"{cid}\t{'complete' if complete else 'partial'}\t{length}\t{ntm.get(cid, 0)}\n")
sys.stderr.write(f"txome candidates: {len(cand_ids)} (aalen missing for {missing})\n")
PYEOF

# --- Build the RefSeq-model + RefSeq-loci module inputs from the RefSeq GFF. ---
# --refseq-models: query, chrom, start, end, strand, complete, length, n_tm
#   (query = versioned XP_ accession; coords from the candidate's gene feature).
# --refseq-loci:   ALL gene intervals (chrom,start,end,strand) from the SAME GFF,
#   so a candidate overlapping a NON-candidate RefSeq gene is still chimeric-detectable.
log "Assembling RefSeq gene-model + loci tables from ${GENOME_GFF}..."
REFSEQ_MODELS_TSV="${WORK}/refseq_models.tsv"
REFSEQ_LOCI_TSV="${WORK}/refseq_loci.tsv"
# The GFF walk lives in scripts/build_refseq_model_table.py (unit-tested): it
# resolves each candidate XP_ accession to its gene coords via the CDS gene=
# attribute or the Parent=rna- -> mRNA -> gene chain, and FAILS LOUD (nonzero,
# aborting the stage under `set -e`) if any candidate can't be resolved so a
# genome candidate is never silently dropped (the module's never-drop contract).
python3 "${SCRIPTS_DIR}/build_refseq_model_table.py" \
    --ids "$REFSEQ_CAND_IDS" \
    --gff "$GENOME_GFF" \
    --proteins "$REFSEQ_CAND_FA" \
    --prediction "$REFSEQ_PRED" \
    --models-out "$REFSEQ_MODELS_TSV" \
    --loci-out "$REFSEQ_LOCI_TSV"

# --- (c) Build the GMAP index once. ---
GMAP_DB_DIR="${WORK}/gmap_index"
GMAP_DB_NAME="berghia_genome"
mkdir -p "$GMAP_DB_DIR"
run_command "gmap_build" ${GMAP_BUILD} -D "$GMAP_DB_DIR" -d "$GMAP_DB_NAME" -t "$THREADS" "$GENOME"

# --- (d) Subset the candidate .mrna and run the nucleotide spliced aligners. ---
# Design decision (confirmed): run minimap2 + gmap + miniprot + RBH on ALL
# candidates and let the module's cascade arbitrate (it prefers concordant
# minimap2+gmap > miniprot > rbh, so extra hits for already-placed transcripts
# are inert). No bash "unplaced" pre-filter — that would duplicate module gating.
CAND_MRNA="${WORK}/txome_candidates.mrna"
run_command "subset_cand_mrna" --stdout="$CAND_MRNA" \
    ${SEQTK} subseq "$BERGHIA_TRANSCRIPTOME_MRNA" "$TXOME_CAND_IDS"

MINIMAP_PAF="${WORK}/minimap2.paf"
GMAP_GFF="${WORK}/gmap.gff3"
: > "$MINIMAP_PAF"; : > "$GMAP_GFF"
if [ -s "$CAND_MRNA" ]; then
    run_command "minimap2_splice" --stdout="$MINIMAP_PAF" \
        ${MINIMAP2} -x splice -t "$THREADS" "$GENOME" "$CAND_MRNA"
    run_command "gmap_align" --stdout="$GMAP_GFF" \
        ${GMAP} -D "$GMAP_DB_DIR" -d "$GMAP_DB_NAME" -f gff3_gene -t "$THREADS" "$CAND_MRNA"
else
    log --level=WARN "No candidate .mrna sequences — skipping minimap2/gmap."
fi

# --- (e) Protein-spliced fallback: miniprot on the candidate proteins. ---
MINIPROT_GFF="${WORK}/miniprot.gff3"
: > "$MINIPROT_GFF"
if [ -s "$TXOME_CANDS_FA" ]; then
    run_command "miniprot_align" --stdout="$MINIPROT_GFF" \
        ${MINIPROT} -t "$THREADS" --gff "$GENOME" "$TXOME_CANDS_FA"
fi

# --- (f) Homology fallback: BLASTp reciprocal-best-hit (both directions). ---
# The module labels every gate-passing --blastp-tab hit 'rbh' WITHOUT re-checking
# reciprocity — establishing it is THIS stage's job (design §3/§6). We BLAST both
# directions and keep only the reciprocal-best forward rows.
RBH_TAB="${WORK}/blastp_rbh.tab"
: > "$RBH_TAB"
if [ -s "$TXOME_CANDS_FA" ] && [ -s "$REFSEQ_CAND_FA" ]; then
    REF_DB="${WORK}/blastdb_refseq"
    TXO_DB="${WORK}/blastdb_txome"
    FWD="${WORK}/blastp_fwd.tab"
    REV="${WORK}/blastp_rev.tab"
    # outfmt 6 + qcovs (13 fields, qcovs last) — exactly what parse_blastp_tab reads.
    BLAST_OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
    run_command "makeblastdb_refseq" ${MAKEBLASTDB} -in "$REFSEQ_CAND_FA" -dbtype prot -out "$REF_DB"
    run_command "makeblastdb_txome"  ${MAKEBLASTDB} -in "$TXOME_CANDS_FA" -dbtype prot -out "$TXO_DB"
    run_command "blastp_fwd" --stdout="$FWD" ${BLASTP} -query "$TXOME_CANDS_FA" -db "$REF_DB" \
        -outfmt "$BLAST_OUTFMT" -evalue "${RECON_BLAST_EVALUE:-1e-5}" -max_target_seqs 5 -num_threads "$THREADS"
    run_command "blastp_rev" --stdout="$REV" ${BLASTP} -query "$REFSEQ_CAND_FA" -db "$TXO_DB" \
        -outfmt "$BLAST_OUTFMT" -evalue "${RECON_BLAST_EVALUE:-1e-5}" -max_target_seqs 5 -num_threads "$THREADS"
    export FWD REV RBH_TAB
    python3 <<'PYEOF'
import os
fwd, rev, out = os.environ["FWD"], os.environ["REV"], os.environ["RBH_TAB"]

def best_by_query(path):
    """query -> (best_bitscore, best_subject) over an outfmt6+qcovs file."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            q, s = f[0], f[1]
            try:
                bit = float(f[11])
            except ValueError:
                continue
            if q not in best or bit > best[q][0]:
                best[q] = (bit, s)
    return best

best_fwd = best_by_query(fwd)   # txome query -> refseq subject
best_rev = best_by_query(rev)   # refseq query -> txome subject
rbh = {(q, s) for q, (_, s) in best_fwd.items()
       if s in best_rev and best_rev[s][1] == q}

# Emit one forward row per reciprocal-best pair (highest-bitscore row).
kept = {}
with open(fwd) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            continue
        pair = (f[0], f[1])
        if pair not in rbh:
            continue
        try:
            bit = float(f[11])
        except ValueError:
            continue
        if pair not in kept or bit > kept[pair][0]:
            kept[pair] = (bit, line.rstrip("\n"))

with open(out, "w") as of:
    for _, row in kept.values():
        of.write(row + "\n")
print(f"RBH pairs kept: {len(kept)}")
PYEOF
else
    log --level=WARN "Empty candidate protein set on one side — skipping BLASTp RBH."
fi

# --- (g) --proteins must cover EVERY representative id (placed/unplaced/genome/txome). ---
RECON_PROTEINS="${WORK}/reconcile_proteins.faa"
cat "$REFSEQ_CAND_FA" "$TXOME_CANDS_FA" > "$RECON_PROTEINS"

# --- (h) Reconcile. Thresholds come from config.sh. GENOME_TRACK_MIN_ID/_COV=95/90
#         (id/cov calibration found no clean separation → kept the working hypothesis).
#         GENOME_TRACK_MIN_MARGIN/_LOW_MARGIN=1/8 CALIBRATED 2026-07-06 (bead x89): the
#         best-vs-second margin is non-discriminating for LSE paralogs, so a low HARD
#         gate + higher ADVISORY flag (keep-and-flag); real gating = %id/%cov/RBH. ---
log "Reconciling genome + transcriptome candidate tracks -> ${RECON_DIR}"
python3 "${SCRIPTS_DIR}/reconcile_candidates.py" \
    --minimap2-paf "$MINIMAP_PAF" \
    --gmap-gff "$GMAP_GFF" \
    --miniprot-gff "$MINIPROT_GFF" \
    --blastp-tab "$RBH_TAB" \
    --txome-candidates "$TXOME_CAND_TSV" \
    --refseq-models "$REFSEQ_MODELS_TSV" \
    --refseq-loci "$REFSEQ_LOCI_TSV" \
    --proteins "$RECON_PROTEINS" \
    --min-id "${GENOME_TRACK_MIN_ID}" \
    --min-cov "${GENOME_TRACK_MIN_COV}" \
    --min-margin "${GENOME_TRACK_MIN_MARGIN}" \
    --low-margin-threshold "${GENOME_TRACK_LOW_MARGIN}" \
    --out-dir "$RECON_DIR"

check_file "$OUT_FAA" "${RECON_DIR}/reconciled_candidates.tsv"
touch "${RESULTS_DIR}/step_completed_02c.txt"
log "Genome-track reconciliation completed -> ${RECON_DIR}/reconciled_candidates.{tsv,faa} + reconciliation_report.md"
