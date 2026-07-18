#!/bin/bash
# build_molluscan_corpus.sh — Stage B0 of the MolluscaGenes domain-adapted-PLM
# program (roadmap docs/plans/2026-07-18-plm-integration-roadmap.md, bead v4bs.7).
#
# Builds the training corpus for LoRA continued-MLM domain adaptation (DAPT) of an
# ESM-2 protein language model on molluscan proteins, emitting TWO products:
#
#   1. FULL production corpus  (${CORPUS_OUT_DIR}/full_corpus.fasta)
#      — all species INCLUDING Berghia, whole-proteome, redundancy-reduced and
#        per-species-capped. This is the reusable / general model's training data
#        (the cross-mollusc, cross-receptor asset).
#
#   2. FIREWALLED validation corpus  (${CORPUS_OUT_DIR}/firewalled_corpus.fasta)
#      — the full corpus with every sequence homologous to the EVALUATION set
#        removed (>= FIREWALL_MINID identity / >= FIREWALL_MINCOV coverage), so a
#        leave-clade-out / held-out-family generalization measurement is honest
#        (no train/eval leakage). VALIDATION-ONLY: never train the shipped model
#        on this — train on full_corpus, measure generalization with the model
#        trained on firewalled_corpus.
#
# Pipeline: concat (MolluscaGenes + in-house proteomes + Berghia, species-tagged,
# unique IDs)  ->  MMseqs2 easy-linclust dedup (50% id / 80% cov)  ->  per-species
# cap (clade-balance step 1)  ->  full_corpus  ->  MMseqs2 easy-search firewall of
# the eval union {Berghia, candidates, references, gold} vs the corpus  ->  drop
# hits  ->  firewalled_corpus  ->  manifest (counts + per-clade balance + firewall
# audit).
#
# COMPUTE-BEARING (MMseqs2 over tens of millions of sequences) — Unity sbatch ONLY,
# never a login node. This script WRITES nothing outside ${CORPUS_OUT_DIR} and the
# repo logs/ dir.
#
# -----------------------------------------------------------------------------
# TWO OPEN USER DECISIONS (flagged; not silently defaulted — confirm before the
# production run, see the roadmap B0 note "the two firewall knobs are user
# decisions"):
#
#   (i)  CLASS-A-ONLY TAPT ABLATION CORPUS — should we ALSO build a class-A-GPCR-
#        only "targeted DAPT" (TAPT) ablation corpus alongside the whole-proteome
#        corpus, to test whether narrow in-domain adaptation beats broad DAPT?
#        Not built here (whole-proteome DAPT is the settled primary method). Env
#        hook BUILD_CLASSA_ABLATION is reserved (see the guarded stub near the end)
#        but requires a decision on the class-A membership source first.
#
#   (ii) EXACT FIREWALL CUTOFF — FIREWALL_MINID (default 0.30 identity) and
#        FIREWALL_MINCOV (default 0.50 coverage). Stricter (lower id) removes more
#        corpus and is a more conservative firewall but shrinks training data;
#        looser risks eval leakage. 30%/50% is the recipe default; confirm before
#        the production build.
# -----------------------------------------------------------------------------
#
# ENVIRONMENT VARIABLES (override at submit time, e.g.
#   sbatch --export=ALL,MOLLUSCAGENES_AA=/path/mollusca_aa.fa,... build_molluscan_corpus.sh
# — the exact Unity paths are set by the operator; NONE are hardcoded here):
#
#   REPO_ROOT              Repo checkout root (default: $SLURM_SUBMIT_DIR, else $PWD).
#   CORPUS_CONDA_ENV       Conda env providing mmseqs + awk (default: berghia-gpcr).
#   CORPUS_OUT_DIR         REQUIRED. Output dir for both corpora + manifest + tmp.
#
#   MOLLUSCAGENES_AA       REQUIRED. MolluscaGenes mollusca_aa proteins as FASTA
#                          (200+ spp, EvidentialGene). If you only have the BLAST v5
#                          DB, extract first:  blastdbcmd -entry all -db mollusca_aa
#                          -out mollusca_aa.fa  (env ncbi_processing_env). NOTE: its
#                          deflines all start with 'type=protein;' (non-unique, no
#                          per-defline species) — this script rewrites them to unique
#                          '>MolluscaGenes<SP_DELIM>mgN' IDs.
#   INHOUSE_PROTEOMES_DIR  REQUIRED. Directory of the ~557 in-house proteomes
#                          (*.fa|*.faa|*.fasta|*.fna|*.pep, one species per file;
#                          species tag := sanitized basename).
#   BERGHIA_PROTEOME       REQUIRED. Berghia stephanieae whole proteome FASTA.
#   CANDIDATES_FA          REQUIRED. 785 class-A candidates (eval set member).
#   REFERENCES_FA          REQUIRED. 1094 characterized anchors (eval set member).
#   GOLD_SET_FA            OPTIONAL. Extra gold/eval sequences (eval set member).
#
#   MINID                  linclust min-seq-id       (default 0.5).
#   COV                    linclust coverage         (default 0.8).
#   COVMODE                linclust --cov-mode       (default 1 = cov of target).
#   SPLIT_MEM              mmseqs --split-memory-limit (default 300G).
#
#   CAP_PER_SPECIES        Max sequences kept per species after linclust
#                          (clade-balance step 1; 0 disables). Default 5000.
#   CAP_EXEMPT_REGEX       Species tags matching this ERE are NOT capped
#                          (default '^MolluscaGenes$': MolluscaGenes has no per-
#                          defline species, so capping it as one bulk "species"
#                          would wrongly decimate the evidence — it is exempt until
#                          a species map is supplied; see the reweight hook).
#   CLADE_MAP_TSV          OPTIONAL species<TAB>clade map (gastropod/bivalve/
#                          cephalopod/scaphopod/polyplacophoran/...). If given, the
#                          manifest reports per-clade balance. FULL taxonomy-aware
#                          clade REWEIGHTING is a documented follow-on (see the
#                          CLADE REWEIGHT HOOK below) — B0 implements per-species
#                          capping only.
#
#   FIREWALL_MINID         Drop corpus hits with fraction-identity >= this (0.30).
#   FIREWALL_MINCOV        ... AND max(query-cov, target-cov) >= this      (0.50).
#   FIREWALL_EVAL          easy-search E-value ceiling                     (1e-3).
#   FIREWALL_SENS          easy-search sensitivity -s (7.5 = max, for firewall recall).
#   FIREWALL_MAXSEQS       easy-search --max-seqs per query (4000, high recall).
#
#   SP_DELIM               Species/id delimiter embedded in headers (default '___').
#   BERGHIA_SPECIES        Species tag for the Berghia proteome (default Berghia_stephanieae).
#   MOLLUSCA_SPECIES_MAP_TSV  OPTIONAL, RESERVED hook (not wired): a map that would
#                          let MolluscaGenes be split into real species for capping.
#   BUILD_CLASSA_ABLATION  Open decision (i); default 0. Setting 1 currently exits
#                          with a "decision pending" message (see stub).
#
# Submit:
#   sbatch --export=ALL,CORPUS_OUT_DIR=...,MOLLUSCAGENES_AA=...,INHOUSE_PROTEOMES_DIR=...,\
#     BERGHIA_PROTEOME=...,CANDIDATES_FA=...,REFERENCES_FA=... \
#     scripts/unity/build_molluscan_corpus.sh
#
# -----------------------------------------------------------------------------
#SBATCH --job-name=build_molluscan_corpus
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360G
#SBATCH --output=logs/build_molluscan_corpus-%j.out
#SBATCH --error=logs/build_molluscan_corpus-%j.err
#SBATCH --account=pi_pkatz_umass_edu
# NOTE: linclust over ~30-40M sequences is the heavy step; 48h is ~2-3x the
# estimate. For a larger corpus that risks the 48h cpu-partition wall, resubmit
# with:  sbatch --partition=cpu --qos=long --time=<up to 14d> ...

set -eo pipefail
# Activate BEFORE `set -u` (CONDA_BACKUP_CXX unbound-variable trap in the gxx
# deactivate hook fires under nounset otherwise).
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate "${CORPUS_CONDA_ENV:-berghia-gpcr}"
set -u

# ---- inputs / knobs ---------------------------------------------------------
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "${REPO_ROOT}"

MOLLUSCAGENES_AA="${MOLLUSCAGENES_AA:-}"
INHOUSE_PROTEOMES_DIR="${INHOUSE_PROTEOMES_DIR:-}"
BERGHIA_PROTEOME="${BERGHIA_PROTEOME:-}"
CANDIDATES_FA="${CANDIDATES_FA:-}"
REFERENCES_FA="${REFERENCES_FA:-}"
GOLD_SET_FA="${GOLD_SET_FA:-}"
CORPUS_OUT_DIR="${CORPUS_OUT_DIR:-}"

MINID="${MINID:-0.5}"
COV="${COV:-0.8}"
COVMODE="${COVMODE:-1}"
SPLIT_MEM="${SPLIT_MEM:-300G}"

CAP_PER_SPECIES="${CAP_PER_SPECIES:-5000}"
CAP_EXEMPT_REGEX="${CAP_EXEMPT_REGEX:-^MolluscaGenes$}"
CLADE_MAP_TSV="${CLADE_MAP_TSV:-}"

FIREWALL_MINID="${FIREWALL_MINID:-0.30}"
FIREWALL_MINCOV="${FIREWALL_MINCOV:-0.50}"
FIREWALL_EVAL="${FIREWALL_EVAL:-1e-3}"
FIREWALL_SENS="${FIREWALL_SENS:-7.5}"
FIREWALL_MAXSEQS="${FIREWALL_MAXSEQS:-4000}"

SP_DELIM="${SP_DELIM:-___}"
BERGHIA_SPECIES="${BERGHIA_SPECIES:-Berghia_stephanieae}"
MOLLUSCA_SPECIES_MAP_TSV="${MOLLUSCA_SPECIES_MAP_TSV:-}"
BUILD_CLASSA_ABLATION="${BUILD_CLASSA_ABLATION:-0}"

CPUS="${SLURM_CPUS_PER_TASK:-48}"

# ---- validation -------------------------------------------------------------
die() { echo "ERROR: $*" >&2; exit 1; }
[ -n "$CORPUS_OUT_DIR" ]        || die "CORPUS_OUT_DIR is required (refusing to root outputs at /)"
[ -n "$MOLLUSCAGENES_AA" ]      || die "MOLLUSCAGENES_AA is required"
[ -n "$INHOUSE_PROTEOMES_DIR" ] || die "INHOUSE_PROTEOMES_DIR is required"
[ -n "$BERGHIA_PROTEOME" ]      || die "BERGHIA_PROTEOME is required"
[ -n "$CANDIDATES_FA" ]         || die "CANDIDATES_FA is required"
[ -n "$REFERENCES_FA" ]         || die "REFERENCES_FA is required"
[ -s "$MOLLUSCAGENES_AA" ]      || die "MolluscaGenes FASTA missing/empty: $MOLLUSCAGENES_AA"
[ -d "$INHOUSE_PROTEOMES_DIR" ] || die "proteome dir missing: $INHOUSE_PROTEOMES_DIR"
[ -s "$BERGHIA_PROTEOME" ]      || die "Berghia proteome missing/empty: $BERGHIA_PROTEOME"
[ -s "$CANDIDATES_FA" ]         || die "candidates FASTA missing/empty: $CANDIDATES_FA"
[ -s "$REFERENCES_FA" ]         || die "references FASTA missing/empty: $REFERENCES_FA"
[ -z "$GOLD_SET_FA" ] || [ -s "$GOLD_SET_FA" ] || die "GOLD_SET_FA set but missing/empty: $GOLD_SET_FA"
[ -z "$CLADE_MAP_TSV" ] || [ -s "$CLADE_MAP_TSV" ] || die "CLADE_MAP_TSV set but missing/empty: $CLADE_MAP_TSV"
command -v mmseqs >/dev/null 2>&1 || die "mmseqs not on PATH (env ${CORPUS_CONDA_ENV:-berghia-gpcr})"

mkdir -p "$CORPUS_OUT_DIR" logs
# script-local mmseqs scratch under CORPUS_OUT_DIR (never honor $TMP=/tmp; guards
# the cleanup rm to inside CORPUS_OUT_DIR only).
MMTMP="${CORPUS_OUT_DIR}/mmseqs_tmp_${SLURM_JOB_ID:-local}"
mkdir -p "$MMTMP"

RAW="${CORPUS_OUT_DIR}/raw_corpus.fasta"
REPTAG="${CORPUS_OUT_DIR}/corpus_linclust${MINID}"
REP="${REPTAG}_rep_seq.fasta"
FULL="${CORPUS_OUT_DIR}/full_corpus.fasta"
FIRE="${CORPUS_OUT_DIR}/firewalled_corpus.fasta"
EVAL_UNION="${CORPUS_OUT_DIR}/eval_union.fasta"
HITS="${CORPUS_OUT_DIR}/firewall_hits.m8"
LEAKED="${CORPUS_OUT_DIR}/firewall_leaked_target_ids.txt"
SPCOUNTS="${CORPUS_OUT_DIR}/species_counts.tsv"
MANIFEST="${CORPUS_OUT_DIR}/manifest.txt"

log() { echo "[$(date +%T)] $*"; }

# ---- open decision (i): class-A-only TAPT ablation corpus -------------------
if [ "$BUILD_CLASSA_ABLATION" = "1" ]; then
    die "BUILD_CLASSA_ABLATION=1 is an OPEN USER DECISION (see header, decision i).
         A class-A-only TAPT ablation corpus is not built by B0. It needs a decided
         class-A membership source first. Resolve the decision, then wire it here."
fi

log "B0 corpus build starting"
log "  out=$CORPUS_OUT_DIR  cpus=$CPUS  linclust(id=$MINID cov=$COV covmode=$COVMODE)"
log "  firewall(minid=$FIREWALL_MINID mincov=$FIREWALL_MINCOV eval=$FIREWALL_EVAL sens=$FIREWALL_SENS)"

# =============================================================================
# STEP 1 — species-tagged, unique-ID concatenation into $RAW.
# Header rewrite:  >${species}${SP_DELIM}${n} ${original_defline}
# The species tag (everything before the FIRST SP_DELIM) survives linclust so the
# per-species cap can read it off the representative IDs. n is unique within a
# species; species tags are globally unique (basenames unique in a dir; MolluscaGenes
# and Berghia distinct) => IDs are globally unique (mmseqs requires unique IDs).
# =============================================================================
sanitize() { echo "$1" | sed -E 's/[^A-Za-z0-9]+/_/g; s/^_+//; s/_+$//'; }

log "STEP 1: concatenate + species-tag -> $RAW"
: > "$RAW"

# 1a. MolluscaGenes (bulk source tag; deflines carry no per-defline species).
#     If MOLLUSCA_SPECIES_MAP_TSV is ever wired, split this tag into real species.
if [ -n "$MOLLUSCA_SPECIES_MAP_TSV" ]; then
    log "  NOTE: MOLLUSCA_SPECIES_MAP_TSV is a reserved hook and is NOT applied in B0;"
    log "        MolluscaGenes is tagged as one bulk species and cap-exempted."
fi
awk -v sp="MolluscaGenes" -v d="$SP_DELIM" '
    /^>/ { printf(">%s%s%d %s\n", sp, d, ++n, substr($0,2)); next }
    { print }
' "$MOLLUSCAGENES_AA" >> "$RAW"
log "  + MolluscaGenes"

# 1b. In-house proteomes: one species per file, tag := sanitized basename.
shopt -s nullglob
n_prot_files=0
for f in "$INHOUSE_PROTEOMES_DIR"/*.fa "$INHOUSE_PROTEOMES_DIR"/*.faa \
         "$INHOUSE_PROTEOMES_DIR"/*.fasta "$INHOUSE_PROTEOMES_DIR"/*.fna \
         "$INHOUSE_PROTEOMES_DIR"/*.pep; do
    base="$(basename "$f")"
    sp="$(sanitize "${base%.*}")"
    [ -n "$sp" ] || sp="unknown_$(sanitize "$base")"
    # file index in the ID part (after the delimiter) guarantees globally-unique
    # IDs even if two files sanitize to the same species tag; species parse (before
    # the first SP_DELIM) is unaffected.
    awk -v sp="$sp" -v d="$SP_DELIM" -v fi="$n_prot_files" '
        /^>/ { printf(">%s%sf%d_%d %s\n", sp, d, fi, ++n, substr($0,2)); next }
        { print }
    ' "$f" >> "$RAW"
    n_prot_files=$((n_prot_files + 1))
done
shopt -u nullglob
[ "$n_prot_files" -gt 0 ] || die "no proteome files matched in $INHOUSE_PROTEOMES_DIR (*.fa|*.faa|*.fasta|*.fna|*.pep)"
log "  + in-house proteomes: $n_prot_files files"

# 1c. Berghia (included in the FULL/production corpus by design).
awk -v sp="$(sanitize "$BERGHIA_SPECIES")" -v d="$SP_DELIM" '
    /^>/ { printf(">%s%s%d %s\n", sp, d, ++n, substr($0,2)); next }
    { print }
' "$BERGHIA_PROTEOME" >> "$RAW"
log "  + Berghia ($BERGHIA_SPECIES)"

n_raw=$(grep -c '^>' "$RAW" || true)
log "STEP 1 done: raw input sequences = $n_raw"

# =============================================================================
# STEP 2 — MMseqs2 easy-linclust redundancy reduction (keep cluster reps).
# 50% id / 80% cov / cov-mode 1 — the recipe default: collapses isoforms +
# near-identical orthologs while keeping genuine in-clade divergence.
# =============================================================================
log "STEP 2: linclust -> $REP"
mmseqs easy-linclust "$RAW" "$REPTAG" "$MMTMP" \
    --min-seq-id "$MINID" -c "$COV" --cov-mode "$COVMODE" \
    --split-memory-limit "$SPLIT_MEM" --threads "$CPUS"
[ -s "$REP" ] || die "linclust produced no representatives: $REP"
n_rep=$(grep -c '^>' "$REP" || true)
log "STEP 2 done: representatives = $n_rep"

# =============================================================================
# STEP 3 — per-species cap (clade-balance step 1) -> FULL production corpus.
# Keeps the first CAP_PER_SPECIES reps per species tag (deterministic; species
# read as the substring before the first SP_DELIM in each rep ID). Species tags
# matching CAP_EXEMPT_REGEX pass uncapped. Emits per-species seen/kept counts.
#
# --- CLADE REWEIGHT HOOK (follow-on, NOT B0) ---------------------------------
# A full taxonomy-aware reweight (inverse-cluster-size + clade-frequency balance
# across gastropod/bivalve/cephalopod/scaphopod/polyplacophoran) slots in RIGHT
# HERE: read CLADE_MAP_TSV (species->clade), compute a per-clade target budget,
# and sample reps toward it instead of the flat per-species cap. B0 ships the
# flat per-species cap; the manifest already reports per-clade balance so the
# reweight can be calibrated. Keep it a sampling step over reps (never re-cluster).
# =============================================================================
log "STEP 3: per-species cap (cap=$CAP_PER_SPECIES exempt=/$CAP_EXEMPT_REGEX/) -> $FULL"
awk -v cap="$CAP_PER_SPECIES" -v d="$SP_DELIM" -v exempt="$CAP_EXEMPT_REGEX" \
    -v scfile="$SPCOUNTS" '
    function species(h,   id, i) {
        id = substr(h, 2)                    # strip ">"
        sub(/[ \t].*$/, "", id)              # first whitespace token = mmseqs id
        i = index(id, d)
        return (i > 0) ? substr(id, 1, i - 1) : id
    }
    /^>/ {
        sp = species($0)
        seen[sp]++
        if (cap > 0 && exempt != "" && sp ~ exempt)      keep = 1
        else if (cap <= 0)                               keep = 1
        else if (kept[sp] < cap)                         keep = 1
        else                                             keep = 0
        if (keep) { kept[sp]++; print }
        next
    }
    { if (keep) print }
    END {
        print "species\tseen\tkept" > scfile
        for (s in seen) printf("%s\t%d\t%d\n", s, seen[s], kept[s]) >> scfile
    }
' "$REP" > "$FULL"
[ -s "$FULL" ] || die "per-species cap produced empty FULL corpus: $FULL"
n_full=$(grep -c '^>' "$FULL" || true)
n_species=$(( $(wc -l < "$SPCOUNTS") - 1 ))
log "STEP 3 done: FULL production corpus = $n_full seqs across $n_species species"

# =============================================================================
# STEP 4 — homology firewall (VALIDATION corpus only).
# Build the eval union {Berghia, candidates, references, gold?}, MMseqs2 easy-search
# it (max sensitivity, high --max-seqs for firewall recall) against the FULL corpus,
# and DROP any corpus (target) sequence hit by an eval query at
#   fident >= FIREWALL_MINID  AND  max(qcov,tcov) >= FIREWALL_MINCOV.
# max(qcov,tcov) is the conservative rule (removes more => stronger firewall).
# =============================================================================
log "STEP 4: homology firewall"
cat "$BERGHIA_PROTEOME" "$CANDIDATES_FA" "$REFERENCES_FA" > "$EVAL_UNION"
[ -z "$GOLD_SET_FA" ] || cat "$GOLD_SET_FA" >> "$EVAL_UNION"
n_eval=$(grep -c '^>' "$EVAL_UNION" || true)
log "  eval union = $n_eval sequences (Berghia + candidates + references${GOLD_SET_FA:+ + gold})"

# Alternative firewall backend: hmmsearch with per-family HMMs (see roadmap B0).
# easy-search chosen for a single sensitive all-vs-corpus pass with explicit
# identity+coverage columns.
FWTMP="${MMTMP}/firewall"
mkdir -p "$FWTMP"
mmseqs easy-search "$EVAL_UNION" "$FULL" "$HITS" "$FWTMP" \
    --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits" \
    -s "$FIREWALL_SENS" --max-seqs "$FIREWALL_MAXSEQS" -e "$FIREWALL_EVAL" \
    --threads "$CPUS"

# leaked target IDs = corpus seqs to drop.
awk -F'\t' -v minid="$FIREWALL_MINID" -v mincov="$FIREWALL_MINCOV" '
    { cov = ($5 > $6 ? $5 : $6)                 # max(qcov, tcov)
      if ($3 >= minid && cov >= mincov) print $2 }
' "$HITS" | sort -u > "$LEAKED"
n_leaked=$(wc -l < "$LEAKED" || echo 0)
log "  firewall: $n_leaked corpus sequences flagged as eval-homologous -> dropped"

# remove leaked targets from FULL -> FIREWALLED (multi-line-FASTA safe stream).
awk -v leakedfile="$LEAKED" '
    BEGIN { while ((getline l < leakedfile) > 0) drop[l] = 1 }
    /^>/ {
        id = substr($0, 2); sub(/[ \t].*$/, "", id)
        keep = (id in drop) ? 0 : 1
        if (keep) print
        next
    }
    { if (keep) print }
' "$FULL" > "$FIRE"
[ -s "$FIRE" ] || die "firewalled corpus is empty (firewall removed everything?): $FIRE"
n_fire=$(grep -c '^>' "$FIRE" || true)
log "STEP 4 done: FIREWALLED validation corpus = $n_fire seqs"

# firewall audit: max identity between any eval seq and the NEAREST SURVIVING
# corpus seq (target not in the leaked/dropped set). Must be below FIREWALL_MINID;
# a surviving hit CAN still show id >= minid if its coverage was < mincov (allowed
# by the id-AND-cov drop rule) — reported honestly if so.
audit_maxid=$(awk -F'\t' -v leakedfile="$LEAKED" '
    BEGIN { while ((getline l < leakedfile) > 0) drop[l] = 1; m = 0 }
    { if (!($2 in drop) && $3 + 0 > m) m = $3 + 0 }
    END { printf("%.4f", m) }
' "$HITS")
log "STEP 4 audit: max eval->surviving-corpus identity = $audit_maxid (cutoff $FIREWALL_MINID)"
audit_flag="OK (below cutoff)"
awk -v a="$audit_maxid" -v c="$FIREWALL_MINID" 'BEGIN{ exit !(a+0 >= c+0) }' \
    && audit_flag="NOTE: >= cutoff — survived on sub-threshold coverage (id-AND-cov rule); acceptable"

# =============================================================================
# STEP 5 — manifest (counts + per-clade balance + firewall audit).
# =============================================================================
log "STEP 5: manifest -> $MANIFEST"
{
    echo "# MolluscaGenes DAPT corpus — build manifest (B0, v4bs.7)"
    echo "# generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)  job: ${SLURM_JOB_ID:-local}"
    echo
    echo "## inputs"
    echo "MOLLUSCAGENES_AA       $MOLLUSCAGENES_AA"
    echo "INHOUSE_PROTEOMES_DIR  $INHOUSE_PROTEOMES_DIR ($n_prot_files files)"
    echo "BERGHIA_PROTEOME       $BERGHIA_PROTEOME ($BERGHIA_SPECIES)"
    echo "CANDIDATES_FA          $CANDIDATES_FA"
    echo "REFERENCES_FA          $REFERENCES_FA"
    echo "GOLD_SET_FA            ${GOLD_SET_FA:-<none>}"
    echo
    echo "## parameters"
    echo "linclust   min-seq-id=$MINID  cov=$COV  cov-mode=$COVMODE"
    echo "cap        per-species=$CAP_PER_SPECIES  exempt=/$CAP_EXEMPT_REGEX/"
    echo "firewall   minid=$FIREWALL_MINID  mincov=$FIREWALL_MINCOV  eval=$FIREWALL_EVAL  sens=$FIREWALL_SENS  max-seqs=$FIREWALL_MAXSEQS"
    echo
    echo "## counts"
    printf "%-34s %s\n" "raw input sequences"            "$n_raw"
    printf "%-34s %s\n" "linclust representatives"       "$n_rep"
    printf "%-34s %s\n" "FULL production corpus"         "$n_full"
    printf "%-34s %s\n" "  distinct species"             "$n_species"
    printf "%-34s %s\n" "eval union sequences"           "$n_eval"
    printf "%-34s %s\n" "firewall-dropped (leaked)"      "$n_leaked"
    printf "%-34s %s\n" "FIREWALLED validation corpus"   "$n_fire"
    echo
    echo "## firewall audit"
    printf "%-34s %s\n" "max eval->surviving id"         "$audit_maxid"
    printf "%-34s %s\n" "cutoff (FIREWALL_MINID)"        "$FIREWALL_MINID"
    printf "%-34s %s\n" "status"                         "$audit_flag"
    echo
    echo "## per-clade balance (FULL corpus, kept counts)"
    if [ -n "$CLADE_MAP_TSV" ]; then
        # species_counts.tsv (species seen kept) x CLADE_MAP_TSV (species clade)
        awk -F'\t' -v cm="$CLADE_MAP_TSV" '
            BEGIN { while ((getline l < cm) > 0) { split(l, a, "\t"); clade[a[1]] = a[2] } }
            NR == 1 { next }
            { c = (clade[$1] != "") ? clade[$1] : "UNMAPPED"; tot[c] += $3 }
            END { for (c in tot) printf("%-24s %d\n", c, tot[c]) }
        ' "$SPCOUNTS" | sort
    else
        echo "CLADE_MAP_TSV not provided — per-species counts in: $SPCOUNTS"
        echo "(top 15 species by kept count)"
        awk -F'\t' 'NR>1' "$SPCOUNTS" | sort -t$'\t' -k3,3nr | head -15 \
            | awk -F'\t' '{ printf("%-30s seen=%-8s kept=%s\n", $1, $2, $3) }'
    fi
    echo
    echo "## products"
    echo "FULL production corpus  : $FULL"
    echo "FIREWALLED validation   : $FIRE"
    echo "per-species counts      : $SPCOUNTS"
    echo "firewall hits (m8)      : $HITS"
    echo "firewall leaked ids     : $LEAKED"
} > "$MANIFEST"

cat "$MANIFEST"

# ---- cleanup mmseqs scratch (never rm outside CORPUS_OUT_DIR) ----------------
case "$MMTMP" in
    "${CORPUS_OUT_DIR}"/*) rm -rf -- "$MMTMP" ;;
    *) echo "WARN: refusing to remove unexpected temp dir '$MMTMP'" >&2 ;;
esac

log "B0 DONE. FULL=$FULL  FIREWALLED=$FIRE  manifest=$MANIFEST"
