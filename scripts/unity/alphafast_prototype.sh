#!/bin/bash
# alphafast_prototype.sh — AlphaFast vs vanilla AF3 go/no-go validation (bead q0o.1).
#
# Compares the AlphaFast backend against vanilla AlphaFold 3 on 3 correctness
# sequences (structural agreement) and measures batch throughput on ~50 Berghia
# candidates. Writes a report ending in a GO/NO-GO recommendation (the VERDICT is
# the user's call — this script only gathers the evidence).
#
# PREREQUISITES (submit only after both are COMPLETE):
#   - DB download  (job 61080110) -> $ALPHAFAST_DB_DIR/mmseqs/uniref90_padded.dbtype
#   - container pull (job 61080189) -> $ALPHAFAST_REPO/alphafast.sif
#   - AF3 weights af3.bin.zst already present in $ALPHAFAST_WEIGHTS_DIR (verified)
#
# NOTE: first execution is a shakedown — AF3-convert JSON naming, AlphaFast output
# .cif paths, and foldseek output columns may need a one-line tweak on first run.
#
#SBATCH --job-name=alphafast_proto
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=logs/alphafast_proto-%j.out
#SBATCH --error=logs/alphafast_proto-%j.err
#SBATCH --account=pi_pkatz_umass_edu

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate base                              # foldseek + python3 + curl
module load alphafold3/latest 2>/dev/null || true  # provides `af3` (convert)
set -u

# --- Paths (override via env) ---
WS="/scratch3/workspace/$USER-jorge"
REPO="${REPO:-$WS/chemogpcrs_2026-05}"
P="${PROTO_DIR:-$WS/alphafast_proto}"
ALPHAFAST_REPO="${ALPHAFAST_REPO:-$WS/alphafast}"
ALPHAFAST_DB_DIR="${ALPHAFAST_DB_DIR:-$WS/alphafast_db}"
ALPHAFAST_WEIGHTS_DIR="${ALPHAFAST_WEIGHTS_DIR:-$WS/alphafold3_models}"   # has af3.bin.zst
SIF="$ALPHAFAST_REPO/alphafast.sif"
CANDS="$REPO/results/chemogpcrs/chemogpcrs_berghia.fa"
REPORT_DIR="$REPO/docs/plans"; REPORT="$REPORT_DIR/2026-06-23-alphafast-validation-report.md"
mkdir -p "$P" "$REPORT_DIR" logs

# --- Preflight: prerequisites must be ready ---
preflight_fail=0
[ -f "$ALPHAFAST_DB_DIR/mmseqs/uniref90_padded.dbtype" ] || { echo "PREFLIGHT: protein DB not ready ($ALPHAFAST_DB_DIR/mmseqs)"; preflight_fail=1; }
[ -s "$SIF" ] || { echo "PREFLIGHT: container missing ($SIF)"; preflight_fail=1; }
[ -s "$ALPHAFAST_WEIGHTS_DIR/af3.bin.zst" ] || { echo "PREFLIGHT: weights af3.bin.zst missing ($ALPHAFAST_WEIGHTS_DIR)"; preflight_fail=1; }
[ -s "$P/ref_P02699.fasta" ] || { echo "PREFLIGHT: reference FASTA missing ($P/ref_P02699.fasta)"; preflight_fail=1; }
[ -s "$CANDS" ] || { echo "PREFLIGHT: candidate FASTA missing ($CANDS)"; preflight_fail=1; }
[ "$preflight_fail" = 0 ] || { echo "Preflight failed — not ready to run."; exit 2; }

# --- Phase 0: assemble correctness + batch FASTAs (clean single-token headers) ---
extract_one() {  # $1=id  $2=outfile  $3=clean_header
  awk -v id="$1" -v hdr="$3" '
    $0 ~ "^>"id"([ \t]|$)" {print ">"hdr; p=1; next}
    /^>/{p=0}
    p' "$CANDS" > "$2"
  [ -s "$2" ] || { echo "ERROR: could not extract $1 from $CANDS"; exit 3; }
}
# correctness #1: reference GPCR (well-represented class-A anchor, known structure)
awk '/^>/{print ">ref_rhodopsin"; next}{print}' "$P/ref_P02699.fasta" > "$P/c_ref.fasta"
# correctness #2 & #3: typical + divergent Berghia candidates. Prefer the
# principled picks from alphafast_select_divergent.sh (divsel/PICKS.txt: typical =
# median best-hit %id to curated chemoreceptor refs, divergent = lowest). The
# hardcoded TYPICAL_ID is only a fallback if that job hasn't produced PICKS yet;
# there is NO divergent fallback (the old Noncode-codepot guess was wrong).
TYPICAL_ID="BersteEVm009179t1"; DIVERGENT_ID=""
if [ -s "$P/divsel/PICKS.txt" ]; then
    TYPICAL_ID=$(awk -F= '/^typical_id=/{print $2}'    "$P/divsel/PICKS.txt")
    DIVERGENT_ID=$(awk -F= '/^divergent_id=/{print $2}' "$P/divsel/PICKS.txt")
fi
[ -n "$DIVERGENT_ID" ] || { echo "ERROR: no divergent pick — run alphafast_select_divergent.sh first (need $P/divsel/PICKS.txt)"; exit 3; }
echo "correctness candidates: typical=$TYPICAL_ID divergent=$DIVERGENT_ID"
extract_one "$TYPICAL_ID"   "$P/c_typical.fasta"   "berste_typical"
extract_one "$DIVERGENT_ID" "$P/c_divergent.fasta" "berste_divergent"
# throughput batch: first 50 candidate records
awk '/^>/{n++} n<=50{print} n>50{exit}' "$CANDS" > "$P/batch50.fasta"
echo "Batch sequences: $(grep -c '^>' "$P/batch50.fasta")"

# --- Phase 1: FASTA -> AF3 JSON (shared input dir for AlphaFast) ---
JIN="$P/alphafast_input_json"; rm -rf "$JIN"; mkdir -p "$JIN"
for f in c_ref c_typical c_divergent; do
    af3 convert --fasta "$P/$f.fasta" --output_dir "$JIN"
done
af3 convert --fasta "$P/batch50.fasta" --output_dir "$JIN"
echo "AlphaFast input JSONs: $(find "$JIN" -maxdepth 1 -name '*.json' | wc -l)"

# --- Phase 2: vanilla AF3 arm on the 3 correctness seqs (timed) ---
declare -A VAN_SEC
for f in c_ref c_typical c_divergent; do
    od="$P/vanilla/$f"; mkdir -p "$od"
    t0=$SECONDS
    bash "$REPO/scripts/run_alphafold.sh" --fasta_paths="$P/$f.fasta" --output_dir="$od" || echo "WARN: vanilla $f failed"
    VAN_SEC[$f]=$((SECONDS - t0))
    echo "VANILLA_${f}_SECONDS=${VAN_SEC[$f]}"
done

# --- Phase 3: AlphaFast batch arm (correctness + 50, single job; timed) ---
AOUT="$P/alphafast_out"; mkdir -p "$AOUT" "$P/tmp" "$P/jax"
t0=$SECONDS
"$ALPHAFAST_REPO/scripts/run_alphafast.sh" \
    --input_dir "$JIN" --output_dir "$AOUT" \
    --db_dir "$ALPHAFAST_DB_DIR" --weights_dir "$ALPHAFAST_WEIGHTS_DIR" \
    --temp_dir "$P/tmp" --container "$SIF" \
    --jax_compilation_cache_dir "$P/jax" --gpu_devices 0 || echo "WARN: alphafast batch returned nonzero"
ALPHAFAST_BATCH_SECONDS=$((SECONDS - t0))
n_af=$(find "$AOUT" -name '*_model.cif' 2>/dev/null | wc -l)
echo "ALPHAFAST_BATCH_SECONDS=$ALPHAFAST_BATCH_SECONDS  structures=$n_af"

# --- Phase 4: compare correctness structures (foldseek TM-score + mean pLDDT) ---
mean_plddt() {  # $1=dir -> mean pLDDT from AF3 confidences json (fallback: cif B-factor)
  python3 - "$1" <<'PY'
import sys, json, glob, os
d = sys.argv[1]
js = [f for f in glob.glob(os.path.join(d, "**", "*confidences*.json"), recursive=True) if "summary" not in f]
if js:
    try:
        o = json.load(open(js[0]))
        a = o.get("atom_plddts") or o.get("plddt") or []
        if a: print(round(sum(a)/len(a), 2)); sys.exit()
    except Exception: pass
cif = glob.glob(os.path.join(d, "**", "*_model.cif"), recursive=True)
if not cif: print("NA"); sys.exit()
vals=[]
for line in open(cif[0]):
    if line.startswith(("ATOM","HETATM")):
        p=line.split()
        # AF3 mmCIF: B-factor (pLDDT) is in the _atom_site B_iso_or_equiv column
        try:
            if p[3]=="CA" or p[5]=="CA": vals.append(float(p[14]))
        except Exception: pass
print(round(sum(vals)/len(vals),2) if vals else "NA")
PY
}
tm_score() {  # $1=cifA $2=cifB -> TM-score via foldseek TMalign mode
  local o="$P/tmaln"; rm -rf "$o" "$P/tmaln_tmp"
  foldseek easy-search "$1" "$2" "$o" "$P/tmaln_tmp" \
    --alignment-type 1 --tmscore-threshold 0.0 -e inf \
    --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt" >/dev/null 2>&1 || true
  awk 'NR==1{print $4}' "$o" 2>/dev/null   # qtmscore
}

{
  echo "# AlphaFast vs vanilla AF3 — validation report (bead q0o.1)"
  echo
  echo "Date: $(date)  |  job: ${SLURM_JOB_ID:-NA}  |  GPU: ${SLURM_JOB_GPUS:-?}"
  echo "DB: $ALPHAFAST_DB_DIR (protein-only)  |  container: $SIF"
  echo
  echo "## Structural agreement (gate: TM-score >= 0.90, esp. the divergent case)"
  echo
  echo "| sequence | vanilla pLDDT | alphafast pLDDT | TM-score (af vs vanilla) | vanilla s |"
  echo "|---|---|---|---|---|"
  for f in c_ref c_typical c_divergent; do
    vdir="$P/vanilla/$f"
    adir=$(find "$AOUT" -maxdepth 1 -type d -iname "*${f#c_}*" | head -1); [ -n "$adir" ] || adir="$AOUT"
    vcif=$(find "$vdir" -name '*_model.cif' | head -1)
    acif=$(find "$adir" -name "*${f#c_}*_model.cif" | head -1); [ -n "$acif" ] || acif=$(find "$AOUT" -name "*${f#c_}*_model.cif" | head -1)
    tm="NA"; [ -n "$vcif" ] && [ -n "$acif" ] && tm=$(tm_score "$acif" "$vcif")
    echo "| $f | $(mean_plddt "$vdir") | $([ -n "$adir" ] && mean_plddt "$adir" || echo NA) | ${tm:-NA} | ${VAN_SEC[$f]:-NA} |"
  done
  echo
  echo "## Throughput"
  echo "- AlphaFast batch (3 correctness + 50): ${ALPHAFAST_BATCH_SECONDS}s for ${n_af} structures (incl. node-local staging)."
  echo "- Vanilla per-seq: ref=${VAN_SEC[c_ref]:-NA}s typical=${VAN_SEC[c_typical]:-NA}s divproxy=${VAN_SEC[c_divergent]:-NA}s."
  echo
  echo "## Caveats"
  echo "- 'divergent' = complete+Code candidate with the LOWEST best-hit %identity to curated"
  echo "  chemoreceptor refs (alphafast_select_divergent.sh) — a principled 'furthest from known'"
  echo "  pick. Once stage-07 ranking exists, optionally re-confirm on the lowest-paralog-identity"
  echo "  / deepest-LSE candidate before final adoption."
  echo "- Reference rhodopsin (P02699) is a well-represented class-A anchor (easy case)."
  echo
  echo "## GO / NO-GO recommendation (verdict = user's call)"
  echo "- Suggested GO if all TM-scores >= 0.90 (esp. divproxy) AND batch throughput beats vanilla."
  echo "- Fill verdict here: ______"
} > "$REPORT"

echo "=== REPORT written: $REPORT ==="
cat "$REPORT"
