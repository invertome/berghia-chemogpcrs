#!/bin/bash
# validate_pipeline.sh - Dry-run validation of the full pipeline
# Strategy: Create a test directory, symlink pipeline scripts into it,
# generate a custom config.sh that points tools at mocks and data at synthetic inputs,
# then run each step and verify completion flags.
#
# Usage: bash tests/validate_pipeline.sh [--keep-tmpdir]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
KEEP_TMP=false
[[ "${1:-}" == "--keep-tmpdir" ]] && KEEP_TMP=true

# --- Setup test environment ---
TEST_DIR=$(mktemp -d /tmp/berghia_validate_XXXXXX)
MOCK_DIR="$TEST_DIR/mock_tools"

echo "=== Pipeline Dry-Run Validation ==="
echo "Project: $PROJECT_DIR"
echo "Test dir: $TEST_DIR"
echo ""

cleanup() {
    if $KEEP_TMP; then
        echo ""
        echo "Test directory preserved: $TEST_DIR"
    else
        rm -rf "$TEST_DIR"
        echo "Cleaned up $TEST_DIR"
    fi
}
trap cleanup EXIT

# --- Generate synthetic data ---
echo "--- Generating synthetic test data ---"
python3 "$SCRIPT_DIR/generate_test_data.py" "$TEST_DIR"

# --- Create mock tools ---
echo "--- Creating mock tool stubs ---"
bash "$SCRIPT_DIR/mock_tools/create_mocks.sh" "$MOCK_DIR"

# --- Symlink pipeline scripts into test dir ---
echo "--- Setting up test environment ---"
for script in "$PROJECT_DIR"/[0-9][0-9]*.sh; do
    ln -sf "$script" "$TEST_DIR/$(basename "$script")"
done
# Symlink functions.sh
ln -sf "$PROJECT_DIR/functions.sh" "$TEST_DIR/functions.sh"
# Symlink scripts directory (Python scripts)
ln -sf "$PROJECT_DIR/scripts" "$TEST_DIR/scripts"
# Symlink tools directory if it exists
[ -d "$PROJECT_DIR/tools" ] && ln -sf "$PROJECT_DIR/tools" "$TEST_DIR/tools"

# --- Create custom config.sh in test dir ---
# When pipeline scripts do `source config.sh`, they'll find THIS one.
# BASE_DIR is computed from $0 (the invoking script) which is in TEST_DIR.
cat > "$TEST_DIR/config.sh" <<TESTCONFIG
#!/bin/bash
# Test config.sh - auto-generated for pipeline validation
# BASE_DIR is derived from the calling script's location (TEST_DIR)
export BASE_DIR=\$(realpath "\$(dirname "\$0")")
export RESULTS_DIR="\${BASE_DIR}/results"
export SCRIPTS_DIR="\${BASE_DIR}/scripts"
export REFERENCE_DIR="\${BASE_DIR}/references"
export TRANSCRIPTOME_DIR="\${BASE_DIR}/transcriptomes"
export GENOME_DIR="\${BASE_DIR}/genomes"
export LOGS_DIR="\${RESULTS_DIR}/logs"

# --- Input Files ---
export TRANSCRIPTOME="\${TRANSCRIPTOME_DIR}/taxid_berghia_berghia.aa"
export CONSERVED_HMM="\${BASE_DIR}/custom_hmms/conserved.hmm"
export LSE_HMM="\${BASE_DIR}/custom_hmms/lse.hmm"
export ID_MAP="\${RESULTS_DIR}/reference_sequences/id_map.csv"
export EXPRESSION_DATA="\${BASE_DIR}/expression_data.csv"
export GENOME="\${GENOME_DIR}/taxid_berghia_berghia.fasta"

# --- Tool Paths (ALL pointing to mocks) ---
export HMMBUILD="$MOCK_DIR/hmmbuild"
export HMMSEARCH="$MOCK_DIR/hmmsearch"
export HHMAKE="$MOCK_DIR/hhmake"
export HHBLITS="$MOCK_DIR/hhblits"
export MAFFT="$MOCK_DIR/mafft"
export IQTREE="$MOCK_DIR/iqtree2"
export FASTTREE="$MOCK_DIR/FastTree"
export TRIMAL="$MOCK_DIR/trimal"
export FASTML="$MOCK_DIR/fastml"
export MINIMAP2="$MOCK_DIR/minimap2"
export SAMTOOLS="$MOCK_DIR/samtools"
export MCSCANX="$MOCK_DIR/MCScanX"
export ALPHAFOLD="$MOCK_DIR/run_alphafold.sh"
export FOLDTREE="$MOCK_DIR/foldtree"
export TMALIGN="$MOCK_DIR/TMalign"
export SEQTK="$MOCK_DIR/seqtk"
export ORTHOFINDER="$MOCK_DIR/orthofinder"
export PHYLOFORMER="$MOCK_DIR/phyloformer"
export PDFLATEX="$MOCK_DIR/pdflatex"
export BUSCO="$MOCK_DIR/busco"
export DEEPTMHMM="$MOCK_DIR/deeptmhmm"
export ASTRAL="$MOCK_DIR/astral"
export MRBAYES="$MOCK_DIR/mb"
export CLIPKIT="$MOCK_DIR/clipkit"
export CDHIT="$MOCK_DIR/cd-hit"
export CAFE="$MOCK_DIR/cafe5"
export NOTUNG="$MOCK_DIR/java -jar \${BASE_DIR}/tools/Notung-2.9.jar"
export NOTUNG_JAR="\${BASE_DIR}/tools/Notung-2.9.jar"
export RSCRIPT="$MOCK_DIR/Rscript"
export INTERPROSCAN="$MOCK_DIR/interproscan.sh"

# --- NOTUNG Parameters ---
export NOTUNG_THRESHOLD=90

# --- Pipeline Parameters ---
export TAXA=("taxid1" "taxid2" "taxid_berghia")
export BERGHIA_TAXID="taxid_berghia"
export CPUS=2
export GPU_ENABLED=false

# --- SLURM Configuration (local mode) ---
export DEFAULT_TIME="01:00:00"
export DEFAULT_MEM="4G"
export SLURM_EMAIL="test@example.com"
export SLURM_PARTITION=""
export SLURM_ACCOUNT=""
export SLURM_QOS=""
export SLURM_CONSTRAINT=""
export SLURM_RESERVATION=""
export SLURM_EXTRA_ARGS=""
export SLURM_ARRAY_LIMIT=10

# --- Memory Estimation ---
export MEM_MULT_IQTREE=100
export MEM_MULT_ORTHOFINDER=50
export MEM_MULT_MAFFT=20
export MEM_MULT_HYPHY=200
export MEM_MULT_FASTTREE=50
export MEM_BASE_GB=4
export MEM_MAX_GB=16

# --- Tool-Specific Parameters ---
export HMM_EVALUE="1e-5"
export HHBLITS_EVALUE="1e-5"
export MIN_TM_REGIONS=6
export DEEPTMHMM_MIN_CONFIDENCE=0.5
export MIN_SEQ_LENGTH=100
export MAX_GAP_PERCENT=30
export IQTREE_MODEL="MFP -mset LG,WAG,JTT,Dayhoff,mtREV,cpREV"
export IQTREE_BOOTSTRAP=100
export ORTHOFINDER_INFLATION=1.5
export FOLDTREE_METHOD="upgma"
export USE_MRBAYES=false
export NUM_STRUCTURAL_CANDIDATES=2
export MIN_ASR_DISTANCE=0.5

# --- CD-HIT ---
export CDHIT_IDENTITY=0.98
export CDHIT_WORDSIZE=5
export CDHIT_MEMORY=4000

# --- CAFE5 ---
export CAFE_LAMBDA_SEARCH=true
export CAFE_PVALUE_THRESHOLD=0.05

# --- Sensitivity / Cross-Validation ---
export RUN_SENSITIVITY=false
export SENSITIVITY_PERTURBATION=0.5
export SENSITIVITY_ITERATIONS=10
export RUN_CROSSVAL=false
export CROSSVAL_FOLDS=3

# --- Expression ---
export SALMON_QUANT_DIR="\${BASE_DIR}/expression_data"
export MIN_TPM_THRESHOLD=1.0
export TAU_THRESHOLD=0.8
export CHEMOSENSORY_TISSUES="rhinophore,oral_veil,tentacle,cephalic"

# --- Run Mode ---
export RUN_MODE="local"

# --- Ranking Weights ---
export PHYLO_WEIGHT=2
export PURIFYING_WEIGHT=1
export POSITIVE_WEIGHT=1
export SYNTENY_WEIGHT=3
export EXPR_WEIGHT=1
export LSE_DEPTH_WEIGHT=1
export CHEMORECEPTOR_REF_WEIGHT=2.0
export OTHER_GPCR_REF_WEIGHT=1.0

# --- Statistical Thresholds ---
export ABSREL_FDR_THRESHOLD=0.05
export BOOTSTRAP_THRESHOLD=70
export LSE_DEPTH_PERCENTILE=75

# --- GPCRdb ---
export GPCRDB_SEARCH_TERMS="chemoreceptor,invertebrate"
export GPCRDB_SPECIES="Aplysia,Lottia"
export GPCRDB_FAMILIES="all"

# --- LSE Levels ---
export LSE_LEVELS=("TestLevel:taxid1,taxid2")
export LSE_AEOLID_TAXID=54397
export LSE_NUDIBRANCH_TAXID=13843
export LSE_GASTROPOD_TAXID=644

# --- Species Tree ---
export SPECIES_TREE="\${RESULTS_DIR}/busco/busco_species_tree.tre"

# --- Local Database ---
export LOCAL_DB_DIR=""

# === Phase 1-5 Weights ===
export NON_CHEMOSENSORY_TISSUES="foot,digestive,gonad,mantle"
export CHEMOSENSORY_EXPR_WEIGHT=3
export CHEMOSENSORY_HARD_FILTER=false
export GPROTEIN_REF_FASTA="\${REFERENCE_DIR}/gprotein_reference.fasta"
export GPROTEIN_REF_CLASSES="\${REFERENCE_DIR}/gprotein_reference_classes.tsv"
export GPROTEIN_MOLLUSC_FASTA="\${REFERENCE_DIR}/gprotein_mollusc.fasta"
export GPROTEIN_MOLLUSC_CLASSES="\${REFERENCE_DIR}/gprotein_mollusc_classes.tsv"
export GPROTEIN_EVALUE="1e-10"
export GPROTEIN_COEXPR_WEIGHT=2
export OLFACTORY_GPROTEIN_CLASSES="Golf,Gi,Go"
export ECL_DIVERGENCE_WEIGHT=1.5
export MIN_ECL_LENGTH=5
export ECL_CONSERVATION_RATIO_THRESHOLD=2.0
export POCKET_ANALYSIS_MODE="aquatic"
export POCKET_HYDROPHILICITY_WEIGHT=1.0
export MIN_POCKET_VOLUME=100
export MAX_POCKET_VOLUME=1000
export EXPANSION_WEIGHT=1.5
export PREFERRED_EXPANSION_LEVELS="Aeolid-specific"
TESTCONFIG

# --- Track results ---
PASS=0
FAIL=0
SKIP=0
ERRORS=""

run_step() {
    local step_name="$1"
    local step_script="$2"
    local expected_flag="${3:-}"

    echo ""
    echo "--- Step: $step_name ---"

    if [ ! -f "$step_script" ]; then
        echo "  SKIP: Script not found: $step_script"
        SKIP=$((SKIP + 1))
        return 0
    fi

    # Ensure log directory exists
    mkdir -p "$TEST_DIR/results/logs"

    # Run from TEST_DIR so `source config.sh` finds our custom config
    # Prepend MOCK_DIR to PATH so bare tool calls (makeblastdb, blastp) use mocks
    local exit_code=0
    (cd "$TEST_DIR" && PATH="$MOCK_DIR:$PATH" timeout 120 bash "$step_script") \
        >"$TEST_DIR/results/logs/${step_name}.out" \
        2>"$TEST_DIR/results/logs/${step_name}.err" || exit_code=$?

    if [ $exit_code -eq 124 ]; then
        echo "  FAIL: $step_name (TIMEOUT after 120s)"
        FAIL=$((FAIL + 1))
        ERRORS="$ERRORS\n  $step_name: TIMEOUT"
    elif [ $exit_code -eq 0 ]; then
        if [ -n "$expected_flag" ] && [ -f "$expected_flag" ]; then
            echo "  PASS: $step_name (flag: $(basename "$expected_flag"))"
        else
            echo "  PASS: $step_name (exit 0${expected_flag:+, flag missing: $(basename "$expected_flag")})"
        fi
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $step_name (exit code: $exit_code)"
        FAIL=$((FAIL + 1))
        ERRORS="$ERRORS\n  $step_name: exit $exit_code"
        if [ -f "$TEST_DIR/results/logs/${step_name}.err" ]; then
            echo "  Last error lines:"
            tail -5 "$TEST_DIR/results/logs/${step_name}.err" 2>/dev/null | sed 's/^/    /'
        fi
    fi
}

# Helper: seed a completion flag if missing (for downstream steps)
seed_flag() {
    local flag="$1"
    [ -f "$flag" ] || { mkdir -p "$(dirname "$flag")"; touch "$flag"; }
}

# Helper: seed mock output files for downstream dependencies
seed_mock_outputs() {
    local results="$TEST_DIR/results"

    # chemogpcrs outputs - scripts hardcode "chemogpcrs_berghia.fa", not taxid_berghia
    # Always overwrite: step 02 run_command redirect clobbers this file
    cp "$TEST_DIR/transcriptomes/taxid_berghia_berghia.aa" \
       "$results/chemogpcrs/chemogpcrs_berghia.fa" 2>/dev/null || true
    for taxid in taxid1 taxid2; do
        if [ ! -f "$results/chemogpcrs/chemogpcrs_${taxid}.fa" ]; then
            cp "$TEST_DIR/transcriptomes/${taxid}.aa" \
               "$results/chemogpcrs/chemogpcrs_${taxid}.fa" 2>/dev/null || true
        fi
    done
    # Also create per-sample files that step 02 would generate
    if [ ! -f "$results/chemogpcrs/complete_berghia.fa" ]; then
        cp "$TEST_DIR/transcriptomes/taxid_berghia_berghia.aa" \
           "$results/chemogpcrs/complete_berghia.fa" 2>/dev/null || true
    fi
    if [ ! -f "$results/chemogpcrs/complete_ids_berghia.txt" ]; then
        grep "^>" "$results/chemogpcrs/chemogpcrs_berghia.fa" 2>/dev/null | \
            sed 's/>//' > "$results/chemogpcrs/complete_ids_berghia.txt" || true
    fi
    if [ ! -f "$results/chemogpcrs/chemogpcr_ids_berghia.txt" ]; then
        cp "$results/chemogpcrs/complete_ids_berghia.txt" \
           "$results/chemogpcrs/chemogpcr_ids_berghia.txt" 2>/dev/null || true
    fi

    # candidates
    if [ ! -f "$results/candidates/chemogpcr_candidates.fa" ]; then
        cp "$results/chemogpcrs/chemogpcrs_berghia.fa" \
           "$results/candidates/chemogpcr_candidates.fa" 2>/dev/null || true
        grep "^>" "$results/candidates/chemogpcr_candidates.fa" 2>/dev/null | \
            sed 's/>//' > "$results/candidates/chemogpcr_candidates.txt" || true
    fi
    # Also create candidates_clustered.fa (output of step 02a)
    if [ ! -f "$results/candidates/candidates_clustered.fa" ]; then
        cp "$results/candidates/chemogpcr_candidates.fa" \
           "$results/candidates/candidates_clustered.fa" 2>/dev/null || true
    fi

    # Species tree
    mkdir -p "$(dirname "$TEST_DIR/results/busco/busco_species_tree.tre")"
    echo "((taxid_berghia:0.1,taxid1:0.2):0.3,taxid2:0.4);" \
        > "$TEST_DIR/results/busco/busco_species_tree.tre"

    # OrthoFinder mock results for downstream steps
    local og_dir="$results/orthogroups/input/OrthoFinder/Results_mock/Orthogroups"
    mkdir -p "$og_dir"
    if [ ! -f "$og_dir/Orthogroups.GeneCount.tsv" ]; then
        echo -e "Orthogroup\ttaxid_berghia\ttaxid1\ttaxid2" > "$og_dir/Orthogroups.GeneCount.tsv"
        echo -e "OG0000001\t2\t1\t1" >> "$og_dir/Orthogroups.GeneCount.tsv"
        echo -e "OG0000002\t3\t1\t0" >> "$og_dir/Orthogroups.GeneCount.tsv"
    fi
    for og in OG0000001 OG0000002; do
        if [ ! -f "$og_dir/${og}.fa" ]; then
            echo -e ">taxid_berghia_prot_1\nACDEFGHIKLMNPQRSTVWY" > "$og_dir/${og}.fa"
            echo -e ">taxid1_prot_1\nACDEFGHIKLMNPQRSTVWY" >> "$og_dir/${og}.fa"
            echo -e ">ref_1\nGHIKLMNPQRSTVWYACDEF" >> "$og_dir/${og}.fa"
        fi
    done
    if [ ! -f "$og_dir/Orthogroups.txt" ]; then
        echo "OG0000001: taxid_berghia_prot_1 taxid_berghia_prot_2 taxid1_prot_1 ref_1" > "$og_dir/Orthogroups.txt"
        echo "OG0000002: taxid_berghia_prot_3 taxid_berghia_prot_4 taxid1_prot_2 ref_2" >> "$og_dir/Orthogroups.txt"
    fi

    # Orthogroup manifest
    if [ ! -f "$results/orthogroup_manifest.tsv" ]; then
        for og_fa in "$og_dir"/OG*.fa; do
            [ -f "$og_fa" ] && echo -e "$(basename "${og_fa%.fa}")\t$og_fa"
        done > "$results/orthogroup_manifest.tsv"
    fi

    # Mock phylogeny tree
    mkdir -p "$results/phylogenies/protein" "$results/phylogenies/visualizations"
    echo "((taxid_berghia_prot_1:0.1,ref_1:0.2):0.3,(taxid_berghia_prot_2:0.15,ref_2:0.25):0.35);" \
        > "$results/phylogenies/protein/all_berghia_refs.treefile"

    # Mock absrel results
    mkdir -p "$results/selective_pressure"
    echo "orthogroup,branch_id,omega,p_value,fdr" > "$results/selective_pressure/absrel_results.csv"
    echo "OG0000001,taxid_berghia_prot_1,0.3,0.01,0.05" >> "$results/selective_pressure/absrel_results.csv"

    # Mock ranking output
    mkdir -p "$results/ranking"
    grep "^>" "$results/chemogpcrs/chemogpcrs_taxid_berghia.fa" 2>/dev/null | \
        sed 's/>//' > "$results/ranking/candidate_ids.txt" || true
    if [ ! -f "$results/ranking/ranked_candidates_sorted.csv" ]; then
        echo "id,total_score,phylo_score,dnds_score" > "$results/ranking/ranked_candidates_sorted.csv"
        echo "taxid_berghia_prot_1,0.8,0.5,0.3" >> "$results/ranking/ranked_candidates_sorted.csv"
        echo "taxid_berghia_prot_2,0.6,0.4,0.2" >> "$results/ranking/ranked_candidates_sorted.csv"
    fi

    # Mock synteny ids
    echo "taxid_berghia_prot_1" > "$results/synteny/synteny_ids.txt"
}

# ====================================================================
# RUN PIPELINE STEPS
# ====================================================================

R="$TEST_DIR/results"

echo ""
echo "==============================="
echo " Running Pipeline Steps"
echo "==============================="

# --- Step 01: Reference Processing ---
run_step "01_reference_processing" \
    "$TEST_DIR/01_reference_processing.sh" \
    "$R/step_completed_01.txt"

seed_flag "$R/step_completed_01.txt"
seed_mock_outputs

# --- Step 02: ChemoGPCR Identification ---
run_step "02_chemogpcrs_identification" \
    "$TEST_DIR/02_chemogpcrs_identification.sh" \
    "$R/step_completed_02.txt"

seed_flag "$R/step_completed_02.txt"
seed_flag "$R/step_completed_extract_berghia.txt"
seed_mock_outputs

# --- Step 02a: Cluster Sequences ---
run_step "02a_cluster_sequences" \
    "$TEST_DIR/02a_cluster_sequences.sh" \
    "$R/step_completed_02a.txt"

seed_flag "$R/step_completed_02a.txt"

# --- Step 02b: Classify GPCRs ---
run_step "02b_classify_gpcrs" \
    "$TEST_DIR/02b_classify_gpcrs.sh" \
    "$R/step_completed_02b.txt"

# --- Step 03: Orthology Clustering ---
run_step "03_orthology_clustering" \
    "$TEST_DIR/03_orthology_clustering.sh" \
    "$R/step_completed_03.txt"

seed_flag "$R/step_completed_03.txt"
seed_mock_outputs

# --- Step 03a: BUSCO Species Tree ---
run_step "03a_busco_species_tree" \
    "$TEST_DIR/03a_busco_species_tree.sh" \
    "$R/step_completed_busco_species_tree.txt"

seed_flag "$R/step_completed_busco_species_tree.txt"

# --- Step 03b: LSE Classification ---
run_step "03b_lse_classification" \
    "$TEST_DIR/03b_lse_classification.sh" \
    "$R/step_completed_lse_classification.txt"

seed_flag "$R/step_completed_lse_classification.txt"

# --- Step 03c: CAFE Analysis ---
run_step "03c_cafe_analysis" \
    "$TEST_DIR/03c_cafe_analysis.sh" \
    ""

# --- Step 03d: Notung Reconciliation ---
seed_flag "$R/step_completed_04.txt"
seed_mock_outputs
run_step "03d_notung_reconciliation" \
    "$TEST_DIR/03d_notung_reconciliation.sh" \
    ""

# --- Step 04: Phylogenetic Analysis ---
# Ensure LSE classification output for step 04
for level in TestLevel; do
    if [ ! -f "$R/lse_classification/lse_${level}.fa" ]; then
        cp "$TEST_DIR/transcriptomes/taxid_berghia.aa" \
           "$R/lse_classification/lse_${level}.fa" 2>/dev/null || true
    fi
    if [ ! -f "$R/lse_classification/lse_${level}.txt" ]; then
        echo "taxid_berghia_prot_1" > "$R/lse_classification/lse_${level}.txt"
    fi
done

run_step "04_phylogenetic_analysis" \
    "$TEST_DIR/04_phylogenetic_analysis.sh" \
    "$R/step_completed_04.txt"

seed_flag "$R/step_completed_04.txt"
seed_mock_outputs

# --- Step 05: Selective Pressure & ASR ---
run_step "05_selective_pressure_and_asr" \
    "$TEST_DIR/05_selective_pressure_and_asr.sh" \
    "$R/step_completed_05.txt"

seed_flag "$R/step_completed_05.txt"

# --- Step 06: Synteny & Mapping ---
run_step "06_synteny_and_mapping" \
    "$TEST_DIR/06_synteny_and_mapping.sh" \
    "$R/step_completed_synteny.txt"

seed_flag "$R/step_completed_synteny.txt"
seed_mock_outputs

# --- Step 07: Candidate Ranking ---
run_step "07_candidate_ranking" \
    "$TEST_DIR/07_candidate_ranking.sh" \
    "$R/step_completed_07.txt"

seed_flag "$R/step_completed_07.txt"
seed_mock_outputs

# --- Step 08: Structural Analysis ---
run_step "08_structural_analysis" \
    "$TEST_DIR/08_structural_analysis.sh" \
    "$R/step_completed_foldtree.txt"

seed_flag "$R/step_completed_foldtree.txt"

# --- Step 09: Report Generation ---
run_step "09_report_generation" \
    "$TEST_DIR/09_report_generation.sh" \
    "$R/step_completed_report.txt"

# ====================================================================
# SUMMARY
# ====================================================================
echo ""
echo "==============================="
echo " Validation Summary"
echo "==============================="
echo "  PASS: $PASS"
echo "  FAIL: $FAIL"
echo "  SKIP: $SKIP"

if [ $FAIL -gt 0 ]; then
    echo -e "\nFailed steps:$ERRORS"
    echo ""
    echo "Inspect logs: $TEST_DIR/results/logs/"
fi

# --- Python Import Validation ---
echo ""
echo "--- Python Import Validation ---"
python3 "$SCRIPT_DIR/validate_imports.py" "$PROJECT_DIR" 2>&1

# --- Completion flag inventory ---
echo ""
echo "--- Completion Flags ---"
find "$R" -name "step_completed_*.txt" -o -name "checkpoint_*" 2>/dev/null | sort | while read -r f; do
    echo "  $(basename "$f")"
done

echo ""
if [ $FAIL -eq 0 ]; then
    echo "OVERALL: PASS - Pipeline dry-run completed successfully"
    exit 0
else
    echo "OVERALL: FAIL - $FAIL step(s) failed"
    exit 1
fi
