#!/bin/bash
# config.sh
# Purpose: Define global variables and paths for the pipeline.
# Usage: Source this file in all pipeline scripts (source config.sh).

# --- Base Directories ---
# BASE_DIR resolves to the directory containing config.sh, regardless of
# how config.sh is sourced. Using BASH_SOURCE[0] (vs $0) makes this work
# whether config.sh is sourced from the project root or from scripts/.
export BASE_DIR="$(realpath "$(dirname "${BASH_SOURCE[0]}")")"
export RESULTS_DIR="${BASE_DIR}/results"
export SCRIPTS_DIR="${BASE_DIR}/scripts"
export REFERENCE_DIR="${BASE_DIR}/references"
export TRANSCRIPTOME_DIR="${BASE_DIR}/transcriptomes"
export GENOME_DIR="${BASE_DIR}/genomes"
export LOGS_DIR="${RESULTS_DIR}/logs"

# --- Input Files ---
# Berghia identification: real NCBI taxid (Berghia stephanieae = 1287507) +
# genus + species, used to build canonical file paths and as the leaf-id
# prefix in headers. Previous "taxid_berghia_berghia" naming repeated genus
# and used the literal string "taxid" instead of a real ID.
export BERGHIA_TAXID="${BERGHIA_TAXID:-1287507}"
export BERGHIA_GENUS="${BERGHIA_GENUS:-berghia}"
export BERGHIA_SPECIES="${BERGHIA_SPECIES:-stephanieae}"
# File-name prefix: e.g. "1287507_berghia_stephanieae"
export BERGHIA_FILE_PREFIX="${BERGHIA_TAXID}_${BERGHIA_GENUS}_${BERGHIA_SPECIES}"
# Header/leaf-id prefix kept identical to file prefix for cross-reference
# (was previously "taxid_berghia"). Existing data files using the legacy name
# are auto-symlinked at the new path by scripts/fetch_berghia_genome.sh.
export BERGHIA_HEADER_PREFIX="${BERGHIA_FILE_PREFIX}"
export TRANSCRIPTOME="${TRANSCRIPTOME_DIR}/${BERGHIA_FILE_PREFIX}.aa"
export CONSERVED_HMM="${BASE_DIR}/custom_hmms/conserved.hmm"  # Optional
export LSE_HMM="${BASE_DIR}/custom_hmms/lse.hmm"              # Optional
export ID_MAP="${RESULTS_DIR}/reference_sequences/id_map.csv"
export EXPRESSION_DATA="${BASE_DIR}/expression_data.csv"
export GENOME="${GENOME_DIR}/${BERGHIA_FILE_PREFIX}.fasta"
# Bead -4nu: Berghia genome with RefSeq annotations (released 2026-01-29).
# NCBI accession: GCF_034508935.2 (UCSD_Bste_1.3 / Goodheart 2024 BMC Biology
# 22:9). Annotation is in the RefSeq GCF_* path, not the GenBank GCA_*.
# 7713 scaffolds, ~1.13 Gb, 24594 genes, 43318 proteins/CDS.
export GENOME_ACCESSION="${GENOME_ACCESSION:-GCF_034508935.2}"
export GENOME_VERSION="${GENOME_VERSION:-UCSD_Bste_1.3}"
export GENOME_GFF="${GENOME_DIR}/${BERGHIA_FILE_PREFIX}.gff3"
export GENOME_PROTEIN="${GENOME_DIR}/${BERGHIA_FILE_PREFIX}.proteins.fa"
export GENOME_CDS="${GENOME_DIR}/${BERGHIA_FILE_PREFIX}.cds.fna"

# --- Tool Paths (update based on Conda environment) ---
export HMMBUILD="hmmbuild"
export HMMSEARCH="hmmsearch"
export HHMAKE="hhmake"
export HHBLITS="hhblits"
export MAFFT="mafft"
# Bead -align: regime-based aligner dispatch. scripts/run_aligner.sh picks:
#   N <  ALIGNER_LINSI_THRESHOLD (200)   -> MAFFT L-INS-i (gold-standard small)
#   N <  ALIGNER_FAMSA_THRESHOLD (1000)  -> MAFFT --auto (MAFFT picks mode)
#   N >= ALIGNER_FAMSA_THRESHOLD         -> FAMSA 2 (best for large)
# Override via ALIGNER_BACKEND={auto,mafft_linsi,mafft_auto,mafft_retree2,famsa}.
# Falls back to MAFFT --retree 2 when preferred backend unavailable.
export FAMSA="${FAMSA:-famsa}"
export ALIGNER_BACKEND="${ALIGNER_BACKEND:-auto}"
export ALIGNER_LINSI_THRESHOLD="${ALIGNER_LINSI_THRESHOLD:-200}"
export ALIGNER_FAMSA_THRESHOLD="${ALIGNER_FAMSA_THRESHOLD:-1000}"
export IQTREE="iqtree3"
export FASTTREE="FastTree"
export TRIMAL="trimal"
export FASTML="fastml"
export MINIMAP2="minimap2"
export SAMTOOLS="samtools"
export MCSCANX="MCScanX"
# Bead -e59: JCVI MCscan (Tang et al. 2024) entry point. Invoked as a Python
# module so that we always use the version installed via pip in the active
# environment (pip install jcvi). Override only if you've put a custom
# wrapper on PATH.
export JCVI_PYTHON="${JCVI_PYTHON:-python3}"
# Bead -iof: TreeShrink (Mai & Mirarab 2018) — single-tree outlier-long-
# branch removal. Distributed via pip as `treeshrink` which installs a
# `run_treeshrink.py` shim on PATH.
export TREESHRINK="${TREESHRINK:-run_treeshrink.py}"
export ALPHAFOLD="run_alphafold.sh"
export FOLDTREE="foldtree"
export TMALIGN="TMalign"
export SEQTK="seqtk"
export ORTHOFINDER="orthofinder"
export PHYLOFORMER="phyloformer"
export PDFLATEX="pdflatex"
export BUSCO="busco"
export DEEPTMHMM="bash ${SCRIPTS_DIR}/run_deeptmhmm.sh"  # Tries: Apptainer SIF, local binary, biolib, Kyte-Doolittle fallback
export DEEPTMHMM_SIF="${BASE_DIR}/containers/deeptmhmm.sif"
export ASTRAL="astral"
export MRBAYES="mb"
export CLIPKIT="clipkit"
# --- Alignment-cleanup stack (replaces HmmCleaner, May 2026) ---
# PREQUAL (Whelan 2018) — bioconda binary; pre-alignment residue mask.
export PREQUAL="${PREQUAL:-prequal}"
# CLOAK (Chatur/Wheeler bioRxiv 2025) — Python single-file consensus mask
# across an alignment ensemble. cloak.py path expected; cloned from
# github.com/phylowheeler/CLOAK by scripts/unity/install_alignment_filters.sh.
export CLOAK="${CLOAK:-${BASE_DIR}/tools/CLOAK/cloak.py}"
# TAPER (Zhang 2021) — Julia script; per-sequence outlier residue masking.
export TAPER="${TAPER:-${BASE_DIR}/tools/TAPER/correction_multi.jl}"
export JULIA="${JULIA:-julia}"
export CDHIT="cd-hit"
export CAFE="cafe5"
export NOTUNG="java -jar ${BASE_DIR}/tools/Notung-2.9.jar"
# Bead -30g: GeneRax (Morel 2020) replaces Notung as default reconciliation;
# auto-falls-back to Notung when the binary is missing.
export GENERAX="${GENERAX:-generax}"
export RECONCILIATION_BACKEND="${RECONCILIATION_BACKEND:-generax}"
# Bead -urk: SELECTION_BACKEND=stack runs GARD → BUSTED-S → BUSTED-MH →
# aBSREL → MEME (modern HyPhy stack); =absrel runs aBSREL only (legacy,
# ~5x faster but no recombination/SRV/multi-hit/site-level inference).
export SELECTION_BACKEND="${SELECTION_BACKEND:-stack}"
# Bead -i61 (May 2026 v2): the alignment-cleanup stack moved from
# HmmCleaner -> PREQUAL + CLOAK + TAPER (HmmCleaner abandoned because
# its Perl XS deps don't build under the conda env's Perl ABI on Unity).
# Each step gated independently so they can be ablated for sensitivity
# analysis. Defaults: all ON.
#   RUN_PREQUAL — pre-alignment residue mask (Whelan 2018).
#   RUN_CLOAK   — alignment-uncertainty mask via consensus across an
#                 ensemble of MAFFT variants + FAMSA (Wheeler 2025).
#   RUN_TAPER   — post-alignment per-sequence residue-outlier mask
#                 (Zhang 2021).
#   RUN_HMMCLEANER — DEPRECATED. Kept for back-compat; if =1 *and* HmmCleaner
#                 is on PATH, will run after CLOAK as a 4th filter pass.
#                 Default 0 going forward.
export RUN_PREQUAL="${RUN_PREQUAL:-1}"
export RUN_CLOAK="${RUN_CLOAK:-1}"
export RUN_TAPER="${RUN_TAPER:-1}"
export RUN_HMMCLEANER="${RUN_HMMCLEANER:-0}"
# Default canonical aligner mode for run_alignment_ensemble.sh
# (auto = regime-based via run_aligner.sh; or one of linsi/ginsi/einsi/famsa)
export ENSEMBLE_CANONICAL_ALIGNER="${ENSEMBLE_CANONICAL_ALIGNER:-auto}"
export RUN_MACSE="${RUN_MACSE:-1}"
export RUN_TREESHRINK="${RUN_TREESHRINK:-1}"
# Bead -j44: ASR backend
export ASR_BACKEND="${ASR_BACKEND:-iqtree}"
# Bead -e59: synteny backend
export SYNTENY_BACKEND="${SYNTENY_BACKEND:-jcvi}"
# Bead -6nh: TM predictor primary
export TM_PREDICTOR_PRIMARY="${TM_PREDICTOR_PRIMARY:-tmbed}"

# --- Local overrides (gitignored) ---
# Allow per-machine config (real SLURM_EMAIL, custom paths, etc.) without
# committing personal info to public git history. Source last so it
# overrides anything defined above.
if [ -f "${BASE_DIR}/config.local.sh" ]; then
    # shellcheck disable=SC1091
    source "${BASE_DIR}/config.local.sh"
fi
export RSCRIPT="Rscript"
export INTERPROSCAN="interproscan.sh"

# --- NOTUNG Parameters ---
export NOTUNG_JAR="${BASE_DIR}/tools/Notung-2.9.jar"
export NOTUNG_THRESHOLD=90  # Bootstrap threshold for rearrangement

# --- Pipeline Parameters ---
# Project taxa (placeholder labels for non-Berghia samples; real Berghia taxid
# is 1287507, defined above as BERGHIA_TAXID). Replace "taxid1"/"taxid2" with
# real taxids if/when adding additional sequenced samples.
export TAXA=("taxid1" "taxid2" "${BERGHIA_FILE_PREFIX}")
# BERGHIA_TAXID is now the real NCBI ID (1287507), set above. The legacy
# "taxid_berghia" string was a project-internal token, not an NCBI id, and is
# replaced by BERGHIA_FILE_PREFIX (e.g. "1287507_berghia_stephanieae").
export CPUS=16
export GPU_ENABLED=false

# --- HMM Build Strategy ---
# per_species     = one HMM per species (fast, good for preliminary analysis)
# per_orthogroup  = DIAMOND + MCL clustering → one HMM per ortholog group
#                   (gold standard, requires >=32GB RAM and DIAMOND/MCL installed)
export HMM_BUILD_STRATEGY="per_orthogroup"
export DIAMOND="diamond"
export MCL="mcl"

# --- SLURM Configuration ---
# Basic job parameters
export DEFAULT_TIME="${DEFAULT_TIME:-24:00:00}"
export DEFAULT_MEM="${DEFAULT_MEM:-64G}"
# Bead -unity: QOS for long-running jobs on Unity. 'long' QOS allows up to
# 14 days on the cpu partition (default 2 days). Always give jobs more
# walltime than estimated to absorb transient slowdowns.
export SLURM_QOS="${SLURM_QOS:-long}"
# Override via config.local.sh (gitignored) or env var to avoid putting
# personal email in public git history.
export SLURM_EMAIL="${SLURM_EMAIL:-your.email@example.com}"

# Cluster-specific options (leave empty to omit from sbatch)
export SLURM_PARTITION=""           # Partition/queue name, e.g., "cpu", "gpu", "highmem"
export SLURM_ACCOUNT=""             # Account/allocation name, e.g., "katzlab", "pi_username"
export SLURM_QOS=""                 # Quality of service, e.g., "normal", "long", "gpu"
export SLURM_CONSTRAINT=""          # Node constraints, e.g., "skylake", "v100"
export SLURM_RESERVATION=""         # Reservation name (if applicable)
export SLURM_EXTRA_ARGS=""          # Any additional sbatch flags, e.g., "--exclusive --nodelist=node01"

# Array job settings
export SLURM_ARRAY_LIMIT=50         # Max concurrent array tasks (e.g., 50 means --array=0-N%50)

# --- Memory Estimation Profiles ---
# These multipliers are used by functions.sh to estimate memory requirements
# for large datasets. Values are bytes per unit (empirically calibrated).
# Override these if you find estimation is too conservative or aggressive.
export MEM_MULT_IQTREE=100       # bytes per site per taxon (ML tree inference)
export MEM_MULT_ORTHOFINDER=50   # bytes per sequence pair (all-vs-all DIAMOND)
export MEM_MULT_MAFFT=20         # bytes per residue pair (multiple alignment)
export MEM_MULT_HYPHY=200        # bytes per site per branch (aBSREL)
export MEM_MULT_FASTTREE=50      # bytes per site per taxon (approximate ML)
export MEM_BASE_GB=4             # base memory overhead in GB
export MEM_MAX_GB=128            # maximum memory to request via SLURM

# --- Tool-Specific Parameters ---
export HMM_EVALUE="1e-5"
export HHBLITS_EVALUE="1e-5"
export MIN_TM_REGIONS=6  # Kept at 6 to capture fragmented transcriptome assemblies
# DeepTMHMM confidence threshold (0-1): Filter predictions below this confidence
# Recommended: 0.5-0.7 (lower = more permissive, higher = more stringent)
export DEEPTMHMM_MIN_CONFIDENCE=0.5
# Drop proteins longer than MAX_AA_LENGTH before feeding them to TMbed.
# Rationale (bead -m1f): ProtT5 has quadratic GPU memory in sequence length;
# transcript-assembly outliers thousands of aa long take 130+ s/seq on a 2080 Ti
# and timed out stage 02 at 99.88% on 2026-05-13. Real chemoreceptor GPCRs are
# 300-500 aa; 1500 is a safe ceiling well above any plausible isoform.
export MAX_AA_LENGTH=1500

# TIAMMAT-revised, mollusc-optimized GPCR HMM library (bead -m1f).
# 17 HMMs: 5x 7tm_N (Pfam Class A-D GPCR families, revised for mollusc
# sensitivity) + 11x 7TM_GPCR_Sr* (invertebrate Sr-style chemoreceptor
# families from C. elegans, revised on molluscs) + ABA_GPCR. Replaces the
# 4-HMM plain Pfam fallback (PF00001/2/3/F) as the broad detection net
# in stage 02's HMM-first identification — better coverage of mollusc
# chemoreceptors AND broad GPCR classes in a single small scan.
# Source: tiammat_mollusca_v2_workspace REVISED_MODELS.
export TIAMMAT_MOLLUSCA_GPCR_HMM="${BASE_DIR}/references/tiammat_mollusca_gpcr.hmm"
export MIN_SEQ_LENGTH=100
export MAX_GAP_PERCENT=30  # 25-30% recommended for robust phylogenetic inference
# Pre-alignment sequence length filtering (removes gene prediction artifacts)
export SEQ_LENGTH_FILTER_MIN=250           # Biological minimum for 6+ TM GPCRs
export SEQ_LENGTH_FILTER_MAX_METHOD="tukey" # Upper bound: Q3 + 1.5*IQR
export SEQ_LENGTH_FILTER_MAX_FLOOR=800      # Floor for computed upper bound
# IQ-TREE model selection: Use standard amino acid models to reduce computation and overfitting
# Options: "TEST" (all 286 models), "MFP" (ModelFinder Plus), or specific models
# Recommended: Standard empirical AA models with gamma rate variation
# Bead -ryr: split IQTREE_MODEL so '-mset' is a separate flag (recent IQ-TREE
# 2.3+ requires this). Old IQTREE_MODEL string kept as a deprecated alias for
# backwards compat — emits a warning if used.
export IQTREE_MODEL_FIND="MFP"
# Bead -5b0: GPCRtm is the empirical Class A GPCR matrix (Rios 2015) — testing
# it is recommended for the global tree. Profile mixtures (LG+C20, LG+C60,
# LG4X) often beat single-matrix on paralog-rich data.
# Add LG+C20+R10 / LG+C60+R10 / LG4X+R10 via -madd in the HPC retest script.
export IQTREE_MODEL_SET="LG,VT,WAG,JTT,Dayhoff,mtREV,cpREV"
# Legacy combined form (deprecated; kept so older invocations don't break)
export IQTREE_MODEL="${IQTREE_MODEL_FIND} -mset ${IQTREE_MODEL_SET}"
# Bead -ryr: deterministic seeds for reproducibility (the 312h tree must be
# re-derivable). Override per-run if needed; log to provenance.
export IQTREE_SEED="${IQTREE_SEED:-12345}"
export FASTTREE_SEED="${FASTTREE_SEED:-12345}"
# Bead -m6k: TBE alongside UFBoot for rogue-taxon-robust support reporting.
# Empty string means do not request TBE; iqtree2 only got --tbe in v2.3+.
export IQTREE_TBE="${IQTREE_TBE:-1}"           # 1 to request TBE; 0 to skip
export IQTREE_BOOTSTRAP=1000
export ORTHOFINDER_INFLATION=1.5
export FOLDTREE_METHOD="upgma"
export USE_MRBAYES=false
export NUM_STRUCTURAL_CANDIDATES=10
export MIN_ASR_DISTANCE=0.5  # Minimum distance for deep node selection in ASR

# --- CD-HIT Clustering Parameters ---
export CDHIT_IDENTITY=0.98       # 98% identity threshold for clustering
export CDHIT_WORDSIZE=5          # Word size for sequence comparison
export CDHIT_MEMORY=16000        # Memory limit in MB

# --- Reference Subsampling for Phylogenetic Trees ---
export MAX_PHYLO_REFS=2000                 # Max reference sequences in phylogenetic tree
export REF_CLUSTER_IDENTITY=0.7            # CD-HIT clustering identity threshold
export REF_LSE_WEIGHT=1.5                  # Multiplier for LSE vs one-to-one orthologs
export REF_TAXONOMY_WEIGHTS="gastropoda:3.0,cephalopoda:1.5,bivalvia:1.5,other_molluscan_classes:1.2,annelida:1.0,platyhelminthes:1.0,other_lophotrochozoan_phyla:1.0"

# --- Phylogenetic Tree ---
# Directory containing the tree file (override for versioned runs, e.g., protein/v2)
export PHYLO_DIR="${RESULTS_DIR}/phylogenies/protein"
# Tree filename (override if different from default IQ-TREE output)
export PHYLO_TREE_FILENAME="all_berghia_refs.treefile"

# --- Synteny Backend (Bead -e59) ---
# `jcvi`     = JCVI MCscan (Tang et al. 2024 iMeta 3:e211) — preferred default.
# `mcscanx`  = legacy MCScanX path (kept as fallback; pre-existing pipeline runs
#              expect this and the AWK-based GFF parsing in 06_synteny_and_mapping.sh).
export SYNTENY_BACKEND="${SYNTENY_BACKEND:-jcvi}"
# Comma-separated list of close-mollusc target species available for JCVI
# synteny. The wrapper expects, for each <species>, files of the form
# ${GENOME_DIR}/<species>.proteins.fa and ${GENOME_DIR}/<species>.gff3.
# Aplysia is the highest-priority reference per Nath et al. 2025; Lottia and
# Biomphalaria are second-tier additions (set up files and add here).
export JCVI_SYNTENY_TARGETS="${JCVI_SYNTENY_TARGETS:-aplysia_californica}"

# --- TreeShrink (Bead -iof) ---
# Run TreeShrink on per-OG IQ-TREE outputs and the global tree before
# downstream selection analysis. Cleaned trees are written next to the
# originals as <prefix>.treefile.shrunk; the original is preserved as
# <prefix>.original.treefile. Set to 0 to skip TreeShrink entirely.
export RUN_TREESHRINK="${RUN_TREESHRINK:-1}"
# False-positive tolerance for the per-clade outlier test. q=0.05 is the
# permissive setting recommended for paralog-rich Class A GPCR families
# (Mai & Mirarab 2018) — keeps legitimate LSE bursts but trims contamination
# / mis-annotation artefacts.
export TREESHRINK_QUANTILE="${TREESHRINK_QUANTILE:-0.05}"

# --- Outgroup for Tree Rooting ---
export OUTGROUP_FASTA="${REFERENCE_DIR}/outgroup.fa"

# --- Orthogroup Confidence Score ---
export OG_CONFIDENCE_WEIGHT=1

# --- InterProScan Classification (production only) ---
export RUN_INTERPROSCAN=false

# --- CAFE5 Parameters ---
export CAFE_LAMBDA_SEARCH=true   # Search for optimal lambda
export CAFE_PVALUE_THRESHOLD=0.05  # P-value threshold for significant expansions

# --- Sensitivity Analysis Parameters ---
export RUN_SENSITIVITY=true              # Run Monte Carlo sensitivity analysis
export SENSITIVITY_PERTURBATION=0.5      # Weight variation range (+/- 50%)
export SENSITIVITY_ITERATIONS=100        # Number of Monte Carlo iterations

# --- Cross-Validation Parameters ---
export RUN_CROSSVAL=false                # Run k-fold cross-validation on reference classification
export CROSSVAL_FOLDS=5                  # Number of CV folds

# --- Expression Analysis Parameters ---
export SALMON_QUANT_DIR="${BASE_DIR}/expression_data"  # Directory with Salmon quant.sf outputs
export MIN_TPM_THRESHOLD=1.0             # Minimum TPM for "expressed"
export TAU_THRESHOLD=0.8                 # Tau index threshold for tissue-specific
export CHEMOSENSORY_TISSUES="rhinophore,oral_veil,tentacle,cephalic"  # Tissues to check for enrichment

# --- Run Mode ---
# Set to "local" to run without SLURM, "slurm" for cluster execution
export RUN_MODE="${RUN_MODE:-slurm}"

# --- Ranking Weights ---
# These weights control how different evidence types contribute to candidate ranking.
# Higher weights give more importance to that criterion.
export PHYLO_WEIGHT=2           # Weight for phylogenetic distance to references
# Bead -ea9: PURIFYING_WEIGHT defaults to 0 because chemoreceptor identification
# rewards diversifying selection on extracellular loops, not whole-gene
# purifying selection. Conserved housekeeping GPCRs (rhodopsin, GABA receptors,
# etc.) should NOT rank highly for the chemoreceptor question. Set this >0
# only when looking for CONSERVED-function GPCRs.
export PURIFYING_WEIGHT=0       # Weight for purifying selection (omega < 1)
export POSITIVE_WEIGHT=1        # Weight for positive selection (omega > 1, primary chemoreceptor signal)
# Bead -e59 / -ar8: with the Berghia genome (GCA_034508935.3) now available,
# real synteny against Aplysia/Lottia/Biomphalaria can run; the previous
# concern about degenerate scores when no genome was available is resolved
# by the log-scale + min_max_anchors=5 normalization in _rank_candidates_lib.
export SYNTENY_WEIGHT=2         # Weight for synteny conservation (intra-genome anchors)
export TANDEM_CLUSTER_WEIGHT=2.5  # Weight for intra-genome tandem clusters (the field's signature chemoreceptor signal)
# Bead -ar8: tandem-cluster sliding-window detection parameters.
export TANDEM_WINDOW_KB="${TANDEM_WINDOW_KB:-100}"   # Max gap between consecutive cluster members (kb)
export TANDEM_MIN_SIZE="${TANDEM_MIN_SIZE:-3}"       # Minimum cluster size to flag
# CSV produced by scripts/compute_tandem_clusters.py
export TANDEM_CLUSTERS_FILE="${TANDEM_CLUSTERS_FILE:-${RESULTS_DIR}/synteny/tandem_clusters.csv}"
# Tandem-cluster plot prefix derived from the CSV path. Single source of truth
# referenced by both run_tandem_detection.sh (producer) and stage 09 (consumer)
# so a TANDEM_CLUSTERS_FILE override doesn't silently break the report-image hook.
export TANDEM_CLUSTERS_PLOT_PREFIX="${TANDEM_CLUSTERS_PLOT_PREFIX:-${TANDEM_CLUSTERS_FILE%.csv}_plot}"
# Bead -qu9: rhinophore RNA-seq has known confounds (low depth + starvation
# artifact, e.g. Galpha_olf absent in seq but expressed by HCR in fed slugs).
# Treat expression as a SOFT signal only; never as a hard gate.
export EXPR_WEIGHT=1            # Weight for expression data (tissue-specific enrichment) — soft signal
export LSE_DEPTH_WEIGHT=1       # Weight for lineage-specific expansion depth

# --- Reference Weighting (for phylogenetic distance) ---
# These weights control how references affect candidate ranking in rank_candidates.py.
# Higher weight = stronger influence on the phylogenetic score component.
#
# If your references are:
#   - Chemoreceptor-focused: Keep CHEMORECEPTOR_REF_WEIGHT > OTHER_GPCR_REF_WEIGHT (default)
#   - General GPCRs: Set both weights equal (e.g., both 1.0) for unbiased ranking
#   - Other focus: Adjust weights according to your biological question
#
export CHEMORECEPTOR_REF_WEIGHT=2.0   # Weight for known chemoreceptor references
export OTHER_GPCR_REF_WEIGHT=1.0      # Weight for other GPCR references

# --- Statistical Thresholds ---
export ABSREL_FDR_THRESHOLD=0.05      # FDR threshold for significant selection
export BOOTSTRAP_THRESHOLD=70         # Minimum bootstrap support for confident nodes
export LSE_DEPTH_PERCENTILE=75        # Percentile for "deep" LSE classification

# --- GPCRdb Fetch Parameters ---
export GPCRDB_SEARCH_TERMS="chemoreceptor,invertebrate"
export GPCRDB_SPECIES="Aplysia,Lottia"
# GPCR families: "all" for all families, or comma-separated: "Class_A,Class_B1,Class_C,Adhesion,Frizzled"
export GPCRDB_FAMILIES="all"

# --- Taxonomic Levels for LSE ---
# Format: "LevelName:taxid1,taxid2,taxid3"
# - LevelName: Human-readable taxonomic level name (NO special characters: avoid & | ; < > $ `)
# - taxids: Comma-separated list of taxonomy IDs belonging to this level
# - Order: Most specific to most general (e.g., Aeolids -> Nudibranchs -> Gastropods)
# - Used by: 03b_lse_classification.sh to classify orthogroups by taxonomic expansion level
# Example with real taxids: "Aeolids:1514845,285658" "Nudibranchs:6524" "Gastropods:6448"
export LSE_LEVELS=("Aeolids:taxid_aeolid1,taxid_aeolid2" "Nudibranchs:taxid_nudi1,taxid_nudi2" "Gastropods:taxid_gastro1,taxid_gastro2")

# --- NCBI Taxonomy IDs for LSE Classification ---
# These are used by lse_refine.py to determine taxonomic levels
export LSE_AEOLID_TAXID=54397      # Aeolidida (aeolid nudibranchs)
export LSE_NUDIBRANCH_TAXID=13843  # Nudibranchia (all nudibranchs)
export LSE_GASTROPOD_TAXID=644     # Gastropoda (all gastropods)

# --- Species Tree ---
export SPECIES_TREE="${RESULTS_DIR}/busco/busco_species_tree.tre"

# --- Local Database Directory (for API-independent mode) ---
# Run setup_databases.py first to populate this directory
# Set to empty string to use online APIs
export LOCAL_DB_DIR="${BASE_DIR}/databases"

# === NEW: Chemosensory Expression (Phase 1) ===
# Tissue-specific expression scoring for olfactory chemoreceptor identification
export NON_CHEMOSENSORY_TISSUES="foot,digestive,gonad,mantle"  # For tau index calculation
export CHEMOSENSORY_EXPR_WEIGHT=3         # Weight for chemosensory-specific expression score
export CHEMOSENSORY_HARD_FILTER=false     # If true, exclude candidates without chemosensory expression
export CHEMOSENSORY_TISSUE_WEIGHTS="rhinophore:2.0,oral-tentacle:1.0"  # Rhinophore prioritized over oral tentacle

# === NEW: G-protein Co-expression (Phase 2) ===
# Identify GPCRs co-expressed with olfactory G-proteins (Golf, Gi, Go)
export GPROTEIN_REF_FASTA="${REFERENCE_DIR}/gprotein_reference.fasta"
export GPROTEIN_REF_CLASSES="${REFERENCE_DIR}/gprotein_reference_classes.tsv"
export GPROTEIN_MOLLUSC_FASTA="${REFERENCE_DIR}/gprotein_mollusc.fasta"      # Optional user-provided
export GPROTEIN_MOLLUSC_CLASSES="${REFERENCE_DIR}/gprotein_mollusc_classes.tsv"
export GPROTEIN_EVALUE="1e-10"
export GPROTEIN_COEXPR_WEIGHT=2           # Weight for G-protein co-expression score
export OLFACTORY_GPROTEIN_CLASSES="Golf,Gs"  # Classes associated with chemosensory signaling

# === NEW: ECL Divergence Analysis (Phase 3) ===
# Analyze extracellular loop divergence vs transmembrane conservation
export ECL_DIVERGENCE_WEIGHT=1.5          # Weight for ECL divergence score
export MIN_ECL_LENGTH=5                   # Minimum ECL length to analyze
export ECL_CONSERVATION_RATIO_THRESHOLD=2.0  # ECL divergence / TM divergence ratio threshold

# === NEW: Binding Pocket Hydrophilicity (Phase 4) ===
# Analyze binding pocket properties for aquatic vs terrestrial ligand detection
export POCKET_ANALYSIS_MODE="aquatic"     # Options: "aquatic", "terrestrial", "both"
export POCKET_HYDROPHILICITY_WEIGHT=1.0   # Weight for pocket hydrophilicity score
export MIN_POCKET_VOLUME=100              # Minimum pocket volume (Å³) for small molecule binding
export MAX_POCKET_VOLUME=1000             # Maximum pocket volume

# === NEW: CAFE Expansion Interpretation (Phase 5) ===
# Interpret CAFE5 gene family expansion with LSE taxonomy
export EXPANSION_WEIGHT=1.5               # Weight for CAFE expansion score
export PREFERRED_EXPANSION_LEVELS="Aeolid-specific"  # Preferred taxonomic expansion levels
