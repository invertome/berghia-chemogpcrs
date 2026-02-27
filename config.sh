#!/bin/bash
# config.sh
# Purpose: Define global variables and paths for the pipeline.
# Usage: Source this file in all pipeline scripts (source config.sh).

# --- Base Directories ---
export BASE_DIR=$(realpath "$(dirname "$0")")
export RESULTS_DIR="${BASE_DIR}/results"
export SCRIPTS_DIR="${BASE_DIR}/scripts"
export REFERENCE_DIR="${BASE_DIR}/references"
export TRANSCRIPTOME_DIR="${BASE_DIR}/transcriptomes"
export GENOME_DIR="${BASE_DIR}/genomes"
export LOGS_DIR="${RESULTS_DIR}/logs"

# --- Input Files ---
export TRANSCRIPTOME="${TRANSCRIPTOME_DIR}/taxid_berghia_berghia.aa"
export CONSERVED_HMM="${BASE_DIR}/custom_hmms/conserved.hmm"  # Optional
export LSE_HMM="${BASE_DIR}/custom_hmms/lse.hmm"              # Optional
export ID_MAP="${RESULTS_DIR}/reference_sequences/id_map.csv"
export EXPRESSION_DATA="${BASE_DIR}/expression_data.csv"
export GENOME="${GENOME_DIR}/taxid_berghia_berghia.fasta"

# --- Tool Paths (update based on Conda environment) ---
export HMMBUILD="hmmbuild"
export HMMSEARCH="hmmsearch"
export HHMAKE="hhmake"
export HHBLITS="hhblits"
export MAFFT="mafft"
export IQTREE="iqtree2"
export FASTTREE="FastTree"
export TRIMAL="trimal"
export FASTML="fastml"
export MINIMAP2="minimap2"
export SAMTOOLS="samtools"
export MCSCANX="MCScanX"
export ALPHAFOLD="run_alphafold.sh"
export FOLDTREE="foldtree"
export TMALIGN="TMalign"
export SEQTK="seqtk"
export ORTHOFINDER="orthofinder"
export PHYLOFORMER="phyloformer"
export PDFLATEX="pdflatex"
export BUSCO="busco"
export DEEPTMHMM="bash ${SCRIPTS_DIR}/run_deeptmhmm.sh"  # Tries: local binary, biolib, Kyte-Doolittle fallback
export ASTRAL="astral"
export MRBAYES="mb"
export CLIPKIT="clipkit"
export CDHIT="cd-hit"
export CAFE="cafe5"
export NOTUNG="java -jar ${BASE_DIR}/tools/Notung-2.9.jar"
export RSCRIPT="Rscript"
export INTERPROSCAN="interproscan.sh"

# --- NOTUNG Parameters ---
export NOTUNG_JAR="${BASE_DIR}/tools/Notung-2.9.jar"
export NOTUNG_THRESHOLD=90  # Bootstrap threshold for rearrangement

# --- Pipeline Parameters ---
export TAXA=("taxid1" "taxid2" "taxid_berghia")
export BERGHIA_TAXID="taxid_berghia"
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
export DEFAULT_TIME="24:00:00"
export DEFAULT_MEM="64G"
export SLURM_EMAIL="your.email@example.com"

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
export MIN_SEQ_LENGTH=100
export MAX_GAP_PERCENT=30  # 25-30% recommended for robust phylogenetic inference
# IQ-TREE model selection: Use standard amino acid models to reduce computation and overfitting
# Options: "TEST" (all 286 models), "MFP" (ModelFinder Plus), or specific models
# Recommended: Standard empirical AA models with gamma rate variation
export IQTREE_MODEL="MFP -mset LG,WAG,JTT,Dayhoff,mtREV,cpREV"
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
export PURIFYING_WEIGHT=1       # Weight for purifying selection (omega < 1)
export POSITIVE_WEIGHT=1        # Weight for positive selection (omega > 1)
export SYNTENY_WEIGHT=3         # Weight for synteny conservation
                                # NOTE: Set to 0 if you don't have multiple genomes for synteny analysis.
                                # Synteny analysis requires genome assemblies + annotations for comparison species.
export EXPR_WEIGHT=1            # Weight for expression data (tissue-specific enrichment)
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
export OLFACTORY_GPROTEIN_CLASSES="Golf,Gi,Go"  # Classes associated with olfactory signaling

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
