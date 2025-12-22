#!/bin/bash
# validate_config.sh
# Purpose: Validate that config.sh definitions match actual files in the filesystem.
# Checks for taxid mismatches, missing files, and configuration consistency.
# Run this before starting the pipeline to catch configuration errors early.
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage: ./validate_config.sh [--fix] [--verbose]
#   --fix      Attempt to suggest fixes for issues found
#   --verbose  Show detailed information about checks

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse arguments
FIX_MODE=false
VERBOSE=false
for arg in "$@"; do
    case $arg in
        --fix) FIX_MODE=true ;;
        --verbose) VERBOSE=true ;;
        --help|-h)
            echo "Usage: $0 [--fix] [--verbose]"
            echo "  --fix      Show suggested fixes for issues"
            echo "  --verbose  Show detailed check information"
            exit 0
            ;;
    esac
done

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Counters
ERRORS=0
WARNINGS=0

# --- Helper Functions ---

print_header() {
    echo -e "\n${BLUE}=== $1 ===${NC}"
}

print_check() {
    if [ "$VERBOSE" = true ]; then
        echo -e "  Checking: $1"
    fi
}

print_ok() {
    echo -e "  ${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "  ${RED}✗${NC} $1"
    ((ERRORS++))
}

print_warning() {
    echo -e "  ${YELLOW}⚠${NC} $1"
    ((WARNINGS++))
}

print_fix() {
    if [ "$FIX_MODE" = true ]; then
        echo -e "    ${BLUE}Fix:${NC} $1"
    fi
}

# --- Validation Functions ---

validate_directories() {
    print_header "Validating Directories"

    local dirs=(
        "BASE_DIR:${BASE_DIR}"
        "RESULTS_DIR:${RESULTS_DIR}"
        "SCRIPTS_DIR:${SCRIPTS_DIR}"
        "REFERENCE_DIR:${REFERENCE_DIR}"
        "TRANSCRIPTOME_DIR:${TRANSCRIPTOME_DIR}"
        "GENOME_DIR:${GENOME_DIR}"
    )

    for dir_entry in "${dirs[@]}"; do
        local name="${dir_entry%%:*}"
        local path="${dir_entry##*:}"

        print_check "$name: $path"

        if [ -d "$path" ]; then
            print_ok "$name exists: $path"
        elif [ "$name" = "RESULTS_DIR" ] || [ "$name" = "LOGS_DIR" ]; then
            print_warning "$name does not exist (will be created): $path"
        else
            print_error "$name does not exist: $path"
            print_fix "mkdir -p '$path'"
        fi
    done
}

validate_taxa() {
    print_header "Validating Taxa Configuration"

    echo "  Configured TAXA: ${TAXA[*]}"
    echo "  BERGHIA_TAXID: ${BERGHIA_TAXID}"

    # Check if BERGHIA_TAXID is in TAXA array
    local berghia_found=false
    for taxid in "${TAXA[@]}"; do
        if [ "$taxid" = "$BERGHIA_TAXID" ]; then
            berghia_found=true
            break
        fi
    done

    if [ "$berghia_found" = true ]; then
        print_ok "BERGHIA_TAXID found in TAXA array"
    else
        print_error "BERGHIA_TAXID '${BERGHIA_TAXID}' not found in TAXA array"
        print_fix "Add '${BERGHIA_TAXID}' to TAXA array in config.sh"
    fi

    # Check for transcriptome files for each taxid
    for taxid in "${TAXA[@]}"; do
        print_check "Transcriptome for $taxid"

        local found=false
        for ext in ".aa" ".fasta" ".fa" ".faa"; do
            local file="${TRANSCRIPTOME_DIR}/${taxid}${ext}"
            local file2="${TRANSCRIPTOME_DIR}/taxid_${taxid}${ext}"
            local file3="${TRANSCRIPTOME_DIR}/${taxid}_${taxid}${ext}"

            if [ -f "$file" ] || [ -f "$file2" ] || [ -f "$file3" ]; then
                found=true
                print_ok "Found transcriptome for $taxid"
                break
            fi
        done

        if [ "$found" = false ]; then
            # Try to find any matching file
            local matches=$(find "${TRANSCRIPTOME_DIR}" -maxdepth 1 -name "*${taxid}*" -type f 2>/dev/null | head -1)
            if [ -n "$matches" ]; then
                print_warning "No exact match for $taxid, but found: $matches"
            else
                print_error "No transcriptome file found for taxid: $taxid"
                print_fix "Add transcriptome file to ${TRANSCRIPTOME_DIR}/ with naming pattern: ${taxid}.aa"
            fi
        fi
    done
}

validate_lse_levels() {
    print_header "Validating LSE Level Configuration"

    if [ ${#LSE_LEVELS[@]} -eq 0 ]; then
        print_warning "LSE_LEVELS is empty - LSE classification will be skipped"
        return
    fi

    for level in "${LSE_LEVELS[@]}"; do
        local level_name=$(echo "$level" | cut -d':' -f1)
        local level_taxids=$(echo "$level" | cut -d':' -f2 | tr ',' ' ')

        print_check "LSE Level: $level_name"

        # Check that each taxid in the level is also in TAXA
        for taxid in $level_taxids; do
            local found=false
            for t in "${TAXA[@]}"; do
                if [ "$t" = "$taxid" ]; then
                    found=true
                    break
                fi
            done

            if [ "$found" = true ]; then
                if [ "$VERBOSE" = true ]; then
                    print_ok "Taxid $taxid is in TAXA array"
                fi
            else
                print_error "LSE level '$level_name' references taxid '$taxid' which is not in TAXA array"
                print_fix "Add '$taxid' to TAXA array or remove from LSE_LEVELS"
            fi
        done

        print_ok "LSE level '$level_name' with taxids: $level_taxids"
    done
}

validate_reference_files() {
    print_header "Validating Reference Files"

    if [ ! -d "${REFERENCE_DIR}" ]; then
        print_warning "Reference directory does not exist: ${REFERENCE_DIR}"
        return
    fi

    local ref_count=$(find "${REFERENCE_DIR}" -name "*.aa" -o -name "*.fa" -o -name "*.fasta" 2>/dev/null | wc -l)

    if [ "$ref_count" -eq 0 ]; then
        print_warning "No reference files found in ${REFERENCE_DIR}"
    else
        print_ok "Found $ref_count reference file(s)"
    fi

    # Check for expected reference patterns
    for taxid in "${TAXA[@]}"; do
        local conserved="${REFERENCE_DIR}/${taxid}_conserved_refs.aa"
        local lse="${REFERENCE_DIR}/${taxid}_lse_refs.aa"

        if [ -f "$conserved" ]; then
            print_ok "Found conserved refs for $taxid"
        elif [ "$VERBOSE" = true ]; then
            print_warning "No conserved refs for $taxid (optional)"
        fi

        if [ -f "$lse" ]; then
            print_ok "Found LSE refs for $taxid"
        elif [ "$VERBOSE" = true ]; then
            print_warning "No LSE refs for $taxid (optional)"
        fi
    done
}

validate_tools() {
    print_header "Validating Tool Paths"

    local tools=(
        "HMMBUILD:${HMMBUILD:-hmmbuild}"
        "HMMSEARCH:${HMMSEARCH:-hmmsearch}"
        "MAFFT:${MAFFT:-mafft}"
        "IQTREE:${IQTREE:-iqtree2}"
        "FASTTREE:${FASTTREE:-FastTree}"
    )

    for tool_entry in "${tools[@]}"; do
        local name="${tool_entry%%:*}"
        local cmd="${tool_entry##*:}"

        print_check "$name: $cmd"

        if command -v "$cmd" &> /dev/null; then
            print_ok "$name found: $(which $cmd)"
        else
            print_warning "$name not found in PATH: $cmd"
        fi
    done
}

validate_input_files() {
    print_header "Validating Input Files"

    # Check main transcriptome
    if [ -n "${TRANSCRIPTOME}" ]; then
        print_check "TRANSCRIPTOME: ${TRANSCRIPTOME}"
        if [ -f "${TRANSCRIPTOME}" ]; then
            local seq_count=$(grep -c "^>" "${TRANSCRIPTOME}" 2>/dev/null || echo 0)
            print_ok "TRANSCRIPTOME exists ($seq_count sequences)"
        else
            print_error "TRANSCRIPTOME not found: ${TRANSCRIPTOME}"
        fi
    fi

    # Check genome (optional)
    if [ -n "${GENOME}" ]; then
        print_check "GENOME: ${GENOME}"
        if [ -f "${GENOME}" ]; then
            print_ok "GENOME exists"
        else
            print_warning "GENOME not found (optional): ${GENOME}"
        fi
    fi

    # Check expression data (optional)
    if [ -n "${EXPRESSION_DATA}" ]; then
        print_check "EXPRESSION_DATA: ${EXPRESSION_DATA}"
        if [ -f "${EXPRESSION_DATA}" ]; then
            print_ok "EXPRESSION_DATA exists"
        else
            print_warning "EXPRESSION_DATA not found (optional): ${EXPRESSION_DATA}"
        fi
    fi

    # Check custom HMMs (optional)
    if [ -n "${CONSERVED_HMM}" ] && [ -f "${CONSERVED_HMM}" ]; then
        print_ok "Custom conserved HMM found"
    fi

    if [ -n "${LSE_HMM}" ] && [ -f "${LSE_HMM}" ]; then
        print_ok "Custom LSE HMM found"
    fi
}

validate_local_databases() {
    print_header "Validating Local Databases"

    if [ -z "${LOCAL_DB_DIR}" ]; then
        print_warning "LOCAL_DB_DIR not set - will use online APIs"
        return
    fi

    print_check "LOCAL_DB_DIR: ${LOCAL_DB_DIR}"

    if [ ! -d "${LOCAL_DB_DIR}" ]; then
        print_warning "LOCAL_DB_DIR does not exist: ${LOCAL_DB_DIR}"
        print_fix "Run: python3 setup_databases.py --db-dir '${LOCAL_DB_DIR}'"
        return
    fi

    # Check for taxonomy database
    if [ -f "${LOCAL_DB_DIR}/taxa.sqlite" ]; then
        print_ok "Taxonomy database found"
    else
        print_warning "Taxonomy database not found"
        print_fix "Run: python3 setup_databases.py --db-dir '${LOCAL_DB_DIR}'"
    fi

    # Check for GPCRdb data
    if [ -d "${LOCAL_DB_DIR}/gpcrdb" ]; then
        local gpcrdb_files=$(find "${LOCAL_DB_DIR}/gpcrdb" -name "*.json" | wc -l)
        if [ "$gpcrdb_files" -gt 0 ]; then
            print_ok "GPCRdb data found ($gpcrdb_files files)"
        else
            print_warning "GPCRdb directory exists but is empty"
        fi
    else
        print_warning "GPCRdb data not found"
    fi
}

check_case_sensitivity() {
    print_header "Checking Case Sensitivity Issues"

    # Check for potential case mismatches in transcriptome directory
    if [ -d "${TRANSCRIPTOME_DIR}" ]; then
        for taxid in "${TAXA[@]}"; do
            # Find files that might match with different cases
            local exact_match=$(find "${TRANSCRIPTOME_DIR}" -maxdepth 1 -name "${taxid}*" -type f 2>/dev/null | wc -l)
            local case_insensitive=$(find "${TRANSCRIPTOME_DIR}" -maxdepth 1 -iname "${taxid}*" -type f 2>/dev/null | wc -l)

            if [ "$exact_match" -eq 0 ] && [ "$case_insensitive" -gt 0 ]; then
                print_error "Case mismatch for taxid '$taxid' - files exist but with different case"
                local found_file=$(find "${TRANSCRIPTOME_DIR}" -maxdepth 1 -iname "${taxid}*" -type f 2>/dev/null | head -1)
                print_fix "Rename file or update config.sh: found '$found_file'"
            fi
        done
    fi
}

validate_parameter_ranges() {
    print_header "Validating Parameter Ranges"

    # Check E-value thresholds
    if [ -n "${HMM_EVALUE}" ]; then
        local evalue_ok=$(echo "${HMM_EVALUE}" | awk '{if ($1 > 0 && $1 <= 1) print "ok"}')
        if [ "$evalue_ok" = "ok" ]; then
            print_ok "HMM_EVALUE is valid: ${HMM_EVALUE}"
        else
            print_warning "HMM_EVALUE may be unusual: ${HMM_EVALUE}"
        fi
    fi

    # Check TM regions threshold
    if [ -n "${MIN_TM_REGIONS}" ]; then
        if [ "${MIN_TM_REGIONS}" -ge 5 ] && [ "${MIN_TM_REGIONS}" -le 8 ]; then
            print_ok "MIN_TM_REGIONS is valid: ${MIN_TM_REGIONS}"
        else
            print_warning "MIN_TM_REGIONS (${MIN_TM_REGIONS}) is outside typical GPCR range (5-8)"
        fi
    fi

    # Check ranking weights sum to something reasonable
    local weight_sum=$((PHYLO_WEIGHT + PURIFYING_WEIGHT + POSITIVE_WEIGHT + SYNTENY_WEIGHT + EXPR_WEIGHT + LSE_DEPTH_WEIGHT))
    if [ "$weight_sum" -gt 0 ]; then
        print_ok "Ranking weights sum to $weight_sum"
    else
        print_error "Ranking weights sum to 0 - no scoring will occur"
    fi

    # Check synteny weight vs genome availability
    if [ "${SYNTENY_WEIGHT:-0}" -gt 0 ]; then
        local genome_count=$(find "${GENOME_DIR}" -name "*.fasta" -o -name "*.fa" 2>/dev/null | wc -l)
        if [ "$genome_count" -lt 2 ]; then
            print_warning "SYNTENY_WEIGHT=${SYNTENY_WEIGHT} but found only ${genome_count} genome(s)"
            print_fix "Synteny analysis requires multiple genomes. Set SYNTENY_WEIGHT=0 if you don't have comparison genomes"
        else
            print_ok "SYNTENY_WEIGHT=${SYNTENY_WEIGHT} with ${genome_count} genomes available"
        fi
    fi

    # Check CPU count
    if [ -n "${CPUS}" ]; then
        local available_cpus=$(nproc 2>/dev/null || echo 1)
        if [ "${CPUS}" -le "$available_cpus" ]; then
            print_ok "CPUS (${CPUS}) <= available CPUs ($available_cpus)"
        else
            print_warning "CPUS (${CPUS}) > available CPUs ($available_cpus)"
        fi
    fi
}

validate_memory_profiles() {
    print_header "Validating Memory Estimation Profiles"

    # Check that memory multipliers are set and reasonable
    local mem_vars=("MEM_MULT_IQTREE" "MEM_MULT_ORTHOFINDER" "MEM_MULT_MAFFT" "MEM_MULT_HYPHY" "MEM_MULT_FASTTREE")

    for var in "${mem_vars[@]}"; do
        local value="${!var}"
        if [ -z "$value" ]; then
            print_warning "$var not set - will use default"
        elif [ "$value" -gt 0 ] 2>/dev/null; then
            if [ "$VERBOSE" = true ]; then
                print_ok "$var = $value"
            fi
        else
            print_warning "$var has invalid value: $value"
        fi
    done

    # Check MEM_BASE_GB
    if [ -n "${MEM_BASE_GB}" ] && [ "${MEM_BASE_GB}" -ge 1 ] 2>/dev/null; then
        print_ok "MEM_BASE_GB = ${MEM_BASE_GB}G"
    else
        print_warning "MEM_BASE_GB not set or invalid - will use default (4G)"
    fi

    # Check MEM_MAX_GB
    if [ -n "${MEM_MAX_GB}" ] && [ "${MEM_MAX_GB}" -ge 8 ] 2>/dev/null; then
        print_ok "MEM_MAX_GB = ${MEM_MAX_GB}G"
    else
        print_warning "MEM_MAX_GB not set or too low - will use default (128G)"
    fi

    # Check if MEM_MAX_GB exceeds available system memory
    local available_mem_gb=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo 0)
    if [ "$available_mem_gb" -gt 0 ] && [ -n "${MEM_MAX_GB}" ]; then
        if [ "${MEM_MAX_GB}" -gt "$available_mem_gb" ]; then
            print_warning "MEM_MAX_GB (${MEM_MAX_GB}G) exceeds available system memory (${available_mem_gb}G)"
        fi
    fi
}

validate_reference_structure() {
    print_header "Validating Reference File Structure"

    # Check for nath_et_al directory structure
    local nath_et_al_dir="${REFERENCE_DIR}/nath_et_al"

    if [ -d "$nath_et_al_dir" ]; then
        print_ok "Detected nath_et_al reference structure"

        # Check for expected subdirectories
        if [ -d "${nath_et_al_dir}/lse" ]; then
            local lse_count=$(find "${nath_et_al_dir}/lse" -name "*.faa" 2>/dev/null | wc -l)
            print_ok "LSE references: $lse_count files"
        else
            print_warning "nath_et_al/lse directory not found"
        fi

        if [ -d "${nath_et_al_dir}/one_to_one_ortholog" ]; then
            local conserved_count=$(find "${nath_et_al_dir}/one_to_one_ortholog" -name "*.faa" 2>/dev/null | wc -l)
            print_ok "Conserved references: $conserved_count files"
        else
            print_warning "nath_et_al/one_to_one_ortholog directory not found"
        fi

        # Check if files have taxid prefixes
        local has_taxid_prefix=$(find "${nath_et_al_dir}" -name "*.faa" 2>/dev/null | head -1 | xargs basename 2>/dev/null | grep -E "^[0-9]+_" || true)
        if [ -n "$has_taxid_prefix" ]; then
            print_ok "Reference files have taxid prefixes"
        else
            print_warning "Reference files may need taxid prefixes"
            print_fix "Run: python3 scripts/rename_references_with_taxid.py ${nath_et_al_dir} --dry-run"
        fi
    else
        # Check for legacy structure
        local legacy_count=$(find "${REFERENCE_DIR}" -maxdepth 1 -name "*_refs.aa" 2>/dev/null | wc -l)
        if [ "$legacy_count" -gt 0 ]; then
            print_ok "Using legacy reference structure ($legacy_count files)"
        else
            print_warning "No reference files found in ${REFERENCE_DIR}"
        fi
    fi
}

validate_api_credentials() {
    print_header "Validating API Credentials"

    # Check for .env file
    local env_file="${SCRIPT_DIR}/.env"

    if [ -f "$env_file" ]; then
        print_ok "Found .env file"

        # Check for NCBI_API_KEY (without revealing the value)
        if grep -q "^NCBI_API_KEY=" "$env_file" && ! grep -q "^NCBI_API_KEY=$" "$env_file"; then
            local key_length=$(grep "^NCBI_API_KEY=" "$env_file" | cut -d'=' -f2 | tr -d '\n' | wc -c)
            if [ "$key_length" -gt 10 ]; then
                print_ok "NCBI_API_KEY is set (${key_length} characters)"
            else
                print_warning "NCBI_API_KEY appears too short"
            fi
        else
            print_warning "NCBI_API_KEY not set in .env file"
            print_fix "Add NCBI_API_KEY=your_key to .env (get key at https://www.ncbi.nlm.nih.gov/account/settings/)"
        fi

        # Check for NCBI_EMAIL
        if grep -q "^NCBI_EMAIL=" "$env_file"; then
            local email=$(grep "^NCBI_EMAIL=" "$env_file" | cut -d'=' -f2)
            if [[ "$email" == *"@"* ]]; then
                print_ok "NCBI_EMAIL is set"
            else
                print_warning "NCBI_EMAIL may be invalid: $email"
            fi
        else
            print_warning "NCBI_EMAIL not set in .env file"
            print_fix "Add NCBI_EMAIL=your.email@example.com to .env"
        fi
    else
        print_warning ".env file not found"
        print_fix "Create .env file with NCBI_API_KEY and NCBI_EMAIL for faster API access"
    fi
}

validate_python_dependencies() {
    print_header "Validating Python Dependencies"

    local required_modules=("Bio" "ete3" "pandas" "numpy")
    local optional_modules=("dotenv" "requests" "matplotlib" "seaborn")

    for module in "${required_modules[@]}"; do
        if python3 -c "import $module" 2>/dev/null; then
            if [ "$VERBOSE" = true ]; then
                print_ok "Python module '$module' is available"
            fi
        else
            print_error "Required Python module '$module' is not installed"
            print_fix "Run: pip install $module"
        fi
    done

    for module in "${optional_modules[@]}"; do
        if python3 -c "import $module" 2>/dev/null; then
            if [ "$VERBOSE" = true ]; then
                print_ok "Python module '$module' is available"
            fi
        else
            print_warning "Optional Python module '$module' is not installed"
        fi
    done

    # Verify critical packages
    if python3 -c "import Bio; import ete3; import pandas" 2>/dev/null; then
        print_ok "Core Python dependencies available (Bio, ete3, pandas)"
    fi
}

validate_lse_taxonomy() {
    print_header "Validating LSE Taxonomy Configuration"

    # Check if LSE taxonomy IDs are valid
    local lse_taxids=("$LSE_AEOLID_TAXID" "$LSE_NUDIBRANCH_TAXID" "$LSE_GASTROPOD_TAXID")
    local lse_names=("Aeolidida" "Nudibranchia" "Gastropoda")

    for i in "${!lse_taxids[@]}"; do
        local taxid="${lse_taxids[$i]}"
        local name="${lse_names[$i]}"

        if [ -n "$taxid" ] && [ "$taxid" -gt 0 ] 2>/dev/null; then
            print_ok "LSE_${name^^}_TAXID = $taxid"
        else
            print_warning "LSE_${name^^}_TAXID not set or invalid"
        fi
    done

    # Check if USE_PYTHON_LSE is set
    if [ "${USE_PYTHON_LSE:-true}" = true ]; then
        print_ok "Python-based LSE classification enabled (recommended)"

        # Check if ete3 taxonomy database is available
        if [ -n "${LOCAL_DB_DIR}" ] && [ -f "${LOCAL_DB_DIR}/taxa.sqlite" ]; then
            print_ok "Local NCBI taxonomy database available for LSE classification"
        else
            print_warning "Local NCBI taxonomy database not found"
            print_fix "Run: python3 setup_databases.py to download taxonomy database"
        fi
    else
        print_ok "Using legacy bash-based LSE classification"
    fi
}

validate_metadata_utilities() {
    print_header "Validating Metadata Utilities"

    # Check for get_metadata.py script
    local metadata_script="${SCRIPTS_DIR}/get_metadata.py"
    if [ -f "$metadata_script" ]; then
        print_ok "Metadata lookup utility found: $metadata_script"

        # Check if it's executable or can be run with python3
        if [ -x "$metadata_script" ] || command -v python3 &>/dev/null; then
            if [ "$VERBOSE" = true ]; then
                print_ok "Metadata utility is executable"
            fi
        else
            print_warning "Metadata utility may not be executable - ensure python3 is available"
        fi
    else
        print_warning "Metadata lookup utility not found: $metadata_script"
        print_fix "The pipeline will fall back to header parsing (less robust)"
    fi

    # Check if ID_MAP exists (only relevant after step 01 has run)
    if [ -n "${ID_MAP}" ]; then
        if [ -f "${ID_MAP}" ]; then
            local entry_count=$(wc -l < "${ID_MAP}" 2>/dev/null || echo 0)
            entry_count=$((entry_count - 1))  # Subtract header row
            print_ok "ID_MAP exists with $entry_count entries: ${ID_MAP}"

            # Validate CSV structure
            local header=$(head -1 "${ID_MAP}" 2>/dev/null)
            if [[ "$header" == *"original_id"* && "$header" == *"short_id"* && "$header" == *"taxid"* ]]; then
                if [ "$VERBOSE" = true ]; then
                    print_ok "ID_MAP has valid CSV structure"
                fi
            else
                print_warning "ID_MAP may have incorrect structure - expected columns: original_id,short_id,taxid,..."
            fi
        else
            print_warning "ID_MAP not found (will be created in step 01): ${ID_MAP}"
        fi
    fi
}

# --- Main ---

echo -e "${BLUE}=======================================${NC}"
echo -e "${BLUE} Berghia ChemoGPCR Config Validator${NC}"
echo -e "${BLUE}=======================================${NC}"
echo ""
echo "Config file: ${SCRIPT_DIR}/config.sh"
echo ""

validate_directories
validate_taxa
validate_lse_levels
validate_lse_taxonomy
validate_reference_files
validate_reference_structure
validate_input_files
validate_local_databases
validate_api_credentials
validate_tools
validate_python_dependencies
check_case_sensitivity
validate_parameter_ranges
validate_memory_profiles
validate_metadata_utilities

# --- Summary ---
print_header "Validation Summary"

if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo -e "${GREEN}All checks passed!${NC}"
    echo "Your configuration is valid. You can run the pipeline."
    exit 0
elif [ $ERRORS -eq 0 ]; then
    echo -e "${YELLOW}$WARNINGS warning(s) found.${NC}"
    echo "The pipeline may still work, but review warnings above."
    exit 0
else
    echo -e "${RED}$ERRORS error(s) and $WARNINGS warning(s) found.${NC}"
    echo "Please fix the errors before running the pipeline."
    if [ "$FIX_MODE" = false ]; then
        echo ""
        echo "Run with --fix to see suggested fixes:"
        echo "  $0 --fix"
    fi
    exit 1
fi
