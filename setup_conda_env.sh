#!/bin/bash
# setup_conda_env.sh
# Purpose: Set up conda environment with all dependencies for the Berghia ChemoGPCR pipeline.
# Usage: ./setup_conda_env.sh [--verify] [--minimal] [--full]
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

ENV_NAME="berghia-gpcr"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Functions ---

print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}\n"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

check_conda() {
    if ! command -v conda &> /dev/null; then
        print_error "Conda not found. Please install Miniconda or Anaconda first."
        echo "  Download from: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    print_success "Conda found: $(conda --version)"
}

verify_installation() {
    print_header "Verifying Installation"

    local errors=0

    # Check if environment exists
    if ! conda env list | grep -q "^${ENV_NAME} "; then
        print_error "Environment '${ENV_NAME}' not found"
        echo "  Run: ./setup_conda_env.sh to create it"
        exit 1
    fi

    # Activate environment
    eval "$(conda shell.bash hook)"
    conda activate "${ENV_NAME}"

    echo "Checking required tools..."

    # Core bioinformatics tools
    local tools=(
        "hmmbuild:HMMER"
        "hmmsearch:HMMER"
        "hhblits:HH-suite"
        "mafft:MAFFT"
        "iqtree2:IQ-TREE"
        "FastTree:FastTree"
        "trimal:TrimAl"
        "clipkit:ClipKit"
        "cd-hit:CD-HIT"
        "orthofinder:OrthoFinder"
        "hyphy:HyPhy"
        "minimap2:minimap2"
        "samtools:samtools"
        "blastp:BLAST+"
        "seqtk:seqtk"
        "busco:BUSCO"
        "Rscript:R"
    )

    for tool_entry in "${tools[@]}"; do
        tool="${tool_entry%%:*}"
        name="${tool_entry##*:}"
        if command -v "$tool" &> /dev/null; then
            print_success "$name ($tool)"
        else
            print_error "$name ($tool) - NOT FOUND"
            ((errors++))
        fi
    done

    echo ""
    echo "Checking Python packages..."

    # Python packages
    local py_packages=(
        "Bio:biopython"
        "ete3:ete3"
        "pandas:pandas"
        "numpy:numpy"
        "matplotlib:matplotlib"
        "seaborn:seaborn"
        "scipy:scipy"
        "requests:requests"
    )

    for pkg_entry in "${py_packages[@]}"; do
        pkg="${pkg_entry%%:*}"
        name="${pkg_entry##*:}"
        if python3 -c "import $pkg" 2>/dev/null; then
            print_success "$name"
        else
            print_error "$name - NOT FOUND"
            ((errors++))
        fi
    done

    echo ""
    echo "Checking R packages..."

    # R packages
    if Rscript -e "library(ape)" 2>/dev/null; then
        print_success "R ape package"
    else
        print_error "R ape package - NOT FOUND"
        ((errors++))
    fi

    echo ""
    echo "Checking optional tools..."

    # Optional tools
    local optional_tools=(
        "cafe5:CAFE5"
        "mb:MrBayes"
        "interproscan.sh:InterProScan"
        "java:Java (for NOTUNG)"
        "pdflatex:pdfLaTeX"
    )

    for tool_entry in "${optional_tools[@]}"; do
        tool="${tool_entry%%:*}"
        name="${tool_entry##*:}"
        if command -v "$tool" &> /dev/null; then
            print_success "$name ($tool) [optional]"
        else
            print_warning "$name ($tool) - not installed [optional]"
        fi
    done

    echo ""
    if [ $errors -eq 0 ]; then
        print_header "All Required Tools Installed Successfully"
        echo "You can now run the pipeline:"
        echo "  conda activate ${ENV_NAME}"
        echo "  ./01_reference_processing.sh"
    else
        print_error "$errors required tool(s) missing. Please install them before running the pipeline."
        exit 1
    fi
}

install_minimal() {
    print_header "Installing Minimal Environment"
    echo "This installs only the core tools needed for basic pipeline functionality."

    # Create environment with core tools
    conda create -n "${ENV_NAME}" -y python=3.10

    eval "$(conda shell.bash hook)"
    conda activate "${ENV_NAME}"

    # Add channels
    conda config --env --add channels conda-forge
    conda config --env --add channels bioconda

    # Install core tools
    echo "Installing core bioinformatics tools..."
    conda install -y \
        hmmer \
        mafft \
        iqtree \
        fasttree \
        trimal \
        cd-hit \
        blast \
        seqtk

    # Install Python packages
    echo "Installing Python packages..."
    pip install biopython ete3 pandas numpy matplotlib seaborn scipy requests

    print_success "Minimal environment installed"
}

install_full() {
    print_header "Installing Full Environment"
    echo "This installs all tools including optional ones."

    # Check if environment.yml exists
    if [ -f "${SCRIPT_DIR}/environment.yml" ]; then
        echo "Using environment.yml..."
        conda env create -f "${SCRIPT_DIR}/environment.yml" --force
    else
        # Create environment from scratch
        conda create -n "${ENV_NAME}" -y python=3.10

        eval "$(conda shell.bash hook)"
        conda activate "${ENV_NAME}"

        # Add channels
        conda config --env --add channels conda-forge
        conda config --env --add channels bioconda
        conda config --env --set channel_priority strict

        # Install all bioinformatics tools
        echo "Installing bioinformatics tools (this may take a while)..."
        conda install -y \
            hmmer=3.3.2 \
            hhsuite=3.3.0 \
            mafft=7.520 \
            iqtree=2.2.2 \
            fasttree=2.1.11 \
            trimal=1.4.1 \
            clipkit=1.3.0 \
            cd-hit=4.8.1 \
            orthofinder=2.5.5 \
            hyphy=2.5.51 \
            minimap2=2.26 \
            samtools=1.17 \
            blast=2.14.0 \
            seqtk=1.4 \
            busco=5.4.7 \
            diamond=2.1.8 \
            r-base=4.3 \
            r-ape=5.7

        # Install Python packages
        echo "Installing Python packages..."
        pip install biopython ete3 pandas numpy matplotlib seaborn scipy requests
    fi

    eval "$(conda shell.bash hook)"
    conda activate "${ENV_NAME}"

    # Install additional optional tools
    echo ""
    echo "Installing optional tools..."

    # Try to install CAFE5 (may not be available on all platforms)
    conda install -y -c bioconda cafe || print_warning "CAFE5 not available via conda"

    # Try to install MrBayes
    conda install -y mrbayes || print_warning "MrBayes not available via conda"

    print_success "Full environment installed"
}

install_external_tools() {
    print_header "Installing External Tools"
    echo "Some tools require manual installation. This section provides guidance."

    local tools_dir="${SCRIPT_DIR}/tools"
    mkdir -p "$tools_dir"

    echo ""
    echo "1. DeepTMHMM (transmembrane prediction)"
    echo "   - Register at: https://services.healthtech.dtu.dk/software.php"
    echo "   - Download and follow installation instructions"
    echo "   - Or use Docker: docker pull dtubioinformatics/deeptmhmm"
    echo ""

    echo "2. NOTUNG (gene tree reconciliation)"
    echo "   - Download from: http://www.cs.cmu.edu/~durand/Notung/"
    echo "   - Place Notung-2.9.jar in: ${tools_dir}/"
    echo ""

    echo "3. MCScanX (synteny detection)"
    echo "   Attempting to install from source..."
    if [ ! -d "${tools_dir}/MCScanX" ]; then
        cd "$tools_dir"
        git clone https://github.com/wyp1125/MCScanX.git 2>/dev/null || true
        if [ -d "MCScanX" ]; then
            cd MCScanX && make 2>/dev/null || print_warning "MCScanX compilation failed"
            print_success "MCScanX installed to ${tools_dir}/MCScanX"
        fi
        cd "$SCRIPT_DIR"
    else
        print_success "MCScanX already installed"
    fi
    echo ""

    echo "4. pal2nal (codon alignment)"
    if [ ! -d "${tools_dir}/pal2nal.v14" ]; then
        echo "   Downloading pal2nal..."
        cd "$tools_dir"
        wget -q http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz 2>/dev/null || true
        if [ -f "pal2nal.v14.tar.gz" ]; then
            tar -xzf pal2nal.v14.tar.gz
            rm pal2nal.v14.tar.gz
            print_success "pal2nal installed to ${tools_dir}/pal2nal.v14"
        else
            print_warning "Could not download pal2nal"
        fi
        cd "$SCRIPT_DIR"
    else
        print_success "pal2nal already installed"
    fi
    echo ""

    echo "5. FastML (ancestral sequence reconstruction)"
    echo "   - Install via conda: conda install -c bioconda fastml"
    echo "   - Or download from: http://fastml.tau.ac.il/"
    echo ""

    echo "6. AlphaFold (structure prediction)"
    echo "   - See: https://github.com/deepmind/alphafold"
    echo "   - Or use ColabFold: pip install colabfold"
    echo ""

    echo "7. InterProScan (domain annotation)"
    echo "   - Download from: https://www.ebi.ac.uk/interpro/download/"
    echo "   - Note: Large download (~15 GB)"
    echo ""

    # Update PATH recommendation
    echo "Add the following to your ~/.bashrc or ~/.zshrc:"
    echo ""
    echo "  export PATH=\"\${PATH}:${tools_dir}/MCScanX\""
    echo "  export PATH=\"\${PATH}:${tools_dir}/pal2nal.v14\""
    echo ""
}

show_usage() {
    echo "Berghia ChemoGPCR Pipeline - Environment Setup"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --verify     Verify installation of all tools"
    echo "  --minimal    Install minimal environment (core tools only)"
    echo "  --full       Install full environment (all tools)"
    echo "  --external   Show instructions for external tools"
    echo "  --help       Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                 # Default: install full environment"
    echo "  $0 --verify        # Check if all tools are installed"
    echo "  $0 --minimal       # Install only core tools"
    echo ""
}

# --- Main ---

case "${1:-}" in
    --verify)
        check_conda
        verify_installation
        ;;
    --minimal)
        check_conda
        install_minimal
        install_external_tools
        echo ""
        print_success "Setup complete. Run '$0 --verify' to check installation."
        ;;
    --full)
        check_conda
        install_full
        install_external_tools
        echo ""
        print_success "Setup complete. Run '$0 --verify' to check installation."
        ;;
    --external)
        install_external_tools
        ;;
    --help|-h)
        show_usage
        ;;
    "")
        # Default: full installation
        check_conda
        install_full
        install_external_tools
        echo ""
        echo "To activate the environment:"
        echo "  conda activate ${ENV_NAME}"
        echo ""
        echo "To verify installation:"
        echo "  $0 --verify"
        ;;
    *)
        print_error "Unknown option: $1"
        show_usage
        exit 1
        ;;
esac
