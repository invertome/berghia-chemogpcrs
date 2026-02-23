#!/bin/bash
# create_mocks.sh - Generate mock tool wrappers for dry-run testing
# Each mock prints what it would do and creates minimal expected outputs.
set -euo pipefail

MOCK_DIR="${1:?Usage: create_mocks.sh <output_dir>}"
mkdir -p "$MOCK_DIR"

# --- Generic mock: just echo and exit 0 ---
create_generic_mock() {
    local name="$1"
    cat > "$MOCK_DIR/$name" <<STUB
#!/bin/bash
echo "[MOCK] $name \$*" >&2
exit 0
STUB
    chmod +x "$MOCK_DIR/$name"
}

# --- HMMBUILD: creates .hmm file from first arg ---
cat > "$MOCK_DIR/hmmbuild" <<'STUB'
#!/bin/bash
echo "[MOCK] hmmbuild $*" >&2
# hmmbuild <hmmfile> <msafile>
touch "$1"
exit 0
STUB
chmod +x "$MOCK_DIR/hmmbuild"

# --- HMMSEARCH: creates domtblout ---
cat > "$MOCK_DIR/hmmsearch" <<'STUB'
#!/bin/bash
echo "[MOCK] hmmsearch $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        --domtblout) echo "# mock domtblout" > "$2"; shift 2 ;;
        -o) touch "$2"; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/hmmsearch"

# --- MAFFT: copies input to stdout (identity alignment) ---
cat > "$MOCK_DIR/mafft" <<'STUB'
#!/bin/bash
echo "[MOCK] mafft $*" >&2
# Last argument is input file; copy it as "aligned"
args=("$@")
input="${args[-1]}"
if [ -f "$input" ]; then
    cat "$input"
else
    echo ">mock_seq"
    echo "ACDEFGHIKLMNPQRSTVWY"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/mafft"

# --- TRIMAL: copies input to output ---
cat > "$MOCK_DIR/trimal" <<'STUB'
#!/bin/bash
echo "[MOCK] trimal $*" >&2
infile=""
outfile=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -in) infile="$2"; shift 2 ;;
        -out) outfile="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$infile" ] && [ -n "$outfile" ]; then
    cp "$infile" "$outfile" 2>/dev/null || touch "$outfile"
elif [ -n "$infile" ]; then
    cat "$infile" 2>/dev/null
fi
exit 0
STUB
chmod +x "$MOCK_DIR/trimal"

# --- ClipKit: creates output file from input ---
cat > "$MOCK_DIR/clipkit" <<'STUB'
#!/bin/bash
echo "[MOCK] clipkit $*" >&2
input="$1"
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) cp "$input" "$2" 2>/dev/null || touch "$2"; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/clipkit"

# --- IQ-TREE: creates .treefile ---
cat > "$MOCK_DIR/iqtree2" <<'STUB'
#!/bin/bash
echo "[MOCK] iqtree2 $*" >&2
prefix=""
alignment=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s) alignment="$2"; shift 2 ;;
        --prefix|-pre) prefix="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$prefix" ]; then
    echo "((A:0.1,B:0.2):0.3,C:0.4);" > "${prefix}.treefile"
    echo "((A:0.1,B:0.2):0.3,C:0.4);" > "${prefix}.contree"
    printf "IQ-TREE mock log\nBest-fit model according to BIC: LG+G4\n" > "${prefix}.log"
    echo "mock iqtree iqtree" > "${prefix}.iqtree"
elif [ -n "$alignment" ]; then
    echo "((A:0.1,B:0.2):0.3,C:0.4);" > "${alignment}.treefile"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/iqtree2"

# --- FastTree: writes Newick to stdout ---
cat > "$MOCK_DIR/FastTree" <<'STUB'
#!/bin/bash
echo "[MOCK] FastTree $*" >&2
echo "((A:0.1,B:0.2):0.3,C:0.4);"
exit 0
STUB
chmod +x "$MOCK_DIR/FastTree"

# --- DeepTMHMM: creates output directory with predicted_topologies ---
cat > "$MOCK_DIR/deeptmhmm" <<'STUB'
#!/bin/bash
echo "[MOCK] deeptmhmm $*" >&2
fasta=""
outdir="deeptmhmm_output"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta|-f) fasta="$2"; shift 2 ;;
        -od|--output-dir|-o) outdir="$2"; shift 2 ;;
        *) shift ;;
    esac
done
mkdir -p "$outdir"
if [ -f "$fasta" ]; then
    # Generate mock 'prediction' file matching DeepTMHMM output format
    # Format: ID prediction_type confidence extra tm_count
    > "$outdir/prediction"
    grep "^>" "$fasta" | sed 's/>//' | while read -r seqid; do
        echo -e "${seqid}\tTMH\t0.95\tSP\t7" >> "$outdir/prediction"
    done
    # Also create other output files
    cp "$outdir/prediction" "$outdir/predicted_topologies.3line" 2>/dev/null
fi
touch "$outdir/TMRs.gff3"
exit 0
STUB
chmod +x "$MOCK_DIR/deeptmhmm"

# --- seqtk: subseq extracts named sequences ---
cat > "$MOCK_DIR/seqtk" <<'STUB'
#!/bin/bash
echo "[MOCK] seqtk $*" >&2
if [ "$1" = "subseq" ] && [ -f "$2" ] && [ -f "$3" ]; then
    python3 -c "
ids = set(line.strip().split()[0] for line in open('$3') if line.strip())
printing = False
for line in open('$2'):
    if line.startswith('>'):
        seqid = line[1:].strip().split()[0]
        printing = seqid in ids
    if printing:
        print(line, end='')
" 2>/dev/null || cat "$2"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/seqtk"

# --- HHblits ---
cat > "$MOCK_DIR/hhblits" <<'STUB'
#!/bin/bash
echo "[MOCK] hhblits $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o)
            cat > "$2" <<'HHR'
Query         mock_query
Match_columns 200

No 1
>ref_1 mock reference
Probab=99.9 E-value=1e-50 Score=200
HHR
            shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/hhblits"

# --- HHmake ---
cat > "$MOCK_DIR/hhmake" <<'STUB'
#!/bin/bash
echo "[MOCK] hhmake $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) touch "$2"; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/hhmake"

# --- OrthoFinder ---
cat > "$MOCK_DIR/orthofinder" <<'STUB'
#!/bin/bash
echo "[MOCK] orthofinder $*" >&2
input_dir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -f) input_dir="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$input_dir" ]; then
    results_dir="$input_dir/OrthoFinder/Results_mock"
    mkdir -p "$results_dir/Orthogroups"
    # Gene count table
    echo -e "Orthogroup\ttaxid_berghia\ttaxid1\ttaxid2" > "$results_dir/Orthogroups/Orthogroups.GeneCount.tsv"
    echo -e "OG0000001\t2\t1\t1" >> "$results_dir/Orthogroups/Orthogroups.GeneCount.tsv"
    echo -e "OG0000002\t3\t1\t0" >> "$results_dir/Orthogroups/Orthogroups.GeneCount.tsv"
    # Orthogroup FASTA files
    for og in OG0000001 OG0000002; do
        echo -e ">taxid_berghia_prot_1\nACDEFGHIKLMNPQRSTVWY" > "$results_dir/Orthogroups/${og}.fa"
        echo -e ">taxid1_prot_1\nACDEFGHIKLMNPQRSTVWY" >> "$results_dir/Orthogroups/${og}.fa"
        echo -e ">ref_1\nGHIKLMNPQRSTVWYACDEF" >> "$results_dir/Orthogroups/${og}.fa"
    done
    # Orthogroups.txt
    echo "OG0000001: taxid_berghia_prot_1 taxid_berghia_prot_2 taxid1_prot_1 taxid2_prot_1 ref_1" > "$results_dir/Orthogroups/Orthogroups.txt"
    echo "OG0000002: taxid_berghia_prot_3 taxid_berghia_prot_4 taxid_berghia_prot_5 taxid1_prot_2 ref_2" >> "$results_dir/Orthogroups/Orthogroups.txt"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/orthofinder"

# --- BUSCO ---
cat > "$MOCK_DIR/busco" <<'STUB'
#!/bin/bash
echo "[MOCK] busco $*" >&2
outdir=""
outpath=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) outdir="$2"; shift 2 ;;
        --out_path) outpath="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$outpath" ] && [ -n "$outdir" ]; then
    busco_dir="$outpath/$outdir/run_metazoa_odb10/busco_sequences/single_copy_busco_sequences"
    mkdir -p "$busco_dir"
    echo -e ">mock_busco_1\nACDEFGHIKLMNPQRSTVWY" > "$busco_dir/100at33208.faa"
    echo -e ">mock_busco_2\nGHIKLMNPQRSTVWYACDEF" > "$busco_dir/200at33208.faa"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/busco"

# --- ASTRAL ---
cat > "$MOCK_DIR/astral" <<'STUB'
#!/bin/bash
echo "[MOCK] astral $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) echo "((taxid_berghia:0.1,taxid1:0.2):0.3,taxid2:0.4);" > "$2"; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/astral"

# --- HyPhy ---
cat > "$MOCK_DIR/hyphy" <<'STUB'
#!/bin/bash
echo "[MOCK] hyphy $*" >&2
output=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) output="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$output" ]; then
    cat > "$output" <<'JSON'
{
  "fits": {"Full model": {"Rate Distributions": {}}},
  "branch attributes": {"0": {
    "taxid_berghia_prot_1": {"Rate classes": 1, "Uncorrected P-value": 0.5, "original name": "taxid_berghia_prot_1"}
  }},
  "test results": {"positive test results": 0}
}
JSON
fi
exit 0
STUB
chmod +x "$MOCK_DIR/hyphy"

# --- FastML ---
cat > "$MOCK_DIR/fastml" <<'STUB'
#!/bin/bash
echo "[MOCK] fastml $*" >&2
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outDir) outdir="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$outdir" ]; then
    mkdir -p "$outdir"
    echo -e ">N1\nACDEFGHIKLMNPQRSTVWY" > "$outdir/seq.marginal_IndelAndChars.txt"
    echo "((A:0.1,B:0.2):0.3,C:0.4);" > "$outdir/tree.newick.txt"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/fastml"

# --- Minimap2 ---
cat > "$MOCK_DIR/minimap2" <<'STUB'
#!/bin/bash
echo "[MOCK] minimap2 $*" >&2
echo "@HD	VN:1.6	SO:coordinate"
echo "@SQ	SN:scaffold_1	LN:3000"
exit 0
STUB
chmod +x "$MOCK_DIR/minimap2"

# --- Samtools ---
cat > "$MOCK_DIR/samtools" <<'STUB'
#!/bin/bash
echo "[MOCK] samtools $*" >&2
case "$1" in
    sort)
        # Find output file from -o flag or last arg
        outfile=""
        shift
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -o) outfile="$2"; shift 2 ;;
                *) shift ;;
            esac
        done
        [ -n "$outfile" ] && touch "$outfile"
        ;;
    index) touch "${2}.bai" 2>/dev/null || true ;;
    view) cat ;;
    *) ;;
esac
exit 0
STUB
chmod +x "$MOCK_DIR/samtools"

# --- BLAST ---
cat > "$MOCK_DIR/makeblastdb" <<'STUB'
#!/bin/bash
echo "[MOCK] makeblastdb $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        -out) touch "$2.pdb" "$2.phr" "$2.pin" "$2.psq" 2>/dev/null; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/makeblastdb"

cat > "$MOCK_DIR/blastp" <<'STUB'
#!/bin/bash
echo "[MOCK] blastp $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        -out) echo "# mock blast output" > "$2"; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/blastp"

# --- MCScanX ---
cat > "$MOCK_DIR/MCScanX" <<'STUB'
#!/bin/bash
echo "[MOCK] MCScanX $*" >&2
prefix="$1"
if [ -n "$prefix" ]; then
    cat > "${prefix}.collinearity" <<'COL'
############### Parameters ###############
# MATCH_SCORE: 50
## Alignment 0: score=100 e_value=1e-10 N=2 taxid_berghia&taxid1
 0-  0:	taxid_berghia_gene_1	taxid1_gene_1	  1e-50
 0-  1:	taxid_berghia_gene_2	taxid1_gene_2	  1e-40
COL
fi
exit 0
STUB
chmod +x "$MOCK_DIR/MCScanX"

# --- AlphaFold ---
cat > "$MOCK_DIR/run_alphafold.sh" <<'STUB'
#!/bin/bash
echo "[MOCK] run_alphafold.sh $*" >&2
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output_dir) outdir="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$outdir" ]; then
    mkdir -p "$outdir"
    cat > "$outdir/ranked_0.pdb" <<'PDB'
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C
END
PDB
fi
exit 0
STUB
chmod +x "$MOCK_DIR/run_alphafold.sh"

# --- FoldTree ---
cat > "$MOCK_DIR/foldtree" <<'STUB'
#!/bin/bash
echo "[MOCK] foldtree $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output) echo "((struct_A:0.1,struct_B:0.2):0.3,struct_C:0.4);" > "$2"; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/foldtree"

# --- cd-hit ---
cat > "$MOCK_DIR/cd-hit" <<'STUB'
#!/bin/bash
echo "[MOCK] cd-hit $*" >&2
input=""
output=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) input="$2"; shift 2 ;;
        -o) output="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$input" ] && [ -n "$output" ]; then
    cp "$input" "$output"
    echo ">Cluster 0" > "${output}.clstr"
    echo "0	200aa, >seq1... *" >> "${output}.clstr"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/cd-hit"

# --- InterProScan ---
cat > "$MOCK_DIR/interproscan.sh" <<'STUB'
#!/bin/bash
echo "[MOCK] interproscan.sh $*" >&2
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) echo -e "seq1\tmd5\t200\tPfam\tPF00001\t7tm_1\t10\t280\t1e-50\tT\t01-01-2024" > "$2"; shift 2 ;;
        *) shift ;;
    esac
done
exit 0
STUB
chmod +x "$MOCK_DIR/interproscan.sh"

# --- cafe5 ---
cat > "$MOCK_DIR/cafe5" <<'STUB'
#!/bin/bash
echo "[MOCK] cafe5 $*" >&2
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) outdir="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "$outdir" ]; then
    mkdir -p "$outdir"
    echo -e "FamilyID\ttaxid_berghia<0>\ttaxid1<1>\ttaxid2<2>" > "$outdir/Base_clade_results.txt"
    echo -e "OG0000001\t2\t1\t1" >> "$outdir/Base_clade_results.txt"
    echo "Lambda: 0.001" > "$outdir/Base_results.txt"
fi
exit 0
STUB
chmod +x "$MOCK_DIR/cafe5"

# --- pal2nal ---
cat > "$MOCK_DIR/pal2nal.pl" <<'STUB'
#!/bin/bash
echo "[MOCK] pal2nal.pl $*" >&2
echo "  2  12"
echo "seq1"
echo "ATGATGATGATG"
echo "seq2"
echo "ATGATGATGATG"
exit 0
STUB
chmod +x "$MOCK_DIR/pal2nal.pl"

# --- pdflatex ---
cat > "$MOCK_DIR/pdflatex" <<'STUB'
#!/bin/bash
echo "[MOCK] pdflatex $*" >&2
for arg in "$@"; do
    if [[ "$arg" == *.tex ]]; then
        touch "${arg%.tex}.pdf"
    fi
done
exit 0
STUB
chmod +x "$MOCK_DIR/pdflatex"

# --- Simple generic mocks ---
create_generic_mock "TMalign"
create_generic_mock "mb"
create_generic_mock "phyloformer"
# --- Rscript: detect output file arg and create it (e.g. ultrametric tree) ---
cat > "$MOCK_DIR/Rscript" <<'STUB'
#!/bin/bash
echo "[MOCK] Rscript $*" >&2
# generate_ultrametric.R <input_tree> <output_tree> [method] [lambda]
# Find output file: second non-.R, non-flag argument
FOUND_SCRIPT=false
ARG_IDX=0
for arg in "$@"; do
    if [[ "$arg" == *.R ]]; then
        FOUND_SCRIPT=true
        continue
    fi
    if $FOUND_SCRIPT; then
        ARG_IDX=$((ARG_IDX + 1))
        if [ "$ARG_IDX" -eq 2 ]; then
            mkdir -p "$(dirname "$arg")"
            echo "(taxid_berghia:1.0,taxid1:1.0,taxid2:1.0);" > "$arg"
            break
        fi
    fi
done
exit 0
STUB
chmod +x "$MOCK_DIR/Rscript"
create_generic_mock "java"

echo "Created $(ls "$MOCK_DIR" | wc -l) mock tools in $MOCK_DIR" >&2
