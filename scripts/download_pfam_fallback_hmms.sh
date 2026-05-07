#!/bin/bash
# download_pfam_fallback_hmms.sh — Phase 2 Task 2.3 of the
# non-chemoreceptor classification pipeline.
#
# Downloads Pfam HMMs for the four canonical 7TM GPCR class signatures.
# These provide a fallback for candidates that don't hit any custom HMM
# (built in Task 2.1 from our curated reference set) — at minimum we'll
# get a coarse Class A/B/C/F call.
#
# Pfam HMMs:
#   PF00001  7tm_1       Class A 7-transmembrane (rhodopsin-like)
#   PF00002  7tm_2       Class B secretin/adhesion 7TM
#   PF00003  7tm_3       Class C glutamate-like 7TM
#   PF01534  Frizzled    Class F frizzled cysteine-rich domain (ligand-binding)
#
# Output:
#   results/classification/hmms/pfam_fallback/<accession>.hmm
#   results/classification/hmms/pfam_fallback/manifest.tsv
#
# Usage:  bash scripts/download_pfam_fallback_hmms.sh

set -euo pipefail

OUTPUT_DIR="${OUTPUT_DIR:-results/classification/hmms/pfam_fallback}"
mkdir -p "$OUTPUT_DIR"

# Pfam HMM download endpoint (InterPro re-exports the Pfam HMMs).
# Each HMM is a single text file at:
#   https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<accession>?annotation=hmm
PFAM_HMMS=(
    "PF00001:7tm_1:class-A-7tm"
    "PF00002:7tm_2:class-B-7tm"
    "PF00003:7tm_3:class-C-7tm"
    "PF01534:Frizzled:class-F-frizzled"
)

MANIFEST="$OUTPUT_DIR/manifest.tsv"
echo -e "pfam_accession\tpfam_name\tclassification_label\thmm_path" > "$MANIFEST"

for entry in "${PFAM_HMMS[@]}"; do
    IFS=":" read -r acc name label <<< "$entry"
    out="$OUTPUT_DIR/$acc.hmm"
    if [[ -s "$out" ]]; then
        echo "[pfam] $acc ($name) — already present, skipping"
    else
        echo "[pfam] Downloading $acc ($name)..."
        url="https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/$acc?annotation=hmm"
        # InterPro returns gzipped HMM
        if curl -sf "$url" -o "$out.gz" 2>/dev/null; then
            gunzip -f "$out.gz"
        else
            echo "  WARN: download failed for $acc; trying Pfam direct"
            # Fallback: try Pfam's direct endpoint (may need different format)
            curl -sf "https://pfam.xfam.org/family/$acc/hmm" -o "$out" || {
                echo "  ERROR: both download paths failed for $acc"
                continue
            }
        fi
    fi
    if [[ -s "$out" ]]; then
        echo -e "${acc}\t${name}\t${label}\t${out}" >> "$MANIFEST"
    fi
done

echo ""
echo "[pfam] DONE. Library at $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"
echo ""
cat "$MANIFEST"
