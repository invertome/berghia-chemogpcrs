#!/bin/bash
#SBATCH --job-name=af3_weights
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/af3_weights-%j.out
#SBATCH --error=logs/af3_weights-%j.err
# download_af3_weights.sh — pull AlphaFold 3 model weights into a Unity
# scratch path so stage 08 can use them. Bead -05o.
#
# AF3 weights workflow (per https://github.com/google-deepmind/alphafold3
# and the WEIGHTS_TERMS_OF_USE):
#
#   1. A human (you) must request access by filling out the form at
#      https://forms.gle/svvpY4u2jsHEwWYS6 and accepting the terms at
#      https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md.
#      Google approves manually, typically within 1-3 business days.
#
#   2. The approval email includes a one-time signed download URL
#      (typically expires in ~24 h). Pass that URL into this script.
#
#   3. This script downloads, extracts (if tarballed), validates the
#      manifest, and reports the path to set as ALPHAFOLD3_MODEL_DIR in
#      config.sh.
#
# Submit:
#   sbatch scripts/unity/download_af3_weights.sh '<google-signed-url>'
# or pre-export AF3_WEIGHTS_URL and submit with no positional arg:
#   AF3_WEIGHTS_URL='<url>' sbatch scripts/unity/download_af3_weights.sh
#
# Idempotent: re-running with an already-present weights tarball/dir is
# a no-op. Cluster-friendly: respects --time of 4 h, default 2 CPU, 8 GB.

set -eo pipefail
mkdir -p logs

URL="${1:-${AF3_WEIGHTS_URL:-}}"
DEST="${AF3_WEIGHTS_DEST:-/scratch3/workspace/${USER}-jorge/alphafold3_models}"

if [ -z "$URL" ]; then
    cat <<EOF >&2
[af3-weights] ERROR: No download URL given.
  Pass as positional arg or AF3_WEIGHTS_URL env var.
  Request access first: https://forms.gle/svvpY4u2jsHEwWYS6
EOF
    exit 2
fi

mkdir -p "$DEST"
cd "$DEST"

# Skip if weights already present.
# AF3 weights ship as a *.tar archive (Google's distribution format).
# Once extracted, the model directory contains *.npz files for each
# parameter set. Probe for either.
if find "$DEST" -maxdepth 3 -name 'af3.bin*' 2>/dev/null | grep -q . \
   || find "$DEST" -maxdepth 3 -iname '*.npz' 2>/dev/null | grep -q .; then
    echo "[af3-weights] Weights already present under $DEST — skipping download."
    echo "[af3-weights] Set in config.sh:"
    echo "    export ALPHAFOLD3_MODEL_DIR=\"$DEST\""
    exit 0
fi

# Pull the weights. Use curl with retry and resume so a transient network
# glitch doesn't waste the queue slot. Single-stream — Google-signed URLs
# don't support range-parallel downloads.
ARCHIVE="af3_params_$(date +%Y%m%d).tar"
echo "[af3-weights] Downloading AF3 weights to $DEST/$ARCHIVE..."
curl -fLC - --retry 5 --retry-delay 30 -o "$ARCHIVE" "$URL"
echo "[af3-weights] Download complete: $(du -h "$ARCHIVE" | cut -f1)"

# Extract if it's a tar archive; otherwise assume Google delivered a
# directory bundle or single file and leave as-is.
if tar -tf "$ARCHIVE" >/dev/null 2>&1; then
    echo "[af3-weights] Extracting tar archive..."
    tar -xf "$ARCHIVE"
    rm -f "$ARCHIVE"
fi

# Surface the resolved model directory.
FOUND="$(find "$DEST" -maxdepth 4 -iname '*.npz' -printf '%h\n' 2>/dev/null | sort -u | head -1)"
if [ -z "$FOUND" ]; then
    FOUND="$(find "$DEST" -maxdepth 4 -name 'af3.bin*' -printf '%h\n' 2>/dev/null | sort -u | head -1)"
fi
if [ -z "$FOUND" ]; then
    echo "[af3-weights] WARN: extracted, but no .npz / af3.bin* files found in $DEST." >&2
    echo "             Inspect the directory manually and set ALPHAFOLD3_MODEL_DIR appropriately." >&2
    ls -la "$DEST" >&2
    exit 1
fi

echo "[af3-weights] Done. AF3 model directory: $FOUND"
echo "[af3-weights] Add to config.sh:"
echo "    export ALPHAFOLD3_MODEL_DIR=\"$FOUND\""
