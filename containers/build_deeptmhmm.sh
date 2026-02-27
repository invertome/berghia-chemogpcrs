#!/bin/bash
# build_deeptmhmm.sh - Build DeepTMHMM Apptainer container
# Requires: DeepTMHMM academic license zip + Apptainer installed
#
# IMPORTANT: Build this on the target architecture (x86_64 for HPC, aarch64 for ARM).
# A container built on one arch won't run on the other.
#
# The ESM1b model is 2.6GB — the build needs ~10GB of working space.
# By default uses /tmp. Set APPTAINER_TMPDIR to override if /tmp is too small.
#
# Usage: bash build_deeptmhmm.sh [/path/to/DeepTMHMM-Academic-License-v1.0.zip]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ZIP="${1:-$HOME/Downloads/DeepTMHMM-Academic-License-v1.0.zip}"
SIF="$SCRIPT_DIR/deeptmhmm.sif"
DEF="$SCRIPT_DIR/deeptmhmm.def"

if [ ! -f "$ZIP" ]; then
    echo "ERROR: DeepTMHMM zip not found: $ZIP" >&2
    echo "Usage: $0 [/path/to/DeepTMHMM-Academic-License-v1.0.zip]" >&2
    exit 1
fi

if ! command -v apptainer &>/dev/null; then
    echo "ERROR: apptainer not found. Install Apptainer first." >&2
    exit 1
fi

if [ ! -f "$DEF" ]; then
    echo "ERROR: Definition file not found: $DEF" >&2
    exit 1
fi

# Use a build directory with sufficient space (ESM model = 2.6GB)
# Prefer APPTAINER_TMPDIR if set, then /tmp, then current dir
BUILD_BASE="${APPTAINER_TMPDIR:-${TMPDIR:-/tmp}}"
BUILD_DIR=$(mktemp -d "$BUILD_BASE/deeptmhmm_build.XXXXXX")
trap "rm -rf '$BUILD_DIR'" EXIT

AVAIL_GB=$(df --output=avail -BG "$BUILD_DIR" | tail -1 | tr -d ' G')
if [ "$AVAIL_GB" -lt 10 ]; then
    echo "WARNING: Only ${AVAIL_GB}GB available in $BUILD_BASE." >&2
    echo "Build needs ~10GB. Set APPTAINER_TMPDIR to a larger filesystem." >&2
fi

echo "Build directory: $BUILD_DIR (${AVAIL_GB}GB available)"
echo "Extracting DeepTMHMM..."
unzip -q "$ZIP" -d "$BUILD_DIR"

echo "Architecture: $(uname -m)"
echo "Building Apptainer container (this may take 15-30 minutes)..."
echo ""

# Build from the extracted directory so %files can find DeepTMHMM-Academic-License-v1.0/
cd "$BUILD_DIR"
apptainer build "$SIF" "$DEF"

echo ""
echo "=== Build complete ==="
echo "Container: $SIF"
echo "Size: $(du -h "$SIF" | cut -f1)"
echo "Architecture: $(uname -m)"
echo ""
echo "Test with:"
echo "  apptainer run $SIF --fasta sample.fasta --output-dir /tmp/test_dtmhmm"
