#!/bin/bash
# Regression: run_deeptmhmm.sh must abort (exit 3) when a GPU is allocated
# but torch.cuda.is_available() returns False — otherwise TMbed silently
# falls back to CPU inference, turning a ~30-minute stage 02 run into a
# 1-3 day wall-time burn (lesson learned bead -m1f, 2026-05-13).
#
# Coverage (mirrors the guard logic at scripts/run_deeptmhmm.sh:97-120):
#   A. GPU env set + cuda unavailable           → exit 3 + ERROR stderr
#   B. GPU env set + cuda unavailable + ALLOW=1 → guard bypassed, no exit 3
#   C. No GPU env vars + cuda unavailable       → guard does NOT fire
#   D. GPU env set + cuda available             → guard does NOT fire
#
# Strategy: PATH-stub `tmbed` with a shebang pointing to a fake python
# that simulates torch.cuda.is_available() via $TEST_CUDA_AVAILABLE.
# The guard extracts the python via `head -1 $(command -v tmbed)`, so
# the shebang is the lever.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
RUN_TMHMM="$PROJECT_DIR/scripts/run_deeptmhmm.sh"

TMPDIR=$(mktemp -d /tmp/test_gpu_guard_XXXX)
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

mkdir -p "$TMPDIR/bin"

# Fake python: when invoked as `fake_python -c 'CODE'` (the guard's torch
# check), exit 0 if TEST_CUDA_AVAILABLE=1 else exit 2. Otherwise (being
# invoked via the tmbed shebang as the script interpreter for `tmbed
# predict ...`), exit 1 — we don't care what happens after the guard,
# the assertions only inspect stderr substrings + the exit-3 signal.
cat > "$TMPDIR/bin/fake_python" <<'EOF'
#!/bin/bash
if [ "${1:-}" = "-c" ]; then
    case "${TEST_CUDA_AVAILABLE:-0}" in
        1) exit 0 ;;
        *) exit 2 ;;
    esac
fi
exit 1
EOF
chmod +x "$TMPDIR/bin/fake_python"

# Fake tmbed: shebang points to the fake python. When the guard runs
# `head -1 "$(command -v tmbed)" | sed 's|^#!||'`, it gets the fake-python
# path; when the script later runs `tmbed predict ...`, the kernel invokes
# fake_python with this script as $0 (which exits 1, fine for our purposes).
cat > "$TMPDIR/bin/tmbed" <<EOF
#!$TMPDIR/bin/fake_python
EOF
chmod +x "$TMPDIR/bin/tmbed"

# Tiny input FASTA (content irrelevant; guard fires before tmbed reads it).
INPUT="$TMPDIR/input.fa"
cat > "$INPUT" <<'EOF'
>seq_test
MASTRSEQUENCE
EOF

export PATH="$TMPDIR/bin:$PATH"
GUARD_MSG='ERROR: GPU allocation present'
BYPASS_MSG='ALLOW_TMBED_CPU=1 set; continuing on CPU'

run_case() {
    local label="$1"
    local out_dir="$TMPDIR/out_${label}"
    mkdir -p "$out_dir"
    local err_log="$TMPDIR/err_${label}.log"
    # Run with the per-case env vars already exported by the caller
    bash "$RUN_TMHMM" -f "$INPUT" -o "$out_dir" >/dev/null 2>"$err_log"
    echo "$?:$err_log"
}

# -------- Case A: GPU env + cuda unavailable → exit 3 --------
(
    export SLURM_JOB_GPUS=1
    export TEST_CUDA_AVAILABLE=0
    unset ALLOW_TMBED_CPU 2>/dev/null || true
    result=$(run_case A)
    rc=${result%%:*}; err=${result#*:}
    [[ "$rc" == "3" ]] \
        || { echo "FAIL[A]: expected exit 3, got $rc" >&2; cat "$err" >&2; exit 1; }
    grep -qF "$GUARD_MSG" "$err" \
        || { echo "FAIL[A]: stderr missing guard error message" >&2; cat "$err" >&2; exit 1; }
    echo "OK[A]: exit 3 + guard ERROR fired on GPU+no-cuda"
)

# -------- Case B: GPU env + cuda unavailable + ALLOW_TMBED_CPU=1 → bypass --------
(
    export SLURM_JOB_GPUS=1
    export TEST_CUDA_AVAILABLE=0
    export ALLOW_TMBED_CPU=1
    result=$(run_case B)
    rc=${result%%:*}; err=${result#*:}
    [[ "$rc" != "3" ]] \
        || { echo "FAIL[B]: guard fired (exit 3) despite ALLOW_TMBED_CPU=1" >&2; cat "$err" >&2; exit 1; }
    grep -qF "$BYPASS_MSG" "$err" \
        || { echo "FAIL[B]: stderr missing bypass message ('$BYPASS_MSG')" >&2; cat "$err" >&2; exit 1; }
    echo "OK[B]: ALLOW_TMBED_CPU=1 bypasses guard with bypass-log message"
)

# -------- Case C: no GPU env vars → guard does not fire --------
(
    unset SLURM_JOB_GPUS SLURM_GPUS SLURM_GPUS_ON_NODE CUDA_VISIBLE_DEVICES 2>/dev/null || true
    export TEST_CUDA_AVAILABLE=0  # cuda unavail, but guard shouldn't even check
    unset ALLOW_TMBED_CPU 2>/dev/null || true
    result=$(run_case C)
    rc=${result%%:*}; err=${result#*:}
    [[ "$rc" != "3" ]] \
        || { echo "FAIL[C]: guard fired without any GPU env var" >&2; cat "$err" >&2; exit 1; }
    if grep -qF "$GUARD_MSG" "$err"; then
        echo "FAIL[C]: guard ERROR emitted without GPU env vars" >&2
        cat "$err" >&2
        exit 1
    fi
    echo "OK[C]: guard silent when no GPU env vars set"
)

# -------- Case D: GPU env + cuda AVAILABLE → guard does not fire --------
(
    export SLURM_JOB_GPUS=1
    export TEST_CUDA_AVAILABLE=1
    unset ALLOW_TMBED_CPU 2>/dev/null || true
    result=$(run_case D)
    rc=${result%%:*}; err=${result#*:}
    [[ "$rc" != "3" ]] \
        || { echo "FAIL[D]: guard fired despite cuda available" >&2; cat "$err" >&2; exit 1; }
    if grep -qF "$GUARD_MSG" "$err"; then
        echo "FAIL[D]: guard ERROR emitted despite cuda available" >&2
        cat "$err" >&2
        exit 1
    fi
    echo "OK[D]: guard silent when cuda is available"
)

# Coverage of the other 3 GPU env vars (SLURM_GPUS, SLURM_GPUS_ON_NODE,
# CUDA_VISIBLE_DEVICES) — any one should trigger the guard the same way
# SLURM_JOB_GPUS does.
for var in SLURM_GPUS SLURM_GPUS_ON_NODE CUDA_VISIBLE_DEVICES; do
    (
        unset SLURM_JOB_GPUS SLURM_GPUS SLURM_GPUS_ON_NODE CUDA_VISIBLE_DEVICES 2>/dev/null || true
        export "$var=1"
        export TEST_CUDA_AVAILABLE=0
        unset ALLOW_TMBED_CPU 2>/dev/null || true
        result=$(run_case "Aalt_${var}")
        rc=${result%%:*}; err=${result#*:}
        [[ "$rc" == "3" ]] \
            || { echo "FAIL[${var}]: expected exit 3, got $rc" >&2; cat "$err" >&2; exit 1; }
        grep -qF "$GUARD_MSG" "$err" \
            || { echo "FAIL[${var}]: stderr missing guard error" >&2; cat "$err" >&2; exit 1; }
    )
done
echo "OK[alt-envs]: guard fires for SLURM_GPUS / SLURM_GPUS_ON_NODE / CUDA_VISIBLE_DEVICES"

echo "PASS test_run_deeptmhmm_gpu_guard.sh"
