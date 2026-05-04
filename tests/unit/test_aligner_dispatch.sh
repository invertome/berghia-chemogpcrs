#!/bin/bash
# Test the regime-based aligner dispatch logic in run_aligner.sh.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
WRAPPER="$PROJECT_DIR/scripts/run_aligner.sh"

if [ ! -x "$WRAPPER" ]; then
    echo "FAIL: wrapper not executable: $WRAPPER" >&2
    exit 1
fi
if ! command -v mafft &>/dev/null; then
    echo "SKIP: mafft not installed in this env" >&2
    exit 0
fi

TMPDIR=$(mktemp -d /tmp/aligner_test_XXXX)
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

gen_fasta() {
    local n="$1" outfile="$2"
    python3 - "$n" "$outfile" <<'PYEOF'
import random, sys
n = int(sys.argv[1]); outfile = sys.argv[2]
random.seed(42)
aa = "ACDEFGHIKLMNPQRSTVWY"
with open(outfile, "w") as f:
    for i in range(n):
        f.write(f">seq_{i:05d}\n")
        f.write("".join(random.choices(aa, k=120)) + "\n")
PYEOF
}

assert_backend() {
    local expected="$1" prov_file="$2"
    local got
    got=$(awk -F= '/^backend=/{print $2; exit}' "$prov_file")
    if [ "$got" != "$expected" ]; then
        echo "FAIL: expected backend=$expected, got backend=$got"
        cat "$prov_file"
        exit 1
    fi
    echo "OK: backend=$expected"
}

# N=10 → mafft_linsi
gen_fasta 10 "$TMPDIR/small.fa"
"$WRAPPER" --input="$TMPDIR/small.fa" --output="$TMPDIR/small_aln.fa" \
    --threads=2 --linsi-threshold=200 --famsa-threshold=1000 \
    >/dev/null 2>"$TMPDIR/small.log"
assert_backend "mafft_linsi" "$TMPDIR/small_aln.fa.aligner_used.txt"

# N=300 → mafft_auto
gen_fasta 300 "$TMPDIR/mid.fa"
"$WRAPPER" --input="$TMPDIR/mid.fa" --output="$TMPDIR/mid_aln.fa" \
    --threads=2 --linsi-threshold=200 --famsa-threshold=1000 \
    >/dev/null 2>"$TMPDIR/mid.log"
assert_backend "mafft_auto" "$TMPDIR/mid_aln.fa.aligner_used.txt"

# N=1500 → famsa preferred but falls back to mafft_auto if FAMSA missing
gen_fasta 1500 "$TMPDIR/large.fa"
"$WRAPPER" --input="$TMPDIR/large.fa" --output="$TMPDIR/large_aln.fa" \
    --threads=2 --linsi-threshold=200 --famsa-threshold=1000 \
    >/dev/null 2>"$TMPDIR/large.log"
backend_got=$(awk -F= '/^backend=/{print $2; exit}' "$TMPDIR/large_aln.fa.aligner_used.txt")
if [ "$backend_got" = "famsa" ] || [ "$backend_got" = "mafft_auto" ]; then
    echo "OK: large case backend=$backend_got (famsa preferred, mafft_auto fallback)"
else
    echo "FAIL: large N expected famsa or mafft_auto; got backend=$backend_got"
    exit 1
fi

# --force-backend override
gen_fasta 50 "$TMPDIR/forced.fa"
"$WRAPPER" --input="$TMPDIR/forced.fa" --output="$TMPDIR/forced_aln.fa" \
    --threads=2 --force-backend=mafft_retree2 \
    >/dev/null 2>"$TMPDIR/forced.log"
assert_backend "mafft_retree2" "$TMPDIR/forced_aln.fa.aligner_used.txt"

# Custom thresholds: N=300 with linsi_threshold=500 → mafft_linsi
gen_fasta 300 "$TMPDIR/custom.fa"
"$WRAPPER" --input="$TMPDIR/custom.fa" --output="$TMPDIR/custom_aln.fa" \
    --threads=2 --linsi-threshold=500 --famsa-threshold=2000 \
    >/dev/null 2>"$TMPDIR/custom.log"
assert_backend "mafft_linsi" "$TMPDIR/custom_aln.fa.aligner_used.txt"

# Output FASTA is aligned (uniform seq length)
seq_lengths=$(awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq $0} END {if (seq) print length(seq)}' \
    "$TMPDIR/small_aln.fa" | sort -u)
n_distinct=$(echo "$seq_lengths" | wc -l)
if [ "$n_distinct" -ne 1 ]; then
    echo "FAIL: aligned FASTA has $n_distinct distinct sequence lengths"
    exit 1
fi
echo "OK: aligned FASTA has uniform sequence length"

echo ""
echo "PASS test_aligner_dispatch.sh"
