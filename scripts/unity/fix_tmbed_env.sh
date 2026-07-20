#!/bin/bash
# fix_tmbed_env.sh — restore the TMbed install in the berghia-gpcr conda env
# after a transformers/tokenizers upgrade breaks it.
#
# Postmortem (2026-05-08, stage 02 chain failure):
#
#   The pipeline wraps TMbed via scripts/run_deeptmhmm.sh. After a routine
#   `pip` operation in the env upgraded transformers to v5, four bugs
#   stacked up and crashed stage 02 in 40 s:
#
#     (1) `-f-format=fasta` in the wrapper — phantom flag, TMbed has no
#         such option. Fixed in scripts/run_deeptmhmm.sh.
#
#     (2) transformers v5 dropped `T5Tokenizer.batch_encode_plus`, which
#         tmbed/embed.py:86 still calls. Pin transformers<5.
#
#     (3) tmbed ships a bundled local copy of the ProtT5 encoder under
#         site-packages/tmbed/models/t5/. The previous install ran
#         `tok.save_pretrained()` under v5 and dropped the slow
#         tokenizer's `spiece.model` (a known transformers gotcha — only
#         the FAST tokenizer is saved by default). Without spiece.model,
#         SentencePieceProcessor.LoadFromFile gets None at load time.
#
#     (4) The bundled tokenizer_config.json was written by v5 with
#         `extra_special_tokens` as a *list*, but v4 expects a *dict*
#         (raises `'list' object has no attribute 'keys'`). Re-saving the
#         tokenizer under v4 fixes it.
#
#   This script applies fixes (2)–(4) idempotently. Fix (1) is in the
#   wrapper script and ships via git.
#
# Usage:
#   bash scripts/unity/fix_tmbed_env.sh
#
# Run interactively on the login node (needs HuggingFace network access).
# A leftover `tmp...` directory inside the bundled t5/ dir from an
# interrupted previous install is also cleaned up.

# `set -u` would trip on conda's deactivate hooks (CONDA_BACKUP_CXX etc.),
# so only enforce -e and -o pipefail. The script doesn't rely on unset-var
# safety beyond what conda needs.
set -eo pipefail

ENV_NAME="${1:-berghia-gpcr}"

# shellcheck disable=SC1091
source ~/.miniconda3/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

PKG_DIR="$(python3 -c 'import tmbed, os; print(os.path.dirname(tmbed.__file__))')"
T5_DIR="$PKG_DIR/models/t5"

echo "[fix-tmbed] tmbed install: $PKG_DIR"
echo "[fix-tmbed] bundled t5: $T5_DIR"

# (2) Pin transformers<5 (the API tmbed depends on).
echo "[fix-tmbed] Pinning transformers<5..."
pip install --quiet 'transformers<5'
TF_VER="$(python3 -c 'import transformers; print(transformers.__version__)')"
echo "[fix-tmbed]   transformers $TF_VER"

# (2b) Ensure torch retains CUDA support (bead -m1f, 2026-05-13).
#
# A prior run of this script reinstalled transformers from the default PyPI
# index, which pulled the CPU-only `torch-*+cpu` wheel as a side dependency.
# That kept TMbed importable but silently broke GPU inference: every stage 02
# run after that point reported `torch.cuda.is_available() = False` and burned
# its GPU allocation on CPU-bound inference (~3 orders of magnitude slower).
#
# Reinstall torch from the cu118 wheel index. cu118 retains Pascal (sm_61)
# support, which is required for the gypsum-gpu* nodes' GTX 1080 Ti — newer
# cu124 wheels have dropped sm_61. Override via TORCH_INDEX_URL if you
# specifically need a different CUDA toolchain.
TORCH_INDEX="${TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu118}"
echo "[fix-tmbed] Reinstalling torch with CUDA support from $TORCH_INDEX..."
pip install --quiet --upgrade --force-reinstall --index-url "$TORCH_INDEX" torch
python3 - <<PY
import torch
build = torch.version.cuda
flavor = '+cpu' if '+cpu' in torch.__version__ else (f'+cu{build.replace(".", "")}' if build else '?')
print(f'[fix-tmbed]   torch {torch.__version__} (cuda build: {build}, flavor: {flavor})')
if build is None or '+cpu' in torch.__version__:
    raise SystemExit('[fix-tmbed] FATAL: torch is still CPU-only after reinstall; '
                     'check TORCH_INDEX_URL and pip resolver output.')
PY

# (3+4) Clear any tokenizer files written by a different transformers
# major and re-save under the pinned version. This regenerates spiece.model
# (a slow-tokenizer file that save_pretrained skips when only the fast
# tokenizer is loaded) and writes a v4-compatible tokenizer_config.json.
echo "[fix-tmbed] Re-saving bundled T5 tokenizer under transformers v$TF_VER..."
python3 - <<PY
import os, shutil
from transformers import T5Tokenizer
DST = "${T5_DIR}"
for f in ('tokenizer_config.json', 'tokenizer.json', 'spiece.model',
          'special_tokens_map.json', 'added_tokens.json'):
    p = os.path.join(DST, f)
    if os.path.exists(p):
        os.remove(p)
        print(f'  rm {f}')
# Clean up any stranded tmp* directory from an interrupted previous install
for d in os.listdir(DST):
    if d.startswith('tmp') and os.path.isdir(os.path.join(DST, d)):
        shutil.rmtree(os.path.join(DST, d))
        print(f'  rm tmp dir {d}')
tok = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc',
                                  do_lower_case=False)
tok.save_pretrained(DST)
print('  saved tokenizer to', DST)
print('  files:', sorted(os.listdir(DST)))
PY

# Smoke test: 3 short Berghia proteins, CPU, no model_path arg (uses bundled).
echo "[fix-tmbed] Smoke test..."
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cat > "$TMP/tiny.fa" <<'FA'
>smoke_1
MSPRSFRDTKDDCLDQLEDAADQLQQGGVAAVHLRTVPRAEAPPLWEQLRSRLLRSPTFQ
>smoke_2
MLKNRIAKYREKKREKKEIAAGKNVEEEEDKVEPDTPSDTFMTIENEINLSKLEVS
>smoke_3
MRWLPEFKGAVKVGGSTKFSSPCFGSNSITATQKDDTTIVLEFQPSQKKSLLCYDWYFI
FA
if tmbed predict -f "$TMP/tiny.fa" -p "$TMP/out.3line" --out-format 3 \
        2>"$TMP/tmbed.log"; then
    n="$(grep -c '^>' "$TMP/out.3line" 2>/dev/null || true)"
    n=${n:-0}
    echo "[fix-tmbed]   OK — predicted $n / 3 sequences"
else
    echo "[fix-tmbed]   FAILED — see $TMP/tmbed.log:" >&2
    tail -20 "$TMP/tmbed.log" >&2
    exit 2
fi

echo "[fix-tmbed] Done. To run stage 02 on Unity, sbatch with a GPU partition"
echo "[fix-tmbed]   for ProtT5 inference at scale (CPU is ~30s/seq, prohibitive"
echo "[fix-tmbed]   for the 86k Berghia transcriptome)."
