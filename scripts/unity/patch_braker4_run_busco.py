#!/usr/bin/env python3
"""Idempotent patch: guard `rule busco_proteins` in run_busco.smk against
runaway gene models (>50,000 aa).

A single >50k-aa "protein" (a fused/runaway BRAKER model) makes BUSCO's
internal hmmsearch abort with SIGABRT, which under `set -euo pipefail` fails
the whole genome — even though the annotation itself is fine (observed on the
masking pilot, Nautilus de-novo-reuse arm: gene g1134 had a 157,627-aa isoform).
Real BUSCO orthologs are never that long, so excluding them does not change the
completeness score.

This patch (berghia-chemogpcrs, bead -pfo):
  A. filters proteins >50,000 aa out of the BUSCO input before running BUSCO;
  B. makes the BUSCO invocation non-fatal (a QC metric must not block the
     annotation), still touching the rule's .done output.

All edits are scoped to the `rule busco_proteins:` block so the sibling
busco_genome rule is untouched. Mirrors scripts/unity/patch_braker4_run_compleasm.py:
copy next to run_busco.smk (rules/quality_control/) and run, or run from that dir.
"""
import sys
from pathlib import Path

path = Path(__file__).parent / "run_busco.smk"
if not path.exists():
    path = Path("run_busco.smk")
if not path.exists():
    print(f"ERROR: run_busco.smk not found next to {__file__} or in CWD", file=sys.stderr)
    sys.exit(1)

src = path.read_text()
MARKER = "berghia-chemogpcrs guard: drop runaway"
if MARKER in src:
    print("Already patched — no-op.")
    sys.exit(0)

# Scope every edit to the busco_proteins rule block.
start = src.find("rule busco_proteins:")
if start < 0:
    print("ERROR: 'rule busco_proteins:' not found", file=sys.stderr)
    sys.exit(1)
nxt = src.find("\nrule ", start + 1)
end = nxt if nxt >= 0 else len(src)
block = src[start:end]

# --- A: filter >50k-aa proteins right after the proteins dir is cleared -----
rm_anchor = 'rm -rf "$OUTDIR_ABS/proteins"\n'
if rm_anchor not in block:
    print("ERROR: 'rm -rf $OUTDIR_ABS/proteins' anchor not found", file=sys.stderr)
    sys.exit(1)
filter_block = (
    "\n"
    "        # --- berghia-chemogpcrs guard: drop runaway gene models (>50k aa) ---\n"
    "        # A single >50,000-aa \"protein\" (a fused/runaway model) makes BUSCO's\n"
    "        # hmmsearch abort with SIGABRT and fails the whole genome. Real BUSCO\n"
    "        # orthologs are never that long, so removing them does not change the\n"
    "        # completeness score; BUSCO runs on the filtered set instead.\n"
    "        FILTERED_AA=\"$OUTDIR_ABS/proteins_busco_input.faa\"\n"
    "        awk -v MAXL=50000 '/^>/ {{ if (n && l<=MAXL) print s; n=1; s=$0; l=0; next }} "
    "{{ s=s ORS $0; l+=length($0) }} END {{ if (n && l<=MAXL) print s }}' "
    "{input.proteins} > \"$FILTERED_AA\"\n"
    "        NDROP=$(( $(grep -c \"^>\" {input.proteins}) - $(grep -c \"^>\" \"$FILTERED_AA\") ))\n"
    "        echo \"[INFO] busco_proteins: dropped $NDROP protein(s) >50000 aa before BUSCO "
    "(runaway models, not orthologs)\" >> {log}\n"
)
block = block.replace(rm_anchor, rm_anchor + filter_block, 1)

# Point BUSCO at the filtered set.
if "-i {input.proteins} \\" not in block:
    print("ERROR: '-i {input.proteins}' anchor not found", file=sys.stderr)
    sys.exit(1)
block = block.replace("-i {input.proteins} \\", '-i "$FILTERED_AA" \\', 1)

# --- B: make the BUSCO call non-fatal (QC metric must not block annotation) -
redirect = ">> {log} 2>&1\n"
if redirect not in block:
    print("ERROR: busco '>> {log} 2>&1' redirect anchor not found", file=sys.stderr)
    sys.exit(1)
block = block.replace(
    redirect,
    '>> {log} 2>&1 || echo "[WARN] busco_proteins: BUSCO exited non-zero even after '
    '>50k filter; continuing (QC metric only, not blocking annotation)" >> {log}\n',
    1,
)

path.write_text(src[:start] + block + src[end:])
print("Patched run_busco.smk: busco_proteins now filters >50k-aa proteins + is non-fatal.")
