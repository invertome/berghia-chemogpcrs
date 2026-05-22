#!/usr/bin/env python3
"""Idempotent patch: append 'exit 0' to the shell block end in run_compleasm.smk
so post-compleasm cleanup (version stamping, citations) can't cascade a
strict-mode failure. Belt-and-braces alongside pre-downloading the BUSCO
metazoa_odb12 lineage (bead -p49).
"""
import re
import sys
from pathlib import Path

path = Path(__file__).parent / "run_compleasm.smk"
if not path.exists():
    path = Path("run_compleasm.smk")

src = path.read_text()

if "# exit 0 added by berghia-chemogpcrs" in src:
    print("Already patched — no-op.")
    sys.exit(0)

marker = 'cite miniprot "$REPORT_DIR"'
closing = '"""'

idx = src.rfind(marker)
if idx < 0:
    print("ERROR: 'cite miniprot' line not found", file=sys.stderr)
    sys.exit(1)

# Find the closing """ after the marker
close_idx = src.find(closing, idx)
if close_idx < 0:
    print("ERROR: closing triple-quote not found after marker", file=sys.stderr)
    sys.exit(1)

# Identify indentation of the closing """
line_start = src.rfind("\n", 0, close_idx) + 1
indent = src[line_start:close_idx]

patch_block = (
    f"\n"
    f"{indent}# exit 0 added by berghia-chemogpcrs (bead -p49): post-compleasm\n"
    f"{indent}# cleanup is best-effort. The rule's actual outputs\n"
    f"{indent}# (compleasm_hints, compleasm_log, compleasm_summary) are independently\n"
    f"{indent}# `touch`-ed above, so a strict-mode failure in version-stamping or\n"
    f"{indent}# citation reporting must not abort the rule and cascade --keep-going.\n"
    f"{indent}exit 0\n"
)

new_src = src[:close_idx] + patch_block + src[close_idx:]
path.write_text(new_src)
print(f"Patched: inserted `exit 0` before closing triple-quote at byte {close_idx}.")
