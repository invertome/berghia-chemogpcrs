#!/usr/bin/env python3
"""patch_braker4_run_genemark_ep.py — make GeneMark-EP non-fatal with a
GeneMark-ES fallback in the vendored BRAKER4 workflow.

GeneMark-EP (`gmes_petap.pl --EP`) hard-crashes (parse_ET.pl div-by-zero) when
protein evidence yields ~0 high-confidence introns — which happens on deeply
isolated lineages (e.g. Nautilus, ~450 My from its nearest proteome relative).
That abort throws away the GeneMark-ES result that already ran for the ProtHint
seeds. This patch wraps the gmes invocation so that on EP failure (or empty
output) it falls back to `GeneMark-ES/genemark.gtf` (which the EP rule only needs
to emit as `GeneMark-EP/genemark.gtf`), letting AUGUSTUS/TSEBRA proceed on the
ab-initio genes. EP is used where it works; ES carries the rest.

Idempotent — safe to re-apply after any `git pull` in external/braker4/.
Mirrors scripts/unity/patch_braker4_run_compleasm.py.

Usage:  python3 scripts/unity/patch_braker4_run_genemark_ep.py [<braker4_dir>]
"""
import sys
from pathlib import Path

SENTINEL = "falling back to GeneMark-ES ab-initio predictions"

OLD = '''        eval $CMD 1>> $LOG_FILE_ABS 2>> $LOG_ABS

        if [ ! -f genemark.gtf ]; then
            echo "ERROR: GeneMark-EP failed to produce genemark.gtf" >> $LOG_ABS
            exit 1
        fi
'''

NEW = '''        # EP can hard-crash (parse_ET.pl div-by-zero) on near-zero intron
        # evidence (deeply isolated lineages). Capture the exit code instead of
        # letting set -e abort, then fall back to the ab-initio GeneMark-ES genes
        # (already computed for the ProtHint seeds) so the run is not lost.
        set +e
        eval $CMD 1>> $LOG_FILE_ABS 2>> $LOG_ABS
        EP_EXIT=$?
        set -e

        if [ $EP_EXIT -ne 0 ] || [ ! -f genemark.gtf ]; then
            ES_GTF="$WORKDIR/output/{wildcards.sample}/GeneMark-ES/genemark.gtf"
            if [ -s "$ES_GTF" ]; then
                echo "WARNING: GeneMark-EP failed (exit=$EP_EXIT) or produced no genemark.gtf; falling back to GeneMark-ES ab-initio predictions ($ES_GTF). Protein evidence was too sparse/divergent for EP on this genome." >> $LOG_ABS
                cp "$ES_GTF" genemark.gtf
            else
                echo "ERROR: GeneMark-EP failed (exit=$EP_EXIT) and no GeneMark-ES fallback at $ES_GTF" >> $LOG_ABS
                exit 1
            fi
        fi
'''


def main() -> int:
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("external/braker4")
    rule = base / "rules" / "genemark" / "run_genemark_ep.smk"
    if not rule.is_file():
        print(f"ERROR: rule not found: {rule}", file=sys.stderr)
        return 1
    text = rule.read_text()
    if SENTINEL in text:
        print(f"[patch_genemark_ep] already patched: {rule}")
        return 0
    if OLD not in text:
        print(f"ERROR: expected EP eval/failure block not found in {rule}; "
              f"upstream may have changed — re-inspect before patching.", file=sys.stderr)
        return 2
    rule.write_text(text.replace(OLD, NEW, 1))
    print(f"[patch_genemark_ep] PATCHED {rule} — EP failure now falls back to GeneMark-ES")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
