"""Regression tests for scripts/unity/patch_braker4_run_busco.py (bead pfo).

The patch idempotently rewrites the `busco_proteins` rule in run_busco.smk to
(A) filter >50,000-aa runaway proteins out of the BUSCO input and (B) make the
BUSCO call non-fatal — scoped to the busco_proteins block so busco_genome is
untouched. The script does its work at import time (no importable function), so
these tests run it as a subprocess against a fixture run_busco.smk, exactly as
it runs on Unity.

Coverage:
  1. the three scoped edits land (filter block, -i swap, non-fatal redirect);
  2. the sibling busco_genome rule is byte-for-byte untouched;
  3. re-running is a no-op (idempotent);
  4. a missing busco_proteins rule fails loudly (exit 1);
  5. the embedded >50k awk filter actually drops runaway proteins.
"""
from __future__ import annotations

import shlex
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PATCH = PROJECT_ROOT / "scripts" / "unity" / "patch_braker4_run_busco.py"

MARKER = "berghia-chemogpcrs guard: drop runaway"
# The exact awk program embedded by the patch (Snakemake-doubled braces).
AWK_DOUBLED = (
    "awk -v MAXL=50000 '/^>/ {{ if (n && l<=MAXL) print s; n=1; s=$0; l=0; next }} "
    "{{ s=s ORS $0; l+=length($0) }} END {{ if (n && l<=MAXL) print s }}'"
)

# Minimal run_busco.smk carrying both rules + the exact anchors the patch keys on.
FIXTURE_SMK = '''rule busco_genome:
    input:
        genome="output/{sample}/genome.fa"
    output:
        done="output/{sample}/busco/genome/.done"
    shell:
        r"""
        set -euo pipefail
        OUTDIR_ABS=$(readlink -f {params.outdir})
        rm -rf "$OUTDIR_ABS/genome"
        busco \\
            -i {input.genome} \\
            -m genome \\
            >> {log} 2>&1
        touch {output.done}
        """


rule busco_proteins:
    input:
        proteins="output/{sample}/braker.longest.aa"
    output:
        done="output/{sample}/busco/proteins/.done"
    shell:
        r"""
        set -euo pipefail
        OUTDIR_ABS=$(readlink -f {params.outdir})
        rm -rf "$OUTDIR_ABS/proteins"

        mkdir -p {params.download_path}
        OFFLINE_FLAG=""
        busco \\
            -i {input.proteins} \\
            -o proteins \\
            -m proteins \\
            $OFFLINE_FLAG \\
            >> {log} 2>&1

        touch {output.done}
        """


rule busco_summary:
    shell:
        "echo done"
'''


def _setup(tmp_path, smk=FIXTURE_SMK):
    (tmp_path / "run_busco.smk").write_text(smk)
    (tmp_path / "patch_braker4_run_busco.py").write_text(PATCH.read_text())


def _run(tmp_path):
    return subprocess.run(
        [sys.executable, str(tmp_path / "patch_braker4_run_busco.py")],
        cwd=str(tmp_path), capture_output=True, text=True,
    )


def test_applies_the_three_scoped_edits(tmp_path):
    _setup(tmp_path)
    res = _run(tmp_path)
    assert res.returncode == 0, res.stderr
    assert "Patched run_busco.smk" in res.stdout
    out = (tmp_path / "run_busco.smk").read_text()

    assert MARKER in out
    # (A) filter block + the exact awk program
    assert 'FILTERED_AA="$OUTDIR_ABS/proteins_busco_input.faa"' in out
    assert AWK_DOUBLED in out
    # BUSCO now reads the filtered set; the raw -i is gone
    assert '-i "$FILTERED_AA"' in out
    assert "-i {input.proteins}" not in out
    # (B) the BUSCO call is now non-fatal
    assert '>> {log} 2>&1 || echo "[WARN] busco_proteins:' in out


def test_busco_genome_rule_untouched(tmp_path):
    _setup(tmp_path)
    _run(tmp_path)
    out = (tmp_path / "run_busco.smk").read_text()
    sep = "rule busco_proteins:"
    # Everything before the busco_proteins rule must be byte-identical.
    assert out.split(sep)[0] == FIXTURE_SMK.split(sep)[0]
    # busco_genome's own -i and redirect are preserved unmodified.
    assert "-i {input.genome}" in out


def test_idempotent_second_run_is_noop(tmp_path):
    _setup(tmp_path)
    _run(tmp_path)
    after_first = (tmp_path / "run_busco.smk").read_text()
    res2 = _run(tmp_path)
    assert res2.returncode == 0
    assert "Already patched — no-op." in res2.stdout
    assert (tmp_path / "run_busco.smk").read_text() == after_first


def test_missing_busco_proteins_rule_fails_loudly(tmp_path):
    smk = "rule busco_genome:\n    shell:\n        'echo hi'\n"
    _setup(tmp_path, smk=smk)
    res = _run(tmp_path)
    assert res.returncode == 1
    assert "rule busco_proteins:" in res.stderr


def test_embedded_awk_drops_runaway_proteins(tmp_path):
    # The awk asserted present above, run with Snakemake braces un-doubled.
    _setup(tmp_path)
    _run(tmp_path)
    assert AWK_DOUBLED in (tmp_path / "run_busco.smk").read_text()
    runnable = AWK_DOUBLED.replace("{{", "{").replace("}}", "}")

    fasta = tmp_path / "in.faa"
    fasta.write_text(">small1\nMMMMM\n>runaway\n" + "A" * 60000
                     + "\n>small2\nKKLLKK\n>multi\nAAAA\nCCCC\n")
    out = tmp_path / "out.faa"
    res = subprocess.run(
        ["bash", "-c", f"{runnable} {shlex.quote(str(fasta))} > {shlex.quote(str(out))}"],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr
    text = out.read_text()
    assert text.count(">") == 3            # small1, small2, multi kept
    assert ">runaway" not in text          # the 60k-aa seq dropped
    assert "AAAA\nCCCC" in text            # multi-line sequence preserved
