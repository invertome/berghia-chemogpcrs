"""Stage 06's MCScanX branch looked for `.gff`, never `.gff3` (audit finding F12.1).

``06:158`` read ``${GENOME_DIR}/${taxid}.gff`` while the JCVI branch in the same
script (``06:57``) reads ``${GENOME_DIR}/${tgt_base}.gff3``. Nothing in the
repository writes a bare ``.gff`` into ``${GENOME_DIR}``:
``scripts/fetch_berghia_genome.sh:32`` writes ``${PREFIX}.gff3`` (and its legacy
alias at :39 is ``.gff3`` too), and every annotation file actually present in
``genomes/`` is ``.gff3``. So the MCScanX branch's
``if [ -f "$gff_file" ]`` was never true and no ``.gff``-derived MCScanX
positional input was ever produced.

Low blast radius -- ``SYNTENY_BACKEND`` defaults to ``jcvi`` and the JCVI branch
was always correct -- but the fallback backend was silently degraded.

The fix probes ``.gff3`` first and keeps ``.gff`` as a second choice, so it is
strictly additive: any hand-staged ``.gff`` that happened to work keeps working.
These tests execute the REAL resolution lifted out of the shipping script.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
STAGE06 = PROJECT_ROOT / "06_synteny_and_mapping.sh"
FETCH_GENOME = PROJECT_ROOT / "scripts" / "fetch_berghia_genome.sh"


def _extract_gff_resolution() -> str:
    """Lift the MCScanX branch's gff_file resolution out of stage 06."""
    lines = STAGE06.read_text().splitlines()
    idxs = [i for i, l in enumerate(lines) if l.lstrip().startswith("gff_file=")]
    assert idxs, "stage 06's MCScanX branch no longer assigns gff_file"
    start = idxs[0]

    # The resolution is either a bare assignment or an assignment followed by a
    # probe loop. Walk forward while block constructs stay open, so the snippet
    # we execute is always syntactically complete.
    depth = 0
    end = start
    for j in range(start, min(start + 20, len(lines))):
        stripped = lines[j].strip()
        depth += len(re.findall(r"\bfor\b|\bwhile\b", stripped))
        depth -= len(re.findall(r"\bdone\b", stripped))
        end = j
        if depth <= 0 and j > start and not stripped.endswith("\\"):
            # stop once we are back at top level and the next line starts
            # something unrelated to gff_file resolution
            nxt = lines[j + 1].strip() if j + 1 < len(lines) else ""
            if stripped.startswith("done") or not nxt or "gff" not in stripped:
                break
    return "\n".join(l.strip() for l in lines[start : end + 1] if l.strip())


def _resolve(tmp_path: Path, filenames: list[str], taxid: str = "6500") -> str:
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir(parents=True, exist_ok=True)
    for name in filenames:
        (genome_dir / name).write_text("##gff-version 3\n")

    snippet = "\n".join(
        [
            f'GENOME_DIR="{genome_dir}"',
            f'taxid="{taxid}"',
            'gff_file=""',
            _extract_gff_resolution(),
            'echo "GFF=${gff_file}"',
        ]
    )
    proc = subprocess.run(["bash", "-c", snippet], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    for line in proc.stdout.splitlines():
        if line.startswith("GFF="):
            return line.split("=", 1)[1]
    pytest.fail("no GFF= line: " + proc.stdout + proc.stderr)


# --------------------------------------------------------------------------
# the regression
# --------------------------------------------------------------------------
def test_gff3_annotation_is_found(tmp_path):
    """The core bug: only .gff3 exists on disk, and it must resolve."""
    resolved = _resolve(tmp_path, ["6500.gff3"])

    assert resolved.endswith("6500.gff3"), (
        "the MCScanX branch must find the .gff3 files the pipeline actually writes"
    )


def test_bare_gff_is_still_honoured(tmp_path):
    """Additive fix: a hand-staged .gff keeps working."""
    resolved = _resolve(tmp_path, ["6500.gff"])

    assert resolved.endswith("6500.gff")


def test_gff3_wins_when_both_are_present(tmp_path):
    """.gff3 is the pipeline's own convention, so it takes priority."""
    resolved = _resolve(tmp_path, ["6500.gff", "6500.gff3"])

    assert resolved.endswith("6500.gff3")


def test_no_annotation_resolves_to_empty(tmp_path):
    """The downstream `if [ -f "$gff_file" ]` guard must still see a falsy value."""
    assert _resolve(tmp_path, ["6500.fasta"]) == ""


def test_another_taxids_annotation_is_not_picked_up(tmp_path):
    """Resolution is keyed on the taxid, not on whatever gff happens to be there."""
    assert _resolve(tmp_path, ["9606.gff3"], taxid="6500") == ""


# --------------------------------------------------------------------------
# the two branches of stage 06 must agree, and match the producer
# --------------------------------------------------------------------------
def test_the_pipeline_only_ever_writes_gff3():
    """Pin the premise: the producer writes .gff3, so readers must accept it."""
    src = FETCH_GENOME.read_text()
    assert '.gff3"' in src, "fetch_berghia_genome.sh no longer writes .gff3"


def test_jcvi_branch_still_reads_gff3():
    """Regression guard on the branch that was already correct."""
    assert '${GENOME_DIR}/${tgt_base}.gff3' in STAGE06.read_text()


def test_mcscanx_branch_no_longer_probes_only_bare_gff():
    block = _extract_gff_resolution()
    assert ".gff3" in block, "the MCScanX branch must probe .gff3"
