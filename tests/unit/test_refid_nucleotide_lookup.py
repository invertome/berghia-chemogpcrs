"""Tests pinning the ID-boundary contract in ``find_nucleotide_sequences``.

The defect
----------
``functions.sh`` gated the CDS lookup on an UNANCHORED grep::

    if grep -q "^>${seq_id}" "$REFERENCE_CDS_FILE"; then

while the extraction immediately below anchors on a record boundary::

    match($0, "^>" id "($|[[:space:]])")

Short IDs are ``ref_<taxid>_<N>`` with an unpadded counter, so ``ref_phau_1``
is a strict prefix of ``ref_phau_12``. With ``ref_phau_12`` present and
``ref_phau_1`` absent, the gate MATCHES, the anchored extraction emits
NOTHING, and the function still counts the sequence as found -- writing zero
bases under a header that downstream codon alignment expects to be real. The
gate and the extraction must agree on where an ID ends.

The same unanchored gate appeared at all four lookup sites (reference CDS,
the two transcriptome branches, and the non-reference CDS fallback); the
transcriptome branch is worse still, because its ``protein_id`` fallback
RENAMES whatever it matched to ``seq_id``, so a prefix collision would relabel
another gene's nucleotides as this one's.

These tests drive the real function out of ``functions.sh`` with the pipeline
helpers stubbed, so they pin observable behavior rather than the text of the
guard.
"""
from __future__ import annotations

import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = REPO_ROOT / "functions.sh"


HARNESS = """
set -uo pipefail
source "{functions}"

# --- stub the pipeline helpers find_nucleotide_sequences depends on ---
log() {{ echo "LOG: $*"; }}
is_reference_seq() {{ [[ "$1" == ref_* ]]; }}
# ref_<taxid>_<N> keeps the taxid in field 2; <taxid>_<accession> in field 1.
get_taxid_for_seq() {{
    if [[ "$1" == ref_* ]]; then echo "$1" | cut -d'_' -f2
    else echo "$1" | cut -d'_' -f1; fi
}}

REFERENCE_CDS_FILE="{cds}"
TRANSCRIPTOME_DIR="{transcriptomes}"
RESULTS_DIR="{results}"

find_nucleotide_sequences "{proteins}" "{out}"
echo "EXIT:$?"
"""


def run_lookup(tmp_path, cds_records, protein_ids, transcriptome_records=None):
    """Run the real ``find_nucleotide_sequences`` against fixture files."""
    cds = tmp_path / "all_references_cds.fna"
    cds.write_text("".join(f">{h}\n{s}\n" for h, s in cds_records))

    proteins = tmp_path / "proteins.fa"
    proteins.write_text("".join(f">{i}\nMAAA\n" for i in protein_ids))

    tdir = tmp_path / "transcriptomes"
    tdir.mkdir()
    if transcriptome_records:
        for fname, records in transcriptome_records.items():
            (tdir / fname).write_text(
                "".join(f">{h}\n{s}\n" for h, s in records)
            )

    out = tmp_path / "nucleotides.fa"
    script = tmp_path / "harness.sh"
    script.write_text(HARNESS.format(
        functions=FUNCTIONS_SH, cds=cds, transcriptomes=tdir,
        proteins=proteins, out=out, results=tmp_path,
    ))
    proc = subprocess.run(["bash", str(script)], capture_output=True, text=True)
    return proc, out


pytestmark = pytest.mark.skipif(
    shutil.which("bash") is None, reason="bash not available"
)


def test_prefix_of_a_longer_id_is_not_treated_as_a_match(tmp_path):
    """``ref_phau_12`` present + ``ref_phau_1`` absent must report MISSING.

    The unanchored gate reported this as found and wrote zero bases.
    """
    proc, out = run_lookup(
        tmp_path,
        cds_records=[("ref_phau_12", "ATGCCCAAA")],
        protein_ids=["ref_phau_1"],
    )
    combined = proc.stdout + proc.stderr
    assert "No nucleotide sequence found for ref_phau_1" in combined, combined
    assert "missing=1" in combined, combined
    assert "ATGCCCAAA" not in out.read_text(), (
        "returned a different record's nucleotides"
    )


def test_every_sequence_counted_as_found_has_bases_written(tmp_path):
    """The found tally must equal the records actually emitted.

    Mixed batch: ``ref_phau_2`` resolves for real, ``ref_phau_1`` only
    prefix-collides with ``ref_phau_12``. The unanchored gate counted both,
    so the tally claimed 2 while only 1 record reached the file.
    """
    proc, out = run_lookup(
        tmp_path,
        cds_records=[("ref_phau_2", "ATGAAAGGG"), ("ref_phau_12", "ATGCCCAAA")],
        protein_ids=["ref_phau_1", "ref_phau_2"],
    )
    combined = proc.stdout + proc.stderr
    emitted = [ln for ln in out.read_text().splitlines() if ln.startswith(">")]
    assert "found=1" in combined and "missing=1" in combined, combined
    assert emitted == [">ref_phau_2"], emitted


def test_exact_id_still_resolves(tmp_path):
    """Positive control: the anchored gate must not break real lookups."""
    proc, out = run_lookup(
        tmp_path,
        cds_records=[("ref_phau_1", "ATGAAAGGG"), ("ref_phau_12", "ATGCCCAAA")],
        protein_ids=["ref_phau_1"],
    )
    combined = proc.stdout + proc.stderr
    text = out.read_text()
    assert ">ref_phau_1\n" in text, text
    assert "ATGAAAGGG" in text, text
    assert "ATGCCCAAA" not in text, "matched the wrong record too"
    assert "found=1" in combined, combined


def test_id_followed_by_description_still_resolves(tmp_path):
    """Record boundary is end-of-ID or whitespace, matching the extractor."""
    proc, out = run_lookup(
        tmp_path,
        cds_records=[("ref_phau_1 some description here", "ATGAAAGGG")],
        protein_ids=["ref_phau_1"],
    )
    assert "ATGAAAGGG" in out.read_text(), out.read_text()


def test_transcriptome_branch_does_not_rename_a_prefix_match(tmp_path):
    """The transcriptome fallback RENAMES what it matches -- worst case.

    ``phau`` is not a reference prefix here, so lookup falls to the
    transcriptome branch. With only ``ACC12`` present, a prefix match on
    ``ACC1`` would relabel ACC12's nucleotides as ``phau_ACC1``.
    """
    proc, out = run_lookup(
        tmp_path,
        cds_records=[],
        protein_ids=["phau_ACC1"],
        transcriptome_records={"phau.cds": [("ACC12", "ATGTTTCCC")]},
    )
    combined = proc.stdout + proc.stderr
    text = out.read_text()
    assert "ATGTTTCCC" not in text, "relabelled another record's nucleotides"
    assert ">phau_ACC1" not in text, text
    assert "No nucleotide sequence found for phau_ACC1" in combined, combined
    assert "missing=1" in combined, combined
