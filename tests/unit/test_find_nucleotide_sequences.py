"""Tests for find_nucleotide_sequences in functions.sh (bead sz6).

The per-OG codon-alignment step looks up each sequence's nucleotides. Target/
Berghia transcriptomes follow the config convention
``<taxid>_<genus>_<species>.<ext>`` (e.g. 1287507_Berghia_stephanieae.mrna), but
the transcriptome search only tried ``<taxid>.<ext>``, ``taxid_<taxid>.<ext>``,
and ``<taxid>_<taxid>.<ext>`` — none of which match. So the focal-species CDS was
silently absent from the codon alignment and dN/dS for the branch of interest was
lost. Production resolves the taxid via id_map.csv; these tests exercise the
header-parse fallback (seq_id ``<taxid>_<rest>``) to isolate the filename match.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"


def _run_lookup(tmp_path: Path, transcriptome_name: str) -> tuple[subprocess.CompletedProcess, Path]:
    tdir = tmp_path / "transcriptomes"
    tdir.mkdir()
    (tdir / transcriptome_name).write_text(">1287507_myseq desc\nATGAAATTTTAG\n")
    protein = tmp_path / "prot.fa"
    protein.write_text(">1287507_myseq\nMKF*\n")
    out = tmp_path / "nuc.fa"
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    env["TRANSCRIPTOME_DIR"] = str(tdir)
    env["ID_MAP"] = str(tmp_path / "no_such_id_map.csv")   # force header-parse fallback
    env["REFERENCE_CDS_FILE"] = str(tmp_path / "no_such_refs.fna")
    script = f'source "{FUNCTIONS_SH}"; find_nucleotide_sequences "{protein}" "{out}"'
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)
    return result, out


def test_finds_taxid_genus_species_named_transcriptome(tmp_path):
    # 1287507_<genus>_<species>.mrna — the real Berghia/target convention.
    result, out = _run_lookup(tmp_path, "1287507_testgenus_sp.mrna")
    assert result.returncode == 0, result.stderr
    text = out.read_text()
    assert ">1287507_myseq" in text
    assert "ATGAAATTTTAG" in text


def test_still_finds_bare_taxid_named_transcriptome(tmp_path):
    # Regression: the pre-existing <taxid>.<ext> pattern still works.
    result, out = _run_lookup(tmp_path, "1287507.mrna")
    assert result.returncode == 0, result.stderr
    assert ">1287507_myseq" in out.read_text()
