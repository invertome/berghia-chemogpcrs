"""Tests for scripts/fetch_reference_cds.py.

Bug context (bead -lfy / 2026-05-05):
    For LSE annelid species, protein headers in Nath et al. references look
    like ``>ampa_OX439002.1_2091693``. The accession parser extracts
    ``OX439002.1`` (a Wellcome Sanger Tree-of-Life genome contig). The
    classifier matches ``[A-Z]{1,2}\\d{5,6}\\.\\d+`` and labels it
    ``ncbi_genbank_nuc``; dispatch routes that to ``fetch_cds_ncbi_tsa``,
    which calls ``efetch nuccore <acc> rettype=fasta`` and gets back the
    ENTIRE 50 MB genome contig. The script then writes that whole contig as
    the "CDS" for the protein. Repeated for ~1500 LSE proteins, this
    produces 60 GB of garbage in ``all_references_cds.fna``.

These tests pin the contract going forward:
    1. Genomic-contig accessions (Sanger ToL prefixes OX/LR/CAJ/CAK,
       short-prefix-with-version chromosomes) are classified separately and
       NOT routed through the "fetch the whole record" path.
    2. ``fetch_cds_ncbi_tsa`` enforces a ``MAX_CDS_LENGTH_BP`` cap. Anything
       larger than the cap is rejected as a fetched-the-wrong-thing failure.
    3. Successful small-CDS fetches still work.
"""
from __future__ import annotations

import io
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# conftest.py adds scripts/ to sys.path
import fetch_reference_cds as frc


class FakeHandle:
    """Minimal stand-in for a Biopython Entrez handle."""
    def __init__(self, payload: str) -> None:
        self._payload = payload

    def read(self) -> str:
        return self._payload

    def close(self) -> None:
        pass


def _make_fasta(header: str, seq: str, line_width: int = 60) -> str:
    body = "\n".join(seq[i:i + line_width] for i in range(0, len(seq), line_width))
    return f">{header}\n{body}\n"


# -------- classify_accession: Sanger ToL prefixes get a dedicated class -----

@pytest.mark.parametrize("acc", ["OX439002.1", "LR761643.1"])
def test_sanger_tol_accessions_classified_as_genome_contig(acc: str) -> None:
    """OX / LR prefix accessions are Wellcome Sanger Tree-of-Life genome
    contigs / chromosomes — unambiguously genomic. The fix routes them to a
    dedicated ``ncbi_genome_contig`` class so the dispatcher does NOT call
    ``efetch fasta`` on them (which would return a 50 MB+ whole contig)."""
    assert frc.classify_accession(acc) == "ncbi_genome_contig"


def test_real_tsa_accession_still_classified_as_tsa() -> None:
    """Sanity check: actual TSA contigs are still routed through the TSA path.
    Multi-letter prefix patterns (CAJxxx, JABxxx, GABxxx) overlap between TSA
    and WGS at the accession-format level; the script can't reliably
    distinguish them from string alone, so it accepts them as TSA and relies
    on the MAX_CDS_LENGTH_BP cap to reject any whole-contig fetches."""
    assert frc.classify_accession("GABX01000123") == "ncbi_tsa"
    assert frc.classify_accession("JABMCL010001063.1") == "ncbi_tsa"


def test_refseq_protein_accession_still_classified_correctly() -> None:
    """Don't break the working RefSeq path."""
    assert frc.classify_accession("XP_005089002.1") == "ncbi_refseq_protein"
    assert frc.classify_accession("NP_001234.1") == "ncbi_refseq_protein"


# -------- fetch_cds_ncbi_tsa: size cap on the returned sequence -------------

def test_fetch_cds_ncbi_tsa_rejects_oversized_response() -> None:
    """If the efetch response is > MAX_CDS_LENGTH_BP (50 kb), reject it.

    This is the safety net: if anything else slips through similar paths,
    we don't poison the output FASTA with whole-genome contigs."""
    huge_seq = "CCCTAA" * 200_000  # 1.2 Mbp of telomere repeat
    payload = _make_fasta("OX439002.1", huge_seq)
    fake_entrez = MagicMock()
    fake_entrez.efetch.return_value = FakeHandle(payload)

    result = frc.fetch_cds_ncbi_tsa("OX439002.1", fake_entrez)
    assert result is None, "oversized response must be rejected"


def test_fetch_cds_ncbi_tsa_accepts_normal_cds() -> None:
    """A normal-sized response (1.5 kb single CDS) is accepted."""
    normal_cds = "ATG" + "GCT" * 499 + "TAA"  # 1.5 kb, valid frame
    payload = _make_fasta("GABXXX01000001.1", normal_cds)
    fake_entrez = MagicMock()
    fake_entrez.efetch.return_value = FakeHandle(payload)

    result = frc.fetch_cds_ncbi_tsa("GABXXX01000001.1", fake_entrez)
    assert result is not None
    acc, seq = result
    assert acc == "GABXXX01000001.1"
    assert seq == normal_cds


def test_fetch_cds_ncbi_tsa_size_cap_constant_is_reasonable() -> None:
    """The cap should be high enough to admit the largest plausible single
    GPCR CDS (the pipeline's actual targets) but well below any genome
    contig. 50 kb is a reasonable middle ground."""
    assert hasattr(frc, "MAX_CDS_LENGTH_BP"), "fix must define MAX_CDS_LENGTH_BP"
    cap = frc.MAX_CDS_LENGTH_BP
    assert 10_000 <= cap <= 200_000, f"MAX_CDS_LENGTH_BP={cap} outside sane range"


# -------- fetch_cds_for_accession: dispatcher does not call TSA for ToL ----

def test_dispatcher_skips_efetch_for_genome_contig_class() -> None:
    """For an accession classified as ncbi_genome_contig, the dispatcher
    must NOT call fetch_cds_ncbi_tsa or fetch_cds_ena. The fetch should be
    marked skipped (no efetch path can produce a real CDS for these — they
    need miniprot recovery against the source assembly)."""
    info = frc.AccessionInfo(
        accession="OX439002.1",
        database=frc.classify_accession("OX439002.1"),
        protein_id="ampa_OX439002.1_2091693",
        original_header=">ampa_OX439002.1_2091693",
    )
    fake_entrez = MagicMock()
    # Any efetch call would return a huge bogus payload; the test asserts
    # the dispatcher does NOT call efetch at all for this class.
    fake_entrez.efetch.return_value = FakeHandle(
        _make_fasta("OX439002.1", "CCCTAA" * 100_000))

    result = frc.fetch_cds_for_accession(info, fake_entrez)
    assert result.status == "skipped", (
        f"expected status=skipped for genome-contig accession, got {result.status!r}"
        f" (cds_sequence length={len(result.cds_sequence) if result.cds_sequence else 0})")
    assert result.cds_sequence is None
    fake_entrez.efetch.assert_not_called()


def test_dispatcher_still_works_for_tsa() -> None:
    """The dispatch for a real TSA accession must continue to call efetch
    and return a CDS (regression guard for the legitimate path)."""
    cds = "ATG" + "AAC" * 199 + "TAA"  # 600 bp
    info = frc.AccessionInfo(
        accession="GABX01000123",
        database=frc.classify_accession("GABX01000123"),
        protein_id="alvmar_GABX01000123_42",
        original_header=">alvmar_GABX01000123_42",
    )
    fake_entrez = MagicMock()
    fake_entrez.efetch.return_value = FakeHandle(_make_fasta("GABX01000123", cds))

    result = frc.fetch_cds_for_accession(info, fake_entrez)
    assert result.status == "success"
    assert result.cds_sequence == cds
