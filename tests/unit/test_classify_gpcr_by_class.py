"""Tests for scripts/classify_gpcr_by_class.py.

P1 of the per-class refactor — sequence-level GPCR class classifier.
Tests cover pure-logic functions only (no real hmmscan / network calls).
Subprocess calls are either monkeypatched or exercised via synthetic
tblout text fed to the parser directly.
"""
from __future__ import annotations

import csv
import gzip
import io
import sys
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest

# Ensure scripts/ is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import classify_gpcr_by_class as cgc


# ---------------------------------------------------------------------------
# 1. Class A from PF00001
# ---------------------------------------------------------------------------

def test_call_class_A_from_PF00001():
    """Best Pfam hit PF00001 -> class A."""
    hits = [("PF00001", "seq1", 1e-50)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"
    assert evidence_pfam == "PF00001"
    assert top_evalue == 1e-50


# ---------------------------------------------------------------------------
# 2. Class B from PF00002
# ---------------------------------------------------------------------------

def test_call_class_B_from_PF00002():
    hits = [("PF00002", "seq1", 1e-40)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "B"
    assert evidence_pfam == "PF00002"


# ---------------------------------------------------------------------------
# 3. Class C from PF00003
# ---------------------------------------------------------------------------

def test_call_class_C_from_PF00003():
    hits = [("PF00003", "seq1", 1e-35)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "C"
    assert evidence_pfam == "PF00003"


# ---------------------------------------------------------------------------
# 4. Class F from PF01534
# ---------------------------------------------------------------------------

def test_call_class_F_from_PF01534():
    hits = [("PF01534", "seq1", 5e-20)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "F"
    assert evidence_pfam == "PF01534"


# ---------------------------------------------------------------------------
# 5. Class A from PF02949 (7tm_4, vertebrate ORs — newly bootstrapped)
# ---------------------------------------------------------------------------

def test_call_class_A_from_PF02949():
    hits = [("PF02949", "seq1", 3e-30)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"
    assert evidence_pfam == "PF02949"


# ---------------------------------------------------------------------------
# 6. PF02949 -> class=A AND evidence_family_hmm=insect_OR_atypical
# ---------------------------------------------------------------------------

def test_PF02949_sets_class_A_and_insect_subfamily():
    """PF02949 is the REAL insect odorant-receptor family (InterPro short name
    `7tm_6`), so the insect_OR_atypical subfamily belongs on it.

    It was previously on PF10324, which InterPro reports as `7TM_GPCR_Srw`
    (a *C. elegans* serpentine chemoreceptor). That mislabelling routed real
    chemoreceptors OUT of the class-A tree."""
    hits = [("PF02949", "seq1", 1e-60)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"
    assert evidence_pfam == "PF02949"

    subfam = cgc.refine_subfamily([], pfam_accession="PF02949")
    assert subfam == "insect_OR_atypical"


def test_PF10324_is_srw_and_carries_no_insect_subfamily():
    """Regression guard: PF10324 must never regain the insect-OR label."""
    assert cgc._PFAM_VERIFIED_NAME["PF10324"] == "7TM_GPCR_Srw"
    assert cgc.refine_subfamily([], pfam_accession="PF10324") != "insect_OR_atypical"


# ---------------------------------------------------------------------------
# 7. Multi-hit: best E-value wins
# ---------------------------------------------------------------------------

def test_call_class_best_evalue_wins():
    """Multiple hits: lowest E-value determines the class call."""
    hits = [("PF00001", "seq1", 1e-50), ("PF00002", "seq1", 1e-10)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"
    assert evidence_pfam == "PF00001"
    assert top_evalue == 1e-50


# ---------------------------------------------------------------------------
# 8. No hit at E<threshold -> unclassified
# ---------------------------------------------------------------------------

def test_call_class_unclassified_when_no_passing_hit():
    hits = [("PF00001", "seq1", 1e-3)]  # above 1e-5 threshold
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "unclassified"
    assert evidence_pfam == ""
    assert top_evalue == ""


def test_call_class_unclassified_when_empty_hits():
    cls, evidence_pfam, top_evalue = cgc.call_class([], cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "unclassified"
    assert evidence_pfam == ""
    assert top_evalue == ""


# ---------------------------------------------------------------------------
# 9. 06c family refines subfamily but doesn't change class
# ---------------------------------------------------------------------------

def test_06c_refines_subfamily_not_class():
    """Pfam sets class=A; 06c hit aminergic_dopamine sets evidence_family_hmm.
    Class must remain A."""
    pfam_hits = [("PF00001", "seq1", 1e-50)]
    cls, evidence_pfam, _ = cgc.call_class(pfam_hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"

    # 06c refinement (target name without extension matches 'aminergic_dopamine')
    family_hits_06c = [("aminergic_dopamine", "seq1", 1e-80)]
    subfam = cgc.refine_subfamily(family_hits_06c)
    assert subfam == "aminergic_dopamine"
    # class remains A — it's determined by Pfam, never by refine_subfamily
    assert cls == "A"


# ---------------------------------------------------------------------------
# 10. 06c class-B/C/F HMMs excluded from refinement
# ---------------------------------------------------------------------------

def test_06c_classBC_F_hmms_excluded_from_refinement():
    """class-B-secretin, class-C, class-F-frizzled hits must NOT appear in
    evidence_family_hmm — they're excluded from the 06c refinement scan."""
    family_hits_06c = [
        ("class-B-secretin", "seq1", 1e-80),
        ("class-C", "seq1", 1e-60),
        ("class-F-frizzled", "seq1", 1e-50),
    ]
    subfam = cgc.refine_subfamily(family_hits_06c)
    assert subfam == ""


def test_06c_keeps_non_class_hits():
    """Non-class-level 06c hits (aminergic, opsin, etc.) ARE returned."""
    family_hits_06c = [
        ("class-B-secretin", "seq1", 1e-80),  # excluded
        ("opsin", "seq1", 1e-60),              # included — best remaining
    ]
    subfam = cgc.refine_subfamily(family_hits_06c)
    assert subfam == "opsin"


# ---------------------------------------------------------------------------
# 11. TIAMMAT preference: if tiammat_*.hmm found, prefer over pfam_fallback
# ---------------------------------------------------------------------------

def test_tiammat_preference(tmp_path):
    """select_hmm_library returns TIAMMAT paths when they exist."""
    pfam_dir = tmp_path / "pfam_fallback"
    pfam_dir.mkdir()
    (pfam_dir / "PF00001.hmm").write_text("NAME PF00001\n")

    tiammat_dir = tmp_path
    (tiammat_dir / "tiammat_PF00001.hmm").write_text("NAME tiammat_PF00001\n")

    result = cgc.select_hmm_library(
        pfam_dir=str(pfam_dir),
        family_hmm_dir=str(tiammat_dir),
        tiammat_glob="tiammat_*.hmm",
    )
    # Should prefer the tiammat file when it matches a required Pfam accession
    assert result["mode"] == "tiammat"
    assert any("tiammat_PF00001" in str(p) for p in result["paths"])


def test_tiammat_fallback_when_absent(tmp_path):
    """select_hmm_library falls back to pfam_fallback when no TIAMMAT HMMs."""
    pfam_dir = tmp_path / "pfam_fallback"
    pfam_dir.mkdir()
    (pfam_dir / "PF00001.hmm").write_text("NAME PF00001\n")

    result = cgc.select_hmm_library(
        pfam_dir=str(pfam_dir),
        family_hmm_dir=str(tmp_path),
        tiammat_glob="tiammat_*.hmm",
    )
    assert result["mode"] == "pfam_fallback"


# ---------------------------------------------------------------------------
# 12. Bootstrap download attempted for missing Pfam (monkeypatched)
# ---------------------------------------------------------------------------

def test_bootstrap_downloads_missing_pfams(tmp_path):
    """bootstrap_pfams() calls the download function for each missing Pfam
    and skips ones already on disk. Idempotent."""
    pfam_dir = tmp_path / "pfam_fallback"
    pfam_dir.mkdir()
    # Pre-seed PF00001 so it should be skipped
    (pfam_dir / "PF00001.hmm").write_text("NAME PF00001\n")

    required = ["PF00001", "PF00002"]  # PF00002 is missing

    downloaded = []

    def fake_download(accession, out_path):
        downloaded.append(accession)
        Path(out_path).write_text(f"NAME {accession}\n")

    cgc.bootstrap_pfams(
        pfam_dir=str(pfam_dir),
        required_pfams=required,
        _download_fn=fake_download,
    )

    # Only PF00002 should have been downloaded (PF00001 already present)
    assert downloaded == ["PF00002"]

    # Idempotency: running again downloads nothing
    downloaded.clear()
    cgc.bootstrap_pfams(
        pfam_dir=str(pfam_dir),
        required_pfams=required,
        _download_fn=fake_download,
    )
    assert downloaded == []


# ---------------------------------------------------------------------------
# 12b. Pfam HMM downloader uses a working endpoint (regression for HTTP 406)
# ---------------------------------------------------------------------------

def test_download_pfam_hmm_uses_working_endpoint(tmp_path, monkeypatch):
    """Regression: `Accept: application/octet-stream` made EBI return HTTP 406,
    and pfam.xfam.org is retired (SSL failure). The downloader must hit the
    public interpro API with no restrictive Accept header, and never touch
    pfam.xfam.org."""
    valid_hmm = b"HMMER3/f [3.3]\nNAME  7tm_4\nLENG  300\n//\n"
    gzipped = gzip.compress(valid_hmm)
    captured = []

    class FakeResp:
        def read(self):
            return gzipped
        def __enter__(self):
            return self
        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        captured.append(req)
        return FakeResp()

    monkeypatch.setattr(cgc.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(cgc, "_hmmpress", lambda p: None)

    out = tmp_path / "PF02949.hmm"
    cgc._download_pfam_hmm("PF02949", str(out))

    # a valid HMM was written (gzip transparently decompressed)
    assert out.read_bytes().startswith(b"HMMER3")
    # primary endpoint is the public interpro API, for the right accession
    assert "ebi.ac.uk/interpro/api/" in captured[0].full_url
    assert "PF02949" in captured[0].full_url
    # the 406-causing Accept header is gone from every request
    for req in captured:
        assert req.get_header("Accept") != "application/octet-stream"
    # the retired host is never contacted
    assert all("pfam.xfam.org" not in r.full_url for r in captured)


# ---------------------------------------------------------------------------
# 13. Output TSV column order
# ---------------------------------------------------------------------------

def test_output_tsv_columns(tmp_path):
    """write_classification_tsv produces the exact required column header
    and at least one data row."""
    rows = [
        {
            "seq_id": "seq1",
            "class": "A",
            "evidence_pfam": "PF00001",
            "evidence_family_hmm": "aminergic_dopamine",
            "top_evalue": "1.2e-50",
        }
    ]
    out_path = tmp_path / "out.tsv"
    cgc.write_classification_tsv(rows, str(out_path))

    with open(out_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        # Check exact column names in exact order
        assert reader.fieldnames == [
            "seq_id", "class", "evidence_pfam", "evidence_family_hmm", "top_evalue"
        ]
        data_rows = list(reader)
    assert len(data_rows) == 1
    assert data_rows[0]["seq_id"] == "seq1"
    assert data_rows[0]["class"] == "A"
    assert data_rows[0]["evidence_pfam"] == "PF00001"
    assert data_rows[0]["evidence_family_hmm"] == "aminergic_dopamine"
    assert data_rows[0]["top_evalue"] == "1.2e-50"


# ---------------------------------------------------------------------------
# 14. CLI argparse
# ---------------------------------------------------------------------------

def test_cli_argparse(tmp_path):
    """build_arg_parser() returns a parser that correctly handles all flags."""
    parser = cgc.build_arg_parser()

    # Minimal required args
    args = parser.parse_args([
        "--input", "test.fa",
        "--out", "out.tsv",
    ])
    assert args.input == "test.fa"
    assert args.out == "out.tsv"
    assert args.evalue == 1e-5           # default
    assert args.bootstrap_pfam is True   # default ON
    assert args.force is False           # default OFF

    # Full args
    args2 = parser.parse_args([
        "--input", "in.fa",
        "--out", "out.tsv",
        "--pfam-dir", "/pfams",
        "--family-hmm-dir", "/hmms",
        "--tiammat-glob", "tiammat_*.hmm",
        "--evalue", "1e-10",
        "--threads", "8",
        "--no-bootstrap",
        "--force",
    ])
    assert args2.pfam_dir == "/pfams"
    assert args2.family_hmm_dir == "/hmms"
    assert args2.tiammat_glob == "tiammat_*.hmm"
    assert args2.evalue == 1e-10
    assert args2.threads == 8
    assert args2.bootstrap_pfam is False
    assert args2.force is True


# ---------------------------------------------------------------------------
# 15. parse_hmmscan_tblout (reuse the same format as classify_via_hmm)
# ---------------------------------------------------------------------------

TBLOUT_SAMPLE = """\
# hmmscan tblout
#                         --- full sequence ---
# target name  accession  query name  accession  E-value  score  bias
# ------------ ---------  ----------- ---------  ------- ------ ----
PF00001        -          seq1        -          1.2e-50  200.0  0.0
PF00002        -          seq1        -          5.0e-10  80.0   0.0
PF00003        -          seq2        -          3.0e-80  350.0  0.0
# [ok]
"""


def test_parse_tblout_groups_and_sorts(tmp_path):
    tbl = tmp_path / "scan.tbl"
    tbl.write_text(TBLOUT_SAMPLE)
    hits = cgc.parse_hmmscan_tblout(str(tbl))
    assert set(hits.keys()) == {"seq1", "seq2"}
    # seq1 hits sorted best E-value first
    assert hits["seq1"][0][0] == "PF00001"
    assert hits["seq1"][0][2] == 1.2e-50
    assert hits["seq1"][1][0] == "PF00002"
    # seq2 single hit
    assert hits["seq2"][0][0] == "PF00003"


# ---------------------------------------------------------------------------
# 16. _PFAM_TO_CLASS completeness — 11 verified GPCR Pfams, 6 blocklisted
# ---------------------------------------------------------------------------

def test_pfam_to_class_holds_only_verified_gpcr_families():
    # Verified against InterPro 2026-07-20. Six accessions previously in this
    # map are NOT GPCRs and were moved to the blocklist; two of the comments
    # they carried (`7tm_8`, `7tm_9`) name Pfam entries that do not exist.
    required = {
        "PF00001", "PF00002", "PF00003", "PF01534",
        "PF02949", "PF05296", "PF10324", "PF08395",
        "PF03402", "PF13853", "PF10326",
    }
    blocklisted = {
        "PF12022",  # COG2_C          - COG complex component
        "PF11399",  # DUF3192
        "PF13863",  # DUF4200
        "PF13886",  # TM7S3_TM198     - 7-TM but NOT a GPCR; cleared the >=6-TM filter
        "PF13887",  # MYRF_ICA
        "PF13889",  # Chromosome_seg  - chromosome segregation during meiosis
    }
    assert required == set(cgc._PFAM_TO_CLASS.keys()), (
        f"Missing: {required - set(cgc._PFAM_TO_CLASS.keys())}\n"
        f"Extra:   {set(cgc._PFAM_TO_CLASS.keys()) - required}"
    )


def test_pfam_to_class_all_valid_class_values():
    valid_classes = {"A", "B", "C", "F"}
    for pfam, cls in cgc._PFAM_TO_CLASS.items():
        assert cls in valid_classes, f"{pfam} has invalid class {cls!r}"


# ---------------------------------------------------------------------------
# 17. refine_subfamily respects evalue threshold
# ---------------------------------------------------------------------------

def test_refine_subfamily_respects_evalue_threshold():
    """Hits above threshold are ignored during 06c refinement."""
    family_hits_06c = [
        ("aminergic", "seq1", 1e-3),  # above default 1e-5 → ignored
    ]
    subfam = cgc.refine_subfamily(family_hits_06c, evalue_threshold=1e-5)
    assert subfam == ""


def test_refine_subfamily_returns_best_hit_below_threshold():
    family_hits_06c = [
        ("opsin", "seq1", 1e-80),
        ("lipid", "seq1", 1e-20),
    ]
    subfam = cgc.refine_subfamily(family_hits_06c, evalue_threshold=1e-5)
    assert subfam == "opsin"


# ---------------------------------------------------------------------------
# 18. NEW: Versioned Pfam accession handling (regression test for P2 bug fix)
# ---------------------------------------------------------------------------

def test_call_class_handles_versioned_pfam_accession():
    """call_class must handle Pfam accessions with version suffix (e.g. PF00001.27).
    The parser strips the suffix, so this test ensures that call_class receives
    the bare accession and correctly identifies the class."""
    # Simulate what the parser produces after stripping "PF00001.27" -> "PF00001"
    hits = [("7tm_1", "PF00001", 1e-50)]
    cls, evidence_pfam, top_evalue = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A", "Versioned Pfam should still resolve to class A"
    assert evidence_pfam in ("7tm_1", "PF00001"), (
        "evidence_pfam should be either target name or bare accession"
    )
    assert top_evalue == 1e-50


def test_parse_tblout_strips_pfam_version_suffix(tmp_path):
    """parse_hmmscan_tblout must strip Pfam version suffixes from accession column.
    Real InterPro HMMs have versioned accessions like PF00001.27; we must strip
    to .split('.')[0] so lookups against _PFAM_TO_CLASS (which uses bare keys) work."""
    tblout_with_versioned = """\
# hmmscan tblout
#                         --- full sequence ---
# target name  accession  query name  accession  E-value  score  bias
# ------------ ---------  ----------- ---------  ------- ------ ----
7tm_1          PF00001.27 seq1        -          1.2e-50  200.0  0.0
"""
    tbl = tmp_path / "versioned.tbl"
    tbl.write_text(tblout_with_versioned)
    hits = cgc.parse_hmmscan_tblout(str(tbl))

    assert "seq1" in hits, "seq1 should be in the parsed results"
    assert len(hits["seq1"]) >= 1, "seq1 should have at least one hit"
    target, accession, evalue = hits["seq1"][0]
    assert accession == "PF00001", (
        f"Accession should be stripped to 'PF00001', got '{accession}'"
    )
    assert evalue == 1.2e-50


def test_pf02949_atypical_with_versioned_accession(tmp_path):
    """PF02949 (the real insect 7tm_6) with a versioned accession should still trigger
    the special subfamily annotation insect_OR_atypical."""
    # Simulate parsing a tblout with versioned PF10324
    tblout_pf10324 = """\
# hmmscan tblout
#                         --- full sequence ---
# target name  accession  query name  accession  E-value  score  bias
# ------------ ---------  ----------- ---------  ------- ------ ----
7tm_6          PF02949.13 seq_insect  -          1.5e-60  250.0  0.0
"""
    tbl = tmp_path / "pf02949_versioned.tbl"
    tbl.write_text(tblout_pf10324)
    hits = cgc.parse_hmmscan_tblout(str(tbl))

    assert "seq_insect" in hits
    target, accession, evalue = hits["seq_insect"][0]
    assert accession == "PF02949", "Version suffix should be stripped"

    # Now call_class and refine_subfamily with the stripped accession
    cls, evidence_pfam, _ = cgc.call_class(hits["seq_insect"], cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"
    assert evidence_pfam == "PF02949"

    # refine_subfamily should return the intrinsic _PFAM_SUBFAMILY annotation
    subfam = cgc.refine_subfamily([], pfam_accession="PF02949")
    assert subfam == "insect_OR_atypical", (
        "PF02949 (InterPro short name 7tm_6, the real insect odorant-receptor "
        "family) should map to insect_OR_atypical in _PFAM_SUBFAMILY"
    )


# ---------------------------------------------------------------------------
# 19. NEW: Bootstrap robustness — fallback URL + integrity check
# ---------------------------------------------------------------------------

def test_download_retries_on_interpro_failure(tmp_path):
    """_download_pfam_hmm should retry the /wwwapi/ fallback when the primary
    /api/ endpoint fails."""
    accession = "PF00001"
    out_path = tmp_path / "PF00001.hmm"

    valid_hmm_content = b"HMMER3/f\nNAME PF00001\n..."

    call_count = [0]

    def mock_urlopen(req, timeout=60):
        call_count[0] += 1
        url = req.full_url

        # First call (public /api/) raises exception
        if "/interpro/api/" in url:
            raise urllib.error.URLError("Connection refused")
        # Second call (/wwwapi/ fallback) succeeds
        elif "/interpro/wwwapi/" in url:
            return io.BytesIO(valid_hmm_content)
        else:
            raise ValueError(f"Unexpected URL: {url}")

    with patch("classify_gpcr_by_class.urllib.request.urlopen", side_effect=mock_urlopen):
        with patch("classify_gpcr_by_class._hmmpress") as mock_hmmpress:
            cgc._download_pfam_hmm(accession, str(out_path))

    # Should have tried both URLs
    assert call_count[0] == 2, "Should attempt both /api/ and /wwwapi/"
    # File should exist and contain the valid content
    assert out_path.exists()
    assert out_path.read_bytes() == valid_hmm_content
    # hmmpress should have been called
    mock_hmmpress.assert_called_once()


def test_download_rejects_html_body(tmp_path):
    """_download_pfam_hmm should reject HTML responses (invalid HMM format)."""
    accession = "PF00001"
    out_path = tmp_path / "PF00001.hmm"

    html_error = b"<html><body>404 Not Found</body></html>"

    def mock_urlopen(req, timeout=60):
        # Both URLs return HTML error
        return io.BytesIO(html_error)

    with patch("classify_gpcr_by_class.urllib.request.urlopen", side_effect=mock_urlopen):
        with patch("classify_gpcr_by_class._hmmpress"):
            # Should raise RuntimeError, not write the file
            with pytest.raises(RuntimeError) as exc_info:
                cgc._download_pfam_hmm(accession, str(out_path))

            error_msg = str(exc_info.value)
            assert "Failed to download valid HMM" in error_msg
            assert "interpro/api" in error_msg     # primary InterPro URL
            assert "interpro/wwwapi" in error_msg  # /wwwapi/ fallback URL

    # File should NOT exist
    assert not out_path.exists(), "Invalid HMM should not be written to disk"


def test_download_accepts_valid_uncompressed_hmm(tmp_path):
    """_download_pfam_hmm should accept valid uncompressed HMMER3 HMM."""
    accession = "PF00002"
    out_path = tmp_path / "PF00002.hmm"

    valid_hmm = b"HMMER3/f\nNAME PF00002\nACC  PF00002.14\n... rest of HMM ...\n"

    def mock_urlopen(req, timeout=60):
        return io.BytesIO(valid_hmm)

    with patch("classify_gpcr_by_class.urllib.request.urlopen", side_effect=mock_urlopen):
        with patch("classify_gpcr_by_class._hmmpress") as mock_hmmpress:
            cgc._download_pfam_hmm(accession, str(out_path))

    # File should exist with correct content
    assert out_path.exists()
    assert out_path.read_bytes() == valid_hmm
    mock_hmmpress.assert_called_once()


def test_download_accepts_valid_gzipped_hmm(tmp_path):
    """_download_pfam_hmm should decompress and accept gzipped HMMER3 HMM."""
    accession = "PF00003"
    out_path = tmp_path / "PF00003.hmm"

    valid_hmm_content = b"HMMER3/f\nNAME PF00003\nACC  PF00003.19\n"
    gzipped = gzip.compress(valid_hmm_content)

    def mock_urlopen(req, timeout=60):
        return io.BytesIO(gzipped)

    with patch("classify_gpcr_by_class.urllib.request.urlopen", side_effect=mock_urlopen):
        with patch("classify_gpcr_by_class._hmmpress") as mock_hmmpress:
            cgc._download_pfam_hmm(accession, str(out_path))

    # File should exist with decompressed content
    assert out_path.exists()
    assert out_path.read_bytes() == valid_hmm_content
    mock_hmmpress.assert_called_once()


# ---------------------------------------------------------------------------
# 20. _is_valid_hmm_format validation
# ---------------------------------------------------------------------------

def test_is_valid_hmm_format_accepts_hmmer3_f():
    """_is_valid_hmm_format should accept HMMER3/f header."""
    assert cgc._is_valid_hmm_format(b"HMMER3/f\nNAME PF00001\n")


def test_is_valid_hmm_format_accepts_hmmer3_slash():
    """_is_valid_hmm_format should accept HMMER3/ header variants."""
    assert cgc._is_valid_hmm_format(b"HMMER3/3\nNAME PF00001\n")
    assert cgc._is_valid_hmm_format(b"HMMER3/\nNAME PF00001\n")


def test_is_valid_hmm_format_accepts_bare_hmmer3():
    """_is_valid_hmm_format should accept bare HMMER3 header."""
    assert cgc._is_valid_hmm_format(b"HMMER3")


def test_is_valid_hmm_format_rejects_html():
    """_is_valid_hmm_format should reject HTML."""
    assert not cgc._is_valid_hmm_format(b"<html><body>404</body></html>")


def test_is_valid_hmm_format_rejects_empty():
    """_is_valid_hmm_format should reject empty or too-short data."""
    assert not cgc._is_valid_hmm_format(b"")
    assert not cgc._is_valid_hmm_format(b"HMM")  # 3 bytes, too short


def test_is_valid_hmm_format_rejects_random_data():
    """_is_valid_hmm_format should reject random binary data."""
    assert not cgc._is_valid_hmm_format(b"\x00\x01\x02\x03\x04\x05")
