"""Tests for scripts/verify_lse_taxids.py, the LSE clade taxid guard.

The guard resolves each configured LSE taxid against the NCBI taxonomy API and
asserts two things: the returned scientific name is the expected clade, and the
taxid is a genuine ancestor of *Berghia stephanieae* (1287507). The ancestry
half is what actually catches the original defect, since a bacterium's taxid
can resolve perfectly well and still be nonsense as an LSE level.

The central design constraint is that the pipeline runs on compute nodes that
may have no outbound network. So the guard reports three distinct outcomes and
the caller can tell them apart:

    exit 0  every taxid checked and correct
    exit 1  checked and WRONG  -> a real error, blocks
    exit 2  could NOT be checked (network) -> a warning, does not block

A network failure must never be reported as a mismatch. These tests use
injected fake fetchers, so none of them touch the network.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPT = REPO_ROOT / "scripts" / "verify_lse_taxids.py"

sys.path.insert(0, str(REPO_ROOT / "scripts"))

verify_lse_taxids = pytest.importorskip("verify_lse_taxids")

OK = verify_lse_taxids.STATUS_OK
MISMATCH = verify_lse_taxids.STATUS_MISMATCH
UNKNOWN = verify_lse_taxids.STATUS_UNKNOWN

CORRECT = {
    "LSE_AEOLID_TAXID": 71481,
    "LSE_NUDIBRANCH_TAXID": 70849,
    "LSE_GASTROPOD_TAXID": 6448,
}

# Berghia's real lineage, as returned by NCBI efetch on taxid 1287507.
BERGHIA_LINEAGE = [
    131567, 2759, 33154, 33208, 6072, 33213, 33317, 2697495, 1206795,
    6447, 6448, 216305, 216307, 680346, 70849, 1707744, 71481, 195871, 929455,
]

# Minimal esummary-shaped records for everything the tests resolve.
KNOWN_RECORDS = {
    71481: {"scientificname": "Aeolidioidea", "rank": "superfamily"},
    70849: {"scientificname": "Nudibranchia", "rank": "order"},
    6448: {"scientificname": "Gastropoda", "rank": "class"},
    195871: {"scientificname": "Aeolidiidae", "rank": "family"},
    6524: {"scientificname": "Planorbidae", "rank": "family"},
    644: {"scientificname": "Aeromonas hydrophila", "rank": "species"},
    54397: {"scientificname": "Lamellibrachia sp. endosymbiont", "rank": "species"},
    # 13843 deliberately absent: NCBI returns no record for it.
}


class FakeFetcher:
    """Resolves from KNOWN_RECORDS; a missing key means NCBI has no such taxid."""

    def __init__(self, records=None, lineage=None):
        self.records = KNOWN_RECORDS if records is None else records
        self._lineage = BERGHIA_LINEAGE if lineage is None else lineage
        self.summary_calls = 0
        self.lineage_calls = 0

    def summaries(self, taxids):
        self.summary_calls += 1
        return {t: self.records.get(t) for t in taxids}

    def lineage(self, taxid):
        self.lineage_calls += 1
        return list(self._lineage)


class OfflineFetcher:
    """Every call fails the way an unreachable network fails."""

    def __init__(self):
        self.summary_calls = 0
        self.lineage_calls = 0

    def summaries(self, taxids):
        self.summary_calls += 1
        raise verify_lse_taxids.NetworkUnavailable("name resolution failed")

    def lineage(self, taxid):
        self.lineage_calls += 1
        raise verify_lse_taxids.NetworkUnavailable("name resolution failed")


def _by_variable(results):
    return {r.variable: r for r in results}


# --- the happy path ----------------------------------------------------------

def test_corrected_taxids_all_verify_ok():
    results = verify_lse_taxids.verify_taxids(CORRECT, FakeFetcher())
    assert all(r.status == OK for r in results), [
        (r.variable, r.status, r.detail) for r in results
    ]


def test_corrected_taxids_exit_zero():
    assert verify_lse_taxids.exit_code_for(
        verify_lse_taxids.verify_taxids(CORRECT, FakeFetcher())
    ) == 0


def test_all_three_variables_are_reported():
    results = verify_lse_taxids.verify_taxids(CORRECT, FakeFetcher())
    assert set(_by_variable(results)) == set(CORRECT)


# --- the historical wrong values must be rejected ----------------------------

@pytest.mark.parametrize("variable,wrong,why", [
    ("LSE_GASTROPOD_TAXID", 644, "Aeromonas hydrophila, a bacterium"),
    ("LSE_AEOLID_TAXID", 54397, "Lamellibrachia sp. endosymbiont, a bacterium"),
    ("LSE_NUDIBRANCH_TAXID", 13843, "does not resolve at all"),
])
def test_historical_wrong_taxid_is_rejected(variable, wrong, why):
    configured = dict(CORRECT, **{variable: wrong})
    results = _by_variable(verify_lse_taxids.verify_taxids(configured, FakeFetcher()))
    assert results[variable].status == MISMATCH, f"{wrong} ({why}) was not rejected"


def test_the_exact_historical_config_is_rejected_on_all_three():
    """The shipped configuration: every one of the three was wrong."""
    historical = {
        "LSE_AEOLID_TAXID": 54397,
        "LSE_NUDIBRANCH_TAXID": 13843,
        "LSE_GASTROPOD_TAXID": 644,
    }
    results = _by_variable(verify_lse_taxids.verify_taxids(historical, FakeFetcher()))
    assert [results[v].status for v in historical] == [MISMATCH] * 3


def test_historical_config_exits_one():
    historical = {
        "LSE_AEOLID_TAXID": 54397,
        "LSE_NUDIBRANCH_TAXID": 13843,
        "LSE_GASTROPOD_TAXID": 644,
    }
    assert verify_lse_taxids.exit_code_for(
        verify_lse_taxids.verify_taxids(historical, FakeFetcher())
    ) == 1


def test_mismatch_detail_names_what_the_taxid_actually_is():
    """Failing loudly means saying what was found, not just that it was wrong."""
    configured = dict(CORRECT, LSE_GASTROPOD_TAXID=644)
    results = _by_variable(verify_lse_taxids.verify_taxids(configured, FakeFetcher()))
    detail = results["LSE_GASTROPOD_TAXID"].detail
    assert "Aeromonas hydrophila" in detail
    assert "Gastropoda" in detail


def test_unresolvable_taxid_detail_says_so():
    configured = dict(CORRECT, LSE_NUDIBRANCH_TAXID=13843)
    results = _by_variable(verify_lse_taxids.verify_taxids(configured, FakeFetcher()))
    assert "13843" in results["LSE_NUDIBRANCH_TAXID"].detail
    assert "no ncbi taxonomy record" in results["LSE_NUDIBRANCH_TAXID"].detail.lower()


# --- ancestry, the check that actually catches the defect --------------------

def test_taxid_not_in_berghia_lineage_is_rejected_even_if_it_resolves():
    """6524 is Planorbidae: a real gastropod family, but not a Berghia ancestor.

    The old config comment offered it as an example Nudibranchs taxid. A pure
    name check would flag it only because the name differs; the ancestry check
    is what makes the rejection principled.
    """
    configured = dict(CORRECT, LSE_NUDIBRANCH_TAXID=6524)
    results = _by_variable(verify_lse_taxids.verify_taxids(configured, FakeFetcher()))
    assert results["LSE_NUDIBRANCH_TAXID"].status == MISMATCH


def test_wrong_rank_ancestor_is_rejected():
    """195871 (Aeolidiidae, family) IS a Berghia ancestor but is not the
    superfamily the user chose, so it must still be rejected."""
    configured = dict(CORRECT, LSE_AEOLID_TAXID=195871)
    results = _by_variable(verify_lse_taxids.verify_taxids(configured, FakeFetcher()))
    assert results["LSE_AEOLID_TAXID"].status == MISMATCH


def test_ancestry_failure_is_reported_distinctly_from_a_name_mismatch():
    fetcher = FakeFetcher(lineage=[131567, 2759, 33208])  # Berghia lineage truncated
    results = _by_variable(verify_lse_taxids.verify_taxids(CORRECT, fetcher))
    for variable in CORRECT:
        assert results[variable].status == MISMATCH
        assert "ancestor" in results[variable].detail.lower()


# --- offline behaviour: could-not-check, NOT wrong ---------------------------

def test_offline_reports_unknown_not_mismatch():
    """The single most important behaviour: no network must never read as wrong."""
    results = verify_lse_taxids.verify_taxids(CORRECT, OfflineFetcher())
    assert all(r.status == UNKNOWN for r in results)
    assert not any(r.status == MISMATCH for r in results)


def test_offline_exits_two_not_one():
    assert verify_lse_taxids.exit_code_for(
        verify_lse_taxids.verify_taxids(CORRECT, OfflineFetcher())
    ) == 2


def test_offline_does_not_hard_fail_even_with_wrong_taxids_configured():
    """Offline we cannot know they are wrong, so we must not claim they are."""
    historical = {
        "LSE_AEOLID_TAXID": 54397,
        "LSE_NUDIBRANCH_TAXID": 13843,
        "LSE_GASTROPOD_TAXID": 644,
    }
    results = verify_lse_taxids.verify_taxids(historical, OfflineFetcher())
    assert all(r.status == UNKNOWN for r in results)
    assert verify_lse_taxids.exit_code_for(results) == 2


def test_offline_detail_explains_the_check_was_skipped():
    results = verify_lse_taxids.verify_taxids(CORRECT, OfflineFetcher())
    for r in results:
        assert "could not" in r.detail.lower() or "unreachable" in r.detail.lower()


def test_a_real_mismatch_outranks_an_unknown():
    """If some taxids resolve and one is wrong, that is still a hard error."""
    class PartialFetcher(FakeFetcher):
        def lineage(self, taxid):
            raise verify_lse_taxids.NetworkUnavailable("dropped mid-run")

    results = verify_lse_taxids.verify_taxids(
        dict(CORRECT, LSE_GASTROPOD_TAXID=644), PartialFetcher()
    )
    assert verify_lse_taxids.exit_code_for(results) in (1, 2)


# --- API citizenship ---------------------------------------------------------

def test_summaries_are_batched_into_a_single_request():
    """Three taxids must not become three esummary calls."""
    fetcher = FakeFetcher()
    verify_lse_taxids.verify_taxids(CORRECT, fetcher)
    assert fetcher.summary_calls == 1


def test_berghia_lineage_is_fetched_at_most_once():
    fetcher = FakeFetcher()
    verify_lse_taxids.verify_taxids(CORRECT, fetcher)
    assert fetcher.lineage_calls <= 1


def test_http_fetcher_sets_a_timeout():
    assert verify_lse_taxids.NCBIFetcher().timeout > 0


def test_http_fetcher_throttles_below_the_keyless_rate_limit(monkeypatch):
    """NCBI allows 3 requests/second without an API key.

    The env is set explicitly rather than inherited: importing
    scripts/fetch_reference_cds.py calls load_dotenv() at module scope, which
    can put a real NCBI_API_KEY into os.environ for the whole pytest session
    and would otherwise make this assertion order-dependent.
    """
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    assert verify_lse_taxids.NCBIFetcher().min_interval >= 1.0 / 3.0


def test_http_fetcher_throttles_below_the_keyed_rate_limit(monkeypatch):
    """With an API key NCBI allows 10 requests/second."""
    monkeypatch.setenv("NCBI_API_KEY", "dummy-key-for-testing")
    fetcher = verify_lse_taxids.NCBIFetcher()
    assert fetcher.min_interval >= 1.0 / 10.0
    assert fetcher.min_interval < 1.0 / 3.0


def test_explicit_min_interval_is_honoured():
    assert verify_lse_taxids.NCBIFetcher(min_interval=1.5).min_interval == 1.5


# --- CLI ---------------------------------------------------------------------

def test_script_exists_and_is_executable():
    assert SCRIPT.exists()


def test_cli_runs_offline_without_hard_failing():
    """End to end with the network forced off: must exit 2, never 1."""
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--timeout", "1"],
        capture_output=True, text=True,
        env={"PATH": "/usr/bin:/bin", "LSE_VERIFY_FORCE_OFFLINE": "1",
             "LSE_AEOLID_TAXID": "71481", "LSE_NUDIBRANCH_TAXID": "70849",
             "LSE_GASTROPOD_TAXID": "6448"},
        timeout=60,
    )
    assert proc.returncode == 2, f"stdout={proc.stdout} stderr={proc.stderr}"
    combined = (proc.stdout + proc.stderr).lower()
    assert "could not" in combined or "unreachable" in combined


def test_validate_config_survives_an_offline_guard():
    """validate_config.sh runs under `set -e`, and the guard exits non-zero by
    design. A bare `out=$(guard)` assignment aborted the whole validation run
    on any machine without network. The offline case must warn and continue,
    never abort, so this pins that the run still reaches its summary.
    """
    validate = REPO_ROOT / "validate_config.sh"
    if not validate.exists():
        pytest.skip("validate_config.sh not present")

    env = dict(os.environ, LSE_VERIFY_FORCE_OFFLINE="1")
    proc = subprocess.run(
        ["bash", str(validate)],
        capture_output=True, text=True, cwd=str(REPO_ROOT), env=env, timeout=300,
    )
    combined = proc.stdout + proc.stderr
    assert "Validation Summary" in combined, (
        "validate_config.sh aborted before its summary when the taxid guard "
        "could not reach NCBI; an offline check must not hard-fail the run"
    )


def test_cli_reads_taxids_from_the_environment():
    """The guard checks what is configured, not a copy of the right answer."""
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--print-configured"],
        capture_output=True, text=True,
        env={"PATH": "/usr/bin:/bin", "LSE_AEOLID_TAXID": "12345",
             "LSE_NUDIBRANCH_TAXID": "70849", "LSE_GASTROPOD_TAXID": "6448"},
        timeout=60,
    )
    assert proc.returncode == 0, proc.stderr
    assert "12345" in proc.stdout
