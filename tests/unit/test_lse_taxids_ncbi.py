"""Live NCBI verification of the LSE clade taxids.

Every test here talks to the NCBI taxonomy API and SKIPS cleanly when there is
no network, so the suite stays green on compute nodes and in CI. Nothing here
may ever fail because of connectivity; a failure means a taxid is genuinely
wrong.

This is the test that would have caught the original defect. The offline tests
in test_lse_taxids_config.py pin the literals; these confirm the literals mean
what the comments next to them claim.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

verify_lse_taxids = pytest.importorskip("verify_lse_taxids")

BERGHIA_TAXID = 1287507

# taxid -> (expected scientific name, expected rank)
EXPECTED = {
    71481: ("Aeolidioidea", "superfamily"),
    70849: ("Nudibranchia", "order"),
    6448: ("Gastropoda", "class"),
}


@pytest.fixture(scope="module")
def fetcher():
    """A live fetcher, or a skip if NCBI is unreachable."""
    f = verify_lse_taxids.NCBIFetcher(timeout=20)
    try:
        f.summaries([BERGHIA_TAXID])
    except verify_lse_taxids.NetworkUnavailable as exc:
        pytest.skip(f"NCBI taxonomy API unreachable: {exc}")
    return f


@pytest.fixture(scope="module")
def berghia_lineage(fetcher):
    try:
        return fetcher.lineage(BERGHIA_TAXID)
    except verify_lse_taxids.NetworkUnavailable as exc:
        pytest.skip(f"NCBI taxonomy API unreachable: {exc}")


@pytest.mark.parametrize("taxid,expected", sorted(EXPECTED.items()))
def test_taxid_resolves_to_expected_clade(fetcher, taxid, expected):
    name, rank = expected
    try:
        record = fetcher.summaries([taxid])[taxid]
    except verify_lse_taxids.NetworkUnavailable as exc:
        pytest.skip(f"NCBI taxonomy API unreachable: {exc}")
    assert record is not None, f"taxid {taxid} does not resolve at NCBI"
    assert record["scientificname"] == name
    assert record["rank"] == rank


@pytest.mark.parametrize("taxid", sorted(EXPECTED))
def test_taxid_is_an_ancestor_of_berghia(berghia_lineage, taxid):
    """The property lse_refine.py actually relies on."""
    assert taxid in berghia_lineage, (
        f"taxid {taxid} is not in the lineage of Berghia stephanieae; "
        f"lse_refine.py's `in common_lineage` test can never fire for it"
    )


def test_clade_taxids_are_properly_nested(berghia_lineage):
    """Aeolidioidea inside Nudibranchia inside Gastropoda.

    lse_refine.py returns the first match of an if/elif chain, so the levels
    must nest most-specific-first or the classification is meaningless.
    """
    order = [berghia_lineage.index(t) for t in (6448, 70849, 71481)]
    assert order == sorted(order), (
        "expected Gastropoda -> Nudibranchia -> Aeolidioidea nesting"
    )


@pytest.mark.parametrize("wrong,claimed", [
    (644, "Gastropoda"),
    (54397, "Aeolidida"),
    (6524, "Nudibranchs (old config comment example)"),
])
def test_historical_wrong_taxid_is_not_a_berghia_ancestor(berghia_lineage, wrong, claimed):
    """Confirms against the live API that the old values were genuinely broken."""
    assert wrong not in berghia_lineage, (
        f"taxid {wrong} was claimed to be {claimed}"
    )


def test_the_unresolvable_historical_taxid_still_does_not_resolve(fetcher):
    """13843 had no NCBI record at all."""
    try:
        record = fetcher.summaries([13843])[13843]
    except verify_lse_taxids.NetworkUnavailable as exc:
        pytest.skip(f"NCBI taxonomy API unreachable: {exc}")
    assert record is None, f"13843 unexpectedly resolves to {record}"


def test_guard_passes_against_the_live_api(fetcher):
    """The whole guard, end to end, against real NCBI."""
    results = verify_lse_taxids.verify_taxids(
        verify_lse_taxids.EXPECTED_CLADES_AS_CONFIG, fetcher
    )
    unknown = [r for r in results if r.status == verify_lse_taxids.STATUS_UNKNOWN]
    if unknown:
        pytest.skip("NCBI became unreachable mid-test")
    assert all(r.status == verify_lse_taxids.STATUS_OK for r in results), [
        (r.variable, r.status, r.detail) for r in results
    ]
