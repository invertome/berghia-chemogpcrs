"""Live NCBI checks for the accessions the 2026-07 audit called into question.

These are the programmatic re-verification of hand-typed accessions required by
this project's standing rule: obtain, copy and verify identifiers
programmatically, never from recall and never hand-typed. Two independent
wrong-organism bugs (a fish and a dolphin standing in for molluscs) reached the
repo because accessions were transcribed rather than resolved.

They hit api.ncbi.nlm.nih.gov, so they are opt-in (the default
unit run stays offline and fast) and additionally SKIP cleanly whenever the API
is unreachable. Requests are throttled and time-limited to stay well inside
NCBI's unauthenticated rate limit.

Run them with:  RUN_NCBI_TESTS=1 pytest tests/unit/test_inventory_ncbi_*.py
"""
from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_NCBI_TESTS") != "1",
    reason="live NCBI checks are opt-in; set RUN_NCBI_TESTS=1 to run",
)

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{}/dataset_report"
TIMEOUT = 30
THROTTLE_S = 1.0


def _report(accession: str) -> dict | None:
    """First dataset report for an accession, or None when NCBI has no such
    assembly. Skips the test when the network/API is unavailable."""
    time.sleep(THROTTLE_S)
    try:
        with urllib.request.urlopen(API.format(accession), timeout=TIMEOUT) as fh:
            payload = json.load(fh)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        pytest.skip(f"NCBI datasets API unreachable: {exc}")
    except json.JSONDecodeError as exc:
        pytest.skip(f"NCBI returned non-JSON: {exc}")
    reports = payload.get("reports") or []
    return reports[0] if reports else None


def _organism(accession: str) -> tuple[str, int] | None:
    rep = _report(accession)
    if rep is None:
        return None
    org = rep.get("organism") or {}
    return org.get("organism_name", ""), int(org.get("tax_id", 0))


# ------------------------------------------------ the live duplicate (#6)

def test_dreissena_accession_belongs_to_the_subspecies_taxid():
    """GCA_055670145.1 is taxid 427924, so the taxid-205083 manifest row that
    also claims it is a mislabelled duplicate of the same genome."""
    got = _organism("GCA_055670145.1")
    assert got is not None, "GCA_055670145.1 should resolve"
    name, taxid = got
    assert taxid == 427924
    assert name == "Dreissena rostriformis bugensis"
    assert taxid != 205083


@pytest.mark.parametrize("accession,taxid,name", [
    ("GCA_034509925.1", 1348078, "Hirudinaria manillensis"),
    ("GCA_015776775.1", 2653900, "Magallana hongkongensis"),
])
def test_synonym_binomial_duplicate_pairs_share_one_assembly(accession, taxid, name):
    """Two further manifest pairs share an accession under synonym binomials;
    NCBI's accepted name settles which row is redundant."""
    got = _organism(accession)
    assert got is not None
    assert got == (name, taxid)


# --------------------------------------- the deleted scaffold accessions (#14)

@pytest.mark.parametrize("accession,intended", [
    ("GCA_019457155.2", "Chrysomallon squamiferum"),
    ("GCA_917563875.2", "Patella vulgata"),
    ("GCA_963853765.1", "Gibbula magus"),
    ("GCA_944038965.1", "Lymnaea stagnalis"),
    ("GCA_011762535.2", "Haliotis discus hannai"),
])
def test_deleted_scaffold_accessions_do_not_match_their_intended_species(
        accession, intended):
    """None of the five resolve to the mollusc they were labelled with --
    the justification for deleting SPECIES_ASSEMBLIES outright."""
    got = _organism(accession)
    if got is None:
        return  # no such assembly at this version: also not the intended species
    name, _taxid = got
    assert intended.split()[0] not in name, (
        f"{accession} unexpectedly matches {intended}; re-open the deletion")


@pytest.mark.parametrize("versionless,name,taxid", [
    ("GCA_019457155", "Streptococcus pneumoniae", 1313),
    ("GCA_011762535", "Tursiops truncatus", 9739),
])
def test_two_scaffold_accessions_are_wrong_kingdom_at_their_real_version(
        versionless, name, taxid):
    """The '.2' versions never existed; the accession NUMBERS belong to a
    bacterium and a dolphin."""
    got = _organism(versionless)
    assert got is not None
    assert got == (name, taxid)


# ------------------------------------------------------ the date field (#17)

def test_assembly_info_exposes_release_date_not_submission_date():
    rep = _report("GCA_015776775.1")
    assert rep is not None
    info = rep.get("assembly_info") or {}
    assert "release_date" in info
    assert "submission_date" not in info
    assert info["release_date"]
