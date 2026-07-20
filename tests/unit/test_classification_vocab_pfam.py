"""Pin the Pfam accession vocabulary used by the GPCR class classifier.

Audit finding #4: `_PFAM_TO_CLASS` shipped 17 accessions all mapped to
class "A", of which SIX are not GPCRs at all (a COG-complex repeat, two
DUFs, a myelin regulatory factor domain, a meiotic chromosome-segregation
family, and a TM7S3/TM198 domain). Twelve of the seventeen inline comments
named a different family than the accession actually encodes.

The dangerous case is PF13886 (TM7S3_TM198): it is genuinely a 7-TM
domain, so a hit clears the downstream >=6-TM filter and enters the
candidate set as a "class A GPCR".

These tests pin every accession to the name InterPro actually returns, so
a future edit cannot silently reintroduce a wrong mapping. The offline
tests pin the recorded ground truth; the network test re-verifies that
recorded truth against the live InterPro API and SKIPS (never fails) when
the network is unavailable.

Ground truth recorded 2026-07-20 from
https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/
"""
from __future__ import annotations

import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

import classify_gpcr_by_class as cgc


# ---------------------------------------------------------------------------
# Ground truth: accession -> InterPro Pfam "short name", verified 2026-07-20.
# ---------------------------------------------------------------------------

VERIFIED_PFAM_NAMES = {
    "PF00001": "7tm_1",
    "PF00002": "7tm_2",
    "PF00003": "7tm_3",
    "PF01534": "Frizzled",
    "PF02949": "7tm_6",
    "PF05296": "TAS2R",
    "PF08395": "7tm_7",
    "PF03402": "V1R",
    "PF13853": "7tm_4",
    "PF10324": "7TM_GPCR_Srw",
    "PF10326": "7TM_GPCR_Str",
}

# Accessions that were in the shipped map but are NOT GPCRs. Verified against
# InterPro 2026-07-20. These must never reappear in _PFAM_TO_CLASS.
NON_GPCR_ACCESSIONS = {
    "PF12022": "COG2_C",           # COG complex component, COG2, C-terminal
    "PF11399": "DUF3192",          # Protein of unknown function
    "PF13863": "DUF4200",          # Domain of unknown function
    "PF13886": "TM7S3_TM198",      # TM7S3/TM198-like domain (7-TM, NOT a GPCR)
    "PF13887": "MYRF_ICA",         # Myelin regulatory factor ICA domain
    "PF13889": "Chromosome_seg",   # Chromosome segregation during meiosis
}


# ---------------------------------------------------------------------------
# 1. The map is fully name-pinned
# ---------------------------------------------------------------------------

def test_every_mapped_accession_has_a_verified_name():
    """_PFAM_VERIFIED_NAME must cover exactly the accessions in _PFAM_TO_CLASS.

    Pinning names in a dict (rather than a comment) is what makes the
    mapping machine-checkable — comments drifted from reality on 12 of 17
    entries precisely because nothing tested them.
    """
    assert set(cgc._PFAM_TO_CLASS) == set(cgc._PFAM_VERIFIED_NAME)


@pytest.mark.parametrize("accession,name", sorted(VERIFIED_PFAM_NAMES.items()))
def test_accession_name_mapping_is_pinned(accession, name):
    """Each accession resolves to the family name InterPro actually reports."""
    assert cgc._PFAM_VERIFIED_NAME[accession] == name


def test_pfam_map_contains_exactly_the_verified_accessions():
    assert set(cgc._PFAM_TO_CLASS) == set(VERIFIED_PFAM_NAMES)


# ---------------------------------------------------------------------------
# 2. No non-GPCR accession may reappear
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("accession,name", sorted(NON_GPCR_ACCESSIONS.items()))
def test_non_gpcr_accession_is_absent(accession, name):
    """A non-GPCR Pfam accession must not be classified as a GPCR class."""
    assert accession not in cgc._PFAM_TO_CLASS, (
        f"{accession} is {name}, not a GPCR — it must not map to a GPCR class"
    )
    assert accession not in cgc._PFAM_VERIFIED_NAME


def test_non_gpcr_blocklist_is_recorded_in_the_module():
    """The removed accessions are recorded in-module so the rationale
    travels with the code and a re-add is an obvious conflict."""
    assert set(cgc._NON_GPCR_PFAM_BLOCKLIST) == set(NON_GPCR_ACCESSIONS)


def test_blocklist_and_active_map_are_disjoint():
    assert not (set(cgc._NON_GPCR_PFAM_BLOCKLIST) & set(cgc._PFAM_TO_CLASS))


# ---------------------------------------------------------------------------
# 3. Class assignments
# ---------------------------------------------------------------------------

def test_class_letters_are_correct_for_the_four_canonical_families():
    assert cgc._PFAM_TO_CLASS["PF00001"] == "A"   # rhodopsin-like
    assert cgc._PFAM_TO_CLASS["PF00002"] == "B"   # secretin
    assert cgc._PFAM_TO_CLASS["PF00003"] == "C"   # glutamate-like
    assert cgc._PFAM_TO_CLASS["PF01534"] == "F"   # frizzled/smoothened


def test_all_classes_are_valid_letters():
    assert set(cgc._PFAM_TO_CLASS.values()) <= {"A", "B", "C", "F"}


# ---------------------------------------------------------------------------
# 4. The insect-OR subfamily tag is on the insect-OR accession
# ---------------------------------------------------------------------------

def test_insect_or_tag_is_on_PF02949_not_PF10324():
    """PF02949 (7tm_6) IS the insect odorant receptor family with inverted
    topology. PF10324 is 7TM_GPCR_Srw, a *C. elegans* serpentine
    chemoreceptor — tagging it 'insect_OR_atypical' routed real
    chemoreceptors out of the class-A tree.
    """
    assert cgc._PFAM_SUBFAMILY.get("PF02949") == "insect_OR_atypical"
    assert "PF10324" not in cgc._PFAM_SUBFAMILY


def test_srw_is_class_A_with_no_atypical_subfamily():
    """A PF10324 (Srw) hit must yield a plain class-A call."""
    hits = [("PF10324", "PF10324", 1e-40)]
    cls, evidence, _ = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"
    assert cgc.refine_subfamily([], pfam_accession=evidence) == ""


def test_insect_or_hit_still_gets_the_atypical_subfamily():
    hits = [("PF02949", "PF02949", 1e-40)]
    cls, evidence, _ = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "A"
    assert cgc.refine_subfamily([], pfam_accession=evidence) == "insect_OR_atypical"


# ---------------------------------------------------------------------------
# 5. Regression: the specific mis-calls the audit found
# ---------------------------------------------------------------------------

def test_chromosome_segregation_hit_is_not_called_a_gpcr():
    """PF13889 is 'Chromosome segregation during meiosis'. At the stage-02
    default of --evalue 1e-5 this used to return class 'A'."""
    hits = [("PF13889", "PF13889", 1e-30)]
    cls, evidence, _ = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "unclassified"
    assert evidence == ""


def test_tm7s3_hit_is_not_called_a_gpcr():
    """PF13886 (TM7S3_TM198) is genuinely 7-TM, so it clears the >=6-TM
    filter downstream — the worst case, because nothing else catches it."""
    hits = [("PF13886", "PF13886", 1e-45)]
    cls, _, _ = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert cls == "unclassified"


def test_a_real_rhodopsin_hit_still_classifies_after_the_pruning():
    """Guard against over-pruning: the canonical path must still work."""
    hits = [("PF00001", "PF00001", 1e-80)]
    cls, evidence, ev = cgc.call_class(hits, cgc._PFAM_TO_CLASS, 1e-5)
    assert (cls, evidence, ev) == ("A", "PF00001", 1e-80)


def test_bootstrap_only_requests_verified_accessions(tmp_path):
    """The bootstrap downloader must never fetch a blocklisted accession."""
    requested: list[str] = []

    def fake_download(accession, out_path):
        requested.append(accession)
        Path(out_path).write_bytes(b"HMMER3/f fake\n")

    cgc.bootstrap_pfams(str(tmp_path), list(cgc._PFAM_TO_CLASS),
                        _download_fn=fake_download)
    assert set(requested) == set(VERIFIED_PFAM_NAMES)
    assert not (set(requested) & set(NON_GPCR_ACCESSIONS))


# ---------------------------------------------------------------------------
# 6. Live re-verification against InterPro (SKIPS cleanly offline)
# ---------------------------------------------------------------------------

# The unit suite is hermetic by convention, so the live re-verification is
# opt-in. Re-run it whenever this vocabulary is edited:
#
#     VERIFY_ACCESSIONS_ONLINE=1 pytest tests/unit/test_classification_vocab_pfam.py
#
# Requests are throttled; any network failure SKIPS so an offline machine
# can never see a red suite. Only a genuine name mismatch fails.
_ONLINE = os.environ.get("VERIFY_ACCESSIONS_ONLINE") == "1"
_online_only = pytest.mark.skipif(
    not _ONLINE,
    reason="set VERIFY_ACCESSIONS_ONLINE=1 to re-verify against InterPro",
)


def _interpro_short_name(accession: str, timeout: float = 30.0) -> str:
    time.sleep(1.0)  # be polite to the EBI API
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.load(resp)
    return payload["metadata"]["name"]["short"]


@_online_only
@pytest.mark.parametrize("accession,name", sorted(VERIFIED_PFAM_NAMES.items()))
def test_recorded_name_matches_live_interpro(accession, name):
    """Re-verify the pinned name against the authoritative source.

    Network failure must SKIP, never fail — offline CI must stay green.
    Only a genuine name mismatch is a failure.
    """
    try:
        live = _interpro_short_name(accession)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        pytest.skip(f"InterPro unreachable ({exc}) — offline, skipping")
    except (KeyError, ValueError) as exc:
        pytest.skip(f"Unexpected InterPro payload for {accession}: {exc}")
    assert live == name, (
        f"{accession} is '{live}' at InterPro but pinned as '{name}'"
    )


@_online_only
@pytest.mark.parametrize("accession,name", sorted(NON_GPCR_ACCESSIONS.items()))
def test_blocklisted_accession_is_still_not_a_gpcr_upstream(accession, name):
    """Confirm the removed accessions really are the non-GPCR families we
    recorded, so the removal rationale stays auditable."""
    try:
        live = _interpro_short_name(accession)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        pytest.skip(f"InterPro unreachable ({exc}) — offline, skipping")
    except (KeyError, ValueError) as exc:
        pytest.skip(f"Unexpected InterPro payload for {accession}: {exc}")
    assert live == name
