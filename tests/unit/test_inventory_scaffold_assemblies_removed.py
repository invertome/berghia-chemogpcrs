"""``fetch_cds_from_scaffolds`` must not carry unverified assembly accessions.

2026-07 audit finding: the module defined a ``SPECIES_ASSEMBLIES`` dict of five
hand-typed accessions, none of which resolve to the intended mollusc. Verified
independently against the NCBI datasets v2alpha API on 2026-07-19:

    chsq  Chrysomallon squamiferum  GCA_019457155.2 -> version .2 does not
                                    exist; GCA_019457155.1 = Streptococcus
                                    pneumoniae (taxid 1313), a bacterium
    pavu  Patella vulgata           GCA_917563875.2 = Diabrotica virgifera
                                    virgifera (taxid 50390), a beetle
    gima  Gibbula magus             GCA_963853765.1 = Magallana gigas
                                    (taxid 29159), an oyster
    lyst  Lymnaea stagnalis         GCA_944038965.1 -> no such assembly at
                                    any version
    hadi  Haliotis discus hannai    GCA_011762535.2 -> version .2 does not
                                    exist; GCA_011762535.1 = Tursiops
                                    truncatus (taxid 9739), a dolphin

5/5 wrong. The dict had no caller anywhere in the repo, so it was deleted
rather than corrected -- correcting it would have meant inventing five
accessions for a code path that does not exist. This test keeps it deleted:
a re-introduced hardcoded accession map must arrive with verification.
"""
from __future__ import annotations

import re
from pathlib import Path

import fetch_cds_from_scaffolds

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
MODULE_SRC = PROJECT_ROOT / "scripts" / "fetch_cds_from_scaffolds.py"

# The five debunked accessions, versionless so a re-versioned copy still trips.
DEBUNKED = ("GCA_019457155", "GCA_917563875", "GCA_963853765",
            "GCA_944038965", "GCA_011762535")


def test_species_assemblies_constant_is_gone():
    assert not hasattr(fetch_cds_from_scaffolds, "SPECIES_ASSEMBLIES")


def test_no_debunked_accession_remains_in_the_source():
    src = MODULE_SRC.read_text()
    # Allow them inside comments/docstrings only if someone documents the
    # removal; assert they are not live code by checking for assignment.
    for acc in DEBUNKED:
        assert not re.search(rf"['\"]{acc}", src), f"{acc} reintroduced as a literal"


def test_module_still_exposes_its_actual_entry_points():
    """Deletion must not have taken working code with it."""
    assert callable(fetch_cds_from_scaffolds.fetch_scaffold_region)
    assert callable(fetch_cds_from_scaffolds.batch_fetch_scaffolds)
