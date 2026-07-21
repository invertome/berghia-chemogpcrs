"""Guards for the beyond-Lophotrochozoa reference sweep (refexp3).

Every protein name, organism, accession and lineage in this module was read
programmatically out of a live UniProt response during the 2026-07-20 sweep and
is reproduced verbatim. They are real records, not invented fixtures: a
synthetic fixture encodes whatever the author assumed and therefore cannot
catch the case where the assumption itself is wrong, which is exactly the
failure these tests exist to catch.

Two properties are pinned here.

1. refexp3 must never DECIDE anything. The gate lives in refexp2 and is the
   artifact the user approved; this script contributes a taxonomic sweep and a
   reporting layer. A test below asserts the reporting layer cannot flip a
   verdict.

2. The reporting layer must not quietly re-introduce the bias the sweep exists
   to remove. The extension vocabulary is allowed to LABEL an entry, never to
   admit one, and it must defer to refexp2 wherever refexp2 already has an
   opinion so the two runs stay comparable.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))

import refexp2_evidence_gate as gate  # noqa: E402
import refexp3_beyond_lophotrochozoa as sweep  # noqa: E402


def entry(lineage, **kwargs):
    base = {
        "primaryAccession": "A0A0K0PUF3",
        "entryType": "UniProtKB unreviewed (TrEMBL)",
        "proteinExistence": "2: Evidence at transcript level",
        "proteinDescription": {},
        "organism": {"scientificName": "Hydra vulgaris", "taxonId": 6087,
                     "lineage": lineage},
        "sequence": {"value": "MAAA", "length": 4},
        "uniProtKBCrossReferences": [{"database": "Pfam", "id": "PF00001"}],
        "comments": [], "features": [], "genes": [], "references": [],
    }
    base.update(kwargs)
    return base


class TestVertebrataIsOutOfScope:
    """Adding vertebrates would DEEPEN the bias the sweep exists to correct.

    The class-A anchor set measured 51.3% Vertebrata / 38.1% Ecdysozoa on
    2026-07-20, so a vertebrate leak is not a rounding error, it is the whole
    problem coming back.
    """

    def test_vertebrate_lineage_is_detected(self):
        assert sweep.is_vertebrate(entry(
            ["Eukaryota", "Metazoa", "Chordata", "Craniata", "Vertebrata"]))

    def test_invertebrate_chordate_is_not_a_vertebrate(self):
        # Branchiostoma is Chordata but NOT Vertebrata. A substring test on
        # 'Chordata' would wrongly reject amphioxus, which is one of the
        # clades the sweep most wants.
        assert not sweep.is_vertebrate(entry(
            ["Eukaryota", "Metazoa", "Chordata", "Cephalochordata",
             "Leptocardii", "Amphioxiformes"]))


class TestCladeAttribution:
    def test_entry_is_attributed_to_its_own_lineage(self):
        record = entry(["Eukaryota", "Metazoa", "Cnidaria", "Hydrozoa"])
        assert sweep.clade_of(record, "Cnidaria") == "Cnidaria"

    def test_unresolvable_lineage_falls_back_to_the_queried_clade(self):
        # An entry can never be silently attributed to a clade it was not
        # found under.
        record = entry(["Eukaryota", "Metazoa"])
        assert sweep.clade_of(record, "Hemichordata") == "Hemichordata"

    def test_ecdysozoan_phyla_are_resolved_below_the_clade(self):
        # 'Ecdysozoa: 953' would hide whether the yield is all Drosophila
        # again. Measured: Arthropoda 932, Onychophora 9, Nematoda 7,
        # Tardigrada 5.
        assert sweep.phylum_within_clade(entry(
            ["Eukaryota", "Metazoa", "Ecdysozoa", "Nematoda"])) == "Nematoda"
        assert sweep.phylum_within_clade(entry(
            ["Eukaryota", "Metazoa", "Ecdysozoa", "Panarthropoda",
             "Tardigrada"])) == "Tardigrada"


class TestFamilyExtensionLabelsOnly:
    """The extension exists because refexp2's vocabulary is measurably blind to
    this universe, NOT because the gate should admit more."""

    @pytest.mark.parametrize("name", [
        "Ultraviolet-sensitive visual pigment",      # 24 survivors
        "Blue-sensitive visual pigment",             # 18 survivors
        "Long wavelength-sensitive visual pigment",  # 15 survivors
    ])
    def test_visual_pigment_is_an_opsin(self, name):
        """refexp2's needle is \\bopsin\\b|\\brhodopsin, and the entire
        arthropod visual literature says 'visual pigment' instead. 58 survivors
        were parked in 'unclassified' by that gap."""
        assert gate.family_for_report(name) == "unclassified"
        assert sweep.family_for_report_extended(name) == "opsin"

    @pytest.mark.parametrize("name", [
        "SIFamide receptor", "Myosuppressin receptor", "AstC receptor",
        "Ecdysis triggering hormone receptor", "Inotocin receptor",
        "Sex peptide receptor", "Cardioacceleratory peptide receptor",
    ])
    def test_ecdysozoan_peptide_vocabulary_is_recognised(self, name):
        """refexp2's peptide needle is a vertebrate/lophotrochozoan vocabulary
        that simply does not contain the ecdysozoan peptide names."""
        assert gate.family_for_report(name) == "unclassified"
        assert sweep.family_for_report_extended(name) == "peptide"

    def test_lgr_is_a_glycoprotein_hormone_receptor(self):
        name = "Leucine-rich repeat-containing G-protein-coupled receptor 2"
        assert sweep.family_for_report_extended(name) == "glycoprotein-hormone"

    def test_refexp2_label_always_wins(self):
        """Deferring to refexp2 keeps the two sweeps comparable. If the
        extension could override, the same accession could carry different
        families in refexp2 and refexp3 output and every cross-run join would
        be quietly wrong."""
        name = "Rhodopsin"
        assert gate.family_for_report(name) == "opsin"
        assert sweep.family_for_report_extended(name) == gate.family_for_report(name)

    def test_extension_does_not_invent_a_family_for_a_bare_label(self):
        assert sweep.family_for_report_extended(
            "G-protein coupled receptors family 1 profile domain-containing "
            "protein") == "unclassified"


class TestResidualPlaceholderLeak:
    """10 of 1133 survivors carry a name committing to no ligand and no
    function. refexp2's placeholder patterns were tuned on 'Orphan GPCR N' and
    'GCR###' and never saw a UniProt automatic name.

    They are FLAGGED, never dropped: silently dropping them would be refexp3
    editing the approved gate's verdict from the outside.
    """

    @pytest.mark.parametrize("name", [
        "G-protein coupled receptors family 1 profile domain-containing protein",
        "7 transmembrane protein",
        "Peptide receptor GPCR",
        "GH23382p",
    ])
    def test_real_leaked_names_are_flagged(self, name):
        assert sweep.is_residual_placeholder(name)
        # ...and the approved gate genuinely does admit them. If this ever
        # fails, the gate was tightened and the flag is now dead code.
        assert gate.name_asserts_function(name)

    @pytest.mark.parametrize("name", [
        "Ultraviolet-sensitive visual pigment",
        "SIFamide receptor",
        "Rhodopsin",
        "Long wavelength-sensitive visual pigment",
    ])
    def test_genuine_names_are_not_flagged(self, name):
        assert not sweep.is_residual_placeholder(name)


class TestPaperBurden:
    """The sweep's most important qualifier. One beetle opsin study contributes
    198 entries; a focused characterization contributes one."""

    def test_burden_counts_entries_sharing_a_paper(self):
        survivors = [
            {"accession": "A", "primary_pmids": "28127058"},
            {"accession": "B", "primary_pmids": "28127058"},
            {"accession": "C", "primary_pmids": "30735632"},
        ]
        sweep.annotate_paper_burden(survivors)
        assert survivors[0]["min_paper_entry_count"] == 2
        assert survivors[2]["min_paper_entry_count"] == 1

    def test_minimum_credits_the_focused_study(self):
        """An entry appearing in BOTH a 198-entry comparative survey and a
        single-receptor study is credited to the focused one."""
        survivors = [
            {"accession": "A", "primary_pmids": "bulk;focused"},
            {"accession": "B", "primary_pmids": "bulk"},
            {"accession": "C", "primary_pmids": "bulk"},
        ]
        sweep.annotate_paper_burden(survivors)
        assert survivors[0]["min_paper_entry_count"] == 1

    def test_no_pmid_yields_no_burden_rather_than_zero(self):
        survivors = [{"accession": "A", "primary_pmids": ""}]
        sweep.annotate_paper_burden(survivors)
        assert survivors[0]["min_paper_entry_count"] == ""


class TestUniverseIsWiderButTheGateIsNot:
    """refexp3 searches (family string OR PF00001 OR IPR000276) where refexp2
    searched the family string alone. Measured on 2026-07-20, Lophotrochozoa
    PE1/PE2 returns 381 by family string and 572 by the union: the curated
    string is a SIMILARITY comment that TrEMBL entries often lack.

    Widening what is LOOKED AT must not widen what is ADMITTED.
    """

    def test_universe_query_covers_all_three_routes(self):
        query = sweep.universe_query(6073)
        assert gate.CURATED_CLASS_A_FAMILY in query
        assert gate.PFAM_7TM_1 in query
        assert gate.INTERPRO_RHODOPSIN in query
        assert "existence:1" in query and "existence:2" in query

    def test_diagnostic_query_drops_only_the_existence_filter(self):
        """'no entries at all' and 'entries exist, every one PE3' are different
        findings. Porifera measured 99 entries at any PE and 0 at PE1/PE2."""
        assert "existence" not in sweep.diagnostic_query(6040)
        assert "taxonomy_id:6040" in sweep.diagnostic_query(6040)

    def test_chemosensory_clan_is_still_rejected_by_the_gate(self):
        """PF08395 / PF02949 are clan CL0176 Chemosens_recp, not CL0192. The
        wider search must not let a chemosensory-clan hit through."""
        record = entry(["Eukaryota", "Metazoa", "Ecdysozoa", "Arthropoda"],
                       uniProtKBCrossReferences=[
                           {"database": "Pfam", "id": "PF02949"}])
        is_class_a, basis, evidence = gate.verify_class_a(record)
        assert not is_class_a
        assert "CL0176" in evidence
        assert basis == "chemosensory-clan"

    def test_reporting_layer_cannot_flip_a_verdict(self):
        """The single most important invariant in this file: refexp3 labels,
        refexp2 decides."""
        record = entry(["Eukaryota", "Metazoa", "Cnidaria"],
                       proteinDescription={"submissionNames": [
                           {"fullName": {"value": "7 transmembrane protein"}}]})
        verdict = gate.gate_entry(record)
        before = verdict["failed_at"]
        verdict["family_extended"] = sweep.family_for_report_extended(
            "7 transmembrane protein")
        verdict["residual_placeholder"] = (
            "yes" if sweep.is_residual_placeholder("7 transmembrane protein")
            else "no")
        assert verdict["failed_at"] == before


class TestCladeTaxidsArePinned:
    """A wrong taxid here would silently sweep the wrong clade and report a
    confident number for it. This project has already been bitten by exactly
    that: two 'LSE clade' taxids turned out to be bacteria.

    verify_clade_taxids re-queries UniProt at runtime; this test pins the
    values so a careless edit is caught without a network round-trip.
    """

    def test_expected_names_and_taxids(self):
        pinned = {
            "Ecdysozoa": ("Ecdysozoa", 1206794),
            "Cnidaria": ("Cnidaria", 6073),
            "Placozoa": ("Placozoa", 10226),
            "Porifera": ("Porifera", 6040),
            "Ctenophora": ("Ctenophora", 10197),
            "Echinodermata": ("Echinodermata", 7586),
            "Hemichordata": ("Hemichordata", 10219),
            "Cephalochordata": ("Cephalochordata", 7735),
            "Tunicata": ("Tunicata", 7712),
            "Xenacoelomorpha": ("Xenacoelomorpha", 1312402),
        }
        assert {label: (expected, taxid)
                for label, expected, taxid in sweep.CLADES} == pinned

    def test_vertebrata_is_not_in_the_sweep_list(self):
        assert "Vertebrata" not in {label for label, _, _ in sweep.CLADES}

    def test_control_clade_is_one_that_is_actually_swept(self):
        """A zero from Porifera is only meaningful if the control proves the
        query template still resolves."""
        assert sweep.CONTROL_CLADE in {label for label, _, _ in sweep.CLADES}
