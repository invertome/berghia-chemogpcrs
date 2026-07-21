"""Guards for the evidence-gated reference expansion (refexp2).

Every protein name, family string, ECO code and paper title in this module was
read programmatically out of a live UniProt or PubMed response during the
2026-07-20 expansion run and is reproduced verbatim. They are real records, not
invented fixtures: a synthetic fixture encodes whatever the author assumed and
therefore cannot catch the case where the assumption itself is wrong, which is
exactly the failure mode these tests exist to catch.

The gate's whole value is its strictness. Every test below pins a specific way
the gate could silently loosen and start admitting homology-inferred entries.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))

import refexp2_evidence_gate as gate  # noqa: E402


def entry(**kwargs):
    """A minimal UniProt JSON entry with real-shaped defaults."""
    base = {
        "primaryAccession": "A0A0K0PUF3",
        "entryType": "UniProtKB unreviewed (TrEMBL)",
        "proteinExistence": "2: Evidence at transcript level",
        "proteinDescription": {},
        "organism": {"scientificName": "Platynereis dumerilii", "taxonId": 6359,
                     "lineage": ["Eukaryota", "Metazoa", "Spiralia",
                                 "Lophotrochozoa", "Annelida"]},
        "sequence": {"value": "MAAA", "length": 4},
        "uniProtKBCrossReferences": [
            {"database": "Pfam", "id": "PF00001"},
            {"database": "InterPro", "id": "IPR000276"},
        ],
        "comments": [], "features": [], "genes": [], "references": [],
    }
    base.update(kwargs)
    return base


def named(name, **kwargs):
    return entry(proteinDescription={
        "submissionNames": [{"fullName": {"value": name}}]}, **kwargs)


def citation(title, pmid=None, citation_type="journal article"):
    reference = {"citation": {"citationType": citation_type, "title": title}}
    if pmid:
        reference["citation"]["citationCrossReferences"] = [
            {"database": "PubMed", "id": pmid}]
    return reference


class TestProteinExistence:
    """PE3 (homology) and PE4 (predicted) ARE the automated inference the gate
    exists to exclude. 5554 of 5700 molluscan class-A entries are PE3."""

    @pytest.mark.parametrize("raw,expected", [
        ("1: Evidence at protein level", 1),
        ("2: Evidence at transcript level", 2),
        ("3: Inferred from homology", 3),
        ("4: Predicted", 4),
        ("5: Uncertain", 5),
    ])
    def test_parses_every_level(self, raw, expected):
        assert gate.protein_existence_level(entry(proteinExistence=raw)) == expected

    def test_unparseable_raises_rather_than_defaulting(self):
        # Defaulting to a passing level would admit the entire PE3 bulk.
        with pytest.raises(ValueError):
            gate.protein_existence_level(entry(proteinExistence="unknown"))

    @pytest.mark.parametrize("raw", ["3: Inferred from homology", "4: Predicted",
                                     "5: Uncertain"])
    def test_pe3_pe4_pe5_are_rejected_by_the_gate(self, raw):
        verdict = gate.gate_entry(named("Melatonin receptor", proteinExistence=raw))
        assert verdict["failed_at"] == "existence"
        assert not verdict["pass_existence"]


class TestEvidenceCodes:
    def test_experimental_code_detected(self):
        record = entry(comments=[{"commentType": "FUNCTION",
                                  "evidences": [{"evidenceCode": "ECO:0000269"}]}])
        assert gate.has_experimental_evidence(record)
        assert not gate.rests_only_on_automated_evidence(record)

    @pytest.mark.parametrize("code", ["ECO:0000256", "ECO:0000259", "ECO:0000313"])
    def test_automated_codes_are_not_experimental(self, code):
        record = entry(comments=[{"commentType": "SIMILARITY",
                                  "evidences": [{"evidenceCode": code}]}])
        assert not gate.has_experimental_evidence(record)
        assert gate.rests_only_on_automated_evidence(record)

    def test_reviewed_entry_without_experimental_code_fails(self):
        """A reviewed entry HAS a curator, so the curator code is required."""
        record = named("Melatonin receptor",
                       entryType="UniProtKB reviewed (Swiss-Prot)",
                       comments=[{"commentType": "SIMILARITY",
                                  "evidences": [{"evidenceCode": "ECO:0000256"}]}],
                       references=[citation("Melatonin signaling controls circadian "
                                            "swimming behavior", "25259919")])
        assert gate.gate_entry(record)["failed_at"] == "evidence_code"

    def test_trembl_cannot_supply_eco_269_so_it_defers_to_literature(self):
        """MEASURED: 0 of 363 unreviewed lophotrochozoan PE1/PE2 class-A entries
        carry ECO:0000269. Requiring it of TrEMBL would make the whole
        relaxation a no-op, so S2 must not be the stage that kills them."""
        record = named("Melatonin receptor",
                       comments=[{"commentType": "SIMILARITY",
                                  "evidences": [{"evidenceCode": "ECO:0000256"}]}],
                       references=[citation("Melatonin signaling controls circadian "
                                            "swimming behavior", "25259919")])
        verdict = gate.gate_entry(record)
        assert verdict["pass_evidence_code"]
        assert verdict["failed_at"] == ""


class TestClassAVerification:
    """Class membership is verified from the record, never inferred from a name."""

    def test_rhodopsin_clan_domain_verifies(self):
        ok, evidence = gate.verify_class_a(entry())
        assert ok
        assert "PF00001" in evidence

    def test_chemosensory_clan_is_not_class_a(self):
        """PF08395 / PF02949 are Pfam clan CL0176 Chemosens_recp, NOT the
        rhodopsin clan CL0192. A Pfam hit alone does not establish class A."""
        for pfam in ("PF08395", "PF02949"):
            record = entry(uniProtKBCrossReferences=[
                {"database": "Pfam", "id": pfam}])
            ok, reason = gate.verify_class_a(record)
            assert not ok
            assert "CL0176" in reason

    def test_curated_family_string_alone_does_not_establish_class_a(self):
        """The string is an annotation; the domain is the evidence."""
        record = entry(
            uniProtKBCrossReferences=[],
            comments=[{"commentType": "SIMILARITY", "texts": [{"value":
                "Belongs to the G-protein coupled receptor 1 family"}]}])
        ok, _ = gate.verify_class_a(record)
        assert not ok

    @pytest.mark.parametrize("name", [
        # Real anchors currently declared class A that carry no PF00001. Diuretic
        # hormone, calcitonin and PDF receptors are secretin-family CLASS B.
        "Diuretic hormone receptor",
        "Calcitonin receptor 2",
        "PDF receptor",
        "Gamma-aminobutyric acid type B receptor subunit 1",
    ])
    def test_class_b_and_c_receptors_without_pf00001_are_rejected(self, name):
        record = named(name, uniProtKBCrossReferences=[])
        assert gate.gate_entry(record)["failed_at"] == "class_a"


class TestNameGate:
    """The submitters' own naming does the discrimination."""

    @pytest.mark.parametrize("name", [
        "Orphan G-protein coupled receptor 1",
        "Orphan G-protein coupled receptor 34",
        "Orphan G-protein coupled receptor 65",
    ])
    def test_orphan_named_receptors_are_rejected(self, name):
        """All 47 new Platynereis entries from the large-scale deorphanization
        study are named this way BY THE AUTHORS -- they are the receptors that
        were tested and not deorphanized."""
        assert not gate.name_asserts_function(name)

    @pytest.mark.parametrize("name", ["GCR002", "GCR484", "GCR107"])
    def test_planarian_screen_placeholders_are_rejected(self, name):
        assert not gate.name_asserts_function(name)

    @pytest.mark.parametrize("name", ["NPYR-10", "NPYR-2", "NPYR-16"])
    def test_abbreviation_plus_index_is_rejected(self, name):
        """A family label assigned by homology and then numbered. These 14
        Schmidtea entries cleared every other stage; admitting them would have
        been the 'lower the bar to produce a number' failure."""
        assert not gate.name_asserts_function(name)

    @pytest.mark.parametrize("name", [
        "Putative 5-Hydroxytryptamine receptor 3",
        "Probable G protein-coupled receptor 160",
        "Serotonin receptor-like planarian receptor 1",
    ])
    def test_hedged_and_homology_named_entries_are_rejected(self, name):
        assert not gate.name_asserts_function(name)

    @pytest.mark.parametrize("name", [
        "G-protein coupled receptor", "GPCR", "G protein coupled receptor",
        "",
    ])
    def test_bare_labels_carry_no_ligand(self, name):
        assert not gate.name_asserts_function(name)

    @pytest.mark.parametrize("name", [
        "Adipokinetic hormone receptor 1A",
        "Allatotropin receptor 1",
        "FLRFamide receptor",
        "Melatonin receptor",
        "5-HT1 receptor",
        "Serotonin receptor 5-HT7",
        "Rhabdomeric opsin 3",
        "Xenopsin",
    ])
    def test_real_ligand_bearing_names_survive(self, name):
        """Over-tightening is a failure too: these are the genuine additions."""
        assert gate.name_asserts_function(name)


class TestLiteratureGate:
    def test_embl_submission_is_not_primary_literature(self):
        """Real trap: submissions carry paper-shaped titles. Q6HA06's reads
        'Identification and characterization of a glycoprotein hormone receptor
        from the oyster Crassostrea gigas' but it is a database deposit with no
        PubMed ID, no journal and no peer review -- and PubMed confirms the work
        was never published."""
        record = named("Glycoprotein hormone receptor", references=[
            citation("Identification and characterization of a glycoprotein "
                     "hormone receptor from the oyster Crassostrea gigas.",
                     citation_type="submission")])
        assert gate.gate_entry(record)["failed_at"] == "literature"

    def test_journal_article_without_pubmed_id_is_rejected(self):
        record = named("Melatonin receptor",
                       references=[citation("Some article", pmid=None)])
        assert gate.gate_entry(record)["failed_at"] == "literature"

    @pytest.mark.parametrize("title", [
        "The Schistosoma japonicum genome reveals features of host-parasite interplay.",
        "Hirudo verbana central nervous system transcriptome analysis of ion channel "
        "and receptor content.",
        "Identification of expressed genes in cDNA library of hemocytes from the "
        "RLO-challenged oyster, Crassostrea ariakensis Gould.",
        "The Gastric Ganglion of Octopus vulgaris: Preliminary Characterization of "
        "Gene- and Putative Neurochemical-Complexity.",
    ])
    def test_sequence_survey_papers_are_rejected(self, title):
        assert gate.is_survey_paper(title)

    @pytest.mark.parametrize("title", [
        "An ancient FMRFamide-related peptide-receptor pair induces defence "
        "behaviour in a brachiopod larva.",
        "Large-Scale Combinatorial Deorphanization of Platynereis Neuropeptide GPCRs.",
        "Unique pharmacological properties of serotoninergic G-protein coupled "
        "receptors from cestodes.",
    ])
    def test_characterization_papers_are_not_survey_papers(self, title):
        assert not gate.is_survey_paper(title)

    def test_survey_detection_reads_titles_not_abstracts(self):
        """Scanning abstracts over-rejected badly: a functional study routinely
        mentions 'genome' in passing. Titles name what a paper IS."""
        record = named("FLRFamide receptor", references=[
            citation("An ancient FMRFamide-related peptide-receptor pair induces "
                     "defence behaviour in a brachiopod larva.", "28835571")])
        abstracts = {"28835571": "We searched the genome and sequenced cDNA ..."}
        assert gate.gate_entry(record, abstracts)["pass_literature"]


class TestEvidenceTier:
    def test_functional_assay_language_promotes_to_higher_tier(self):
        record = named("FLRFamide receptor", references=[
            citation("An ancient FMRFamide-related peptide-receptor pair induces "
                     "defence behaviour in a brachiopod larva.", "28835571")])
        verdict = gate.gate_entry(record, {"28835571":
            "We deorphanized the receptor and measured a dose-response."})
        assert verdict["deorphanized"]
        assert verdict["evidence_tier"] == "functionally-characterized"

    def test_named_but_unassayed_stays_in_the_lower_tier(self):
        """Clearing the gate does NOT by itself mean a ligand was applied. Most
        surviving opsins are in this tier and must not be reported as
        deorphanized."""
        record = named("Xenopsin", references=[
            citation("Co-expression of xenopsin and rhabdomeric opsin in "
                     "photoreceptors bearing microvilli and cilia.", "28876222")])
        verdict = gate.gate_entry(record, {"28876222":
            "We describe expression patterns in photoreceptor cells."})
        assert not verdict["deorphanized"]
        assert verdict["evidence_tier"] == "published-not-deorphanized"
        assert verdict["failed_at"] == ""  # still a legitimate addition


class TestFunnelMonotonicity:
    def test_survivors_never_increase_along_the_funnel(self):
        records = [
            named("Melatonin receptor", references=[
                citation("Melatonin signaling controls circadian swimming "
                         "behavior in marine zooplankton.", "25259919")]),
            named("Orphan G-protein coupled receptor 34", references=[
                citation("Large-Scale Combinatorial Deorphanization.", "26190115")]),
            named("Xenopsin", proteinExistence="3: Inferred from homology"),
        ]
        counts = gate.funnel([gate.gate_entry(r) for r in records])
        ordered = [counts["input"]] + [counts[s] for s in gate.GATE_STAGES]
        assert ordered == sorted(ordered, reverse=True)
        assert counts["input"] == 3
        assert counts["characterized"] == 1


class TestEmptyAndTruncatedResults:
    """A renamed field returns HTTP 200 with an empty body, which already
    overwrote this project's reference set with headers only, exit 0."""

    def test_write_tsv_refuses_empty_rows(self, tmp_path):
        target = tmp_path / "out.tsv"
        with pytest.raises(RuntimeError, match="refusing to write an empty"):
            gate.write_tsv([], ["accession"], str(target))
        assert not target.exists()

    def test_write_tsv_does_not_clobber_on_empty(self, tmp_path):
        target = tmp_path / "out.tsv"
        target.write_text("accession\tpre-existing\n")
        with pytest.raises(RuntimeError):
            gate.write_tsv([], ["accession"], str(target))
        assert target.read_text() == "accession\tpre-existing\n"

    def test_fetch_refuses_empty_response(self, monkeypatch):
        monkeypatch.setattr(gate, "_get", lambda url, retries=4:
                            ('{"results": []}', "", 0))
        with pytest.raises(RuntimeError, match="zero rows"):
            gate.fetch_json_paginated("some query")

    def test_fetch_detects_lost_pagination(self, monkeypatch):
        body = '{"results": [{"primaryAccession": "P12345"}]}'
        monkeypatch.setattr(gate, "_get", lambda url, retries=4: (body, "", 500))
        with pytest.raises(RuntimeError, match="pagination lost rows"):
            gate.fetch_json_paginated("some query")

    def test_unresolvable_accession_raises_never_silently_drops(self, monkeypatch):
        monkeypatch.setattr(gate, "_get", lambda url, retries=4:
                            ('{"results": []}', "", 0))
        with pytest.raises(RuntimeError, match="did not resolve"):
            gate.fetch_entries_by_accession(["P99999"])


class TestSequenceIntegrity:
    """Accumulating by accession has previously DOUBLED sequences in this project
    while passing every count check."""

    def test_accepts_matching_sequences(self):
        entries = {"P1": {"sequence": {"value": "MAAA", "length": 4}}}
        gate.verify_sequence_integrity(entries, [{"accession": "P1", "sequence": "MAAA"}])

    def test_rejects_a_doubled_sequence(self):
        entries = {"P1": {"sequence": {"value": "MAAA", "length": 4}}}
        with pytest.raises(RuntimeError, match="differs from source"):
            gate.verify_sequence_integrity(
                entries, [{"accession": "P1", "sequence": "MAAAMAAA"}])

    def test_rejects_duplicate_accessions(self):
        entries = {"P1": {"sequence": {"value": "MAAA", "length": 4}}}
        with pytest.raises(RuntimeError, match="duplicate accession"):
            gate.verify_sequence_integrity(entries, [
                {"accession": "P1", "sequence": "MAAA"},
                {"accession": "P1", "sequence": "MAAA"}])

    def test_rejects_length_disagreeing_with_declared(self):
        entries = {"P1": {"sequence": {"value": "MAAA", "length": 99}}}
        with pytest.raises(RuntimeError, match="declared"):
            gate.verify_sequence_integrity(entries, [{"accession": "P1",
                                                      "sequence": "MAAA"}])


class TestWriteOnceSafety:
    def test_module_never_writes_the_anchor_set(self):
        """Identifiers are write-once. Expansion is additive: new entries get new
        identifiers, existing ones keep theirs, gaps included."""
        source = open(gate.__file__).read()
        for forbidden in ('anchor_set_PROD.tsv"', "anchor_set_PROD.tsv'"):
            assert forbidden not in source
        # The only writer is write_tsv, targeting --output alone.
        assert source.count("os.replace(") == 2  # write_tsv + the FASTA writer

    def test_anchor_tsv_is_only_ever_opened_for_reading(self):
        source = open(gate.__file__).read()
        assert 'open(args.anchor_tsv, newline="")' in source
        assert 'open(args.anchor_tsv, "w"' not in source
