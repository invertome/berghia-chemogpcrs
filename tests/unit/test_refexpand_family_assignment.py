"""Guards for the reference-expansion family census.

Every protein name and gene-name string in this module was read programmatically
out of a live UniProt response during the 2026-07-20 reference-expansion census
and is reproduced verbatim. They are real records, not invented fixtures — a
synthetic fixture encodes whatever key the author assumed and therefore cannot
catch the case where the assumption itself is wrong, which is precisely the
failure mode these tests exist to catch.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))

import census_reviewed_classA_gpcr_families as census  # noqa: E402


class TestGeneAliasTraps:
    """Legacy gene aliases name a family the receptor does not belong to.

    Weighting gene names equally with the curated protein name silently
    misfiles every one of these. Each case below was confirmed against the
    live UniProt record.
    """

    def test_ltb4r_is_lipid_despite_cmkrl1_chemokine_alias(self):
        # Q15722, Homo sapiens. Alias CMKRL1 = "chemokine receptor-like 1".
        assert census.classify_family(
            "Leukotriene B4 receptor 1",
            "LTB4R BLT BLT1 BLTR CMKRL1 GPR16",
        ) == "lipid"

    def test_lpar6_is_lipid_despite_p2ry5_nucleotide_alias(self):
        # P43657, Homo sapiens. Alias P2RY5 reads as a purinergic receptor.
        assert census.classify_family(
            "Lysophosphatidic acid receptor 6", "LPAR6 P2RY5",
        ) == "lipid"

    def test_lpar4_is_lipid_despite_p2ry9_nucleotide_alias(self):
        # Q99677, Homo sapiens.
        assert census.classify_family(
            "Lysophosphatidic acid receptor 4", "LPAR4 GPR23 LPA4 P2RY9",
        ) == "lipid"

    def test_mouse_lpar6_ortholog_also_resists_the_alias(self):
        # Q8BMC0, Mus musculus — lowercase alias spelling must not leak through.
        assert census.classify_family(
            "Lysophosphatidic acid receptor 6", "Lpar6 P2ry5 P2y5",
        ) == "lipid"


class TestGonadotropinSubstringTrap:
    """'Gonadotropin-releasing hormone receptor' is a PEPTIDE receptor.

    The glycoprotein hormones are FSH, LH, TSH and CG. A substring match on
    'gonadotropin' pulls GnRH receptors into glycoprotein-hormone. Two
    Platynereis entries are misfiled this way in the current anchor set.
    """

    def test_gnrh_receptor_is_not_glycoprotein_hormone(self):
        # Q2V2K5, Octopus vulgaris — the only molluscan entry at risk here.
        assert census.classify_family(
            "Gonadotropin-releasing hormone receptor", "",
        ) != "glycoprotein-hormone"

    def test_platynereis_gnrh_akh_receptor_is_not_glycoprotein_hormone(self):
        # A0A6B9MSD4 — currently filed as glycoprotein-hormone in the anchor set.
        assert census.classify_family(
            "Gonadotropin releasing hormone receptor 1/adipokinetic hormone "
            "receptor 1", "",
        ) != "glycoprotein-hormone"

    @pytest.mark.parametrize("name,genes", [
        ("Thyrotropin receptor", "TSHR LGR3"),
        ("Lutropin-choriogonadotropic hormone receptor", "LHCGR LGR2 LHR"),
        ("Follicle-stimulating hormone receptor", "FSHR LGR1"),
        ("Probable glycoprotein hormone G protein-coupled receptor", ""),
    ])
    def test_true_glycoprotein_hormone_receptors_still_match(self, name, genes):
        assert census.classify_family(name, genes) == "glycoprotein-hormone"


class TestFamilyAssignment:
    @pytest.mark.parametrize("name", [
        "C-C chemokine receptor type 5",
        "C-X-C chemokine receptor type 4",
        "Atypical chemokine receptor 1",
        "CX3C chemokine receptor 1",
    ])
    def test_chemokine(self, name):
        assert census.classify_family(name, "") == "chemokine"

    @pytest.mark.parametrize("name", [
        "Adenosine receptor A1",
        "Adenosine receptor A2b",
        "P2Y purinoceptor 1",
        "Putative P2Y purinoceptor 10",
    ])
    def test_nucleotide(self, name):
        assert census.classify_family(name, "") == "nucleotide"

    @pytest.mark.parametrize("name", [
        "Prostaglandin E2 receptor EP1 subtype",
        "Cannabinoid receptor 1",
        "Sphingosine 1-phosphate receptor 1",
        "Thromboxane A2 receptor",
    ])
    def test_lipid(self, name):
        assert census.classify_family(name, "") == "lipid"

    def test_unrelated_receptor_matches_nothing(self):
        # Q16950, Aplysia californica — aminergic, outside the four families.
        assert census.classify_family(
            "5-hydroxytryptamine receptor 1", "") == ""

    def test_grl101_is_not_silently_called_glycoprotein_hormone(self):
        """P46023, Lymnaea stagnalis.

        UniProt assigns GRL101 the BARE 'G protein-coupled receptor 1 family'
        with no FSH/LSH/TSH subfamily, so promoting it to glycoprotein-hormone
        would be an inference rather than a curated fact. It is the only
        molluscan entry anywhere near the four target families, so the
        temptation to claim it is real — and it must stay a deliberate,
        user-approved call rather than a silent regex side-effect.
        """
        assert census.classify_family(
            "G protein-coupled receptor GRL101", "") == ""


class TestPhylumAssignment:
    """Lineage strings below are verbatim prefixes of real UniProt lineages."""

    def test_mollusca_detected(self):
        lineage = ("Eukaryota, Metazoa, Spiralia, Lophotrochozoa, Mollusca, "
                   "Gastropoda, Heterobranchia")
        assert census.phylum_of(lineage) == "Mollusca"
        assert census.is_lophotrochozoan(lineage)

    def test_annelida_is_lophotrochozoan(self):
        lineage = ("Eukaryota, Metazoa, Spiralia, Lophotrochozoa, Annelida, "
                   "Clitellata")
        assert census.phylum_of(lineage) == "Annelida"
        assert census.is_lophotrochozoan(lineage)

    def test_vertebrate_is_not_lophotrochozoan(self):
        lineage = "Eukaryota, Metazoa, Chordata, Craniata, Vertebrata, Mammalia"
        assert census.phylum_of(lineage) == "Vertebrata"
        assert not census.is_lophotrochozoan(lineage)

    def test_cnidaria_is_neither_vertebrate_nor_lophotrochozoan(self):
        # P35409, Anthopleura elegantissima.
        lineage = "Eukaryota, Metazoa, Cnidaria, Anthozoa, Hexacorallia"
        assert census.phylum_of(lineage) == "Cnidaria"
        assert not census.is_lophotrochozoan(lineage)

    def test_arthropod_is_not_lophotrochozoan(self):
        lineage = "Eukaryota, Metazoa, Ecdysozoa, Arthropoda, Hexapoda, Insecta"
        assert census.phylum_of(lineage) == "Arthropoda"
        assert not census.is_lophotrochozoan(lineage)


class TestEmptyResultRefusal:
    """An empty UniProt result must never reach an output artifact.

    A renamed keyword returns HTTP 200 with an empty body, which already
    overwrote this project's reference set with headers only, exit 0.
    """

    def test_write_tsv_refuses_empty_rows(self, tmp_path):
        target = tmp_path / "census.tsv"
        with pytest.raises(RuntimeError, match="refusing to write an empty"):
            census.write_tsv([], str(target))
        assert not target.exists()

    def test_write_tsv_does_not_clobber_on_empty(self, tmp_path):
        target = tmp_path / "census.tsv"
        target.write_text("family\tpre-existing\n")
        with pytest.raises(RuntimeError):
            census.write_tsv([], str(target))
        assert target.read_text() == "family\tpre-existing\n"

    def test_fetch_paginated_refuses_empty_response(self, monkeypatch):
        header = "\t".join([
            "Entry", "Reviewed", "Protein names", "Organism",
            "Organism (ID)", "Gene Names", "Length", "Protein families",
            "Taxonomic lineage",
        ])
        monkeypatch.setattr(census, "_get", lambda url: (header + "\n", "", 0))
        with pytest.raises(RuntimeError, match="zero rows"):
            census.fetch_paginated("some query")

    def test_fetch_paginated_detects_lost_pagination(self, monkeypatch):
        header = "\t".join([
            "Entry", "Reviewed", "Protein names", "Organism",
            "Organism (ID)", "Gene Names", "Length", "Protein families",
            "Taxonomic lineage",
        ])
        row = "\t".join([
            "P12345", "reviewed", "C-C chemokine receptor type 5", "Homo sapiens",
            "9606", "CCR5", "352", "G protein-coupled receptor 1 family",
            "Eukaryota, Metazoa, Chordata, Vertebrata",
        ])
        # Server claims 500 results but only one row is returned.
        monkeypatch.setattr(
            census, "_get", lambda url: (header + "\n" + row + "\n", "", 500))
        with pytest.raises(RuntimeError, match="pagination lost rows"):
            census.fetch_paginated("some query")


class TestWriteOnceSafety:
    def test_census_module_never_writes_the_anchor_set(self):
        """The census is read-only w.r.t. every existing reference artifact.

        Identifiers are write-once in this project; a census that edited the
        anchor set could renumber or re-mint one. The only writer in the module
        is write_tsv, and it targets the --output path alone.
        """
        source = open(census.__file__).read()
        for forbidden in ("anchor_set_PROD.tsv\"", "anchor_set_PROD.tsv'"):
            assert forbidden not in source
        assert source.count("os.replace(") == 1
