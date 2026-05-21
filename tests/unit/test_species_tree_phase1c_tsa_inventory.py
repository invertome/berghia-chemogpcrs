"""Tests for scripts/build_species_tree_phase1c_tsa_inventory.py.

Bead -2x3 (Phase 1c of -dnk epic). Phase 1a found only 91/239 species
have public proteomes in NCBI Datasets. This phase queries NCBI nuccore
for TSA (Transcriptome Shotgun Assembly) entries on the 148 dropped
species, to decide whether Phase 1d (EvidentialGene processing of
transcriptomes into proteome proxies) is worth committing to.

Confirmed working query format (manual test 2026-05-21):
  esearch -db nuccore -query "txid6359[Organism] AND srcdb_genbank[Properties] AND TSA[Keyword]"
  -> 10928 TSA records for Platynereis dumerilii, master GBZT00000000
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

import build_species_tree_phase1c_tsa_inventory as tsa


# ----------------------------------------------------------------------
# read_dropped_taxa_from_manifest
# ----------------------------------------------------------------------

class TestReadDroppedTaxaFromManifest:
    """Phase 1c input is the Phase 1a proteome_manifest.tsv. We filter
    to rows where drop_reason is non-empty (no proteome found) and feed
    those taxa to the TSA queries.
    """

    def _write_manifest(self, path: Path, rows: list[dict]) -> None:
        cols = (
            "taxid", "binomial", "clade", "source", "accession",
            "assembly_level", "annotation_status", "est_protein_count",
            "submission_date", "drop_reason",
        )
        with path.open("w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(r.get(c, "") for c in cols) + "\n")

    def test_returns_only_dropped_rows(self, tmp_path: Path) -> None:
        m = tmp_path / "proteome_manifest.tsv"
        self._write_manifest(m, [
            {"taxid": "6500", "binomial": "Aplysia californica",
             "clade": "gastropoda", "source": "RefSeq",
             "accession": "GCF_000002075.1", "drop_reason": ""},
            {"taxid": "6359", "binomial": "Platynereis dumerilii",
             "clade": "annelida", "drop_reason": "no_proteome_in_ncbi"},
            {"taxid": "6426", "binomial": "Riftia pachyptila",
             "clade": "annelida", "drop_reason": "no_proteome_in_ncbi"},
        ])
        dropped = tsa.read_dropped_taxa_from_manifest(m)
        assert len(dropped) == 2
        taxids = {t.taxid for t in dropped}
        assert taxids == {6359, 6426}
        # Should preserve binomial + clade
        p = next(t for t in dropped if t.taxid == 6359)
        assert p.binomial == "Platynereis dumerilii"
        assert p.clade == "annelida"

    def test_empty_manifest_returns_empty_list(self, tmp_path: Path) -> None:
        m = tmp_path / "proteome_manifest.tsv"
        self._write_manifest(m, [])
        assert tsa.read_dropped_taxa_from_manifest(m) == []

    def test_all_successful_returns_empty_list(self, tmp_path: Path) -> None:
        m = tmp_path / "proteome_manifest.tsv"
        self._write_manifest(m, [
            {"taxid": "6500", "binomial": "Aplysia californica",
             "clade": "gastropoda", "source": "RefSeq",
             "accession": "GCF_000002075.1", "drop_reason": ""},
        ])
        assert tsa.read_dropped_taxa_from_manifest(m) == []

    def test_missing_manifest_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            tsa.read_dropped_taxa_from_manifest(tmp_path / "nope.tsv")


# ----------------------------------------------------------------------
# derive_master_prefix
# ----------------------------------------------------------------------

class TestDeriveMasterPrefix:
    """TSA contig accessions follow WGS format: 4-char-prefix + 2 digits
    + 6 digits (e.g., GBZT01000001). The master record uses the same
    4-char prefix + 8 zeros (GBZT00000000). For the inventory we just
    want the 4-char prefix as a stable group identifier.
    """

    def test_standard_tsa_contig_accession(self) -> None:
        assert tsa.derive_master_prefix("GBZT01000001") == "GBZT"

    def test_master_accession(self) -> None:
        assert tsa.derive_master_prefix("GBZT00000000") == "GBZT"

    def test_versioned_accession(self) -> None:
        # NCBI sometimes returns versioned accessions
        assert tsa.derive_master_prefix("GBZT01000001.1") == "GBZT"

    def test_short_or_malformed_returns_input(self) -> None:
        # Defensive: if accession doesn't look like WGS/TSA, return it
        # verbatim so the inventory still surfaces something useful
        # rather than erroring.
        assert tsa.derive_master_prefix("FOO") == "FOO"
        assert tsa.derive_master_prefix("") == ""


# ----------------------------------------------------------------------
# query_tsa_for_taxon (Bio.Entrez wrapper)
# ----------------------------------------------------------------------

class TestQueryTsaForTaxon:
    """Tests for the Bio.Entrez query wrapper. Real NCBI calls happen
    only on Unity; we test the call construction + result parsing with
    Entrez.esearch and Entrez.read mocked.
    """

    def _patch_entrez(self, monkeypatch, count: int, idlist: list[str]) -> dict:
        """Install fake Entrez.esearch + Entrez.read; return captured args."""
        captured: dict = {}

        def fake_esearch(**kwargs):
            captured["esearch_kwargs"] = kwargs
            return "FAKE_HANDLE"

        def fake_read(handle):
            assert handle == "FAKE_HANDLE"
            return {"Count": str(count), "IdList": idlist}

        # Also need a noop esummary for accession lookup
        def fake_esummary(**kwargs):
            captured["esummary_kwargs"] = kwargs
            return "FAKE_SUMMARY_HANDLE"

        def fake_parse(handle):
            assert handle == "FAKE_SUMMARY_HANDLE"
            return [
                {"Id": gid, "AccessionVersion": f"GBZT0100000{i+1}"}
                for i, gid in enumerate(idlist)
            ]

        monkeypatch.setattr(tsa.Entrez, "esearch", fake_esearch)
        monkeypatch.setattr(tsa.Entrez, "esummary", fake_esummary)
        monkeypatch.setattr(tsa.Entrez, "read", fake_read)
        monkeypatch.setattr(tsa.Entrez, "parse", fake_parse)
        return captured

    def test_query_format(self, monkeypatch) -> None:
        captured = self._patch_entrez(monkeypatch, count=100, idlist=["123"])
        tsa.query_tsa_for_taxon(6359)
        term = captured["esearch_kwargs"]["term"]
        assert "txid6359[Organism]" in term
        assert "TSA[Keyword]" in term
        assert "srcdb_genbank[Properties]" in term
        assert captured["esearch_kwargs"]["db"] == "nuccore"

    def test_zero_hits(self, monkeypatch) -> None:
        self._patch_entrez(monkeypatch, count=0, idlist=[])
        result = tsa.query_tsa_for_taxon(6426)
        assert result.has_tsa is False
        assert result.n_records == 0
        assert result.master_prefix == ""
        assert result.first_accession == ""

    def test_many_hits_derives_master_prefix(self, monkeypatch) -> None:
        self._patch_entrez(monkeypatch, count=10928, idlist=["123", "124"])
        result = tsa.query_tsa_for_taxon(6359)
        assert result.has_tsa is True
        assert result.n_records == 10928
        assert result.master_prefix == "GBZT"
        assert result.first_accession.startswith("GBZT")

    def test_entrez_email_and_api_key_set(self, monkeypatch) -> None:
        # The wrapper should set Entrez.email and Entrez.api_key if env
        # vars are provided. Setting these correctly is what lets NCBI
        # allow >3 req/s.
        self._patch_entrez(monkeypatch, count=0, idlist=[])
        tsa.query_tsa_for_taxon(
            6359, email="me@example.com", api_key="abc123"
        )
        assert tsa.Entrez.email == "me@example.com"
        assert tsa.Entrez.api_key == "abc123"

    def test_handles_entrez_failure(self, monkeypatch) -> None:
        # If esearch raises, the wrapper should propagate so the
        # orchestrator can mark the row with a query_error drop_reason.
        def boom(**kwargs):
            raise RuntimeError("NCBI 503")
        monkeypatch.setattr(tsa.Entrez, "esearch", boom)
        with pytest.raises(RuntimeError, match="NCBI"):
            tsa.query_tsa_for_taxon(6359)


# ----------------------------------------------------------------------
# write_tsa_inventory_tsv
# ----------------------------------------------------------------------

class TestWriteTsaInventoryTsv:
    """Output schema (locked):
      taxid | binomial | clade | has_tsa | n_tsa_records |
      tsa_master_prefix | tsa_first_accession | query_error
    """

    def test_columns_match(self, tmp_path: Path) -> None:
        out = tmp_path / "tsa_inventory.tsv"
        tsa.write_tsa_inventory_tsv(out, [])
        header = out.read_text().splitlines()[0].split("\t")
        assert header == [
            "taxid",
            "binomial",
            "clade",
            "has_tsa",
            "n_tsa_records",
            "tsa_master_prefix",
            "tsa_first_accession",
            "query_error",
        ]

    def test_successful_row(self, tmp_path: Path) -> None:
        entry = tsa.TsaInventoryEntry(
            taxon=tsa.DroppedTaxon(taxid=6359, binomial="Platynereis dumerilii", clade="annelida"),
            result=tsa.TsaQueryResult(
                has_tsa=True, n_records=10928,
                master_prefix="GBZT", first_accession="GBZT01000001",
            ),
            query_error="",
        )
        out = tmp_path / "tsa.tsv"
        tsa.write_tsa_inventory_tsv(out, [entry])
        row = out.read_text().splitlines()[1].split("\t")
        assert row[0] == "6359"
        assert row[3] == "yes"
        assert row[4] == "10928"
        assert row[5] == "GBZT"
        assert row[6] == "GBZT01000001"
        assert row[7] == ""

    def test_no_tsa_row(self, tmp_path: Path) -> None:
        entry = tsa.TsaInventoryEntry(
            taxon=tsa.DroppedTaxon(taxid=6426, binomial="Riftia pachyptila", clade="annelida"),
            result=tsa.TsaQueryResult(
                has_tsa=False, n_records=0,
                master_prefix="", first_accession="",
            ),
            query_error="",
        )
        out = tmp_path / "tsa.tsv"
        tsa.write_tsa_inventory_tsv(out, [entry])
        row = out.read_text().splitlines()[1].split("\t")
        assert row[0] == "6426"
        assert row[3] == "no"
        assert row[4] == "0"

    def test_query_error_row(self, tmp_path: Path) -> None:
        entry = tsa.TsaInventoryEntry(
            taxon=tsa.DroppedTaxon(taxid=9999, binomial="Foo bar", clade="annelida"),
            result=None,
            query_error="NCBI 503",
        )
        out = tmp_path / "tsa.tsv"
        tsa.write_tsa_inventory_tsv(out, [entry])
        row = out.read_text().splitlines()[1].split("\t")
        assert row[0] == "9999"
        assert row[3] == ""           # has_tsa empty when query failed
        assert row[7] == "NCBI 503"


# ----------------------------------------------------------------------
# build_tsa_inventory orchestrator
# ----------------------------------------------------------------------

class TestBuildTsaInventory:
    """End-to-end: read manifest -> filter to dropped -> query each ->
    write tsa_inventory.tsv -> return summary.
    """

    def _write_manifest(self, path: Path) -> None:
        cols = ("taxid", "binomial", "clade", "source", "accession",
                "assembly_level", "annotation_status", "est_protein_count",
                "submission_date", "drop_reason")
        rows = [
            # successful entry - skipped by Phase 1c
            ("6500", "Aplysia californica", "gastropoda", "RefSeq",
             "GCF_000002075.1", "Scaffold", "Current", "17654", "2024-01-01", ""),
            # dropped - has TSA
            ("6359", "Platynereis dumerilii", "annelida", "", "", "", "", "", "",
             "no_proteome_in_ncbi"),
            # dropped - no TSA
            ("6426", "Riftia pachyptila", "annelida", "", "", "", "", "", "",
             "no_proteome_in_ncbi"),
        ]
        with path.open("w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(r) + "\n")

    def test_happy_path(self, tmp_path: Path) -> None:
        manifest = tmp_path / "proteome_manifest.tsv"
        self._write_manifest(manifest)
        out = tmp_path / "tsa_inventory.tsv"

        def fake_query(taxid: int) -> tsa.TsaQueryResult:
            if taxid == 6359:
                return tsa.TsaQueryResult(
                    has_tsa=True, n_records=10928,
                    master_prefix="GBZT", first_accession="GBZT01000001",
                )
            return tsa.TsaQueryResult(
                has_tsa=False, n_records=0, master_prefix="", first_accession="",
            )

        summary = tsa.build_tsa_inventory(
            manifest_path=manifest, out_path=out, query_fn=fake_query,
        )
        # 6500 skipped (had proteome); 6359 + 6426 queried
        assert summary.n_dropped_total == 2
        assert summary.n_with_tsa == 1
        assert summary.n_without_tsa == 1
        # TSV has both rows
        rows = out.read_text().splitlines()[1:]
        assert len(rows) == 2

    def test_query_error_recorded_not_raised(self, tmp_path: Path) -> None:
        manifest = tmp_path / "proteome_manifest.tsv"
        self._write_manifest(manifest)
        out = tmp_path / "tsa_inventory.tsv"

        def flaky_query(taxid: int) -> tsa.TsaQueryResult:
            if taxid == 6426:
                raise RuntimeError("simulated NCBI 503")
            return tsa.TsaQueryResult(has_tsa=False, n_records=0,
                                       master_prefix="", first_accession="")

        summary = tsa.build_tsa_inventory(
            manifest_path=manifest, out_path=out, query_fn=flaky_query,
        )
        assert summary.n_query_errors == 1
        # Row appears with query_error populated
        rows = out.read_text().splitlines()[1:]
        riftia_row = next(r for r in rows if r.startswith("6426\t"))
        assert "simulated NCBI 503" in riftia_row
