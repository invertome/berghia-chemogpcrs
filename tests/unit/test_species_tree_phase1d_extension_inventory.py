"""Tests for scripts/build_species_tree_phase1d_extension_inventory.py.

Bead -pcc (Audit C of -v1c epic, under -dnk umbrella).

Phase 1d "extension inventory" widens the species tree beyond the 239
Nath-et-al references using tight-scope policies (decision 2026-05-21):
  - Heterobranchia ALL (Berghia's subclass) — any assembly quality
  - rare basals (Polyplacophora / Scaphopoda / Caudofoveata /
    Monoplacophora) — any quality
  - other Mollusca subclasses (Bivalvia / Cephalopoda /
    Caenogastropoda / Vetigastropoda / Neritimorpha /
    Patellogastropoda) — CHROMOSOME-LEVEL only
  - outgroup phyla (Annelida / Platyhelminthes / Nemertea /
    Brachiopoda / Phoronida / Bryozoa) — CHROMOSOME-LEVEL with
    per-phylum cap (sample, not exhaustive)

The script does NOT download anything; it queries NCBI Datasets for
each clade and writes a manifest TSV at
references/species_tree/extension_inventory.tsv. Phase 1f BRAKER4 picks
that up via the same samples.csv builder.
"""
from __future__ import annotations

from pathlib import Path

import pytest

import build_species_tree_phase1d_extension_inventory as ext
from build_species_tree_phase1a_inventory import AssemblyChoice


# ----------------------------------------------------------------------
# Fake datasets-summary records (one species, one assembly each unless noted)
# ----------------------------------------------------------------------

def _record(
    accession: str,
    *,
    taxid: int,
    organism_name: str,
    source: str = "SOURCE_DATABASE_GENBANK",
    assembly_level: str = "Scaffold",
    submission_date: str = "2020-01-01",
    contig_n50: int = 50000,
    total_length: int = 1_000_000_000,
    annotated: bool = False,
    protein_count: int = 0,
) -> dict:
    rec: dict = {
        "accession": accession,
        "source_database": source,
        "organism": {"tax_id": taxid, "organism_name": organism_name},
        "assembly_info": {
            "assembly_level": assembly_level,
            "submission_date": submission_date,
        },
        "assembly_stats": {
            "contig_n50": contig_n50,
            "total_sequence_length": total_length,
        },
    }
    if annotated:
        rec["annotation_info"] = {
            "status": "Current",
            "stats": {"gene_counts": {"protein_coding": protein_count}},
        }
    return rec


# ----------------------------------------------------------------------
# load_existing_taxids
# ----------------------------------------------------------------------

class TestLoadExistingTaxids:
    """Reads existing taxids from union of:
      - nath_et_al/<clade>/*.faa filename taxid prefix
      - any number of Phase-1a-style manifest TSVs (taxid column)
    """

    def _make_refs(self, root: Path, layout: dict) -> Path:
        refs = root / "one_to_one_ortholog"
        for clade, files in layout.items():
            d = refs / clade
            d.mkdir(parents=True, exist_ok=True)
            for fn in files:
                (d / fn).write_text(">x\nA\n")
        return refs

    def _write_manifest(self, path: Path, rows: list[tuple[int, str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        lines = ["taxid\tbinomial\tdrop_reason"]
        for taxid, reason in rows:
            lines.append(f"{taxid}\tSpecies x\t{reason}")
        path.write_text("\n".join(lines) + "\n")

    def test_reads_faa_filenames(self, tmp_path: Path) -> None:
        refs = self._make_refs(
            tmp_path,
            {"gastropoda": ["100_Foo_bar.faa", "200_Baz_qux.faa"]},
        )
        got = ext.load_existing_taxids(refs)
        assert got == {100, 200}

    def test_reads_manifest_tsv(self, tmp_path: Path) -> None:
        refs = self._make_refs(tmp_path, {"gastropoda": ["100_Foo_bar.faa"]})
        m = tmp_path / "m.tsv"
        self._write_manifest(m, [(300, ""), (400, "no_proteome_in_ncbi")])
        got = ext.load_existing_taxids(refs, m)
        assert got == {100, 300, 400}  # taxids merged regardless of drop_reason

    def test_reads_multiple_manifests(self, tmp_path: Path) -> None:
        refs = self._make_refs(tmp_path, {"gastropoda": ["100_Foo_bar.faa"]})
        m1 = tmp_path / "m1.tsv"
        m2 = tmp_path / "m2.tsv"
        self._write_manifest(m1, [(300, "")])
        self._write_manifest(m2, [(400, ""), (500, "")])
        got = ext.load_existing_taxids(refs, m1, m2)
        assert got == {100, 300, 400, 500}

    def test_missing_manifest_silently_skipped(self, tmp_path: Path) -> None:
        refs = self._make_refs(tmp_path, {"gastropoda": ["100_Foo_bar.faa"]})
        got = ext.load_existing_taxids(refs, tmp_path / "does_not_exist.tsv")
        assert got == {100}

    def test_missing_refs_root_returns_manifest_taxids(self, tmp_path: Path) -> None:
        # When refs_root doesn't exist, fall back to manifests only
        # (refs_root may live on a workstation that isn't running this script)
        m = tmp_path / "m.tsv"
        self._write_manifest(m, [(300, "")])
        got = ext.load_existing_taxids(tmp_path / "no_such_dir", m)
        assert got == {300}

    def test_garbled_taxid_in_filename_skipped(self, tmp_path: Path) -> None:
        refs = self._make_refs(
            tmp_path,
            {"gastropoda": ["100_Foo_bar.faa", "garbage_name.faa"]},
        )
        got = ext.load_existing_taxids(refs)
        assert got == {100}


# ----------------------------------------------------------------------
# group_records_by_taxid + extract_binomial
# ----------------------------------------------------------------------

class TestGrouping:

    def test_groups_by_taxid(self) -> None:
        records = [
            _record("GCA_1", taxid=10, organism_name="Aaa bbb"),
            _record("GCA_2", taxid=10, organism_name="Aaa bbb"),
            _record("GCA_3", taxid=20, organism_name="Ccc ddd"),
        ]
        groups = ext.group_records_by_taxid(records)
        assert set(groups.keys()) == {10, 20}
        assert len(groups[10]) == 2
        assert len(groups[20]) == 1

    def test_drops_records_missing_taxid(self) -> None:
        records = [
            _record("GCA_1", taxid=10, organism_name="Aaa bbb"),
            {"accession": "GCA_X"},  # no organism / tax_id
        ]
        groups = ext.group_records_by_taxid(records)
        assert set(groups.keys()) == {10}

    def test_extract_binomial_returns_first_nonempty(self) -> None:
        records = [
            {"organism": {"organism_name": ""}},
            _record("GCA_1", taxid=10, organism_name="Genus species"),
        ]
        assert ext.extract_binomial(records) == "Genus species"

    def test_extract_binomial_empty_when_all_missing(self) -> None:
        assert ext.extract_binomial([{"organism": {}}]) == ""


# ----------------------------------------------------------------------
# apply_assembly_level_filter
# ----------------------------------------------------------------------

class TestAssemblyLevelFilter:

    def _choice(self, level: str) -> AssemblyChoice:
        return AssemblyChoice(
            accession="GCA_1", source="GenBank",
            assembly_level=level, annotation_status="",
            est_protein_count=0, submission_date="2020-01-01",
        )

    def test_pass_through_when_min_empty(self) -> None:
        c = self._choice("Contig")
        assert ext.apply_assembly_level_filter(c, "") is c

    def test_chromosome_threshold_keeps_chromosome(self) -> None:
        c = self._choice("Chromosome")
        assert ext.apply_assembly_level_filter(c, "Chromosome") is c

    def test_chromosome_threshold_keeps_complete_genome(self) -> None:
        c = self._choice("Complete Genome")
        assert ext.apply_assembly_level_filter(c, "Chromosome") is c

    def test_chromosome_threshold_drops_scaffold(self) -> None:
        c = self._choice("Scaffold")
        assert ext.apply_assembly_level_filter(c, "Chromosome") is None

    def test_chromosome_threshold_drops_contig(self) -> None:
        c = self._choice("Contig")
        assert ext.apply_assembly_level_filter(c, "Chromosome") is None

    def test_none_passes_through(self) -> None:
        assert ext.apply_assembly_level_filter(None, "Chromosome") is None


# ----------------------------------------------------------------------
# select_for_clade  (per-species best-pick + filter + dedup + cap)
# ----------------------------------------------------------------------

class TestSelectForClade:

    def _policy(self, **over) -> ext.ClaadePolicy:
        defaults = dict(
            clade_name="TestClade",
            policy_class="test",
            min_assembly_level="",
            require_annotation=False,
            max_count=None,
        )
        defaults.update(over)
        return ext.ClaadePolicy(**defaults)

    def test_returns_one_entry_per_species(self) -> None:
        records = [
            _record("GCA_1", taxid=10, organism_name="A a"),
            _record("GCA_2", taxid=10, organism_name="A a"),  # 2nd assembly same sp
            _record("GCA_3", taxid=20, organism_name="B b"),
        ]
        out = ext.select_for_clade(records, self._policy(), exclude_taxids=set())
        assert {e.taxid for e in out} == {10, 20}

    def test_excludes_already_seen_taxids(self) -> None:
        records = [
            _record("GCA_1", taxid=10, organism_name="A a"),
            _record("GCA_2", taxid=20, organism_name="B b"),
        ]
        out = ext.select_for_clade(records, self._policy(), exclude_taxids={20})
        assert {e.taxid for e in out} == {10}

    def test_drops_below_min_assembly_level(self) -> None:
        records = [
            _record("GCA_1", taxid=10, organism_name="A a",
                    assembly_level="Scaffold"),
            _record("GCA_2", taxid=20, organism_name="B b",
                    assembly_level="Chromosome"),
        ]
        out = ext.select_for_clade(
            records,
            self._policy(min_assembly_level="Chromosome"),
            exclude_taxids=set(),
        )
        assert {e.taxid for e in out} == {20}

    def test_max_count_caps_with_quality_ranking(self) -> None:
        records = [
            _record("GCA_1", taxid=10, organism_name="A a",
                    assembly_level="Scaffold"),
            _record("GCA_2", taxid=20, organism_name="B b",
                    assembly_level="Chromosome"),
            _record("GCA_3", taxid=30, organism_name="C c",
                    assembly_level="Complete Genome"),
            _record("GCA_4", taxid=40, organism_name="D d",
                    assembly_level="Chromosome"),
        ]
        out = ext.select_for_clade(
            records,
            self._policy(max_count=2),
            exclude_taxids=set(),
        )
        # Top 2 by quality should be Complete Genome (30) + Chromosome (20 or 40)
        assert len(out) == 2
        top = {e.taxid for e in out}
        assert 30 in top  # best quality always kept
        # Among the two Chromosomes, lower taxid wins by stable secondary sort
        assert 20 in top

    def test_policy_class_propagated(self) -> None:
        records = [_record("GCA_1", taxid=10, organism_name="A a")]
        out = ext.select_for_clade(
            records,
            self._policy(policy_class="heterobranchia",
                         clade_name="Heterobranchia"),
            exclude_taxids=set(),
        )
        assert len(out) == 1
        assert out[0].policy_class == "heterobranchia"
        assert out[0].clade_name == "Heterobranchia"

    def test_binomial_populated_from_record(self) -> None:
        records = [_record("GCA_1", taxid=10,
                           organism_name="Genus species")]
        out = ext.select_for_clade(records, self._policy(), set())
        assert out[0].binomial == "Genus species"


# ----------------------------------------------------------------------
# build_extension_inventory (orchestrator)
# ----------------------------------------------------------------------

class TestBuildExtensionInventory:

    def test_dedups_across_policies(self) -> None:
        # Same taxid shows up under two clade queries (e.g. Heterobranchia
        # and Mollusca both return Aplysia). First policy that finds it wins.
        def fake_query(clade: str) -> list[dict]:
            if clade == "Heterobranchia":
                return [_record("GCA_1", taxid=10, organism_name="Aplysia x",
                                assembly_level="Chromosome")]
            if clade == "Bivalvia":
                return [
                    _record("GCA_2", taxid=10, organism_name="Aplysia x",
                            assembly_level="Chromosome"),  # dup of 10
                    _record("GCA_3", taxid=20, organism_name="Crassostrea y",
                            assembly_level="Chromosome"),
                ]
            return []

        policies = [
            ext.ClaadePolicy("Heterobranchia", "heterobranchia",
                             "", False, None),
            ext.ClaadePolicy("Bivalvia", "other_mollusca",
                             "Chromosome", False, None),
        ]
        entries = ext.build_extension_inventory(
            policies, exclude_taxids=set(), query_fn=fake_query,
        )
        # taxid 10 keeps heterobranchia label (first wins), 20 added under bivalvia
        seen = {(e.taxid, e.policy_class) for e in entries}
        assert seen == {(10, "heterobranchia"), (20, "other_mollusca")}

    def test_respects_exclude_taxids(self) -> None:
        def fake_query(clade: str) -> list[dict]:
            return [_record("GCA_1", taxid=999, organism_name="X y")]

        policies = [ext.ClaadePolicy("Whatever", "test", "", False, None)]
        entries = ext.build_extension_inventory(
            policies, exclude_taxids={999}, query_fn=fake_query,
        )
        assert entries == []


# ----------------------------------------------------------------------
# write_extension_manifest
# ----------------------------------------------------------------------

class TestWriteManifest:

    def _entry(self, *, taxid: int, policy_class: str, level: str = "Chromosome",
               binomial: str = "Genus species") -> ext.ExtensionEntry:
        return ext.ExtensionEntry(
            taxid=taxid, binomial=binomial,
            policy_class=policy_class, clade_name="X",
            choice=AssemblyChoice(
                accession=f"GCA_{taxid}", source="GenBank",
                assembly_level=level, annotation_status="",
                est_protein_count=0, submission_date="2020-01-01",
                contig_n50=50000, total_length_bp=1_000_000_000,
            ),
        )

    def test_writes_header_and_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "extension.tsv"
        ext.write_extension_manifest(out, [
            self._entry(taxid=10, policy_class="heterobranchia"),
            self._entry(taxid=20, policy_class="other_mollusca"),
        ])
        lines = out.read_text().splitlines()
        assert lines[0].split("\t")[:4] == [
            "taxid", "binomial", "policy_class", "clade_name",
        ]
        assert len(lines) == 3  # header + 2

    def test_sorted_by_policy_class_then_taxid(self, tmp_path: Path) -> None:
        out = tmp_path / "extension.tsv"
        ext.write_extension_manifest(out, [
            self._entry(taxid=20, policy_class="other_mollusca"),
            self._entry(taxid=10, policy_class="other_mollusca"),
            self._entry(taxid=5, policy_class="heterobranchia"),
        ])
        rows = out.read_text().splitlines()[1:]
        taxids = [int(r.split("\t")[0]) for r in rows]
        policy_classes = [r.split("\t")[2] for r in rows]
        # heterobranchia comes before other_mollusca alphabetically
        assert policy_classes == [
            "heterobranchia", "other_mollusca", "other_mollusca",
        ]
        assert taxids == [5, 10, 20]

    def test_includes_assembly_stats_columns(self, tmp_path: Path) -> None:
        out = tmp_path / "extension.tsv"
        ext.write_extension_manifest(out, [
            self._entry(taxid=10, policy_class="heterobranchia"),
        ])
        header = out.read_text().splitlines()[0].split("\t")
        assert "contig_n50" in header
        assert "total_length_bp" in header


# ----------------------------------------------------------------------
# query_datasets_for_clade (subprocess wrapper)
# ----------------------------------------------------------------------

class TestQueryDatasetsForClade:
    """Mirrors Phase 1a's query_datasets_for_taxon contract — must
    tolerate the same set of NCBI "no hits" error phrases without
    raising, since some clades may genuinely have no chromosome-level
    assemblies (e.g. Caudofoveata only has 2 species globally).
    """

    def test_parses_jsonl_stdout(self, monkeypatch) -> None:
        import subprocess
        captured: dict = {}

        class FakeResult:
            stdout = (
                '{"accession":"GCA_1","organism":{"tax_id":10,"organism_name":"A a"}}\n'
                '{"accession":"GCA_2","organism":{"tax_id":20,"organism_name":"B b"}}\n'
            )
            stderr = ""
            returncode = 0

        def fake_run(cmd, **kw):
            captured["cmd"] = cmd
            return FakeResult()

        monkeypatch.setattr(subprocess, "run", fake_run)
        out = ext.query_datasets_for_clade("Heterobranchia",
                                            datasets_bin="datasets")
        assert len(out) == 2
        assert captured["cmd"][0] == "datasets"
        assert "Heterobranchia" in captured["cmd"]
        assert "--as-json-lines" in captured["cmd"]

    def test_no_hit_phrase_returns_empty(self, monkeypatch) -> None:
        import subprocess

        class FakeResult:
            stdout = ""
            stderr = "Error: did not match any genomes"
            returncode = 1

        monkeypatch.setattr(subprocess, "run", lambda *a, **k: FakeResult())
        assert ext.query_datasets_for_clade("NoSuchClade") == []

    def test_other_error_raises(self, monkeypatch) -> None:
        import subprocess

        class FakeResult:
            stdout = ""
            stderr = "Error: server unreachable"
            returncode = 1

        monkeypatch.setattr(subprocess, "run", lambda *a, **k: FakeResult())
        with pytest.raises(RuntimeError, match="server unreachable"):
            ext.query_datasets_for_clade("Heterobranchia")
