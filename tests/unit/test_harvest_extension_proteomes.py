import os
import sys
import json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
import harvest_extension_proteomes as hep  # noqa: E402

EXT = (
    "taxid\tbinomial\tpolicy_class\tclade_name\tsource\taccession\tassembly_level\t"
    "annotation_status\test_protein_count\tsubmission_date\tcontig_n50\ttotal_length_bp\n"
    "231223\tElysia crispata\t\tgastropoda\tGenBank\tGCA_033675545.1\tContig\t\t67429\t2023-01-01\t1000\t500000\n"
    "6161\tDugesia japonica\t\tplatyhelminthes\tGenBank\tGCA_001938525.1\tScaffold\t\t0\t2016-01-01\t1166\t854187087\n"
)


def _jsonl(*recs):
    return "\n".join(json.dumps(r) for r in recs) + "\n"


JSONL = _jsonl(
    {"accession": "GCA_033675545.1", "annotation_info": {"name": "x", "stats": {}}},
    {"accession": "GCA_001938525.1"},  # no annotation_info
)


def test_annotated_accessions(tmp_path):
    p = tmp_path / "s.jsonl"
    p.write_text(JSONL)
    assert hep.annotated_accessions(str(p)) == {"GCA_033675545.1"}


def test_select_annotated_extension(tmp_path):
    e = tmp_path / "ext.tsv"
    e.write_text(EXT)
    j = tmp_path / "s.jsonl"
    j.write_text(JSONL)
    got = hep.select_annotated_extension(str(e), str(j))
    assert len(got) == 1
    t = got[0]
    assert t["taxid"] == "231223" and t["binomial"] == "Elysia crispata"
    assert t["clade"] == "gastropoda" and t["accession"] == "GCA_033675545.1"
    assert t["est_protein_count"] == "67429"


def test_write_download_manifest(tmp_path):
    e = tmp_path / "ext.tsv"
    e.write_text(EXT)
    j = tmp_path / "s.jsonl"
    j.write_text(JSONL)
    targets = hep.select_annotated_extension(str(e), str(j))
    out = tmp_path / "dl.tsv"
    hep.write_download_manifest(targets, str(out))
    lines = out.read_text().splitlines()
    assert lines[0] == "taxid\tbinomial\tclade\taccession"
    assert lines[1] == "231223\tElysia crispata\tgastropoda\tGCA_033675545.1"


def test_write_staged_manifest(tmp_path):
    e = tmp_path / "ext.tsv"
    e.write_text(EXT)
    j = tmp_path / "s.jsonl"
    j.write_text(JSONL)
    targets = hep.select_annotated_extension(str(e), str(j))
    out = tmp_path / "staged.tsv"
    hep.write_staged_manifest(targets, str(out))
    lines = out.read_text().splitlines()
    expected_header = ("taxid\tbinomial\tclade\tsource\taccession\tassembly_level\t"
                       "annotation_status\test_protein_count\tsubmission_date\tdrop_reason")
    assert lines[0] == expected_header
    assert len(lines) == 2
    row_fields = lines[1].split("\t")
    col_list = list(hep.PROTEOME_MANIFEST_COLUMNS)
    assert row_fields[col_list.index("taxid")] == "231223"
    assert row_fields[col_list.index("binomial")] == "Elysia crispata"
    assert row_fields[col_list.index("drop_reason")] == ""


def test_select_empty_extension_returns_empty(tmp_path):
    e = tmp_path / "ext.tsv"
    e.write_text("taxid\tbinomial\tclade_name\taccession\n")  # header only
    j = tmp_path / "s.jsonl"
    j.write_text(JSONL)
    assert hep.select_annotated_extension(str(e), str(j)) == []


def test_annotated_accessions_malformed_json_raises(tmp_path):
    import pytest
    p = tmp_path / "bad.jsonl"
    p.write_text('{"accession": "GCA_1"}\nNOT JSON\n')
    with pytest.raises(json.JSONDecodeError):
        hep.annotated_accessions(str(p))


PM = ("taxid\tbinomial\tclade\tsource\taccession\tassembly_level\t"
      "annotation_status\test_protein_count\tsubmission_date\tdrop_reason\n"
      "6161\tDugesia japonica\tplatyhelminthes\t\t\t\t\t\t\tno_proteome_in_ncbi\n")
STAGED = ("taxid\tbinomial\tclade\tsource\taccession\tassembly_level\t"
          "annotation_status\test_protein_count\tsubmission_date\tdrop_reason\n"
          "231223\tElysia crispata\tgastropoda\tGenBank\tGCA_033675545.1\tContig\t\t67429\t2023-01-01\t\n"
          "999\tFoo bar\tx\tGenBank\tGCA_000.1\tContig\t\t10\t2020-01-01\t\n")
RESULTS = "taxid\tstatus\n231223\tok\n999\tdownload_failed\n"


def test_append_only_ok_idempotent(tmp_path):
    import csv as _csv
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM)
    st = tmp_path / "staged.tsv"
    st.write_text(STAGED)
    rs = tmp_path / "res.tsv"
    rs.write_text(RESULTS)
    n = hep.append_to_proteome_manifest(str(st), str(rs), str(pm))
    assert n == 1                                   # only 231223 (ok); 999 failed
    with open(pm) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    taxids = {r["taxid"] for r in rows}
    assert taxids == {"6161", "231223"}             # existing preserved + Elysia added
    elysia = [r for r in rows if r["taxid"] == "231223"][0]
    assert elysia["est_protein_count"] == "67429" and elysia["drop_reason"] == ""
    # idempotent re-run appends nothing
    assert hep.append_to_proteome_manifest(str(st), str(rs), str(pm)) == 0


def test_append_uses_lf_line_endings(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM)
    st = tmp_path / "staged.tsv"
    st.write_text(STAGED)
    rs = tmp_path / "res.tsv"
    rs.write_text(RESULTS)
    hep.append_to_proteome_manifest(str(st), str(rs), str(pm))
    assert "\r" not in pm.read_text()


def test_append_includes_ok_no_cds(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM)
    st = tmp_path / "staged.tsv"
    st.write_text("taxid\tbinomial\tclade\tsource\taccession\tassembly_level\t"
                  "annotation_status\test_protein_count\tsubmission_date\tdrop_reason\n"
                  "555\tBar baz\ty\tGenBank\tGCA_5.1\tContig\t\t5\t2021-01-01\t\n")
    rs = tmp_path / "res.tsv"
    rs.write_text("taxid\tstatus\n555\tok_no_cds\n")
    assert hep.append_to_proteome_manifest(str(st), str(rs), str(pm)) == 1


# ---------------------------------------------------------------------------
# CLI tests (Task 2)
# ---------------------------------------------------------------------------
import pytest


def test_main_help_exits_zero():
    with pytest.raises(SystemExit) as e:
        hep.main(["--help"])
    assert e.value.code == 0


def test_main_select_help_exits_zero():
    with pytest.raises(SystemExit) as e:
        hep.main(["select", "--help"])
    assert e.value.code == 0


def test_main_select_missing_args_exits_2():
    with pytest.raises(SystemExit) as e:
        hep.main(["select"])
    assert e.value.code == 2


def test_main_append_missing_args_exits_2():
    with pytest.raises(SystemExit) as e:
        hep.main(["append"])
    assert e.value.code == 2


def test_main_select_happy_path(tmp_path):
    e = tmp_path / "ext.tsv"
    e.write_text(EXT)
    j = tmp_path / "s.jsonl"
    j.write_text(JSONL)
    dl = tmp_path / "dl.tsv"
    staged = tmp_path / "staged.tsv"
    hep.main(["select", "--extension-tsv", str(e), "--datasets-jsonl", str(j),
              "--download-manifest-out", str(dl), "--staged-out", str(staged)])
    assert dl.exists() and staged.exists()
    assert "GCA_033675545.1" in dl.read_text()


def test_main_select_missing_extension_tsv_exits_2(tmp_path):
    j = tmp_path / "s.jsonl"
    j.write_text(JSONL)
    dl = tmp_path / "dl.tsv"
    staged = tmp_path / "staged.tsv"
    with pytest.raises(SystemExit) as e:
        hep.main(["select",
                  "--extension-tsv", str(tmp_path / "no_such.tsv"),
                  "--datasets-jsonl", str(j),
                  "--download-manifest-out", str(dl),
                  "--staged-out", str(staged)])
    assert e.value.code == 2


def test_main_select_missing_jsonl_exits_2(tmp_path):
    ext = tmp_path / "ext.tsv"
    ext.write_text(EXT)
    dl = tmp_path / "dl.tsv"
    staged = tmp_path / "staged.tsv"
    with pytest.raises(SystemExit) as e:
        hep.main(["select",
                  "--extension-tsv", str(ext),
                  "--datasets-jsonl", str(tmp_path / "no_such.jsonl"),
                  "--download-manifest-out", str(dl),
                  "--staged-out", str(staged)])
    assert e.value.code == 2


def test_main_append_happy_path(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM)
    st = tmp_path / "staged.tsv"
    st.write_text(STAGED)
    rs = tmp_path / "res.tsv"
    rs.write_text(RESULTS)
    hep.main(["append",
              "--staged", str(st),
              "--download-results", str(rs),
              "--proteome-manifest", str(pm)])
    import csv as _csv
    with open(pm) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    assert any(r["taxid"] == "231223" for r in rows)


def test_main_append_missing_staged_exits_2(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM)
    rs = tmp_path / "res.tsv"
    rs.write_text(RESULTS)
    with pytest.raises(SystemExit) as e:
        hep.main(["append",
                  "--staged", str(tmp_path / "no_staged.tsv"),
                  "--download-results", str(rs),
                  "--proteome-manifest", str(pm)])
    assert e.value.code == 2


# ---------------------------------------------------------------------------
# mark subcommand tests
# ---------------------------------------------------------------------------

# 14-column genome_inventory header
GI_HEADER = (
    "taxid\tbinomial\tclade\tpolicy_class\tsource\taccession\t"
    "assembly_level\tannotation_status\test_protein_count\tsubmission_date\t"
    "contig_n50\ttotal_length_bp\tdrop_reason\tsource_batch"
)

# proteome_manifest rows for taxids 231223 and 999
PM_MARK = (
    "taxid\tbinomial\tclade\tsource\taccession\tassembly_level\t"
    "annotation_status\test_protein_count\tsubmission_date\tdrop_reason\n"
    "231223\tElysia crispata\tgastropoda\tGenBank\tGCA_033675545.1\tContig\t\t67429\t2023-01-01\t\n"
    "999\tFoo bar\tx\tGenBank\tGCA_000.1\tContig\t\t10\t2020-01-01\t\n"
)


def _make_gi(tmp_path, rows, linesep="\n"):
    r"""Write a minimal genome_inventory.tsv with the standard 14-col header.

    Bytes are written verbatim so `linesep` controls the on-disk terminator
    (default LF; pass "\r\n" for the canonical CRLF genome_inventory.tsv).
    """
    p = tmp_path / "genome_inventory.tsv"
    lines = [GI_HEADER]
    for r in rows:
        lines.append("\t".join(r.get(c, "") for c in GI_HEADER.split("\t")))
    p.write_bytes((linesep.join(lines) + linesep).encode("utf-8"))
    return p


def test_mark_sets_drop_reason_for_proteome_taxids(tmp_path):
    """Rows whose taxid is in the proteome_manifest get drop_reason='harvested_annotated'."""
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "binomial": "Elysia crispata", "accession": "GCA_033675545.1"},
        {"taxid": "6161", "binomial": "Dugesia japonica", "accession": "GCA_001938525.1"},
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    n = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    assert n == 1  # only 231223; 999 is in pm but not in gi
    import csv as _csv
    with open(gi) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    by_taxid = {r["taxid"]: r for r in rows}
    assert by_taxid["231223"]["drop_reason"] == "harvested_annotated"
    assert by_taxid["6161"]["drop_reason"] == ""


def test_mark_skips_no_proteome_manifest_entries(tmp_path):
    """A taxid in proteome_manifest flagged drop_reason='no_proteome_in_ncbi' must
    NOT be marked: those species have no usable proteome and still need BRAKER
    de-novo annotation. Regression for the real-data over-exclusion — the manifest
    tracks BOTH species-with-proteomes and confirmed-no-proteome species, so the
    mark must key only off usable (empty-drop_reason) entries."""
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "binomial": "Elysia crispata",
         "accession": "GCA_033675545.1"},                                  # usable in PM
        {"taxid": "6359", "binomial": "Platynereis dumerilii",
         "accession": "GCA_043381215.1"},                                  # no_proteome in PM
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(
        "taxid\tbinomial\tclade\tsource\taccession\tassembly_level\t"
        "annotation_status\test_protein_count\tsubmission_date\tdrop_reason\n"
        "231223\tElysia crispata\tgastropoda\tGenBank\tGCA_033675545.1\tContig\t\t67429\t2023-01-01\t\n"
        "6359\tPlatynereis dumerilii\tannelida\tGenBank\tGCA_043381215.1\tContig\t\t0\t2024-01-01\tno_proteome_in_ncbi\n"
    )
    n = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    assert n == 1  # only 231223 (usable); 6359 (no_proteome_in_ncbi) must NOT be marked
    import csv as _csv
    with open(gi) as fh:
        by_taxid = {r["taxid"]: r for r in _csv.DictReader(fh, delimiter="\t")}
    assert by_taxid["231223"]["drop_reason"] == "harvested_annotated"
    assert by_taxid["6359"]["drop_reason"] == ""   # no-proteome species stays a BRAKER target


def test_mark_preserves_accession_and_all_columns(tmp_path):
    """Every non-drop_reason column must be byte-identical after marking."""
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "binomial": "Elysia crispata",
         "accession": "GCA_033675545.1", "clade": "gastropoda",
         "source": "GenBank", "assembly_level": "Contig",
         "annotation_status": "Full", "est_protein_count": "67429",
         "submission_date": "2023-01-01", "contig_n50": "1000",
         "total_length_bp": "500000", "source_batch": "phase1d"},
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)

    import csv as _csv

    def _read(p):
        with open(p) as fh:
            return list(_csv.DictReader(fh, delimiter="\t"))

    before = _read(gi)[0]
    hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    after = _read(gi)[0]

    for col in ("taxid", "binomial", "clade", "policy_class", "source",
                "accession", "assembly_level", "annotation_status",
                "est_protein_count", "submission_date", "contig_n50",
                "total_length_bp", "source_batch"):
        assert after[col] == before[col], f"column {col!r} changed"
    assert after["drop_reason"] == "harvested_annotated"


def test_mark_preserves_header_and_column_order(tmp_path):
    """The rewritten file must have the same header in the same order."""
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "accession": "GCA_033675545.1"},
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    first_line = open(gi).readline().rstrip("\n")
    assert first_line == GI_HEADER


def test_mark_idempotent(tmp_path):
    """A second call returns 0 and produces bit-identical file bytes."""
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "accession": "GCA_033675545.1"},
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    n1 = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    bytes_after_first = open(gi, "rb").read()
    n2 = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    assert n1 == 1
    assert n2 == 0
    assert open(gi, "rb").read() == bytes_after_first


def test_mark_preserves_crlf_line_endings(tmp_path):
    """A CRLF-input file stays CRLF after a mark that changes a row, and the
    marked row's field values carry no stray carriage returns.
    """
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "binomial": "Elysia crispata", "accession": "GCA_033675545.1"},
    ], linesep="\r\n")
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    n = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    assert n == 1
    raw = open(gi, "rb").read()
    assert b"\r\n" in raw                            # still CRLF
    assert b"\n" not in raw.replace(b"\r\n", b"")   # ... no bare LF leaked
    # The csv reader (newline="") strips \r, so field values must be clean.
    import csv as _csv
    with open(gi, newline="") as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    row = rows[0]
    assert row["drop_reason"] == "harvested_annotated"
    for v in row.values():
        assert "\r" not in (v or "")


def test_mark_does_not_touch_rows_not_in_manifest(tmp_path):
    """A no-match run leaves the file byte-identical (CRLF endings preserved)."""
    gi = _make_gi(tmp_path, [
        {"taxid": "6161", "binomial": "Dugesia japonica", "accession": "GCA_001938525.1"},
        {"taxid": "6162", "binomial": "Girardia tigrina", "accession": "GCA_001938485.1"},
    ], linesep="\r\n")
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)  # contains 231223 and 999, not 6161/6162
    before = open(gi, "rb").read()
    n = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    assert n == 0
    assert open(gi, "rb").read() == before  # bytes unchanged, CRLF intact


def test_mark_does_not_overwrite_existing_drop_reason(tmp_path):
    """A row that already carries a non-empty drop_reason is left unchanged."""
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "accession": "GCA_033675545.1",
         "drop_reason": "no_proteome_in_ncbi"},
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    n = hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    assert n == 0
    import csv as _csv
    with open(gi) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["drop_reason"] == "no_proteome_in_ncbi"


def test_mark_uses_lf_line_endings(tmp_path):
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "accession": "GCA_033675545.1"},
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    hep.mark_harvested_in_genome_inventory(str(gi), str(pm))
    assert b"\r" not in open(gi, "rb").read()


def test_mark_cli_help_exits_zero():
    with pytest.raises(SystemExit) as e:
        hep.main(["mark", "--help"])
    assert e.value.code == 0


def test_mark_cli_missing_genome_inventory_exits_2(tmp_path):
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    with pytest.raises(SystemExit) as e:
        hep.main(["mark",
                  "--genome-inventory", str(tmp_path / "no_such.tsv"),
                  "--proteome-manifest", str(pm)])
    assert e.value.code == 2


def test_mark_cli_missing_proteome_manifest_exits_2(tmp_path):
    gi = _make_gi(tmp_path, [])
    with pytest.raises(SystemExit) as e:
        hep.main(["mark",
                  "--genome-inventory", str(gi),
                  "--proteome-manifest", str(tmp_path / "no_such.tsv")])
    assert e.value.code == 2


def test_mark_cli_happy_path(tmp_path):
    gi = _make_gi(tmp_path, [
        {"taxid": "231223", "accession": "GCA_033675545.1"},
    ])
    pm = tmp_path / "pm.tsv"
    pm.write_text(PM_MARK)
    hep.main(["mark", "--genome-inventory", str(gi), "--proteome-manifest", str(pm)])
    import csv as _csv
    with open(gi) as fh:
        rows = list(_csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["drop_reason"] == "harvested_annotated"
