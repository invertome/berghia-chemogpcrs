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
