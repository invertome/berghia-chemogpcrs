"""The assembly date tie-break must read the field NCBI actually returns.

2026-07 audit finding: ``pick_best_assembly`` read
``assembly_info["submission_date"]``, which does not exist in the NCBI
datasets v2alpha ``dataset_report`` payload. Verified 2026-07-19 against
GCA_015776775.1 -- ``assembly_info`` keys are::

    assembly_level, assembly_method, assembly_name, assembly_status,
    assembly_type, bioproject_accession, bioproject_lineage, biosample,
    blast_url, comments, refseq_category, release_date, sequencing_tech,
    submitter

so the real field is ``release_date`` ('2020-12-08'). The date tie-break was
therefore a permanent no-op and WHICH genome got selected could flip between
re-queries whenever two candidates tied on source/annotation/assembly level.
466 of 467 manifest rows carry an empty date, consistent with this.
"""
from __future__ import annotations

from build_species_tree_phase1a_inventory import pick_best_assembly


def _rec(accession, level="Chromosome", release_date=None, annotated=True):
    info = {"assembly_level": level}
    if release_date is not None:
        info["release_date"] = release_date
    rec = {"accession": accession, "assembly_info": info}
    if annotated:
        rec["annotation_info"] = {"status": "Full annotation",
                                  "stats": {"gene_counts": {"protein_coding": 20000}}}
    return rec


def test_release_date_breaks_a_tie_toward_the_newer_assembly():
    older = _rec("GCA_000000001.1", release_date="2018-01-01")
    newer = _rec("GCA_000000002.1", release_date="2023-06-15")
    assert pick_best_assembly([older, newer]).accession == "GCA_000000002.1"
    # order-independent
    assert pick_best_assembly([newer, older]).accession == "GCA_000000002.1"


def test_release_date_is_captured_on_the_choice():
    choice = pick_best_assembly([_rec("GCA_000000001.1", release_date="2021-09-30")])
    assert choice.submission_date == "2021-09-30"


def test_legacy_submission_date_key_is_not_used():
    """A payload carrying only the non-existent key yields no date, and must
    not be mistaken for a real one."""
    rec = {"accession": "GCA_000000001.1",
           "assembly_info": {"assembly_level": "Chromosome",
                             "submission_date": "2099-01-01"},
           "annotation_info": {"status": "Full annotation",
                               "stats": {"gene_counts": {"protein_coding": 1}}}}
    assert pick_best_assembly([rec]).submission_date == ""


def test_missing_release_date_still_selects_without_error():
    a = _rec("GCA_000000001.1", release_date=None)
    b = _rec("GCA_000000002.1", release_date="2020-01-01")
    # The dated candidate wins the tie-break; the undated one must not crash.
    assert pick_best_assembly([a, b]).accession == "GCA_000000002.1"


def test_higher_priority_class_still_outranks_a_newer_date():
    """Date is the LAST tie-break -- it must not override source/annotation."""
    refseq_old = _rec("GCF_000000001.1", release_date="2010-01-01")
    refseq_old["source_database"] = "SOURCE_DATABASE_REFSEQ"
    genbank_new = _rec("GCA_000000002.1", release_date="2025-01-01")
    assert pick_best_assembly([genbank_new, refseq_old]).accession == "GCF_000000001.1"


def test_assembly_level_still_outranks_a_newer_date():
    contig_new = _rec("GCA_000000002.1", level="Contig", release_date="2025-01-01")
    chrom_old = _rec("GCA_000000001.1", level="Chromosome", release_date="2010-01-01")
    assert pick_best_assembly([contig_new, chrom_old]).accession == "GCA_000000001.1"
