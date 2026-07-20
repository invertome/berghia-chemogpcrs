"""Tests for draw_outgroup_fasta in functions.sh (bead 444).

Stage 04 roots each per-class tree on an outgroup drawn from the sister class's
reference pool (swap-map A<-C, B<-A, C<-A, F<-A; locked decision 2026-05-28).
The draw was::

    awk -v n="${OUTGROUP_BUDGET_PER_CLASS:-10}" \\
        '/^>/{c++} c>n{exit} {print}' "${OUTGROUP_SOURCE_FA}" > "${OUTGROUP_FA}"

i.e. the first N records of the pool file. But
``build_per_class_reference_pools.py`` writes::

    selected = berghia_pairs + anchor_pairs + selected_refs

so *Berghia records come FIRST in every pool file*. The head-N draw therefore
took exactly Berghia's own class-B/C/F candidates, and every per-class tree was
rooted on unannotated Berghia paralogs whose class assignment came from this
pipeline's own classifier. That is circular (the classifier's output defines the
root that tests the classifier's output), single-species (no outgroup diversity),
and not the characterized reference GPCRs the design intended.

``draw_outgroup_fasta`` draws instead from records whose provenance is not
Berghia, using the ``pool_members_class_<C>.tsv`` manifest the pool builder
already writes (``seq_id, taxid, source``; source in {anchor, berghia, ref}).
Anchors and refs are both valid outgroup material; only ``berghia`` is excluded.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FUNCTIONS_SH = PROJECT_ROOT / "functions.sh"
STAGE04 = PROJECT_ROOT / "04_phylogenetic_analysis.sh"
POOL_BUILDER = PROJECT_ROOT / "scripts" / "build_per_class_reference_pools.py"


def _draw(tmp_path: Path, source_fa, members_tsv, n, out_fa) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["LOGS_DIR"] = str(tmp_path)
    script = (
        f'source "{FUNCTIONS_SH}"; trap - EXIT; '
        f'draw_outgroup_fasta "{source_fa}" "{members_tsv}" "{n}" "{out_fa}"'
    )
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env)


def _write_pool(tmp_path: Path, records: list[tuple[str, str, str]], name: str = "refs_class_C") -> tuple[Path, Path]:
    """Write a pool FASTA + its members manifest in pool-builder order.

    *records* is a list of (seq_id, source, sequence). The FASTA and the TSV are
    written in the given order, exactly as the pool builder co-writes them.
    """
    fa = tmp_path / f"{name}.fa"
    fa.write_text("".join(f">{sid}\n{seq}\n" for sid, _src, seq in records))
    tsv = tmp_path / f"pool_members_{name.replace('refs_', '')}.tsv"
    lines = ["seq_id\ttaxid\tsource"]
    lines += [f"{sid}\t9999\t{src}" for sid, src, _seq in records]
    tsv.write_text("\n".join(lines) + "\n")
    return fa, tsv


def _ids(fa: Path) -> list[str]:
    if not fa.exists():
        return []
    return [ln[1:].split()[0] for ln in fa.read_text().splitlines() if ln.startswith(">")]


def _seq_of(fa: Path, seq_id: str) -> str:
    out, capture = [], False
    for ln in fa.read_text().splitlines():
        if ln.startswith(">"):
            capture = ln[1:].split()[0] == seq_id
            continue
        if capture:
            out.append(ln)
    return "\n".join(out)


# --- provenance filtering -------------------------------------------------

def test_berghia_records_are_excluded(tmp_path):
    """The defect: Berghia is first in the file, so head-N took only Berghia."""
    fa, tsv = _write_pool(tmp_path, [
        ("Berste_0001", "berghia", "MKB1"),
        ("Berste_0002", "berghia", "MKB2"),
        ("ANCHOR_P08908", "anchor", "MKA1"),
        ("ref_6500_7", "ref", "MKR1"),
    ])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 2, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ANCHOR_P08908", "ref_6500_7"]


def test_anchors_and_refs_are_both_eligible(tmp_path):
    fa, tsv = _write_pool(tmp_path, [
        ("Berste_0001", "berghia", "MKB1"),
        ("ANCHOR_P08908", "anchor", "MKA1"),
        ("ref_6500_7", "ref", "MKR1"),
        ("ANCHOR_P41143", "anchor", "MKA2"),
    ])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 10, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ANCHOR_P08908", "ref_6500_7", "ANCHOR_P41143"]


def test_berghia_interleaved_anywhere_is_excluded(tmp_path):
    """Order-independence: don't rely on Berghia being a leading prefix block."""
    fa, tsv = _write_pool(tmp_path, [
        ("ref_a", "ref", "MK1"),
        ("Berste_x", "berghia", "MKB"),
        ("ref_b", "ref", "MK2"),
        ("Berste_y", "berghia", "MKB"),
        ("ref_c", "ref", "MK3"),
    ])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 3, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ref_a", "ref_b", "ref_c"]


# --- record integrity -----------------------------------------------------

def test_multiline_sequence_blocks_are_copied_intact(tmp_path):
    """Real pool FASTAs are wrapped; a header-only copy would corrupt the input."""
    fa = tmp_path / "refs_class_C.fa"
    fa.write_text(
        ">Berste_0001\nMKBBBB\nMKBBBB\n"
        ">ref_6500_7\nMKAAAA\nMKCCCC\nMKDDDD\n"
        ">ANCHOR_P08908\nMKEEEE\nMKFFFF\n"
    )
    tsv = tmp_path / "pool_members_class_C.tsv"
    tsv.write_text(
        "seq_id\ttaxid\tsource\n"
        "Berste_0001\t1287507\tberghia\n"
        "ref_6500_7\t6500\tref\n"
        "ANCHOR_P08908\t9606\tanchor\n"
    )
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 2, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ref_6500_7", "ANCHOR_P08908"]
    assert _seq_of(out, "ref_6500_7") == "MKAAAA\nMKCCCC\nMKDDDD"
    assert _seq_of(out, "ANCHOR_P08908") == "MKEEEE\nMKFFFF"
    assert "MKBBBB" not in out.read_text(), "no Berghia residues may leak into the outgroup"


def test_header_description_after_whitespace_still_matches_the_manifest(tmp_path):
    """Manifest ids are the bare rec.id; FASTA headers may carry a description."""
    fa = tmp_path / "refs_class_C.fa"
    fa.write_text(
        ">Berste_0001 hypothetical protein\nMKB\n"
        ">ref_6500_7 Aplysia opsin-like receptor\nMKR\n"
    )
    tsv = tmp_path / "pool_members_class_C.tsv"
    tsv.write_text(
        "seq_id\ttaxid\tsource\n"
        "Berste_0001\t1287507\tberghia\n"
        "ref_6500_7\t6500\tref\n"
    )
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 5, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ref_6500_7"]
    assert "Berste_0001" not in out.read_text()


# --- N handling -----------------------------------------------------------

def test_n_is_respected_when_more_non_berghia_exist(tmp_path):
    fa, tsv = _write_pool(tmp_path, [
        ("Berste_0001", "berghia", "MKB"),
    ] + [(f"ref_{i}", "ref", f"MK{i}") for i in range(10)])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 3, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ref_0", "ref_1", "ref_2"]


def test_fewer_non_berghia_than_n_returns_all_available(tmp_path):
    fa, tsv = _write_pool(tmp_path, [
        ("Berste_0001", "berghia", "MKB"),
        ("Berste_0002", "berghia", "MKB"),
        ("ref_a", "ref", "MK1"),
        ("ref_b", "ref", "MK2"),
    ])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 10, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ref_a", "ref_b"]


def test_exactly_n_non_berghia_records(tmp_path):
    fa, tsv = _write_pool(tmp_path, [
        ("Berste_0001", "berghia", "MKB"),
        ("ref_a", "ref", "MK1"),
        ("ref_b", "ref", "MK2"),
    ])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 2, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ref_a", "ref_b"]
    assert _seq_of(out, "ref_b") == "MK2", "the last drawn record must keep its sequence"


# --- degenerate inputs ----------------------------------------------------

def test_all_berghia_source_yields_empty_output_and_warns(tmp_path):
    """Not an error (the tree is simply unrooted), but it must not be silent."""
    fa, tsv = _write_pool(tmp_path, [
        ("Berste_0001", "berghia", "MKB1"),
        ("Berste_0002", "berghia", "MKB2"),
    ])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 5, out)
    assert result.returncode == 0, result.stderr
    assert out.exists(), "an empty outgroup file must still be created for the downstream cat"
    assert _ids(out) == []
    combined = result.stdout + result.stderr
    assert "0 non-berghia" in combined.lower(), \
        f"an empty draw must warn; got:\n{combined}"


def test_missing_members_tsv_falls_back_unfiltered_with_a_loud_warning(tmp_path):
    """Provenance is undeterminable without the manifest: warn, do not pretend."""
    fa, _tsv = _write_pool(tmp_path, [
        ("Berste_0001", "berghia", "MKB1"),
        ("ref_a", "ref", "MK1"),
    ])
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tmp_path / "nope.tsv", 1, out)
    assert result.returncode == 0, result.stderr
    # Documented fallback = the previous unfiltered head-N behaviour.
    assert _ids(out) == ["Berste_0001"]
    combined = result.stdout + result.stderr
    assert "berghia" in combined.lower() and "manifest" in combined.lower(), \
        f"the fallback must name the risk (manifest missing, Berghia may be present); got:\n{combined}"
    assert "[WARN]" in combined, "the fallback must log at WARN level"


def test_missing_source_fasta_is_an_error(tmp_path):
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, tmp_path / "nope.fa", tmp_path / "nope.tsv", 5, out)
    assert result.returncode == 1
    assert "not found" in result.stderr.lower()


def test_manifest_without_a_header_row_still_works(tmp_path):
    """Defensive: only a literal `seq_id` first field is treated as a header."""
    fa = tmp_path / "refs_class_C.fa"
    fa.write_text(">Berste_0001\nMKB\n>ref_a\nMK1\n")
    tsv = tmp_path / "pool_members_class_C.tsv"
    tsv.write_text("Berste_0001\t1287507\tberghia\nref_a\t6500\tref\n")
    out = tmp_path / "outgroup.fa"
    result = _draw(tmp_path, fa, tsv, 5, out)
    assert result.returncode == 0, result.stderr
    assert _ids(out) == ["ref_a"]


# --- stage 04 wiring ------------------------------------------------------

def test_stage04_uses_the_helper_not_the_head_n_awk(tmp_path):
    text = STAGE04.read_text()
    assert "draw_outgroup_fasta" in text, \
        "stage 04 must draw the outgroup through draw_outgroup_fasta"
    assert "'/^>/{c++} c>n{exit} {print}'" not in text, \
        "the unfiltered head-N awk draw must be gone from stage 04"


def test_stage04_reads_the_sister_class_members_tsv():
    """The outgroup source is the SISTER class, so its manifest must be used."""
    text = STAGE04.read_text()
    assert "pool_members_class_${OUTGROUP_SOURCE_CLASS}.tsv" in text, \
        "stage 04 must read the outgroup SOURCE class's members TSV, not the current class's"


def test_stage04_keeps_the_og_count_log_and_missing_source_warn():
    text = STAGE04.read_text()
    assert "OG_COUNT" in text, "the outgroup-size logging must be preserved"
    assert "tree will be unrooted" in text, \
        "the missing-outgroup-source WARN branch must be preserved"


def test_pool_builder_still_writes_the_provenance_manifest():
    """The helper depends on this contract; pin it so a rename is caught here."""
    text = POOL_BUILDER.read_text()
    assert 'f"pool_members_class_{cls}.tsv"' in text
    assert '["seq_id", "taxid", "source"]' in text
    assert 'source = "berghia"' in text


def test_stage04_still_parses():
    result = subprocess.run(["bash", "-n", str(STAGE04)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_functions_sh_still_parses():
    result = subprocess.run(["bash", "-n", str(FUNCTIONS_SH)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# --- Fail-loud on FASTA/manifest disagreement (user decision 2026-07-19) -----
#
# The pool FASTA and its members manifest are written back-to-back from the same
# `selected` list in the same loop of build_per_class_reference_pools.py, so
# under normal operation they cannot disagree. A record present in the FASTA but
# carrying no manifest row therefore means the two files are OUT OF SYNC (a run
# interrupted between the two writes, or a stale manifest left beside a freshly
# written FASTA).
#
# That case must NOT be guessed at. Assuming an unmapped record is non-Berghia
# (fail-open) silently readmits the exact circular-rooting bug this helper
# exists to prevent; quietly dropping it (fail-closed) silently shrinks the
# outgroup. Both failure modes are invisible. So the helper refuses: it reports
# the disagreement and returns non-zero, and stage 04 skips that class.


def _write_desynced_pool(tmp_path: Path, fasta_records, manifest_rows):
    """Write a pool FASTA and manifest that deliberately disagree."""
    fa = tmp_path / "refs_class_C.fa"
    fa.write_text("".join(f">{sid}\n{seq}\n" for sid, seq in fasta_records))
    tsv = tmp_path / "pool_members_class_C.tsv"
    lines = ["seq_id\ttaxid\tsource"]
    lines += [f"{sid}\t9999\t{src}" for sid, src in manifest_rows]
    tsv.write_text("\n".join(lines) + "\n")
    return fa, tsv


def test_id_present_in_fasta_but_absent_from_manifest_is_an_error(tmp_path):
    fa, tsv = _write_desynced_pool(
        tmp_path,
        fasta_records=[("ref_1", "AAAA"), ("ghost_2", "CCCC"), ("ref_3", "GGGG")],
        manifest_rows=[("ref_1", "ref"), ("ref_3", "ref")],  # ghost_2 missing
    )
    out = tmp_path / "outgroup.fa"
    res = _draw(tmp_path, fa, tsv, 2, out)
    assert res.returncode != 0, (
        "an unmapped pool record means the FASTA and manifest are out of sync; "
        "the draw must refuse rather than guess its provenance"
    )


def test_manifest_mismatch_error_names_the_offending_id(tmp_path):
    fa, tsv = _write_desynced_pool(
        tmp_path,
        fasta_records=[("ref_1", "AAAA"), ("ghost_2", "CCCC")],
        manifest_rows=[("ref_1", "ref")],
    )
    out = tmp_path / "outgroup.fa"
    res = _draw(tmp_path, fa, tsv, 2, out)
    combined = res.stdout + res.stderr
    assert "ghost_2" in combined, (
        "the operator has to be able to see WHICH record is unmapped; "
        f"got: {combined!r}"
    )


def test_manifest_mismatch_does_not_leave_a_usable_outgroup(tmp_path):
    fa, tsv = _write_desynced_pool(
        tmp_path,
        fasta_records=[("ref_1", "AAAA"), ("ghost_2", "CCCC")],
        manifest_rows=[("ref_1", "ref")],
    )
    out = tmp_path / "outgroup.fa"
    _draw(tmp_path, fa, tsv, 2, out)
    assert _ids(out) == [], (
        "a refused draw must not leave a partially-filtered outgroup behind for "
        "a caller that ignores the return code"
    )


def test_berghia_only_missing_from_manifest_is_still_an_error(tmp_path):
    """The dangerous case: the unmapped record IS Berghia. Fail-open would root
    the tree on it silently. Named explicitly so the intent cannot be lost."""
    fa, tsv = _write_desynced_pool(
        tmp_path,
        fasta_records=[("berghia_x", "AAAA"), ("ref_1", "CCCC")],
        manifest_rows=[("ref_1", "ref")],  # the Berghia record has no row
    )
    out = tmp_path / "outgroup.fa"
    res = _draw(tmp_path, fa, tsv, 1, out)
    assert res.returncode != 0
    assert "berghia_x" in (res.stdout + res.stderr)


def test_fully_consistent_manifest_still_succeeds(tmp_path):
    """Regression guard: the strict check must not break the normal path."""
    fa, tsv = _write_pool(
        tmp_path,
        [("berghia_1", "berghia", "AAAA"),
         ("ANCHOR_1", "anchor", "CCCC"),
         ("ref_1", "ref", "GGGG")],
    )
    out = tmp_path / "outgroup.fa"
    res = _draw(tmp_path, fa, tsv, 2, out)
    assert res.returncode == 0, res.stderr
    assert _ids(out) == ["ANCHOR_1", "ref_1"]


def test_stage04_skips_the_class_when_the_draw_refuses(tmp_path):
    """Consistent with the empty-pool decision: skip the class, don't kill the
    stage, so one desynced sister pool cannot cost the class-A tree."""
    src = STAGE04.read_text()
    assert "draw_outgroup_fasta" in src
    idx = src.index("draw_outgroup_fasta")
    window = src[idx:idx + 700]
    assert "continue" in window, (
        "stage 04 must skip the class when draw_outgroup_fasta returns non-zero"
    )
