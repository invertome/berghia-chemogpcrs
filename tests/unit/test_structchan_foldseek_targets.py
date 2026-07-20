"""The structural-exclusion channel must join on the target ids Foldseek
ACTUALLY emits, and must not substring-match its family vocabulary.

Two defects, one file (scripts/structural_evidence.py):

1. **The join key was never emitted.** ``classify_hit`` looked the raw
   Foldseek target straight up in a ``family_map`` keyed on bare UniProt
   accessions (built from ``references/anchors/anchor_set.tsv``). But
   Foldseek targets come from structure databases and are *structure*
   identifiers -- ``AF-P31356-F1-model_v4`` (AFDB/AFDB50), ``1f88_A``
   (PDB), a GPCRdb file basename -- never ``P31356``. Every lookup missed,
   so a confident 0.92-TM hit to a known opsin reported ``known_other``
   and ``known_non_chemoreceptor`` -- the entire point of the channel --
   was never emitted. Downstream, ``struct_nonchemo_corrob`` was a
   constant 0 that rank aggregation still admitted as a voter.

   The pre-existing fixture in tests/unit/test_build_structural_channel.py
   hid this by writing the target as a bare accession -- a value real
   Foldseek never produces -- so the test encoded the assumed key rather
   than the real one. Every target id below is in a format a real Foldseek
   database emits.

2. **Latent twin: the family vocabulary was substring-matched.**
   ``"opsin" in "Class A (Rhodopsin)"`` is True (``Rh-odopsin``). Today's
   anchor-set vocabulary happens not to collide, but wiring the documented
   ``--gpcrdb-meta`` (GPCRdb labels its families "Class A (Rhodopsin)")
   would classify every class-A hit -- our chemoreceptors -- as a
   non-chemoreceptor and silently gut the shortlist. The match is anchored
   to the ``<coarse>`` / ``<coarse>_<subfamily>`` label grammar instead.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

# conftest.py adds scripts/ to sys.path
import build_structural_channel as bsc
import structural_evidence as se

ANCHOR_FIELDS = ["accession", "tier", "taxid", "species", "family", "class", "evidence"]


def _write_foldseek(tmp_path: Path, name: str, lines: list) -> str:
    p = tmp_path / name
    p.write_text("\n".join(lines) + ("\n" if lines else ""))
    return str(p)


def _write_anchor_set(tmp_path: Path, rows: list, name: str = "anchor_set.tsv") -> str:
    p = tmp_path / name
    with open(p, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ANCHOR_FIELDS, delimiter="\t")
        writer.writeheader()
        for accession, family in rows:
            writer.writerow({
                "accession": accession, "tier": "1", "taxid": "9913",
                "species": "Bos taurus", "family": family,
                "class": "A", "evidence": "reviewed",
            })
    return str(p)


def _hit(target, alntmscore=0.92):
    return {"target": target, "fident": 0.61, "alntmscore": alntmscore,
            "evalue": 1e-30}


# --------------------------------------------------------------------------- #
# target_keys: the parse of a real Foldseek target id into lookup keys
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("target,expected_accession", [
    # AFDB / AFDB50: `foldseek databases Alphafold/UniProt50` entry names.
    ("AF-P31356-F1-model_v4", "P31356"),
    ("AF-P08100-F1-model_v4", "P08100"),
    # Older model version, still in circulation in prebuilt DBs.
    ("AF-P08100-F1-model_v3", "P08100"),
    # Long-protein fragment entries carry a non-1 fragment number.
    ("AF-Q8WZ42-F5-model_v4", "Q8WZ42"),
    # 10-character accessions (the current UniProt format).
    ("AF-A0A0A0MRZ7-F1-model_v4", "A0A0A0MRZ7"),
    # easy-search against a createdb'd directory keeps the file extension.
    ("AF-P31356-F1-model_v4.cif.gz", "P31356"),
    ("AF-P31356-F1-model_v4.pdb", "P31356"),
    # ...and may append the chain foldseek assigned.
    ("AF-P31356-F1-model_v4_A", "P31356"),
])
def test_target_keys_recovers_uniprot_accession_from_afdb_target(target, expected_accession):
    assert expected_accession in se.target_keys(target)


@pytest.mark.parametrize("target,expected", [
    # Foldseek's prebuilt PDB database names entries `<pdbid>_<chain>`.
    ("1f88_A", ["1f88_A", "1f88"]),
    ("6cmo_R", ["6cmo_R", "6cmo"]),
    ("7e2y_A.cif.gz", ["7e2y_A", "7e2y"]),
])
def test_target_keys_recovers_pdb_id_and_chain(target, expected):
    keys = se.target_keys(target)
    for want in expected:
        assert want in keys, f"{want!r} missing from {keys!r}"


def test_target_keys_always_offers_the_raw_target_first():
    """A GPCRdb DB built with `foldseek createdb <structures/>` names targets
    by file basename, so an explicit --gpcrdb-meta mapping keyed on the
    literal target must win before any derived key."""
    keys = se.target_keys("ClassA_opsin_5dys_A.pdb")
    assert keys[0] == "ClassA_opsin_5dys_A.pdb"
    assert "ClassA_opsin_5dys_A" in keys


def test_target_keys_strips_a_leading_path():
    keys = se.target_keys("/refs/foldseek/pdb/1f88_A.pdb")
    assert "1f88_A" in keys
    assert "1f88" in keys


def test_target_keys_is_deduplicated_and_order_stable():
    keys = se.target_keys("AF-P31356-F1-model_v4")
    assert len(keys) == len(set(keys))
    assert keys == se.target_keys("AF-P31356-F1-model_v4")


@pytest.mark.parametrize("target", ["", "   ", None])
def test_target_keys_empty_target_yields_no_keys(target):
    assert se.target_keys(target) == []


# --------------------------------------------------------------------------- #
# classify_hit: the exclusion signal must actually fire on real target ids
# --------------------------------------------------------------------------- #

def test_confident_afdb_hit_to_an_opsin_is_an_exclusion_not_known_other():
    """THE finding. A 0.92-TM hit to bovine rhodopsin's AlphaFold model is
    the exact case the exclusion channel exists to catch."""
    family_map = {"P31356": "opsin"}  # anchor_set.tsv is keyed on accession

    state = se.classify_hit(_hit("AF-P31356-F1-model_v4"), family_map)

    assert state == "known_non_chemoreceptor"


def test_confident_pdb_hit_resolves_via_a_pdb_keyed_family_map():
    """PDB targets carry no accession, so they resolve through a
    target-keyed meta map -- on the id foldseek emits, chain included."""
    family_map = {"1f88": "opsin"}

    assert se.classify_hit(_hit("1f88_A"), family_map) == "known_non_chemoreceptor"


def test_confident_gpcrdb_hit_resolves_on_the_literal_basename():
    family_map = {"Homology_model_5ht2a_human_active": "aminergic_5HT"}

    state = se.classify_hit(
        _hit("Homology_model_5ht2a_human_active.pdb"), family_map)

    assert state == "known_non_chemoreceptor"


def test_unmapped_target_still_falls_back_to_known_other():
    assert se.classify_hit(_hit("AF-Q9Y5N1-F1-model_v4"), {}) == "known_other"


def test_low_tm_hit_is_still_novel_regardless_of_target_format():
    family_map = {"P31356": "opsin"}
    assert se.classify_hit(_hit("AF-P31356-F1-model_v4", alntmscore=0.31),
                           family_map) == "novel"


# --------------------------------------------------------------------------- #
# Latent twin: anchored family matching
# --------------------------------------------------------------------------- #

def test_gpcrdb_class_a_rhodopsin_label_is_not_an_opsin_exclusion():
    """`"opsin" in "Class A (Rhodopsin)"` is True. GPCRdb labels EVERY
    class-A receptor -- including our chemoreceptor candidates -- that way,
    so a substring match would exclude the entire shortlist."""
    family_map = {"1f88": "Class A (Rhodopsin)"}

    assert se.classify_hit(_hit("1f88_A"), family_map) == "known_other"


@pytest.mark.parametrize("label", [
    "Class A (Rhodopsin)",
    "class A rhodopsin-like",
    "chemoreceptor-like_opsin-adjacent-fold",
    "unclassified-gpcr",
    "",
])
def test_family_is_non_chemoreceptor_rejects_unanchored_labels(label):
    assert se.family_is_non_chemoreceptor(label) is False


@pytest.mark.parametrize("label", [
    "opsin",
    "aminergic",
    "aminergic_5HT",
    "peptide_NPY",
    "class-B-secretin",
    "class-C",
    "class-C_mGluR",
    "class-F-frizzled",
    "glycoprotein-hormone",
    "metabotropic-neurotransmitter",
    "lipid",
    "nucleotide",
])
def test_family_is_non_chemoreceptor_accepts_the_real_vocabulary(label):
    assert se.family_is_non_chemoreceptor(label) is True


# --------------------------------------------------------------------------- #
# End-to-end through the PRODUCER with realistic per-DB target ids
# --------------------------------------------------------------------------- #

def test_producer_end_to_end_with_realistic_targets_from_every_db(tmp_path):
    """The regression the original fixture could not catch: three DBs, three
    real target-id formats, and the exclusion signal must fire."""
    pdb = _write_foldseek(tmp_path, "PDB.tsv", [
        "BersteEVm000001t1\t1f88_A\t0.58\t0.71\t1e-18",
        "BersteEVm000002t1\t6cmo_R\t0.12\t0.19\t2e-01",
    ])
    afdb = _write_foldseek(tmp_path, "AFDB50.tsv", [
        "BersteEVm000001t1\tAF-P31356-F1-model_v4\t0.61\t0.92\t1e-30",
        "BersteEVm000003t1\tAF-A0A0A0MRZ7-F1-model_v4\t0.55\t0.88\t1e-24",
    ])
    gpcrdb = _write_foldseek(tmp_path, "GPCRdb.tsv", [
        "BersteEVm000004t1\tHomology_model_5ht2a_human_active.pdb\t0.49\t0.83\t1e-20",
    ])
    anchor_tsv = _write_anchor_set(tmp_path, [("P31356", "opsin")])
    meta = tmp_path / "gpcrdb_meta.tsv"
    meta.write_text(
        "target\tfamily\n"
        "Homology_model_5ht2a_human_active\taminergic_5HT\n")
    family_map = bsc.build_family_map(anchor_tsv, gpcrdb_meta=str(meta))

    channel = bsc.build_structural_channel([pdb, afdb, gpcrdb], family_map)

    # AFDB hit (0.92) beats the PDB hit (0.71) and maps to opsin -> exclusion.
    assert channel["BersteEVm000001t1"]["struct_state"] == "known_non_chemoreceptor"
    assert channel["BersteEVm000001t1"]["struct_nonchemo_corrob"] == 1
    # No confident hit anywhere -> the recall/novelty signal.
    assert channel["BersteEVm000002t1"]["struct_state"] == "novel"
    # Confident, but the accession isn't in the anchor family map.
    assert channel["BersteEVm000003t1"]["struct_state"] == "known_other"
    # GPCRdb basename resolves through --gpcrdb-meta -> exclusion.
    assert channel["BersteEVm000004t1"]["struct_state"] == "known_non_chemoreceptor"


def test_exclusion_signal_is_not_degenerate_on_realistic_targets(tmp_path):
    """The downstream harm: with the broken join every candidate got
    struct_nonchemo_corrob=0, so the column was a dead constant that rank
    aggregation still admitted as a voter. It must now vary."""
    afdb = _write_foldseek(tmp_path, "AFDB50.tsv", [
        "cand_opsin\tAF-P31356-F1-model_v4\t0.61\t0.92\t1e-30",
        "cand_orphan\tAF-Q9Y5N1-F1-model_v4\t0.44\t0.79\t1e-15",
    ])
    anchor_tsv = _write_anchor_set(tmp_path, [("P31356", "opsin")])

    channel = bsc.build_structural_channel(
        [afdb], bsc.build_family_map(anchor_tsv))

    corrob = {c: r["struct_nonchemo_corrob"] for c, r in channel.items()}
    assert len(set(corrob.values())) > 1, f"still degenerate: {corrob}"


def test_cli_end_to_end_with_realistic_afdb_target(tmp_path):
    import pandas as pd

    afdb = _write_foldseek(tmp_path, "AFDB50.tsv", [
        "cand1\tAF-P31356-F1-model_v4\t0.61\t0.92\t1e-30"])
    anchor_tsv = _write_anchor_set(tmp_path, [("P31356", "opsin")])
    out = tmp_path / "struct_channel.tsv"

    rc = bsc.main(["--foldseek-tsvs", afdb, "--anchor-set", anchor_tsv,
                   "--out", str(out)])

    assert rc == 0
    row = pd.read_csv(out, sep="\t").set_index("id").loc["cand1"]
    assert row["struct_state"] == "known_non_chemoreceptor"
    assert int(row["struct_nonchemo_corrob"]) == 1
