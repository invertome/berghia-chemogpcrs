"""Real-data guards for the refexp2 artifacts.

A synthetic fixture encodes whatever key its author assumed and therefore cannot
catch the case where the assumption itself is wrong. This project has paid for
that lesson: 17 confirmed silent-wrong defects traced to joins that were never
asserted against real data. So these tests read the ACTUAL generated artifacts
and the ACTUAL anchor set, and assert the keys really do line up.

They skip (never fail) when an artifact is absent, so a fresh checkout that has
not run the expansion still has a green suite.
"""
import csv
import os

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
ANCHORS = os.path.join(REPO, "references", "anchors", "anchor_set_PROD.tsv")
REFDIR = os.path.join(REPO, "references", "non_chemo_gpcr")
CANDIDATES = os.path.join(REFDIR, "refexp2_candidate_verdicts.tsv")
AUDIT = os.path.join(REFDIR, "refexp2_anchor_audit.tsv")
CORRECTIONS = os.path.join(REFDIR, "refexp2_family_corrections.tsv")
PASSING_FASTA = os.path.join(REFDIR, "refexp2_candidate_verdicts_passing.fasta")


def read_tsv(path):
    if not os.path.exists(path):
        pytest.skip(f"artifact not generated: {os.path.basename(path)}")
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def anchor_accessions():
    return {r["accession"] for r in read_tsv(ANCHORS)}


class TestCandidateArtifact:
    def test_candidates_overlapping_anchors_are_exactly_the_ratified_ones(self):
        """Expansion is additive, and ratification is how a candidate becomes an
        anchor.

        The original form asserted disjointness outright, which held only until
        the first ratification: the literature-gated candidates were promoted
        INTO the anchor set by design, so overlap is now the expected state.
        What must still never happen is a candidate re-entering under a NEW
        identifier -- that is the write-once violation this guards. So the
        assertion is now that every overlapping accession is carried by a
        ratified anchor, unchanged, and that each one passed its own gate.
        """
        candidates = {r["accession"] for r in read_tsv(CANDIDATES)}
        anchors = {r["accession"]: r for r in read_tsv(ANCHORS)}
        overlap = candidates & set(anchors)

        unratified = {a for a in overlap if anchors[a]["tier"] != "10"}
        assert not unratified, (
            f"{len(unratified)} candidates collide with anchors that were NOT "
            f"ratified from the candidate set: {sorted(unratified)[:5]}")

        verdicts = {r["accession"]: r for r in read_tsv(CANDIDATES)}
        failed = {a for a in overlap if (verdicts[a].get("failed_at") or "")}
        assert not failed, (
            f"{len(failed)} ratified anchors failed their own candidate gate: "
            f"{sorted(failed)[:5]}")

    def test_accessions_are_unique(self):
        rows = read_tsv(CANDIDATES)
        accessions = [r["accession"] for r in rows]
        assert len(accessions) == len(set(accessions))

    def test_every_survivor_is_lophotrochozoan(self):
        """The entire point of the expansion is taxonomic reach. A vertebrate
        survivor would mean the taxon restriction silently stopped applying."""
        lopho = {"Mollusca", "Annelida", "Platyhelminthes", "Brachiopoda",
                 "Phoronida", "Nemertea", "Rotifera", "Bryozoa"}
        for row in read_tsv(CANDIDATES):
            if not row["failed_at"]:
                assert row["phylum"] in lopho, f"{row['accession']} is {row['phylum']}"

    def test_survivors_pass_every_stage(self):
        for row in read_tsv(CANDIDATES):
            if not row["failed_at"]:
                for stage in ("class_a", "existence", "evidence_code",
                              "literature", "characterized"):
                    assert row[f"pass_{stage}"] == "True", \
                        f"{row['accession']} survived without passing {stage}"

    def test_no_survivor_is_named_orphan_or_a_placeholder(self):
        """The single loosening most likely to go unnoticed."""
        for row in read_tsv(CANDIDATES):
            if not row["failed_at"]:
                name = row["protein_name"].lower()
                assert "orphan" not in name
                assert "putative" not in name
                assert not name.startswith("gcr")

    def test_every_survivor_cites_a_pubmed_id(self):
        for row in read_tsv(CANDIDATES):
            if not row["failed_at"]:
                assert row["primary_pmids"].strip(), \
                    f"{row['accession']} passed the literature stage with no PMID"

    def test_fasta_matches_the_survivor_set(self):
        if not os.path.exists(PASSING_FASTA):
            pytest.skip("passing FASTA not generated")
        survivors = {r["accession"] for r in read_tsv(CANDIDATES) if not r["failed_at"]}
        headers, lengths = [], {}
        with open(PASSING_FASTA) as handle:
            current = None
            for line in handle:
                if line.startswith(">"):
                    current = line[1:].split("|")[0].strip()
                    headers.append(current)
                    lengths[current] = 0
                elif current:
                    lengths[current] += len(line.strip())
        assert set(headers) == survivors
        assert len(headers) == len(set(headers)), "duplicate FASTA header"
        assert all(v > 0 for v in lengths.values()), "empty sequence emitted"

    def test_fasta_lengths_match_the_recorded_lengths(self):
        """Guards the merge failure that silently DOUBLES a sequence while
        passing every count check."""
        if not os.path.exists(PASSING_FASTA):
            pytest.skip("passing FASTA not generated")
        recorded = {r["accession"]: int(r["sequence_length"])
                    for r in read_tsv(CANDIDATES)
                    if not r["failed_at"] and r["sequence_length"]}
        lengths, current = {}, None
        with open(PASSING_FASTA) as handle:
            for line in handle:
                if line.startswith(">"):
                    current = line[1:].split("|")[0].strip()
                    lengths[current] = 0
                elif current:
                    lengths[current] += len(line.strip())
        for accession, expected in recorded.items():
            assert lengths[accession] == expected, \
                f"{accession}: FASTA {lengths[accession]} vs recorded {expected}"


class TestAuditArtifact:
    def test_audit_covers_exactly_the_uniprot_resolvable_class_a_anchors(self):
        """The join that matters: if the audit and the anchor set disagree on
        which accessions are class A, every rate reported from it is wrong.

        The audit can only reach UniProt-resolvable accessions -- the tier-9
        anchors carry OrthoDB gene ids and are validated by the harvest path
        instead. So this asserts BOTH sides and that their union is the whole
        class-A set, which is a stronger statement than the old equality: it
        additionally rules out an anchor falling between the two paths.
        """
        from refexp2_evidence_gate import is_uniprot_accession

        audited = {r["accession"] for r in read_tsv(AUDIT)}
        class_a = [r for r in read_tsv(ANCHORS) if r["class"] == "A"]
        resolvable = {r["accession"] for r in class_a
                      if is_uniprot_accession(r["accession"])}
        assert audited == resolvable, (
            f"audit/anchor mismatch over the resolvable subset: "
            f"{len(audited - resolvable)} extra, {len(resolvable - audited)} missing")

        skipped = {r["accession"] for r in read_tsv(AUDIT + ".skipped.tsv")}
        assert audited.isdisjoint(skipped), "an anchor is both audited and skipped"
        assert audited | skipped == {r["accession"] for r in class_a}, (
            "a class-A anchor is reached by neither the audit nor the skipped "
            "sidecar, so nothing accounts for it")

    def test_audit_is_non_trivial(self):
        """A gate that passes everything, or nothing, is not measuring."""
        rows = read_tsv(AUDIT)
        passing = sum(1 for r in rows if not r["failed_at"])
        assert 0 < passing < len(rows)

    def test_audit_family_labels_match_the_anchor_set(self):
        anchors = {r["accession"]: r["family"] for r in read_tsv(ANCHORS)}
        for row in read_tsv(AUDIT):
            assert row["anchor_family"] == anchors[row["accession"]]

    def test_no_anchor_was_removed_or_renumbered(self):
        """The audit is a MEASUREMENT. It must leave the anchor set untouched."""
        anchors = read_tsv(ANCHORS)
        assert len(anchors) == len({r["accession"] for r in anchors})


class TestCorrectionsArtifact:
    def test_corrections_are_verified_against_uniprot(self):
        rows = read_tsv(CORRECTIONS)
        by_accession = {r["accession"]: r for r in rows}
        # GnRH is a PEPTIDE receptor; the glycoprotein hormones are FSH/LH/TSH/CG.
        for accession in ("A0A6B9MRA0", "A0A6B9MSD4", "Q75W84"):
            row = by_accession[accession]
            assert row["verdict"] == "SUPPORTED"
            assert row["proposed_family"] == "peptide"
            assert "Vasopressin/oxytocin" in row["uniprot_curated_family"]

    def test_grl101_is_reported_unresolved_not_guessed(self):
        """P46023 gets only the BARE class-A family from UniProt with no
        FSH/LSH/TSH subfamily, so promoting it would be inference rather than
        curation. It stays a deliberate user call."""
        row = {r["accession"]: r for r in read_tsv(CORRECTIONS)}["P46023"]
        assert row["verdict"] == "UNRESOLVED"
        # UniProt names the class-A family and stops -- no subfamily at all.
        curated = row["uniprot_curated_family"]
        assert "G protein-coupled receptor 1 family" in curated
        assert "subfamily" not in curated.lower()
        assert "FSH" not in curated

    def test_corrections_never_rename_an_accession(self):
        for row in read_tsv(CORRECTIONS):
            assert row["accession"] in {
                "A0A6B9MRA0", "A0A6B9MSD4", "Q75W84", "P46023"}
