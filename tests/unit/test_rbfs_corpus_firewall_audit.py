"""Bead rbfs: the train/test homology firewall audit must be able to FAIL.

`audit_maxid` is computed from $HITS. If mmseqs returns zero hits -- wrong DB,
over-strict -e, output-format drift -- then n_leaked=0, nothing is dropped,
audit_maxid=0.0000, and the script prints "OK (below cutoff)" and exits 0. A
completely non-functional firewall produces the most reassuring possible output,
and audit_maxid=0 is indistinguishable from a perfect firewall. The consequence
is homology leakage into the DAPT training corpus, which inflates every
downstream validation metric.

These tests drive the REAL scripts/unity/build_molluscan_corpus.sh with mmseqs
stubbed, exercising the three states the audit must now distinguish:
search failed / search returned nothing / search returned hits.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CORPUS_SH = REPO_ROOT / "scripts" / "unity" / "build_molluscan_corpus.sh"

MMSEQS_STUB = r"""#!/bin/bash
set -eo pipefail
cmd="$1"; shift
case "$cmd" in
  easy-linclust)
    # <raw> <reptag> <tmp> ... -- pass everything through as "representatives".
    cp "$1" "${2}_rep_seq.fasta"
    ;;
  easy-search)
    # <query> <target> <out> <tmp> ...
    FULLF="$2"; OUTF="$3"
    ids=$(grep '^>' "$FULLF" | sed 's/^>//; s/[ \t].*//')
    t1=$(echo "$ids" | sed -n '1p')
    t2=$(echo "$ids" | sed -n '2p')
    case "${STUB_SCENARIO:-normal}" in
      nofile)  : ;;                                   # exit 0, write nothing
      nohits)  : > "$OUTF" ;;                         # exit 0, empty file
      badcols) printf 'q\t%s\t0.9\t1e-5\t100\n' "$t1" > "$OUTF" ;;
      belowcutoff)
        printf 'evalq\t%s\t0.10\t50\t0.20\t0.20\t1e-5\t60\n' "$t1" > "$OUTF" ;;
      normal)
        { printf 'evalq\t%s\t1.00\t100\t1.00\t1.00\t1e-99\t500\n' "$t1"
          printf 'evalq\t%s\t0.12\t40\t0.15\t0.15\t1e-4\t60\n'   "$t2"
        } > "$OUTF" ;;
      highid_lowcov)
        # Survives legitimately: identity over the cutoff, coverage under it.
        { printf 'evalq\t%s\t1.00\t100\t1.00\t1.00\t1e-99\t500\n' "$t1"
          printf 'evalq\t%s\t0.95\t40\t0.10\t0.10\t1e-4\t60\n'    "$t2"
        } > "$OUTF" ;;
    esac
    ;;
esac
"""


def fasta(prefix, n, length=60):
    out = []
    for i in range(1, n + 1):
        out.append(f">{prefix}_{i} some description")
        out.append("M" + "A" * (length - 1))
    return "\n".join(out) + "\n"


@pytest.fixture
def sandbox(tmp_path):
    home = tmp_path / "home"
    (home / ".miniconda3" / "etc" / "profile.d").mkdir(parents=True)
    (home / ".miniconda3" / "etc" / "profile.d" / "conda.sh").write_text(
        "conda() { :; }\n"
    )

    repo = tmp_path / "repo"
    repo.mkdir()

    data = tmp_path / "data"
    data.mkdir()
    (data / "mollusca_aa.fa").write_text(fasta("MG", 3))
    prot = data / "proteomes"
    prot.mkdir()
    (prot / "Aplysia_californica.faa").write_text(fasta("APL", 3))
    (data / "berghia.faa").write_text(fasta("BST", 3))
    (data / "candidates.fa").write_text(fasta("CAND", 2))
    (data / "references.fa").write_text(fasta("REF", 2))

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "mmseqs"
    stub.write_text(MMSEQS_STUB)
    stub.chmod(0o755)

    out_dir = tmp_path / "corpus_out"

    return {"tmp": tmp_path, "home": home, "repo": repo, "data": data,
            "bin": bin_dir, "out": out_dir}


def run_corpus(sandbox, scenario="normal", extra_env=None):
    env = dict(os.environ)
    env["PATH"] = f"{sandbox['bin']}{os.pathsep}{env['PATH']}"
    env["HOME"] = str(sandbox["home"])
    env["REPO_ROOT"] = str(sandbox["repo"])
    env["STUB_SCENARIO"] = scenario
    data = sandbox["data"]
    env.update({
        "CORPUS_OUT_DIR": str(sandbox["out"]),
        "MOLLUSCAGENES_AA": str(data / "mollusca_aa.fa"),
        "INHOUSE_PROTEOMES_DIR": str(data / "proteomes"),
        "BERGHIA_PROTEOME": str(data / "berghia.faa"),
        "CANDIDATES_FA": str(data / "candidates.fa"),
        "REFERENCES_FA": str(data / "references.fa"),
        "CAP_EXEMPT_REGEX": "^NOTHING$",
    })
    env.update(extra_env or {})
    return subprocess.run(
        ["bash", str(CORPUS_SH)], cwd=sandbox["repo"], env=env,
        capture_output=True, text=True,
    )


# --------------------------------------------------------------------------- #
# STATE 1 -- the search failed (no hits file at all).
# --------------------------------------------------------------------------- #

def test_missing_hits_file_is_fatal(sandbox):
    proc = run_corpus(sandbox, scenario="nofile")
    assert proc.returncode != 0, "a search that produced no output must not exit 0"
    assert "produced no hits file" in proc.stderr
    assert "NOT verified" in proc.stderr
    assert not (sandbox["out"] / "firewalled_corpus.fasta").exists()


# --------------------------------------------------------------------------- #
# STATE 2 -- the search returned nothing.
# --------------------------------------------------------------------------- #

def test_zero_hits_is_fatal_not_a_clean_bill_of_health(sandbox):
    """The core regression: empty $HITS must never read as a perfect firewall."""
    proc = run_corpus(sandbox, scenario="nohits")
    assert proc.returncode != 0
    assert "returned ZERO hits" in proc.stderr
    # Must NOT have reported the reassuring verdict.
    assert "OK (below cutoff)" not in proc.stdout
    assert "0.0000" not in proc.stdout
    assert not (sandbox["out"] / "firewalled_corpus.fasta").exists()


def test_zero_hits_explains_why_it_is_impossible(sandbox):
    """The message must state the invariant, not just 'no hits'."""
    proc = run_corpus(sandbox, scenario="nohits")
    assert "Berghia is in BOTH" in proc.stderr


# --------------------------------------------------------------------------- #
# STATE 3 -- the search returned hits.
# --------------------------------------------------------------------------- #

def test_normal_run_drops_leaked_sequences_and_reports_evidence(sandbox):
    proc = run_corpus(sandbox, scenario="normal")
    assert proc.returncode == 0, proc.stderr

    full = (sandbox["out"] / "full_corpus.fasta").read_text().count(">")
    fire = (sandbox["out"] / "firewalled_corpus.fasta").read_text().count(">")
    assert full == 9, full
    assert fire == full - 1, "exactly the one above-threshold target is dropped"

    manifest = (sandbox["out"] / "manifest.txt").read_text()
    assert "OK (below cutoff)" in manifest
    # The verdict now carries the evidence it rests on.
    assert re.search(r"\[verified over \d+ hits, \d+ dropped, 0 breaches\]", manifest)


def test_column_count_drift_is_fatal(sandbox):
    """Positional column indexing must not silently read the wrong fields."""
    proc = run_corpus(sandbox, scenario="badcols")
    assert proc.returncode != 0
    assert "expected 8 tab-separated columns" in proc.stderr
    assert "wrong columns" in proc.stderr


def test_hits_that_drop_nothing_are_fatal(sandbox):
    """Hits present but none clearing the rule means the firewall did nothing."""
    proc = run_corpus(sandbox, scenario="belowcutoff")
    assert proc.returncode != 0
    assert "NOT ONE cleared the drop rule" in proc.stderr
    assert not (sandbox["out"] / "firewalled_corpus.fasta").exists()


def test_high_identity_low_coverage_survivor_is_reported_not_failed(sandbox):
    """The legitimate case still passes, and is still flagged honestly."""
    proc = run_corpus(sandbox, scenario="highid_lowcov")
    assert proc.returncode == 0, proc.stderr

    manifest = (sandbox["out"] / "manifest.txt").read_text()
    assert "survived on sub-threshold coverage" in manifest
    assert "0 breaches" in manifest


# --------------------------------------------------------------------------- #
# The breach branch is an invariant guard: the drop rule selects exactly what
# the audit forbids surviving, so it is unreachable through a correct drop step.
# Exercise the audit's awk program directly, extracted from the script source.
# --------------------------------------------------------------------------- #

def extract_audit_awk():
    """The audit's awk program, read out of the real script.

    Anchored on the `read -r audit_maxid audit_breaches` line: the drop-set awk
    earlier in the script opens with the same `-v mincov=...` prefix.
    """
    src = CORPUS_SH.read_text()
    anchor = src.index("read -r audit_maxid audit_breaches")
    marker = '-v mincov="$FIREWALL_MINCOV" \''
    start = src.index(marker, anchor) + len(marker)
    end = src.index("\n' \"$HITS\"", start)
    program = src[start:end]
    assert "b++" in program and 'printf("%.4f %d"' in program, program
    return program


def run_audit_awk(tmp_path, hits_rows, leaked_ids, minid="0.30", mincov="0.50"):
    hits = tmp_path / "hits.m8"
    hits.write_text("".join("\t".join(map(str, r)) + "\n" for r in hits_rows))
    leaked = tmp_path / "leaked.txt"
    leaked.write_text("".join(f"{i}\n" for i in leaked_ids))
    proc = subprocess.run(
        ["awk", "-F", "\t", "-v", f"leakedfile={leaked}",
         "-v", f"minid={minid}", "-v", f"mincov={mincov}",
         extract_audit_awk(), str(hits)],
        capture_output=True, text=True, check=True,
    )
    maxid, breaches = proc.stdout.split()
    return float(maxid), int(breaches)


def test_audit_awk_counts_a_survivor_clearing_both_thresholds(tmp_path):
    """A target that should have been dropped but wasn't is a breach."""
    rows = [
        # query, target, fident, alnlen, qcov, tcov, evalue, bits
        ("q1", "corpus_A", 0.90, 100, 0.90, 0.90, "1e-50", 300),
        ("q2", "corpus_B", 0.10, 40, 0.10, 0.10, "1e-2", 40),
    ]
    # corpus_A qualifies for dropping but is absent from the leaked set.
    maxid, breaches = run_audit_awk(tmp_path, rows, leaked_ids=[])
    assert breaches == 1
    assert maxid == pytest.approx(0.90)


def test_audit_awk_reports_no_breach_when_the_drop_step_worked(tmp_path):
    rows = [
        ("q1", "corpus_A", 0.90, 100, 0.90, 0.90, "1e-50", 300),
        ("q2", "corpus_B", 0.10, 40, 0.10, 0.10, "1e-2", 40),
    ]
    maxid, breaches = run_audit_awk(tmp_path, rows, leaked_ids=["corpus_A"])
    assert breaches == 0
    assert maxid == pytest.approx(0.10), "max is over SURVIVORS only"


def test_audit_awk_high_identity_low_coverage_is_not_a_breach(tmp_path):
    """The id-AND-cov drop rule legitimately leaves these in place."""
    rows = [("q1", "corpus_A", 0.95, 40, 0.10, 0.10, "1e-4", 60)]
    maxid, breaches = run_audit_awk(tmp_path, rows, leaked_ids=[])
    assert breaches == 0
    assert maxid == pytest.approx(0.95)


def test_audit_awk_uses_max_of_qcov_and_tcov(tmp_path):
    """max(qcov,tcov) is the conservative rule; either side can trigger."""
    rows = [("q1", "corpus_A", 0.90, 100, 0.10, 0.80, "1e-50", 300)]
    _, breaches = run_audit_awk(tmp_path, rows, leaked_ids=[])
    assert breaches == 1, "tcov alone clearing mincov must count"
