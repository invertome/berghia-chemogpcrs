"""Bead 5ubd: the structural consumers must find stage 08's nested models.

Stage 08 creates cand_dir=alphafold/<candidate_id> and lets AlphaFold 3 nest its
own output below it, so models land at

    alphafold/<candidate_id>/<af3_job_name>/<af3_job_name>_model.cif

(08_structural_analysis.sh:298). Both consumers pointed AF_DIR at `alphafold/`
and iterated it at depth 1, where every entry is a DIRECTORY -- an `isfile`
guard skipped all of them, so there were ZERO iterations, no output, and exit 0.
Stage 07 then logged the or_microswitch and structural channels as "stays
dormant" while two SLURM jobs reported success.

The second defect on the same path: the candidate id is the directory name, NOT
the file stem (which is AF3's internal job name), so fixing only the traversal
would produce a join that matches nothing downstream.

These tests drive the REAL scripts with TMalign/foldseek/conda stubbed.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
UNITY = REPO_ROOT / "scripts" / "unity"
FOLDSEEK_SH = UNITY / "run_foldseek_candidates.sh"
MICROSWITCH_SH = UNITY / "run_microswitch_bw_transfer.sh"

# The candidate ids are the directory names; the AF3 job name inside each is
# deliberately different, so a filename-derived id is distinguishable from a
# directory-derived one.
CANDIDATES = {
    "BersteEVm001t1": "berste_evm001t1_af3job",
    "BersteEVm042t3": "berste_evm042t3_af3job",
    "BersteEVm117t1": "berste_evm117t1_af3job",
}


N_RESIDUES = 8


def minimal_pdb(n_residues=N_RESIDUES, offset=0.0):
    """A tiny but parseable PDB: one ALA CA per residue along the x axis."""
    lines = []
    for i in range(1, n_residues + 1):
        x = i * 3.8 + offset
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


_CIF_FIELDS = (
    "group_PDB id type_symbol label_atom_id label_alt_id label_comp_id "
    "label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code "
    "Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv "
    "auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num"
).split()


def minimal_cif(n_residues=N_RESIDUES, offset=0.0):
    """Same geometry as minimal_pdb, in the mmCIF stage 08 actually writes."""
    lines = ["data_stub", "#", "loop_"]
    lines += [f"_atom_site.{f}" for f in _CIF_FIELDS]
    for i in range(1, n_residues + 1):
        x = i * 3.8 + offset
        lines.append(
            f"ATOM {i} C CA . ALA A 1 {i} ? "
            f"{x:.3f} 0.000 0.000 1.00 20.00 {i} ALA A CA 1"
        )
    return "\n".join(lines) + "\n#\n"


def make_af_tree(root: Path, nested=True):
    """Build stage 08's alphafold/ output layout."""
    af_dir = root / "structural_analysis" / "alphafold"
    af_dir.mkdir(parents=True)
    for cand_id, job_name in CANDIDATES.items():
        if nested:
            d = af_dir / cand_id / job_name
            d.mkdir(parents=True)
            (d / f"{job_name}_model.cif").write_text(minimal_cif())
        else:
            (af_dir / f"{cand_id}_model.cif").write_text(minimal_cif())
    return af_dir


@pytest.fixture
def sandbox(tmp_path):
    """A fake REPO_ROOT + HOME so the real scripts can run unmodified."""
    home = tmp_path / "home"
    (home / ".miniconda3" / "etc" / "profile.d").mkdir(parents=True)
    (home / ".miniconda3" / "etc" / "profile.d" / "conda.sh").write_text(
        "conda() { :; }\n"
    )

    repo = tmp_path / "repo"
    results = repo / "results"
    results.mkdir(parents=True)
    (repo / "config.sh").write_text(
        f'RESULTS_DIR="{results}"\n'
        f'SCRIPTS_DIR="{repo}/scripts"\n'
    )

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()

    return {
        "tmp": tmp_path, "home": home, "repo": repo,
        "results": results, "bin": bin_dir,
    }


def run_script(script: Path, sandbox, extra_env=None):
    env = dict(os.environ)
    env["PATH"] = f"{sandbox['bin']}{os.pathsep}{env['PATH']}"
    env["HOME"] = str(sandbox["home"])
    env["REPO_ROOT"] = str(sandbox["repo"])
    env.update(extra_env or {})
    return subprocess.run(
        ["bash", str(script)], cwd=sandbox["repo"], env=env,
        capture_output=True, text=True,
    )


# --------------------------------------------------------------------------- #
# Foldseek
# --------------------------------------------------------------------------- #

def install_foldseek_stub(sandbox, query_log: Path):
    stub = sandbox["bin"] / "foldseek"
    stub.write_text(
        '#!/bin/bash\n'
        # easy-search <query> <targetDB> <out> <tmp> ...
        'echo "$2" >> "' + str(query_log) + '"\n'
        'ls -1 "$2" >> "' + str(query_log) + '"\n'
        'printf "q\\tt\\t0.9\\t0.9\\t1e-30\\n" > "$4"\n'
    )
    stub.chmod(0o755)


def test_foldseek_stages_nested_models_under_candidate_ids(sandbox):
    """The flat query dir is named by candidate id, not by AF3 job name."""
    make_af_tree(sandbox["results"])
    query_log = sandbox["tmp"] / "foldseek_calls.txt"
    install_foldseek_stub(sandbox, query_log)

    db = sandbox["repo"] / "references" / "foldseek"
    db.mkdir(parents=True)
    (db / "pdb").write_text("stub db")

    proc = run_script(FOLDSEEK_SH, sandbox, {"PDB_FOLDSEEK_DB": str(db / "pdb")})
    assert proc.returncode == 0, proc.stderr

    query_dir = sandbox["results"] / "foldseek" / "candidates" / "tmp" / "query_structures"
    staged = sorted(p.name for p in query_dir.iterdir())
    assert staged == sorted(f"{c}.cif" for c in CANDIDATES), staged

    # Each link resolves to the real nested model.
    for cand_id, job_name in CANDIDATES.items():
        target = (query_dir / f"{cand_id}.cif").resolve()
        assert target.name == f"{job_name}_model.cif"
        assert target.is_file()

    # The AF3 job name must never become the query id.
    for job_name in CANDIDATES.values():
        assert not (query_dir / f"{job_name}_model.cif").exists()

    assert "staged 3 candidate structures" in proc.stdout


def test_foldseek_searches_the_staged_dir_not_the_nested_root(sandbox):
    make_af_tree(sandbox["results"])
    query_log = sandbox["tmp"] / "foldseek_calls.txt"
    install_foldseek_stub(sandbox, query_log)

    db = sandbox["repo"] / "references" / "foldseek"
    db.mkdir(parents=True)
    (db / "pdb").write_text("stub db")

    proc = run_script(FOLDSEEK_SH, sandbox, {"PDB_FOLDSEEK_DB": str(db / "pdb")})
    assert proc.returncode == 0, proc.stderr

    logged = query_log.read_text()
    assert "query_structures" in logged
    # foldseek listed the query dir: it must contain files, not per-candidate dirs.
    for cand_id in CANDIDATES:
        assert f"{cand_id}.cif" in logged


def test_foldseek_fails_loudly_when_no_structures_exist(sandbox):
    """Empty AF tree: exit non-zero rather than 'succeed' with no output."""
    (sandbox["results"] / "structural_analysis" / "alphafold").mkdir(parents=True)
    install_foldseek_stub(sandbox, sandbox["tmp"] / "unused.txt")

    proc = run_script(FOLDSEEK_SH, sandbox)
    assert proc.returncode != 0
    assert "no candidate structures found" in proc.stderr


def test_foldseek_fails_loudly_when_no_database_is_available(sandbox):
    """All DBs missing means no search ran -- not a success."""
    make_af_tree(sandbox["results"])
    install_foldseek_stub(sandbox, sandbox["tmp"] / "unused.txt")

    proc = run_script(FOLDSEEK_SH, sandbox)
    assert proc.returncode != 0
    assert "none of the PDB/AFDB50/GPCRdb foldseek DBs were found" in proc.stderr


# --------------------------------------------------------------------------- #
# Microswitch BW transfer
# --------------------------------------------------------------------------- #

def install_tmalign_stub(sandbox):
    """TMalign <cand> <ref> -m <matrix>: identity rotation, plausible scores."""
    stub = sandbox["bin"] / "TMalign"
    stub.write_text(
        '#!/bin/bash\n'
        'cat > "$4" <<EOF\n'
        ' ------ The rotation matrix to rotate Chain_1 to Chain_2 ------\n'
        'm               t[m]        u[m,0]        u[m,1]        u[m,2]\n'
        '1     0.0000000000  1.0000000000  0.0000000000  0.0000000000\n'
        '2     0.0000000000  0.0000000000  1.0000000000  0.0000000000\n'
        '3     0.0000000000  0.0000000000  0.0000000000  1.0000000000\n'
        'EOF\n'
        'echo "TM-score= 0.85321 (if normalized by length of Chain_1)"\n'
        'echo "TM-score= 0.79210 (if normalized by length of Chain_2)"\n'
    )
    stub.chmod(0o755)


def make_reference(sandbox):
    ref = sandbox["repo"] / "references" / "structural"
    ref.mkdir(parents=True)
    path = ref / "rhodopsin_bovine_reference.pdb"
    path.write_text(minimal_pdb())
    return path


def test_microswitch_writes_one_report_per_candidate_id(sandbox):
    """Reports are keyed on the candidate id build_microswitch_channel expects."""
    make_af_tree(sandbox["results"])
    make_reference(sandbox)
    install_tmalign_stub(sandbox)

    proc = run_script(MICROSWITCH_SH, sandbox)
    assert proc.returncode == 0, proc.stderr + proc.stdout

    out_dir = sandbox["results"] / "ranking" / "microswitch" / "alignments"
    reports = sorted(p.name for p in out_dir.glob("*.bw_report.txt"))
    assert reports == sorted(f"{c}.bw_report.txt" for c in CANDIDATES), reports

    # The AF3 job name must never leak into the id.
    for job_name in CANDIDATES.values():
        assert not (out_dir / f"{job_name}_model.bw_report.txt").exists()

    assert "resolved 3 candidate structures" in proc.stdout


def test_microswitch_reports_have_the_documented_format(sandbox):
    """Sanity: the reports parse as build_microswitch_channel.py expects."""
    make_af_tree(sandbox["results"])
    make_reference(sandbox)
    install_tmalign_stub(sandbox)

    proc = run_script(MICROSWITCH_SH, sandbox)
    assert proc.returncode == 0, proc.stderr

    out_dir = sandbox["results"] / "ranking" / "microswitch" / "alignments"
    lines = (out_dir / "BersteEVm001t1.bw_report.txt").read_text().splitlines()
    assert lines[0] == "TM-score\t0.7921"
    assert lines[1] == "ref_resnum\tcand_resnum\tcand_residue"
    # Identity rotation over identical coordinates: every reference residue pairs.
    assert len(lines) == 2 + N_RESIDUES
    for line in lines[2:]:
        ref_resnum, cand_resnum, cand_residue = line.split("\t")
        assert cand_resnum == ref_resnum
        assert cand_residue == "A"


def test_microswitch_fails_loudly_when_no_structures_exist(sandbox):
    (sandbox["results"] / "structural_analysis" / "alphafold").mkdir(parents=True)
    make_reference(sandbox)
    install_tmalign_stub(sandbox)

    proc = run_script(MICROSWITCH_SH, sandbox)
    assert proc.returncode != 0
    assert "no candidate structures found" in proc.stderr


def test_microswitch_fails_loudly_when_every_alignment_fails(sandbox):
    """Structures found but zero reports written must not exit 0."""
    make_af_tree(sandbox["results"])
    make_reference(sandbox)
    failing = sandbox["bin"] / "TMalign"
    failing.write_text('#!/bin/bash\nexit 1\n')
    failing.chmod(0o755)

    proc = run_script(MICROSWITCH_SH, sandbox)
    assert proc.returncode != 0
    assert "NOT ONE produced a BW report" in proc.stderr


def test_microswitch_handles_an_already_flat_layout(sandbox):
    """A flattened AF_DIR still resolves ids from the filename stem."""
    make_af_tree(sandbox["results"], nested=False)
    make_reference(sandbox)
    install_tmalign_stub(sandbox)

    proc = run_script(MICROSWITCH_SH, sandbox)
    assert proc.returncode == 0, proc.stderr

    out_dir = sandbox["results"] / "ranking" / "microswitch" / "alignments"
    reports = sorted(p.name for p in out_dir.glob("*.bw_report.txt"))
    assert reports == sorted(f"{c}_model.bw_report.txt" for c in CANDIDATES)


# --------------------------------------------------------------------------- #
# The regression itself: depth-1 iteration finds nothing.
# --------------------------------------------------------------------------- #

def test_nested_layout_has_no_files_at_depth_one(sandbox):
    """Documents why the old `os.listdir` + isfile loop ran zero times."""
    af_dir = make_af_tree(sandbox["results"])
    entries = sorted(af_dir.iterdir())
    assert len(entries) == len(CANDIDATES)
    assert all(e.is_dir() for e in entries), \
        "stage 08 writes per-candidate DIRECTORIES, never files, at this level"
    assert not any(e.is_file() for e in entries)
