#!/bin/bash
# run_microswitch_bw_transfer.sh — structural Ballesteros-Weinstein (BW)
# number-transfer for the OR-microswitch evidence channel (Glue task G5,
# docs/plans/2026-07-01-ml-plm-chemoreceptor-ranking.md).
#
# For each candidate AlphaFold model (stage 08 output), superposes it onto
# a reference class-A GPCR structure with known BW numbering (bovine
# rhodopsin: BW 6.48 = Trp265, BW 6.50 = Pro267 -- see
# scripts/build_microswitch_channel.py's REFERENCE_BW) using TM-align, then
# derives the per-residue structural correspondence via a nearest-Calpha-
# neighbor pass (Biopython) over TM-align's OWN rotation matrix (its `-m`
# output) -- NOT from TM-align's raw text alignment block, which renumbers
# residues positionally and does not reliably preserve either structure's
# own author-assigned residue numbers. Applying TM-align's rotation
# ourselves and nearest-neighbor-pairing Calpha atoms keeps
# ref_resnum/cand_resnum in the written report exactly equal to each
# structure file's own residue numbering, with no offset-guessing -- this
# is the "Biopython superposition + Calpha nearest-neighbor mapping"
# option documented as acceptable alongside a raw TM-align/US-align
# residue alignment; here it is TM-align's rotation specifically, applied
# by Biopython, that anchors the pairing.
#
# Writes one alignment-report file per candidate, in the format
# scripts/build_microswitch_channel.py's parse_alignment() consumes (see
# that module's docstring for the full spec):
#   TM-score\t<float>
#   ref_resnum\tcand_resnum\tcand_residue
#   <ref_resnum>\t<cand_resnum>\t<cand_residue>
#   ...
#   (cand_resnum/cand_residue are "-"/"-" for a reference residue with no
#   candidate Calpha within NEIGHBOR_CUTOFF_A -- a gap.)
#
# Aligner binary: TM-align ($TMALIGN in config.sh, default "TMalign" on
# PATH) -- the same tool already used elsewhere in this pipeline for
# structural comparison (e.g. scripts/cluster_structures.py /
# scripts/plot_pca.py's tmalign_*.txt convention). If your installed
# TM-align build predates mmCIF support, convert AF3 *_model.cif candidates
# to PDB first (e.g. via Biopython or the `maxit`/`gemmi` CLI) -- this
# script does not do that conversion itself.
#
# Reference structure: $REFERENCE_STRUCTURE (env-overridable below). No
# reference structure ships in this repo -- fetch a bovine rhodopsin
# structure (e.g. RCSB PDB 1U19) yourself and VERIFY residue 265=Trp /
# residue 267=Pro in the actual file before trusting this pipeline stage --
# obtain and verify this programmatically against the file itself, never
# assume any specific deposition's numbering from memory. If you use a
# differently-numbered reference,
# override REFERENCE_BW to match it when calling
# scripts/build_microswitch_channel.py downstream (its `reference_bw`
# parameter / this script has no override for it, since REFERENCE_BW lives
# in the Python producer, not here).
#
# Downstream:
#   python3 scripts/build_microswitch_channel.py \
#       --alignments <OUT_DIR> --out results/ranking/microswitch/microswitch_channel.tsv
#
# CPU-only: structural alignment (not structure prediction) does not need
# a GPU, same reasoning as run_foldseek_candidates.sh. The per-candidate
# loop below is serial (one TM-align subprocess + one Biopython pass at a
# time) -- --time is sized generously (2-3x a rough few-hundred-candidate
# estimate) rather than tightly, per HPC etiquette; override with
# `sbatch --time=... run_microswitch_bw_transfer.sh` if your candidate set
# is much larger.
#
# Submit from the repo root:
#   sbatch scripts/unity/run_microswitch_bw_transfer.sh
#
#SBATCH --job-name=microswitch_bw_transfer
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --account=pi_pkatz_umass_edu
#SBATCH --output=logs/microswitch_bw_transfer-%j.out
#SBATCH --error=logs/microswitch_bw_transfer-%j.err

set -eo pipefail
source "$HOME/.miniconda3/etc/profile.d/conda.sh"
conda activate berghia-gpcr
set -u

REPO_ROOT="${REPO_ROOT:-/scratch3/workspace/jperezmoreno_umass_edu-jorge/chemogpcrs_2026-05}"
cd "${REPO_ROOT}"
source "${REPO_ROOT}/config.sh"

# Candidate AlphaFold models (stage 08 output). Same default as
# run_foldseek_candidates.sh.
#
# Bead 5ubd: stage 08 does NOT write one structure file per candidate at the
# top level. It creates cand_dir=alphafold/<candidate_id> and lets AF3 nest its
# own output below that, so models land at
#     alphafold/<candidate_id>/<af3_job_name>/<af3_job_name>_model.cif
# (08_structural_analysis.sh:298). Iterating this directory at depth 1 yields
# only DIRECTORIES, every one of which an `isfile` guard skips -- zero
# iterations, no reports, exit 0, and stage 07 then logs the microswitch
# channel as "stays dormant" while both jobs report success. Walk the tree.
AF_DIR="${AF_DIR:-${RESULTS_DIR}/structural_analysis/alphafold}"

# See the REFERENCE_STRUCTURE comment block above.
REFERENCE_STRUCTURE="${REFERENCE_STRUCTURE:-${REPO_ROOT}/references/structural/rhodopsin_bovine_reference.pdb}"

OUT_DIR="${OUT_DIR:-${RESULTS_DIR}/ranking/microswitch/alignments}"
TMP_DIR="${OUT_DIR}/tmp"
# Matches TM-align's own ":" (d < 5.0 Angstrom) "confidently aligned"
# pairing convention.
NEIGHBOR_CUTOFF_A="${NEIGHBOR_CUTOFF_A:-5.0}"
mkdir -p "${OUT_DIR}" "${TMP_DIR}" logs

TMALIGN_BIN="$(command -v "${TMALIGN:-TMalign}" || true)"
if [ -z "${TMALIGN_BIN}" ]; then
    echo "ERROR: TM-align not found (TMALIGN=${TMALIGN:-TMalign}). " \
         "Install via scripts/unity/install_structural_tools.sh." >&2
    exit 2
fi
if [ ! -s "${REFERENCE_STRUCTURE}" ]; then
    echo "ERROR: reference structure not found: ${REFERENCE_STRUCTURE}" >&2
    echo "       (see the REFERENCE_STRUCTURE comment near the top of this script)" >&2
    exit 2
fi
if [ ! -d "${AF_DIR}" ]; then
    echo "ERROR: candidate AlphaFold model dir not found: ${AF_DIR}" >&2
    echo "       (run stage 08_structural_analysis.sh first)" >&2
    exit 2
fi

echo "[microswitch-bw-transfer] using TM-align: ${TMALIGN_BIN}"
echo "[microswitch-bw-transfer] reference structure: ${REFERENCE_STRUCTURE}"
echo "[microswitch-bw-transfer] candidate models: ${AF_DIR} -> ${OUT_DIR}"

# --- Bead 5ubd: resolve <candidate_id> -> <model path> ------------------------
# The candidate id is the FIRST path component under AF_DIR -- that is the
# `cand_dir` stage 08 itself created (08:222, `alphafold/${candidate_id}`), and
# it is what stage 08 reads back when rendering figures (08:303,
# `basename $(dirname $(dirname $cif))`). It is NOT the file stem: the file is
# named after AF3's internal job name, so deriving the id from the filename
# produces `<af3_job_name>_model`, which joins against nothing downstream --
# build_microswitch_channel.py keys candidates on the part before the
# ".bw_report.txt" suffix and matches them against the ranking table's ids.
# Taking the first component (rather than literally the grandparent) also stays
# correct if AF3 changes how deeply it nests inside cand_dir.
#
# Deliberately NOT results/structural_analysis/all_pdb/, even though that
# directory is flat: stage 08 pools GPCRdb REFERENCE structures into it
# alongside the predictions (08:274), so it would emit bw_reports keyed on
# reference PDB accessions and contaminate the ranking join with non-candidates.
# It also flattens to the AF3 job name, losing the candidate id entirely.
MODEL_MANIFEST="${TMP_DIR}/candidate_models.tsv"
: > "${MODEL_MANIFEST}"
while IFS= read -r _model; do
    [ -n "${_model}" ] || continue
    _rel="${_model#"${AF_DIR%/}"/}"
    case "${_rel}" in
        */*) _cand_id="${_rel%%/*}" ;;                    # <cand_id>/.../<file>
        *)   _cand_id="$(basename "${_rel%.*}")" ;;       # already flat
    esac
    printf '%s\t%s\n' "${_cand_id}" "${_model}" >> "${MODEL_MANIFEST}"
done < <(find "${AF_DIR}" \( -name '*_model.cif' -o -name '*_model.pdb' \) | sort)

N_MODELS=$(wc -l < "${MODEL_MANIFEST}")
if [ "${N_MODELS}" -eq 0 ]; then
    echo "ERROR: no candidate structures found under ${AF_DIR}" >&2
    echo "       (searched recursively for *_model.cif / *_model.pdb)" >&2
    echo "       Stage 08 writes alphafold/<candidate_id>/<af3_job>/<af3_job>_model.cif;" >&2
    echo "       run 08_structural_analysis.sh first, or point AF_DIR at the tree that holds them." >&2
    echo "       Refusing to exit 0 with zero reports -- that is what made the" >&2
    echo "       or_microswitch ranking axis silently dormant (bead 5ubd)." >&2
    exit 3
fi
echo "[microswitch-bw-transfer] resolved ${N_MODELS} candidate structures"

# One TM-align run per candidate (rotation matrix + TM-score) followed by a
# Biopython nearest-Calpha-neighbor pass under that SAME rotation, so the
# residue pairing reflects TM-align's own optimal superposition rather
# than an independently re-derived one. This inline driver is the ONLY
# place TM-align or structure parsing is invoked -- everything downstream
# (scripts/build_microswitch_channel.py) is pure and unit-tested without
# needing either.
MODEL_MANIFEST="${MODEL_MANIFEST}" REFERENCE_STRUCTURE="${REFERENCE_STRUCTURE}" OUT_DIR="${OUT_DIR}" \
TMP_DIR="${TMP_DIR}" TMALIGN_BIN="${TMALIGN_BIN}" NEIGHBOR_CUTOFF_A="${NEIGHBOR_CUTOFF_A}" \
python3 - <<'PYEOF'
import os
import subprocess
import sys

from Bio.PDB import MMCIFParser, PDBParser

MODEL_MANIFEST = os.environ["MODEL_MANIFEST"]
REFERENCE_STRUCTURE = os.environ["REFERENCE_STRUCTURE"]
OUT_DIR = os.environ["OUT_DIR"]
TMP_DIR = os.environ["TMP_DIR"]
TMALIGN_BIN = os.environ["TMALIGN_BIN"]
CUTOFF = float(os.environ["NEIGHBOR_CUTOFF_A"])

# Hardcoded rather than Bio.PDB.Polypeptide.three_to_one/protein_letters_3to1
# -- that helper's location has moved across Biopython versions -- so this
# driver doesn't depend on which Biopython version happens to be installed
# in the conda env. Standard 20 amino acids only; anything else (modified
# residues, ligands, waters) is skipped, same as AlphaFold models only ever
# predicting the standard 20 anyway.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load_ca_atoms(path):
    """{resnum: (Atom, one_letter_residue)} for every standard-amino-acid
    Calpha in the first model of `path` (PDB or mmCIF, by extension)."""
    parser = (
        MMCIFParser(QUIET=True)
        if path.lower().endswith((".cif", ".mmcif"))
        else PDBParser(QUIET=True)
    )
    structure = parser.get_structure("s", path)
    out = {}
    for residue in structure[0].get_residues():
        one_letter = THREE_TO_ONE.get(residue.get_resname())
        if one_letter is None or "CA" not in residue:
            continue
        out[residue.get_id()[1]] = (residue["CA"], one_letter)
    return out


def run_tmalign(candidate_path, matrix_path):
    """TM-align (Chain_1=candidate, Chain_2=reference); returns the
    TM-score normalized by the REFERENCE (Chain_2) length -- the second
    "TM-score=" occurrence in TM-align's stdout, matching TM-align's own
    printed recommendation for structures of unequal length. None on any
    failure (missing binary error, unparseable output, etc.)."""
    result = subprocess.run(
        [TMALIGN_BIN, candidate_path, REFERENCE_STRUCTURE, "-m", matrix_path],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        return None
    scores = [
        line.split("TM-score=")[1].split()[0]
        for line in result.stdout.splitlines() if "TM-score=" in line
    ]
    if len(scores) < 2:
        return None
    try:
        return float(scores[1])
    except ValueError:
        return None


def read_rotation_matrix(matrix_path):
    """(t, u) from TM-align's `-m` rotation-matrix file: three data rows
    "m  t(m)  u(m,1)  u(m,2)  u(m,3)" such that a candidate point (x,y,z)
    maps into the reference frame via
        X(m) = t(m) + u(m,1)*x + u(m,2)*y + u(m,3)*z.
    None if the file is missing/malformed."""
    if not os.path.exists(matrix_path):
        return None
    rows = []
    with open(matrix_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 5 and parts[0] in ("1", "2", "3"):
                rows.append(parts)
    if len(rows) != 3:
        return None
    t = [0.0, 0.0, 0.0]
    u = [[0.0] * 3 for _ in range(3)]
    for parts in rows:
        m = int(parts[0]) - 1
        t[m] = float(parts[1])
        u[m] = [float(parts[2]), float(parts[3]), float(parts[4])]
    return t, u


def apply_rotation(coord, t, u):
    x, y, z = coord
    return tuple(
        t[m] + u[m][0] * x + u[m][1] * y + u[m][2] * z for m in range(3)
    )


def distance(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


reference_ca = load_ca_atoms(REFERENCE_STRUCTURE)
print(f"[microswitch-bw-transfer] reference: {len(reference_ca)} CA atoms", file=sys.stderr)

# (candidate_id, model_path) pairs, resolved by the shell above from stage 08's
# nested alphafold/<candidate_id>/<af3_job>/<af3_job>_model.cif layout. The id
# comes from the directory tree, never from the filename -- see the
# MODEL_MANIFEST comment block in the shell section.
models = []
with open(MODEL_MANIFEST) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        cand_id, _, model_path = line.partition("\t")
        if cand_id and model_path:
            models.append((cand_id, model_path))

if not models:
    sys.exit("ERROR: candidate model manifest is empty: " + MODEL_MANIFEST)

n_ok, n_fail = 0, 0
for cand_id, model_path in models:
    if not os.path.isfile(model_path):
        print(f"WARNING: model file vanished for {cand_id}: {model_path}", file=sys.stderr)
        n_fail += 1
        continue
    matrix_path = os.path.join(TMP_DIR, f"{cand_id}.tmalign_matrix.txt")
    report_path = os.path.join(OUT_DIR, f"{cand_id}.bw_report.txt")

    tmscore = run_tmalign(model_path, matrix_path)
    if tmscore is None:
        print(f"WARNING: TM-align failed for {cand_id}", file=sys.stderr)
        n_fail += 1
        continue

    rotation = read_rotation_matrix(matrix_path)
    if rotation is None:
        print(f"WARNING: could not read rotation matrix for {cand_id}", file=sys.stderr)
        n_fail += 1
        continue
    t, u = rotation

    try:
        candidate_ca = load_ca_atoms(model_path)
    except Exception as exc:  # noqa: BLE001 -- one bad structure must not kill the batch
        print(f"WARNING: could not parse structure for {cand_id}: {exc}", file=sys.stderr)
        n_fail += 1
        continue

    rotated = {
        resnum: (apply_rotation(atom.get_coord(), t, u), residue)
        for resnum, (atom, residue) in candidate_ca.items()
    }

    with open(report_path, "w") as out:
        out.write(f"TM-score\t{tmscore}\n")
        out.write("ref_resnum\tcand_resnum\tcand_residue\n")
        for ref_resnum in sorted(reference_ca):
            ref_coord = reference_ca[ref_resnum][0].get_coord()
            best_resnum, best_dist = None, None
            for cand_resnum, (cand_coord, _cand_residue) in rotated.items():
                d = distance(ref_coord, cand_coord)
                if best_dist is None or d < best_dist:
                    best_resnum, best_dist = cand_resnum, d
            if best_resnum is not None and best_dist <= CUTOFF:
                out.write(f"{ref_resnum}\t{best_resnum}\t{rotated[best_resnum][1]}\n")
            else:
                out.write(f"{ref_resnum}\t-\t-\n")
    n_ok += 1

print(f"[microswitch-bw-transfer] DONE: {n_ok} reports written, {n_fail} failed -> {OUT_DIR}",
      file=sys.stderr)

# Bead 5ubd: zero reports must never be a success. Stage 07 gates the
# or_microswitch channel on this directory being non-empty and silently logs
# "stays dormant" otherwise, so an exit-0-with-nothing-written run removes a
# scored ranking axis without any failure surfacing anywhere.
if n_ok == 0:
    sys.exit(
        f"ERROR: {len(models)} candidate structures were found but NOT ONE "
        f"produced a BW report ({n_fail} failed). Refusing to exit 0 -- the "
        f"or_microswitch ranking axis would silently stay dormant."
    )
PYEOF

echo "[microswitch-bw-transfer] next: python3 ${SCRIPTS_DIR}/build_microswitch_channel.py" \
     " --alignments ${OUT_DIR} --out ${RESULTS_DIR}/ranking/microswitch/microswitch_channel.tsv"
