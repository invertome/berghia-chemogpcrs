#!/usr/bin/env python3
"""
Recover CDS for reference proteins by downloading genome assemblies
and using miniprot (protein-to-genome alignment) to extract coding sequences.

Pipeline per species:
  1. Download genome assembly from NCBI (datasets CLI)
  2. Extract protein sequences for that species from OG FASTAs
  3. Run miniprot (protein → genome) to get GFF3 with CDS coordinates
  4. Extract CDS nucleotide sequences from genome using GFF3 coordinates
  5. Validate CDS (length divisible by 3, no internal stops, matches protein)

Usage:
  # Test with one species
  python3 recover_cds_from_assemblies.py --species losc --threads 2

  # Run all species
  python3 recover_cds_from_assemblies.py --threads 4

  # Run specific list
  python3 recover_cds_from_assemblies.py --species alvmar,cobe,dipe --threads 4

  # HPC mode (skip download, assemblies already in genome_dir)
  python3 recover_cds_from_assemblies.py --threads 8 --genome-dir /path/to/genomes
"""

import argparse
import glob
import gzip
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Species → (taxid, species_name, assembly_accession)
SPECIES_MAP = {
    "alvmar": ("1491186", "Alviniconcha_marisindica", "GCA_018857735.1"),
    "amba": ("582868", "Ampullaceana_balthica", "GCA_944989445.1"),
    "anvo": ("271030", "Anisus_vortex", "GCA_949126835.1"),
    "aplcal": ("6500", "Aplysia_californica", "GCA_041379995.1"),
    "arvu": ("1028688", "Arion_vulgaris", "GCA_020796225.1"),
    "baar": ("370345", "Batillaria_attramentaria", "GCA_018292915.2"),
    "baare": ("304850", "Babylonia_areolata", "GCF_041734735.1"),
    "bigl": ("6526", "Biomphalaria_glabrata", "GCA_025434175.2"),
    "bipf": ("112525", "Biomphalaria_pfeifferi", "GCA_030265305.1"),
    "bist": ("112526", "Biomphalaria_straminea", "GCA_021533235.1"),
    "caun": ("100452", "Candidula_unifasciata", "GCA_905116865.2"),
    "chori": ("499925", "Chromodoris_orientalis", "GCA_028571245.1"),
    "chsq": ("216257", "Chrysomallon_squamiferum", "GCA_012295275.1"),
    "cobe": ("89764", "Conus_betulinus", "GCA_016801955.1"),
    "coco": ("101297", "Conus_consors", "GCA_004193615.1"),
    "cotr": ("101761", "Conus_tribblei", "GCA_001262575.1"),
    "dipe": ("2172688", "Dirona_pellucida", "GCA_030265205.1"),
    "drsu": ("2038759", "Dracogyra_subfuscus", "GCA_016106625.1"),
    "elch": ("188477", "Elysia_chlorotica", "GCA_003991915.1"),
    "elma": ("1093978", "Elysia_marginata", "GCA_019649035.1"),
    "giae": ("1735272", "Gigantopelta_aegis", "GCF_016097555.1"),
    "gima": ("703304", "Gibbula_magus", "GCA_936450465.1"),
    "goco": ("2723957", "Goniobranchus_coi", "GCA_025762795.1"),
    "gofi": ("508130", "Goniobranchus_fidelis", "GCA_028565935.1"),
    "goge": ("262603", "Goniobranchus_geometricus", "GCA_028565895.1"),
    "goku": ("262604", "Goniobranchus_kuniei", "GCA_025770095.1"),
    "gole": ("262605", "Goniobranchus_leopardus", "GCA_028566475.1"),
    "hacra": ("6455", "Haliotis_cracherodii", "GCA_022045225.1"),
    "hadi": ("42344", "Haliotis_discus_hannai", "GCA_044707095.1"),
    "hala": ("36097", "Haliotis_laevigata", "GCA_008038995.1"),
    "haru": ("6454", "Haliotis_rufescens", "GCF_003343065.1"),
    "logi": ("225164", "Lottia_gigantea", "GCF_000327385.1"),
    "losc": ("225169", "Lottia_scabra", "GCA_029955385.1"),
    "lyst": ("6523", "Lymnaea_stagnalis", "GCA_964033795.1"),
    "metu": ("55729", "Melanoides_tuberculata", "GCA_028565955.2"),
    "orid": ("2584915", "Oreohelix_idahoensis", "GCA_024509875.1"),
    "pade": ("87960", "Patella_depressa", "GCA_948474765.1"),
    "pape": ("88005", "Patella_pellucida", "GCA_917208275.1"),
    "pavu": ("6465", "Patella_vulgata", "GCF_932274485.2"),
    "phau": ("109671", "Physella_acuta", "GCF_028476545.2"),
    "phli": ("1620919", "Phorcus_lineatus", "GCA_921293015.1"),
    "ploc": ("259542", "Plakobranchus_ocellatus", "GCA_019648995.1"),
    "poca": ("400727", "Pomacea_canaliculata", "GCA_004794335.1"),
    "raau": ("52793", "Radix_auricularia", "GCA_002072015.1"),
    "rave": ("55521", "Rapana_venosa", "GCA_028751875.1"),
    "trte": ("2780533", "Tritonia_tetraquetra", "GCA_030265355.1"),
    # --- LSE annelid + non-mollusc lophotrochozoan models (added 2026-05-05) ---
    # Curated keep-list after bead -lfy (Sanger ToL OX/LR-prefix annelid
    # genomes don't have efetchable protein records and were producing
    # whole-contig garbage). These have downloadable reference assemblies
    # whose Nath et al. proteins can be miniprot-recovered to CDS.
    "pldu": ("6359",    "Platynereis_dumerilii",    "GCA_043381215.1"),
    "cate": ("283909",  "Capitella_teleta",         "GCA_000328365.1"),
    "hero": ("6412",    "Helobdella_robusta",       "GCF_000326865.2"),
    "owfu": ("6347",    "Owenia_fusiformis",        "GCA_903813345.2"),
    "digy": ("2664684", "Dimorphilus_gyrociliatus", "GCA_904063045.1"),
    "lian":   ("7574",    "Lingula_anatina",         "GCF_001039355.2"),
    "lilo":   ("88925",   "Lineus_longissimus",      "GCF_910592395.1"),
    "noge":   ("416868",  "Notospermus_geniculatus", "GCA_002633025.1"),
    # Bryozoa (added 2026-05-05 to cover the previously-missing phylum;
    # both have NCBI chromosome-level reference assemblies).
    "bune":   ("10212",   "Bugula_neritina",         "GCA_964035545.1"),
    "memmem": ("95170",   "Membranipora_membranacea", "GCA_914767715.1"),
    # NOTE — Phoronis australis (taxid 115415, GCA_055505105.1) is deliberately
    # ABSENT. It used to sit here under the code "phaust", the target of a
    # documented one-time rename:
    #     sed -i 's/^>phau_/>phaust_/' \
    #       references/nath_et_al/**/115415_Phoronis_australis.faa
    # That rename was never applied. Measured 2026-07-20: ZERO sequences
    # anywhere under references/ carry a "phaust_" prefix, while 428 phoronid
    # sequences still carry "phau_" (colliding with 1,078 Physella acuta ones).
    # The entry was therefore a PHANTOM key — species_code_lookup injected it
    # into the code map, so species_for('ref_phaust_...') confidently named a
    # species for an identifier that cannot exist.
    # The collision is registered in AMBIGUOUS_SPECIES_CODES below instead.
    # Re-adding "phaust" is only correct AFTER the rename is applied to the
    # data; validate_species_map() fails loudly on a code matching no sequence.
}

# Header-prefix codes carried by proteomes of MORE THAN ONE species.
# SPECIES_MAP maps a code to exactly one assembly, so for these codes that
# mapping is a guess: it silently picks one species' genome for both species'
# sequences. Declaring the collision here lets validate_species_map() assert
# that the register still matches the data, and lets callers refuse rather
# than resolve by curation fiat.
#
# 'phau' (measured 2026-07-20 across references/): 1,078 headers in
# 109671_Physella_acuta.faa (Gastropoda) and 428 in
# 115415_Phoronis_australis.faa (Phoronida). Identifiers here are write-once,
# so the register stays honest instead of the data being re-keyed.
AMBIGUOUS_SPECIES_CODES = {
    "phau": ("Physella_acuta", "Phoronis_australis"),
}


class SpeciesCodeCollisionError(RuntimeError):
    """A species code is claimed by >1 species; CDS recovery refuses to run.

    Distinct from ``species_code_lookup.AmbiguousSpeciesCodeError`` (which
    guards leaf-name -> species lookups) because this module cannot import
    that one: ``species_code_lookup`` imports ``SPECIES_MAP`` from here.
    Both express the same rule at different seams.
    """


def assert_species_code_unambiguous(prefix, ambiguous=None, species_map=None):
    """Refuse `prefix` if it is a DECLARED multi-species header code.

    ``get_missing_proteins`` collects every sequence whose id starts with
    ``<prefix>_``, and ``SPECIES_MAP`` resolves ``prefix`` to exactly ONE
    assembly. For a code carried by two species that combination silently
    aligns one species' proteins against the other species' genome: miniprot
    returns worse hits rather than no hits, so nothing in the run fails and
    the wrong-species CDS propagates into the dN/dS axis.

    Reads :data:`AMBIGUOUS_SPECIES_CODES` rather than testing any literal
    code, so a collision is refused the first time it is declared.

    Raises :class:`SpeciesCodeCollisionError`; returns None when clean.
    """
    amb = AMBIGUOUS_SPECIES_CODES if ambiguous is None else ambiguous
    if prefix not in amb:
        return
    smap = SPECIES_MAP if species_map is None else species_map
    species = sorted(amb[prefix])
    entry = smap.get(prefix)
    if entry:
        resolved = (f"assembly {entry[2]} ({entry[1]}, taxid {entry[0]})")
    else:
        resolved = "no assembly (the code is not in SPECIES_MAP)"
    raise SpeciesCodeCollisionError(
        f"species code {prefix!r} is claimed by {len(species)} species "
        f"({', '.join(species)}), but SPECIES_MAP resolves it to a single "
        f"{resolved}.\n"
        f"Every sequence with a {prefix!r} header prefix would be aligned "
        f"against that one genome, so the other species' proteins would have "
        f"their CDS recovered from the WRONG species' genome. miniprot returns "
        f"degraded hits rather than no hits, so this corruption is silent: no "
        f"stage fails, and the bad CDS flows on into the codon alignments and "
        f"the dN/dS axis. No CDS recovered under this code can be trusted.\n"
        f"REFUSING to run. Fix the DATA -- give one species a distinct header "
        f"prefix -- then update AMBIGUOUS_SPECIES_CODES and SPECIES_MAP. "
        f"Identifiers here are write-once, so that rename needs explicit "
        f"sign-off; it is tracked separately and must NOT be applied as a "
        f"side effect of this run.")


def assert_species_list_unambiguous(prefixes, ambiguous=None, species_map=None):
    """Pre-flight a whole run: refuse if ANY requested code is ambiguous.

    Reports every offending code at once so a data problem is fixed in one
    pass, and fires before the first download so a multi-hour run does not
    abort halfway through.
    """
    amb = AMBIGUOUS_SPECIES_CODES if ambiguous is None else ambiguous
    offenders = [p for p in prefixes if p in amb]
    if not offenders:
        return
    if len(offenders) == 1:
        assert_species_code_unambiguous(offenders[0], amb, species_map)
    raise SpeciesCodeCollisionError(
        f"{len(offenders)} requested species code(s) are claimed by more than "
        f"one species: "
        + '; '.join(f"{p} ({', '.join(sorted(amb[p]))})" for p in offenders)
        + ".\nSPECIES_MAP resolves each to ONE genome, so every one of these "
          "would recover CDS for at least one species from the wrong species' "
          "assembly -- silently. REFUSING to run. Fix the header prefixes in "
          "the data (write-once identifiers: needs explicit sign-off), then "
          "update AMBIGUOUS_SPECIES_CODES.")

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}


# ---------------------------------------------------------------------------
# Paralog-discrimination thresholds (bead -8st, 2026-05-08).
# ---------------------------------------------------------------------------
# The pipeline targets chemoreceptor LSE expansions where paralogs commonly
# share 80–95 % protein identity. The previous 0.5 translation-vs-protein
# identity threshold was useless for paralog discrimination — a wrong-paralog
# miniprot hit easily clears 50 %, attaches the *other* paralog's CDS to
# this protein's identity, and silently poisons the dN/dS axis.
#
# CDS_PROTEIN_IDENTITY_MIN
#   Minimum identity of the back-translated CDS to the query protein, after
#   global pairwise alignment, for the CDS to be accepted. 0.95 is high
#   enough to reject paralog matches at typical LSE-expansion divergence
#   while tolerating routine back-translation artifacts (rare wobble codons,
#   miniprot frameshift handling). Override only when intentionally working
#   on more divergent gene families.
#
# CDS_AMBIGUITY_MARGIN
#   When a query has multiple miniprot hits, the best hit's identity must
#   exceed the second-best's by at least this margin. If the top two hits
#   are closer than this, the assignment is ambiguous and the CDS is
#   rejected — we'd rather have no CDS than the wrong CDS.
#
# CDS_OVERLAP_FRACTION
#   When two query proteins map to overlapping genomic intervals, only the
#   query with higher identity at that locus keeps the CDS (the other was
#   matching its paralog's locus). Two intervals "overlap" when the smaller
#   is at least this fraction contained in the larger.
CDS_PROTEIN_IDENTITY_MIN = float(
    os.environ.get('CDS_PROTEIN_IDENTITY_MIN', '0.95'))
CDS_AMBIGUITY_MARGIN = float(
    os.environ.get('CDS_AMBIGUITY_MARGIN', '0.05'))
CDS_OVERLAP_FRACTION = float(
    os.environ.get('CDS_OVERLAP_FRACTION', '0.5'))


def _alignment_identity(seq1: str, seq2: str) -> float:
    """Return identity = matches / aligned-length from a global pairwise
    alignment of `seq1` vs `seq2`.

    Uses Bio.Align.PairwiseAligner (Needleman-Wunsch) with a BLOSUM-style
    scoring tilt so that gaps and mismatches both reduce identity. The
    earlier code zip-compared `seq1[:N]` vs `seq2[:N]` position-by-position,
    which was wrong for any pair that was not perfectly co-linear from
    position 0 (signal peptide differences, miniprot frame-shift handling,
    etc.) — a problem that became invisible noise in the dN/dS axis.

    Returns 0.0 on either empty input.
    """
    if not seq1 or not seq2:
        return 0.0
    try:
        from Bio.Align import PairwiseAligner
    except ImportError:
        # Last-resort fallback: position-by-position over the shorter length.
        # Strictly worse than alignment-based identity but better than crash.
        n = min(len(seq1), len(seq2))
        if n == 0:
            return 0.0
        return sum(1 for a, b in zip(seq1[:n], seq2[:n]) if a == b) / n
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    try:
        alignments = aligner.align(seq1, seq2)
    except Exception:
        return 0.0
    # Bio.Align.Alignments.__len__ enumerates the internal _paths and
    # overflows int64 for pairs with very many tied optimal alignments
    # (highly divergent paralogs, repetitive low-complexity regions —
    # this killed the bead -8st re-extract chain on Physella acuta).
    # next(iter(...)) yields the first path lazily, without enumeration.
    try:
        aln = next(iter(alignments))
    except StopIteration:
        return 0.0
    # Modern Biopython: aln[0] / aln[1] return the aligned sequences as
    # strings with '-' gaps. Older `str(aln)` text format is too unstable
    # to parse reliably.
    try:
        aligned1 = str(aln[0])
        aligned2 = str(aln[1])
    except (TypeError, IndexError):
        # Fallback for very old Biopython where indexing isn't supported.
        text = str(aln).split('\n')
        if len(text) < 3:
            return 0.0
        aligned1, aligned2 = text[0], text[2]
    if not aligned1:
        return 0.0
    matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-')
    return matches / max(len(aligned1), 1)


def _intervals_overlap(coords_a: list[tuple[int, int]],
                       coords_b: list[tuple[int, int]],
                       overlap_fraction: float = CDS_OVERLAP_FRACTION) -> bool:
    """Test whether the genomic spans covered by `coords_a` and `coords_b`
    overlap by at least `overlap_fraction` of the smaller span.

    Used to detect two query proteins whose miniprot hits land at the same
    genomic locus — a sign that one query is matching the other's paralog.
    """
    if not coords_a or not coords_b:
        return False
    min_a = min(s for s, _ in coords_a)
    max_a = max(e for _, e in coords_a)
    min_b = min(s for s, _ in coords_b)
    max_b = max(e for _, e in coords_b)
    span_a = max_a - min_a + 1
    span_b = max_b - min_b + 1
    overlap = max(0, min(max_a, max_b) - max(min_a, min_b) + 1)
    smaller = min(span_a, span_b)
    if smaller <= 0:
        return False
    return (overlap / smaller) >= overlap_fraction


def parse_fasta(path):
    """Parse FASTA file, return dict of id -> sequence."""
    seqs = {}
    current_id = None
    current_seq = []
    opener = gzip.open if str(path).endswith('.gz') else open
    mode = 'rt' if str(path).endswith('.gz') else 'r'
    with opener(path, mode) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            elif current_id:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = ''.join(current_seq)
    return seqs


def write_fasta(seqs, path):
    """Write dict of id -> sequence to FASTA file."""
    with open(path, 'w') as f:
        for sid, seq in seqs.items():
            f.write(f'>{sid}\n{seq}\n')


def get_missing_proteins(og_dir, cds_file, species_prefix):
    """Get protein IDs for a species that lack CDS.

    The Nath et al. reference tree uses ``.faa`` and stores files in
    nested subdirectories (``lse/<phylum>/<taxid>_<species>.faa``,
    ``one_to_one_ortholog/OG*.faa``). We recurse through ``og_dir`` and
    accept the common protein-FASTA extensions ``.faa`` / ``.fa`` /
    ``.fasta``.

    Refuses a ``species_prefix`` declared in :data:`AMBIGUOUS_SPECIES_CODES`:
    the prefix match below is exactly where two species' sequences get fused
    into one query batch. This function is importable and callable on its own,
    so it carries its own refusal rather than relying on ``process_species``.
    """
    assert_species_code_unambiguous(species_prefix)

    # Load existing CDS IDs
    cds_ids = set()
    if os.path.exists(cds_file):
        with open(cds_file) as f:
            for line in f:
                if line.startswith('>'):
                    cds_ids.add(line[1:].split()[0])

    # Find proteins for this species in OG FASTAs (recursive; multi-extension)
    missing = {}
    seen_paths = set()
    for ext in ('faa', 'fa', 'fasta'):
        for fa_path in glob.glob(os.path.join(og_dir, '**', f'*.{ext}'),
                                 recursive=True):
            if fa_path in seen_paths:
                continue
            seen_paths.add(fa_path)
            seqs = parse_fasta(fa_path)
            for sid, seq in seqs.items():
                if sid.startswith(species_prefix + '_') and sid not in cds_ids:
                    missing[sid] = seq

    return missing


def download_assembly(accession, genome_dir, threads=2):
    """Download genome assembly using NCBI datasets CLI."""
    genome_fa = os.path.join(genome_dir, f'{accession}.fa')
    if os.path.exists(genome_fa) and os.path.getsize(genome_fa) > 0:
        return genome_fa

    zip_path = os.path.join(genome_dir, f'{accession}.zip')

    # Download
    cmd = [
        'datasets', 'download', 'genome', 'accession', accession,
        '--include', 'genome',
        '--filename', zip_path
    ]
    print(f'  Downloading {accession}...')
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        print(f'  ERROR downloading {accession}: {r.stderr[:200]}')
        return None

    # Extract
    import zipfile
    with zipfile.ZipFile(zip_path) as zf:
        # Find the genome FASTA inside the zip
        fasta_files = [n for n in zf.namelist() if n.endswith('.fna')]
        if not fasta_files:
            print(f'  ERROR: no .fna found in {zip_path}')
            return None

        # Concatenate all FASTA files (some assemblies split by chromosome)
        with open(genome_fa, 'w') as out:
            for fn in fasta_files:
                with zf.open(fn) as f:
                    for line in f:
                        out.write(line.decode())

    # Clean up zip
    os.remove(zip_path)

    size_mb = os.path.getsize(genome_fa) / 1e6
    print(f'  Downloaded {accession}: {size_mb:.0f} MB')
    return genome_fa


def run_miniprot(genome_fa, proteins_fa, output_gff, threads=2):
    """Run miniprot: protein → genome alignment."""
    if os.path.exists(output_gff) and os.path.getsize(output_gff) > 0:
        return True

    cmd = [
        'miniprot', '-t', str(threads),
        '--gff', genome_fa, proteins_fa
    ]
    print(f'  Running miniprot ({threads} threads)...')
    with open(output_gff, 'w') as out:
        r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, timeout=7200)

    if r.returncode != 0:
        print(f'  ERROR miniprot: {r.stderr[:200]}')
        return False

    return True


def extract_cds_from_gff(genome_fa, gff_path, protein_seqs):
    """
    Parse miniprot GFF3 output and extract CDS nucleotide sequences
    from the genome using exon coordinates.

    Returns dict of protein_id -> CDS nucleotide sequence.
    """
    # Load genome sequences
    genome = parse_fasta(genome_fa)

    # Parse GFF3 to get CDS exon coordinates per mRNA
    # miniprot GFF3 format:
    #   mRNA line has Target=PROTEIN_ID ...
    #   CDS lines follow with coordinates
    mrnas = {}  # mRNA_id -> {chrom, strand, target, cds_coords: [(start, end), ...]}
    current_mrna = None

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = parts

            if feature == 'mRNA':
                # Extract Target protein ID and mRNA ID
                target_match = re.search(r'Target=(\S+)', attrs)
                id_match = re.search(r'ID=([^;]+)', attrs)
                if target_match and id_match:
                    mrna_id = id_match.group(1)
                    target_id = target_match.group(1)
                    identity_match = re.search(r'Identity=([0-9.]+)', attrs)
                    identity = float(identity_match.group(1)) if identity_match else 0
                    mrnas[mrna_id] = {
                        'chrom': chrom,
                        'strand': strand,
                        'target': target_id,
                        'identity': identity,
                        'score': float(score) if score != '.' else 0,
                        'cds_coords': []
                    }
                    current_mrna = mrna_id

            elif feature == 'CDS' and current_mrna:
                parent_match = re.search(r'Parent=([^;]+)', attrs)
                if parent_match:
                    parent = parent_match.group(1)
                    if parent in mrnas:
                        mrnas[parent]['cds_coords'].append(
                            (int(start), int(end))
                        )

    # For each protein, pick the best miniprot hit and extract CDS
    # Group hits by target protein
    hits_by_protein = defaultdict(list)
    for mrna_id, info in mrnas.items():
        hits_by_protein[info['target']].append((mrna_id, info))

    cds_seqs = {}
    # Track per-protein meta for the second pass (genomic-overlap dedup) and
    # the manifest. accepted_meta[prot_id] = (chrom, sorted_coords, identity)
    accepted_meta = {}
    stats = {
        'found': 0, 'no_hit': 0, 'bad_cds': 0, 'multi_hit': 0,
        'ambiguous_top_hits': 0,    # bead -8st: best/2nd-best within margin
        'low_identity': 0,          # bead -8st: translation-vs-protein < threshold
        'paralog_overlap_rejected': 0,  # bead -8st: another query won this locus
    }

    for prot_id, prot_seq in protein_seqs.items():
        if prot_id not in hits_by_protein:
            stats['no_hit'] += 1
            continue

        hits = hits_by_protein[prot_id]
        if len(hits) > 1:
            stats['multi_hit'] += 1

        # Pick best hit by identity, then score
        hits.sort(key=lambda x: (x[1]['identity'], x[1]['score']), reverse=True)
        best_mrna_id, best_info = hits[0]

        # Bead -8st: paralog-discrimination ambiguity gate. If multiple hits
        # are within CDS_AMBIGUITY_MARGIN of the best by identity, miniprot
        # cannot reliably tell which paralog the query corresponds to; refuse
        # to assign the CDS rather than risk a wrong-paralog attribution.
        if len(hits) > 1:
            second_id = hits[1][1].get('identity', 0.0)
            if (best_info['identity'] - second_id) < CDS_AMBIGUITY_MARGIN:
                stats['ambiguous_top_hits'] += 1
                continue

        chrom = best_info['chrom']
        strand = best_info['strand']
        coords = sorted(best_info['cds_coords'])

        if chrom not in genome:
            stats['bad_cds'] += 1
            continue

        chrom_seq = genome[chrom]

        # Extract and concatenate CDS exons
        cds_parts = []
        for start, end in coords:
            # GFF is 1-based, inclusive
            exon_seq = chrom_seq[start - 1:end]
            cds_parts.append(exon_seq)

        cds = ''.join(cds_parts)

        # Reverse complement if on minus strand
        if strand == '-':
            cds = reverse_complement(cds)

        # Validate CDS
        cds_upper = cds.upper()
        if len(cds_upper) < 9:  # minimum viable CDS
            stats['bad_cds'] += 1
            continue

        if len(cds_upper) % 3 != 0:
            # Try trimming to nearest codon boundary
            cds_upper = cds_upper[:len(cds_upper) - (len(cds_upper) % 3)]
            cds = cds[:len(cds) - (len(cds) % 3)]

        # Check for internal stop codons (excluding final codon)
        codons = [cds_upper[i:i+3] for i in range(0, len(cds_upper) - 3, 3)]
        if any(c in STOP_CODONS for c in codons):
            stats['bad_cds'] += 1
            continue

        # Bead -8st: tightened translation-vs-protein identity check.
        # Previous code zip-compared seq1[:N] vs seq2[:N] from position 0
        # at a 0.5 threshold. For chemoreceptor LSE paralogs sharing
        # 80–95 % identity, both pieces were wrong: position-by-position
        # comparison overweights co-linear regions (signal peptides,
        # conserved TM helices) and the threshold doesn't discriminate
        # paralog-on-paralog confusion. Use a proper global alignment and
        # the higher CDS_PROTEIN_IDENTITY_MIN.
        translated = ''.join(CODON_TABLE.get(c, 'X') for c in
                           [cds_upper[i:i+3] for i in range(0, len(cds_upper), 3)])
        translated = translated.rstrip('*')
        prot_clean = prot_seq.replace('-', '').replace('*', '')

        identity = 0.0
        if len(translated) > 0 and len(prot_clean) > 0:
            identity = _alignment_identity(translated, prot_clean)
            if identity < CDS_PROTEIN_IDENTITY_MIN:
                stats['low_identity'] += 1
                continue

        cds_seqs[prot_id] = cds
        accepted_meta[prot_id] = (chrom, coords, identity)
        stats['found'] += 1

    # Bead -8st: post-process layer for genomic-locus dedup. If two query
    # proteins were both accepted with overlapping coords on the same
    # chromosome, only the higher-identity one keeps the CDS — the other
    # was matching its paralog's locus and the wrong-paralog CDS would
    # poison dN/dS. We process accepted proteins in order of decreasing
    # identity and reject any later one whose locus is already claimed.
    by_chrom = defaultdict(list)  # chrom -> list of (prot_id, coords, identity)
    for prot_id, (chrom, coords, identity) in accepted_meta.items():
        by_chrom[chrom].append((prot_id, coords, identity))
    for chrom, entries in by_chrom.items():
        entries.sort(key=lambda x: x[2], reverse=True)
        kept_intervals: list[tuple[str, list[tuple[int, int]]]] = []
        for prot_id, coords, _ident in entries:
            if any(_intervals_overlap(coords, kept_coords)
                   for _, kept_coords in kept_intervals):
                cds_seqs.pop(prot_id, None)
                accepted_meta.pop(prot_id, None)
                stats['paralog_overlap_rejected'] += 1
                stats['found'] -= 1
            else:
                kept_intervals.append((prot_id, coords))

    return cds_seqs, stats


# --- SPECIES_MAP drift validation ------------------------------------------
#
# SPECIES_MAP is hand-maintained against auto-generated header prefixes, so it
# drifts silently: every consumer looks a code up with a `.get()`-shaped miss
# (get_missing_proteins finds no sequences; species_code_lookup injects the
# entry regardless), which means a stale key and a wrong key both sail through
# without ever failing a run. The functions below derive the truth from the
# proteomes themselves and refuse a map that no longer matches.
#
# Deliberately NOT checked: codes present in the data but absent from
# SPECIES_MAP. SPECIES_MAP is a curated subset of the species that have a
# recoverable genome assembly, not a census of the reference set. Measured
# 2026-07-20: 57 mapped codes vs 238 observed header codes, so 182 unmapped
# codes are correct. This is the one place PREFIX_TO_GROUP's checklist must NOT
# be ported verbatim — its map is meant to be total, this one is not.

# <taxid>_<Genus_species[...]>.faa — the filename is the authoritative species
# label (same convention species_code_lookup.py reads).
_SPECIES_FAA_RE = re.compile(r'(\d+)_([A-Z][A-Za-z]+_[a-z0-9_]+)\.(?:faa|fasta)$')


class SpeciesMapDriftError(AssertionError):
    """SPECIES_MAP no longer matches the reference proteomes on disk."""


class ObservedSpecies:
    """What the proteomes actually say about one header-prefix code."""

    def __init__(self, species, count, ambiguous_over, files):
        self.species = species                  # None when ambiguous
        self.count = count                      # headers carrying this code
        self.ambiguous_over = ambiguous_over    # set of species claiming it
        self.files = files                      # basenames of claiming proteomes

    def __repr__(self):
        return (f'ObservedSpecies(species={self.species!r}, count={self.count}, '
                f'ambiguous_over={sorted(self.ambiguous_over)})')


def scan_reference_species(references_dir=None):
    """Scan reference proteomes -> ``{code: ObservedSpecies}``.

    Only files whose basename matches ``<taxid>_<Genus_species>`` are in scope;
    Swiss-Prot/outgroup references use pipe-delimited headers that carry no
    species code and belong to no assembly.
    """
    root = Path(references_dir) if references_dir else PROJECT_ROOT / 'references'
    per_species = defaultdict(lambda: defaultdict(int))   # code -> {species: n}
    files = defaultdict(set)                              # code -> {basename}
    patterns = ('*.faa', '*.fasta')
    paths = sorted({p for pat in patterns
                    for p in glob.glob(os.path.join(str(root), '**', pat),
                                       recursive=True)})
    for faa in paths:
        m = _SPECIES_FAA_RE.search(os.path.basename(faa))
        if not m:
            continue
        species = m.group(2)
        with open(faa) as fh:
            for line in fh:
                if line.startswith('>'):
                    code = line[1:].split()[0].split('_')[0]
                    per_species[code][species] += 1
                    files[code].add(os.path.basename(faa))
    observed = {}
    for code, counts in per_species.items():
        claimed = set(counts)
        observed[code] = ObservedSpecies(
            species=next(iter(claimed)) if len(claimed) == 1 else None,
            count=sum(counts.values()),
            ambiguous_over=claimed if len(claimed) > 1 else set(),
            files=files[code])
    return observed


def validate_species_map(references_dir=None, species_map=None, ambiguous=None):
    """Raise :class:`SpeciesMapDriftError` unless the map matches the proteomes.

    Checks, in the order they have historically broken:
      1. dead keys    — a mapped code matching ZERO sequences (the 'phaust'
                        phantom: a rename target recorded as though applied)
      2. mis-named    — a mapped code whose proteome names another species
      3. collisions   — a code claimed by two species must be DECLARED
                        ambiguous, since SPECIES_MAP resolves it to one genome
      4. stale/changed register — a declared collision that no longer exists,
                        or whose membership changed (e.g. a third species)

    Returns the ``{code: ObservedSpecies}`` scan on success.

    Not run at import: it needs ``references/`` on disk, which the CDS-recovery
    job on the cluster does not necessarily have. It runs via
    ``--validate-species-map`` and in the test suite against the real
    reference directory.
    """
    smap = SPECIES_MAP if species_map is None else species_map
    amb = dict(AMBIGUOUS_SPECIES_CODES if ambiguous is None else ambiguous)
    observed = scan_reference_species(references_dir)
    real_amb = {c: o.ambiguous_over for c, o in observed.items() if o.ambiguous_over}
    problems = []

    dead = sorted(set(smap) - set(observed))
    if dead:
        problems.append(
            f'{len(dead)} mapped code(s) match ZERO sequences (stale or '
            f'half-applied rename): {dead}')

    misnamed = [(c, smap[c][1], observed[c].species)
                for c in sorted(set(smap) & set(observed))
                if observed[c].species is not None
                and smap[c][1] != observed[c].species]
    if misnamed:
        problems.append('mis-named code(s): ' + ', '.join(
            f'{c} mapped {m!r} but proteome is {o!r}' for c, m, o in misnamed))

    undeclared = sorted(set(real_amb) - set(amb))
    if undeclared:
        problems.append(
            'code(s) claimed by >1 species but not declared ambiguous — '
            'SPECIES_MAP resolves each to ONE genome, so the other species\' '
            'proteins would be recovered against the wrong assembly: '
            + ', '.join(f'{c} {sorted(real_amb[c])}' for c in undeclared))

    stale = sorted(set(amb) - set(real_amb))
    if stale:
        problems.append(
            f'code(s) declared ambiguous but no longer collide in the data '
            f'(rename applied? update the register): {stale}')

    changed = [(c, sorted(amb[c]), sorted(real_amb[c]))
               for c in sorted(set(amb) & set(real_amb))
               if set(amb[c]) != set(real_amb[c])]
    if changed:
        problems.append('code(s) declared ambiguous over the wrong species: '
                        + ', '.join(f'{c} declared {d} but data says {r}'
                                    for c, d, r in changed))

    if problems:
        raise SpeciesMapDriftError(
            'SPECIES_MAP is out of sync with the reference proteomes:\n  - '
            + '\n  - '.join(problems))
    return observed


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]


def process_species(prefix, og_dir, cds_file, output_dir, genome_dir,
                    threads=2, reextract_only=False):
    """Process one species: download genome, run miniprot, extract CDS.

    When `reextract_only` is True, skip miniprot (re-using the cached
    GFF and genome at output_dir + genome_dir) and just re-run
    extract_cds_from_gff. Used to apply tightened paralog-discrimination
    thresholds (bead -8st) to a previously-recovered set without paying
    the multi-hour miniprot cost again.

    Raises :class:`SpeciesCodeCollisionError` for a code declared in
    :data:`AMBIGUOUS_SPECIES_CODES`. This is a hard refusal, NOT a skip:
    every other early exit here returns 0, which a caller cannot distinguish
    from "this species had no missing proteins" -- a refused species that
    looks like a completed one is the same silent-corruption defect wearing a
    different costume. Both the miniprot and the ``--reextract-only`` branch
    resolve the code to a genome, so the guard sits above both.
    """
    assert_species_code_unambiguous(prefix)

    if prefix not in SPECIES_MAP:
        print(f'  SKIP {prefix}: not in species map')
        return 0

    taxid, species_name, accession = SPECIES_MAP[prefix]
    print(f'\n{"="*60}')
    print(f'Processing {prefix} ({species_name})')
    print(f'  TaxID: {taxid}, Assembly: {accession}')

    # Get missing protein IDs for this species
    missing = get_missing_proteins(og_dir, cds_file, prefix)
    if not missing:
        print(f'  No missing proteins for {prefix}, skipping')
        return 0
    print(f'  Missing proteins: {len(missing)}')

    proteins_fa = os.path.join(output_dir, f'{prefix}_query_proteins.fa')
    gff_path = os.path.join(output_dir, f'{prefix}_miniprot.gff')

    if reextract_only:
        if not os.path.exists(gff_path):
            print(f'  SKIP {prefix}: --reextract-only set but no cached GFF '
                  f'at {gff_path}')
            return 0
        # Locate the genome FASTA. It may have been deleted to save disk
        # (default behavior of the orchestration is to clean up genomes
        # post-recovery), so we accept either an existing file or fall
        # back to redownloading just for the extraction.
        genome_fa = os.path.join(genome_dir, f'{accession}.fa')
        if not os.path.exists(genome_fa):
            print(f'  Re-downloading genome for re-extraction…')
            genome_fa = download_assembly(accession, genome_dir, threads)
            if not genome_fa:
                print(f'  FAILED: cannot fetch genome for re-extraction')
                return 0
        print(f'  Re-extracting CDS from cached GFF under tightened '
              f'paralog-discrimination thresholds…')
    else:
        # Write protein queries to temp FASTA
        write_fasta(missing, proteins_fa)

        # Download assembly
        genome_fa = download_assembly(accession, genome_dir, threads)
        if not genome_fa:
            print(f'  FAILED: could not download assembly for {prefix}')
            return 0

        # Run miniprot
        if not run_miniprot(genome_fa, proteins_fa, gff_path, threads):
            print(f'  FAILED: miniprot error for {prefix}')
            return 0

    # Extract CDS
    print(f'  Extracting CDS from miniprot alignments...')
    cds_seqs, stats = extract_cds_from_gff(genome_fa, gff_path, missing)

    print(f'  Results: {stats["found"]} CDS recovered, '
          f'{stats["no_hit"]} no hit, {stats["bad_cds"]} bad CDS, '
          f'{stats["multi_hit"]} multi-hit')
    # Bead -8st: paralog-discrimination rejection breakdown
    print(f'    paralog filters: {stats.get("ambiguous_top_hits", 0)} ambiguous '
          f'(best/2nd-best within {CDS_AMBIGUITY_MARGIN} identity), '
          f'{stats.get("low_identity", 0)} below {CDS_PROTEIN_IDENTITY_MIN} '
          f'translation identity, {stats.get("paralog_overlap_rejected", 0)} '
          f'paralog-locus collisions')

    # Write species CDS
    if cds_seqs:
        species_cds = os.path.join(output_dir, f'{prefix}_cds.fna')
        write_fasta(cds_seqs, species_cds)
        print(f'  Wrote {len(cds_seqs)} CDS to {species_cds}')

    # Clean up query proteins file (no-op in --reextract-only mode where
    # we never wrote it)
    if os.path.exists(proteins_fa):
        os.remove(proteins_fa)

    return len(cds_seqs)


def main():
    # Standalone map audit. Checked before argparse so it can run without the
    # recovery-run arguments (--og-dir/--cds-file/--output-dir), which are
    # required=True and irrelevant to validating the map.
    if '--validate-species-map' in sys.argv[1:]:
        try:
            observed = validate_species_map()
        except SpeciesMapDriftError as exc:
            print(exc, file=sys.stderr)
            return 1
        print(f'SPECIES_MAP OK: {len(SPECIES_MAP)} mapped code(s) all match '
              f'real sequences; {len(observed)} code(s) observed; '
              f'{len(AMBIGUOUS_SPECIES_CODES)} declared ambiguous.')
        return 0

    parser = argparse.ArgumentParser(
        description='Recover CDS from genome assemblies via miniprot')
    parser.add_argument('--validate-species-map', action='store_true',
                       help='Audit SPECIES_MAP against the reference proteomes '
                            'and exit (no recovery is run). Fails if a mapped '
                            'code matches zero sequences, names the wrong '
                            'species, or collides with another species.')
    parser.add_argument('--species', type=str, default=None,
                       help='Comma-separated species prefixes (default: all)')
    parser.add_argument('--threads', type=int, default=2,
                       help='CPU threads for miniprot (default: 2)')
    parser.add_argument('--og-dir', type=str, required=True,
                       help='Directory with orthogroup FASTA files')
    parser.add_argument('--cds-file', type=str, required=True,
                       help='Existing CDS file to check for already-fetched sequences')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory for per-species CDS files')
    parser.add_argument('--genome-dir', type=str, default=None,
                       help='Directory for genome downloads (default: output-dir/genomes)')
    parser.add_argument('--keep-genomes', action='store_true',
                       help='Keep genome files after processing (default: delete to save space)')
    parser.add_argument('--manifest', type=str, default=None,
                       help='Optional CSV manifest of per-sequence CDS provenance '
                            '(seq_id, source=miniprot, species_prefix, assembly_accession). '
                            'Rank_candidates / parse_absrel can join on this to flag '
                            'miniprot-recovered branches in dN/dS results (literature '
                            'review noted potential frameshift/splice-site bias).')
    parser.add_argument('--reextract-only', action='store_true',
                        help='Skip miniprot. Re-run CDS extraction against '
                             'existing cached GFFs in --output-dir. Used to '
                             'apply the bead -8st paralog-discrimination '
                             'thresholds to a previously-recovered set '
                             'without paying the multi-hour miniprot cost. '
                             'Genome FASTAs are re-downloaded if they were '
                             'cleaned up after the original recovery.')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    genome_dir = args.genome_dir or os.path.join(args.output_dir, 'genomes')
    os.makedirs(genome_dir, exist_ok=True)

    # Determine species to process
    if args.species:
        species_list = args.species.split(',')
    else:
        species_list = sorted(SPECIES_MAP.keys())

    # Pre-flight the whole requested set before the first download, so an
    # ambiguous code aborts at second 0 rather than after hours of miniprot on
    # the species that happened to sort earlier.
    assert_species_list_unambiguous(species_list)

    print(f'CDS Recovery Pipeline')
    print(f'  Species: {len(species_list)}')
    print(f'  Threads: {args.threads}')
    print(f'  OG dir: {args.og_dir}')
    print(f'  Output: {args.output_dir}')
    print(f'  Genomes: {genome_dir}')

    total_recovered = 0
    results = []

    for prefix in species_list:
        n = process_species(
            prefix, args.og_dir, args.cds_file,
            args.output_dir, genome_dir, args.threads,
            reextract_only=args.reextract_only,
        )
        total_recovered += n
        results.append((prefix, n))

        # Delete genome to save disk space (unless --keep-genomes)
        if not args.keep_genomes and prefix in SPECIES_MAP:
            acc = SPECIES_MAP[prefix][2]
            genome_fa = os.path.join(genome_dir, f'{acc}.fa')
            if os.path.exists(genome_fa):
                os.remove(genome_fa)
                print(f'  Cleaned up genome file')

    # Summary
    print(f'\n{"="*60}')
    print(f'SUMMARY')
    print(f'{"="*60}')
    for prefix, n in results:
        species = SPECIES_MAP.get(prefix, ("", "unknown", ""))[1]
        print(f'  {prefix:<10} {species:<35} {n:>6} CDS')
    print(f'  {"TOTAL":<10} {"":35} {total_recovered:>6} CDS')

    # Merge all species CDS into one file
    merged = os.path.join(args.output_dir, 'all_recovered_cds.fna')
    with open(merged, 'w') as out:
        for prefix, n in results:
            species_cds = os.path.join(args.output_dir, f'{prefix}_cds.fna')
            if os.path.exists(species_cds):
                with open(species_cds) as f:
                    out.write(f.read())
    merged_count = sum(1 for line in open(merged) if line.startswith('>'))
    print(f'\n  Merged file: {merged} ({merged_count} sequences)')

    # Write report
    report = os.path.join(args.output_dir, 'recovery_report.tsv')
    with open(report, 'w') as f:
        f.write('species_prefix\tspecies_name\tassembly\tcds_recovered\n')
        for prefix, n in results:
            info = SPECIES_MAP.get(prefix, ("", "unknown", "unknown"))
            f.write(f'{prefix}\t{info[1]}\t{info[2]}\t{n}\n')
    print(f'  Report: {report}')

    # Per-sequence provenance manifest (bead -325 follow-up). Rank_candidates
    # / parse_absrel can join on this to flag miniprot-recovered branches
    # in dN/dS results (literature review flagged frameshift / splice-site
    # bias risk for back-translated CDS).
    if args.manifest:
        os.makedirs(os.path.dirname(args.manifest) or '.', exist_ok=True)
        # Append-mode if file exists with the right header (multi-run safe);
        # otherwise create.
        write_header = not (os.path.exists(args.manifest)
                            and open(args.manifest).readline().startswith('seq_id'))
        with open(args.manifest, 'a') as out:
            if write_header:
                out.write('seq_id,source,species_prefix,assembly_accession,recovery_method\n')
            for prefix, n in results:
                if n <= 0:
                    continue
                info = SPECIES_MAP.get(prefix, ("", "", ""))
                species_cds = os.path.join(args.output_dir, f'{prefix}_cds.fna')
                if not os.path.exists(species_cds):
                    continue
                for line in open(species_cds):
                    if line.startswith('>'):
                        sid = line[1:].split()[0]
                        out.write(f'{sid},miniprot,{prefix},{info[2]},miniprot+ncbi-datasets\n')
        print(f'  Manifest: {args.manifest}')


if __name__ == '__main__':
    # main() returns 1 when --validate-species-map fails its audit. Calling it
    # bare discarded that, so the audit reported a drifted map and still exited
    # 0 -- a guard that cannot fail the caller is not a guard.
    sys.exit(main())
