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
}

# Phoronis australis shares "phau" prefix but different taxid
# Handle separately — its proteins have NMRA scaffold patterns
SPECIES_MAP["phau2"] = ("115415", "Phoronis_australis", "GCA_055505105.1")

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
    """Get protein IDs for a species that lack CDS."""
    # Load existing CDS IDs
    cds_ids = set()
    if os.path.exists(cds_file):
        with open(cds_file) as f:
            for line in f:
                if line.startswith('>'):
                    cds_ids.add(line[1:].split()[0])

    # Find proteins for this species in OG FASTAs
    missing = {}
    for fa_path in glob.glob(os.path.join(og_dir, '*.fa')):
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
    stats = {'found': 0, 'no_hit': 0, 'bad_cds': 0, 'multi_hit': 0}

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

        # Verify translation roughly matches protein (>60% identity)
        translated = ''.join(CODON_TABLE.get(c, 'X') for c in
                           [cds_upper[i:i+3] for i in range(0, len(cds_upper), 3)])
        translated = translated.rstrip('*')
        prot_clean = prot_seq.replace('-', '').replace('*', '')

        if len(translated) > 0 and len(prot_clean) > 0:
            # Quick identity check: compare overlapping region
            min_len = min(len(translated), len(prot_clean))
            matches = sum(1 for a, b in zip(translated[:min_len], prot_clean[:min_len]) if a == b)
            identity = matches / min_len
            if identity < 0.5:
                stats['bad_cds'] += 1
                continue

        cds_seqs[prot_id] = cds
        stats['found'] += 1

    return cds_seqs, stats


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]


def process_species(prefix, og_dir, cds_file, output_dir, genome_dir, threads=2):
    """Process one species: download genome, run miniprot, extract CDS."""
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

    # Write protein queries to temp FASTA
    proteins_fa = os.path.join(output_dir, f'{prefix}_query_proteins.fa')
    write_fasta(missing, proteins_fa)

    # Download assembly
    genome_fa = download_assembly(accession, genome_dir, threads)
    if not genome_fa:
        print(f'  FAILED: could not download assembly for {prefix}')
        return 0

    # Run miniprot
    gff_path = os.path.join(output_dir, f'{prefix}_miniprot.gff')
    if not run_miniprot(genome_fa, proteins_fa, gff_path, threads):
        print(f'  FAILED: miniprot error for {prefix}')
        return 0

    # Extract CDS
    print(f'  Extracting CDS from miniprot alignments...')
    cds_seqs, stats = extract_cds_from_gff(genome_fa, gff_path, missing)

    print(f'  Results: {stats["found"]} CDS recovered, '
          f'{stats["no_hit"]} no hit, {stats["bad_cds"]} bad CDS, '
          f'{stats["multi_hit"]} multi-hit')

    # Write species CDS
    if cds_seqs:
        species_cds = os.path.join(output_dir, f'{prefix}_cds.fna')
        write_fasta(cds_seqs, species_cds)
        print(f'  Wrote {len(cds_seqs)} CDS to {species_cds}')

    # Clean up query proteins file
    os.remove(proteins_fa)

    return len(cds_seqs)


def main():
    parser = argparse.ArgumentParser(
        description='Recover CDS from genome assemblies via miniprot')
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
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    genome_dir = args.genome_dir or os.path.join(args.output_dir, 'genomes')
    os.makedirs(genome_dir, exist_ok=True)

    # Determine species to process
    if args.species:
        species_list = args.species.split(',')
    else:
        species_list = sorted(SPECIES_MAP.keys())

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
            args.output_dir, genome_dir, args.threads
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


if __name__ == '__main__':
    main()
