#!/usr/bin/env python3
"""
fetch_cds_from_scaffolds.py
Extract CDS for Nath et al. coordinate-based gene models.

Strategy:
1. Parse scaffold accession from gene IDs (e.g., alvmar_JABMCL010000060.1_7503 -> JABMCL010000060.1)
2. Fetch scaffold sequences from NCBI in batches
3. Run tBLASTn (protein vs scaffold) to find CDS exons
4. Extract and concatenate CDS

Author: Jorge L. Perez-Moreno, Ph.D.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Entrez.email = "jperezmoreno@umass.edu"
Entrez.tool = "berghia_chemogpcrs_pipeline"


def parse_scaffold_from_id(gene_id):
    """
    Extract scaffold accession from a Nath et al. gene ID.

    Formats:
    - alvmar_JABMCL010000060.1_7503 -> JABMCL010000060.1
    - chsq_chr10_391081 -> chr10 (need assembly accession)
    - hadi_HDSC00048_47552 -> HDSC00048 (Haliotis scaffolds)
    - lyst_gLs.1.0.scaf02289_1650 -> gLs.1.0.scaf02289 (custom scaffold)
    - pavu_2_3274693 -> chromosome 2 (need assembly)
    - gima_11_234097 -> chromosome 11 (need assembly)
    - hadi_HDSC00071CG00110 -> HDSC00071 (Haliotis gene model)
    - lyst_g43305.t1 -> gene model ID (no scaffold)
    """
    sp = gene_id.split('_')[0]
    rest = '_'.join(gene_id.split('_')[1:])

    # NCBI-style accessions: NW_NNN.N, NC_NNN.N, CMNNNN.N, OXNNNN.N, etc.
    # NW_/NC_ have underscore between prefix and digits
    m = re.match(r'((?:NW_|NC_)\d+\.\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'ncbi_scaffold', sp

    # CM/OX/OU/OV style (no underscore between prefix and digits)
    m = re.match(r'((?:CM|OX|OU|OV)\d+\.\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'ncbi_scaffold', sp

    # WGS accessions: 4+ uppercase letters followed by digits
    m = re.match(r'([A-Z]{4,}\d+\.\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'wgs_scaffold', sp

    # Haliotis scaffolds: HDSC followed by digits
    m = re.match(r'(HDSC\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'hadi_scaffold', sp

    # Haliotis gene models: HDSCnnnnnCGnnnnn
    m = re.match(r'(HDSC\d+)CG\d+$', rest)
    if m:
        return m.group(1), 'hadi_gene_model', sp

    # Chrysomallon v2 assembly: Csqv2.0_NNN_NN.NN or Csqv2.0S_NNN_NN.NN
    m = re.match(r'(Csqv2\.0S?_\d+)_[\d.]+$', rest)
    if m:
        return m.group(1), 'chsq_scaffold', sp

    # Short scaffold accessions (2-3 uppercase + digits): KB, GF, etc. (Lottia, Elysia)
    m = re.match(r'([A-Z]{2,3}\d+\.\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'ncbi_scaffold', sp

    # BRAKER gene models (no scaffold info)
    m = re.match(r'BRAKERGMAP\d+\.\d+$', rest)
    if m:
        return None, 'braker_gene_model', sp

    # Chromosome-based: chrN_coord or N_coord (Gibbula, Patella)
    m = re.match(r'(chr\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'chromosome', sp

    m = re.match(r'(\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'chromosome_num', sp

    # Lymnaea scaffolds: gLs.1.0.scafNNNNN_coord
    m = re.match(r'(gLs\.\d+\.\d+\.scaf\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'lyst_scaffold', sp

    # Lymnaea gene models: gNNNNN.t1
    m = re.match(r'(g\d+\.t\d+)$', rest)
    if m:
        return None, 'gene_model', sp

    # 3-letter ENA protein accessions with version
    m = re.match(r'([A-Z]{3}\d+\.\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'ena_scaffold', sp

    # VKKT-type (Haliotis laevigata)
    m = re.match(r'(VKKT\d+\.\d+)_(\d+)$', rest)
    if m:
        return m.group(1), 'wgs_scaffold', sp

    return None, 'unknown', sp


# Map species codes to NCBI assembly accessions (for chromosome-based IDs)
SPECIES_ASSEMBLIES = {
    'chsq': 'GCA_019457155.2',   # Chrysomallon squamiferum (chromosomal)
    'pavu': 'GCA_917563875.2',   # Patella vulgata (chromosomal)
    'gima': 'GCA_963853765.1',   # Gibbula magus
    'lyst': 'GCA_944038965.1',   # Lymnaea stagnalis
    'hadi': 'GCA_011762535.2',   # Haliotis discus hannai
}


def fetch_scaffold_region(scaffold_acc, start=None, end=None, retries=3):
    """Fetch a scaffold sequence or region from NCBI."""
    for attempt in range(retries):
        try:
            params = {'db': 'nucleotide', 'id': scaffold_acc, 'rettype': 'fasta', 'retmode': 'text'}
            if start and end:
                params['seq_start'] = str(min(start, end))
                params['seq_stop'] = str(max(start, end))
            handle = Entrez.efetch(**params)
            records = list(SeqIO.parse(handle, 'fasta'))
            handle.close()
            if records:
                return records[0]
            return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"    Failed to fetch {scaffold_acc}: {e}", file=sys.stderr)
                return None


def batch_fetch_scaffolds(scaffold_ids, batch_size=50):
    """Fetch multiple scaffold sequences from NCBI in batches."""
    results = {}
    scaffold_list = list(set(scaffold_ids))

    for i in range(0, len(scaffold_list), batch_size):
        batch = scaffold_list[i:i + batch_size]
        ids_str = ','.join(batch)
        try:
            handle = Entrez.efetch(db='nucleotide', id=ids_str,
                                   rettype='fasta', retmode='text')
            for rec in SeqIO.parse(handle, 'fasta'):
                # Match back to original accession
                for sid in batch:
                    if sid in rec.id or rec.id.startswith(sid.split('.')[0]):
                        results[sid] = rec
                        break
            handle.close()
            time.sleep(0.5)  # Rate limit
        except Exception as e:
            print(f"  Batch fetch failed ({len(batch)} IDs): {e}", file=sys.stderr)
            # Fall back to individual fetch
            for sid in batch:
                rec = fetch_scaffold_region(sid)
                if rec:
                    results[sid] = rec
                time.sleep(0.4)

    return results


def run_tblastn(protein_seq, scaffold_seq, protein_id):
    """Run tBLASTn to find CDS on scaffold. Returns list of (start, end, strand, sseq) tuples."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as pf:
        pf.write(f">{protein_id}\n{protein_seq}\n")
        prot_file = pf.name

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as sf:
        sf.write(f">scaffold\n{scaffold_seq}\n")
        scaffold_file = sf.name

    # Create BLAST database
    db_file = scaffold_file
    subprocess.run(['makeblastdb', '-in', scaffold_file, '-dbtype', 'nucl',
                    '-out', db_file], capture_output=True)

    # Run tBLASTn with tabular output
    # Format: qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq
    result = subprocess.run(
        ['tblastn', '-query', prot_file, '-db', db_file,
         '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq',
         '-evalue', '1e-5', '-max_target_seqs', '1', '-num_threads', '1'],
        capture_output=True, text=True
    )

    # Clean up temp files
    for f in [prot_file, scaffold_file]:
        try:
            os.unlink(f)
        except OSError:
            pass
    for ext in ['.nhr', '.nin', '.nsq', '.ndb', '.not', '.ntf', '.nto']:
        try:
            os.unlink(db_file + ext)
        except OSError:
            pass

    if result.returncode != 0:
        return []

    hits = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            continue
        pident = float(fields[2])
        qstart, qend = int(fields[4]), int(fields[5])
        sstart, send = int(fields[6]), int(fields[7])
        evalue = float(fields[8])
        sseq = fields[10].replace('-', '')  # Remove alignment gaps

        strand = '+' if sstart < send else '-'
        hits.append({
            'pident': pident,
            'qstart': qstart, 'qend': qend,
            'sstart': min(sstart, send), 'send': max(sstart, send),
            'strand': strand,
            'evalue': evalue,
            'sseq': sseq
        })

    return hits


def extract_cds_from_hits(hits, scaffold_seq, protein_len):
    """
    Extract CDS from tBLASTn hits.
    For single-exon genes, use the full hit region.
    For multi-exon, concatenate HSPs sorted by query position.
    """
    if not hits:
        return None

    # Sort hits by query start position (to handle multi-exon)
    hits.sort(key=lambda h: h['qstart'])

    # Check if single hit covers most of the protein
    best_hit = max(hits, key=lambda h: h['qend'] - h['qstart'])
    coverage = (best_hit['qend'] - best_hit['qstart'] + 1) / protein_len

    if coverage >= 0.8:
        # Single-exon or near-complete hit
        start = best_hit['sstart'] - 1  # Convert to 0-based
        end = best_hit['send']
        if best_hit['strand'] == '+':
            cds = str(scaffold_seq[start:end])
        else:
            cds = str(Seq(str(scaffold_seq[start:end])).reverse_complement())
        return cds

    # Multi-exon: concatenate HSPs
    # Filter overlapping HSPs (keep best per query region)
    filtered = []
    for hit in hits:
        overlap = False
        for existing in filtered:
            if hit['qstart'] <= existing['qend'] and hit['qend'] >= existing['qstart']:
                # Overlap - keep better one
                if hit['pident'] > existing['pident']:
                    filtered.remove(existing)
                else:
                    overlap = True
                break
        if not overlap:
            filtered.append(hit)

    filtered.sort(key=lambda h: h['qstart'])

    # Extract and concatenate
    cds_parts = []
    strand = filtered[0]['strand']

    for hit in filtered:
        start = hit['sstart'] - 1
        end = hit['send']
        if hit['strand'] == '+':
            part = str(scaffold_seq[start:end])
        else:
            part = str(Seq(str(scaffold_seq[start:end])).reverse_complement())
        cds_parts.append(part)

    cds = ''.join(cds_parts)

    # Verify: translate and check length
    try:
        translated = str(Seq(cds).translate())
        if abs(len(translated) - protein_len) <= 5:
            return cds
    except Exception:
        pass

    # Fallback: return concatenated anyway
    return cds if len(cds) >= protein_len * 2.5 else None


def main():
    parser = argparse.ArgumentParser(description='Fetch CDS for coordinate-based gene models')
    parser.add_argument('--missing-ids', required=True, help='File with missing gene IDs (one per line)')
    parser.add_argument('--protein-fasta', required=True, help='Protein FASTA with sequences for missing IDs')
    parser.add_argument('--output', required=True, help='Output CDS FASTA file')
    parser.add_argument('--report', help='Output report TSV')
    parser.add_argument('--existing-cds', help='Existing CDS file to append to')
    parser.add_argument('--batch-size', type=int, default=20, help='NCBI batch size')
    args = parser.parse_args()

    # Load missing IDs
    with open(args.missing_ids) as f:
        missing_ids = [line.strip() for line in f if line.strip()]

    # Load protein sequences
    proteins = {}
    for rec in SeqIO.parse(args.protein_fasta, 'fasta'):
        if rec.id in missing_ids:
            proteins[rec.id] = str(rec.seq)

    print(f"Missing IDs: {len(missing_ids)}", file=sys.stderr)
    print(f"Protein sequences loaded: {len(proteins)}", file=sys.stderr)

    # Parse scaffold accessions
    scaffold_to_genes = defaultdict(list)  # scaffold_acc -> [(gene_id, sp)]
    skipped = []

    for gid in missing_ids:
        if gid not in proteins:
            skipped.append((gid, 'no_protein'))
            continue
        scaffold, id_type, sp = parse_scaffold_from_id(gid)
        if scaffold is None:
            skipped.append((gid, f'unparseable_{id_type}'))
            continue
        scaffold_to_genes[scaffold].append((gid, sp, id_type))

    print(f"Unique scaffolds to fetch: {len(scaffold_to_genes)}", file=sys.stderr)
    print(f"Skipped (unparseable): {len(skipped)}", file=sys.stderr)

    # Fetch scaffolds and extract CDS
    results = []
    report_rows = []
    total_fetched = 0
    total_cds = 0

    scaffolds_list = list(scaffold_to_genes.keys())

    for i, scaffold_acc in enumerate(scaffolds_list):
        genes = scaffold_to_genes[scaffold_acc]
        gene_ids = [g[0] for g in genes]

        if (i + 1) % 10 == 0 or i == 0:
            print(f"  [{i + 1}/{len(scaffolds_list)}] Fetching {scaffold_acc} "
                  f"({len(genes)} genes)...", file=sys.stderr)

        # Fetch scaffold
        scaffold_rec = fetch_scaffold_region(scaffold_acc)
        time.sleep(0.4)  # Rate limit

        if scaffold_rec is None:
            for gid, sp, id_type in genes:
                report_rows.append({'gene_id': gid, 'scaffold': scaffold_acc,
                                    'status': 'fetch_failed', 'cds_length': 0})
            continue

        total_fetched += 1
        scaffold_seq = str(scaffold_rec.seq)

        # For each gene on this scaffold, run tBLASTn
        for gid, sp, id_type in genes:
            protein_seq = proteins.get(gid, '')
            if not protein_seq:
                continue

            hits = run_tblastn(protein_seq, scaffold_seq, gid)

            if not hits:
                report_rows.append({'gene_id': gid, 'scaffold': scaffold_acc,
                                    'status': 'no_blast_hit', 'cds_length': 0})
                continue

            cds = extract_cds_from_hits(hits, scaffold_seq, len(protein_seq))

            if cds:
                results.append(SeqRecord(Seq(cds), id=gid, description=''))
                total_cds += 1
                report_rows.append({'gene_id': gid, 'scaffold': scaffold_acc,
                                    'status': 'success', 'cds_length': len(cds)})
            else:
                report_rows.append({'gene_id': gid, 'scaffold': scaffold_acc,
                                    'status': 'extraction_failed', 'cds_length': 0})

    # Write output CDS FASTA
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.existing_cds and os.path.exists(args.existing_cds):
        # Append to existing
        existing = list(SeqIO.parse(args.existing_cds, 'fasta'))
        existing_ids = {r.id for r in existing}
        new_records = [r for r in results if r.id not in existing_ids]
        all_records = existing + new_records
        SeqIO.write(all_records, args.output, 'fasta')
        print(f"Appended {len(new_records)} new CDS to {len(existing)} existing", file=sys.stderr)
    else:
        SeqIO.write(results, args.output, 'fasta')

    # Write report
    if args.report:
        with open(args.report, 'w') as rf:
            writer = csv.DictWriter(rf, fieldnames=['gene_id', 'scaffold', 'status', 'cds_length'],
                                    delimiter='\t')
            writer.writeheader()
            for row in report_rows:
                writer.writerow(row)
            for gid, reason in skipped:
                writer.writerow({'gene_id': gid, 'scaffold': '', 'status': reason, 'cds_length': 0})

    print(f"\n=== CDS Extraction Summary ===", file=sys.stderr)
    print(f"  Scaffolds fetched: {total_fetched}/{len(scaffolds_list)}", file=sys.stderr)
    print(f"  CDS extracted: {total_cds}/{len(missing_ids)}", file=sys.stderr)
    print(f"  Skipped: {len(skipped)}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
