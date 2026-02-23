#!/usr/bin/env python3
"""Generate minimal synthetic test data for pipeline dry-run validation."""
import os
import sys
import random
import string


def random_protein_seq(length=200):
    """Generate random protein sequence."""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(aa) for _ in range(length))


def random_dna_seq(length=600):
    """Generate random DNA sequence (length must be multiple of 3)."""
    codons = ["ATG", "TAC", "GCT", "AAA", "TTT", "GGG", "CCC",
              "GAT", "CAG", "AGC", "TGC", "ACG", "CTG", "GCA"]
    seq = ""
    while len(seq) < length:
        seq += random.choice(codons)
    return seq[:length]


def write_fasta(path, entries):
    """Write FASTA file. entries = list of (id, seq) tuples."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for seq_id, seq in entries:
            f.write(f">{seq_id}\n{seq}\n")


def write_csv(path, header, rows):
    """Write CSV file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")


def generate_all(base_dir):
    """Generate all synthetic test data."""
    random.seed(42)  # Reproducible test data

    results = os.path.join(base_dir, "results")
    transcriptomes = os.path.join(base_dir, "transcriptomes")
    genomes = os.path.join(base_dir, "genomes")
    references = os.path.join(base_dir, "references")
    ref_results = os.path.join(results, "reference_sequences")

    # Taxa for testing
    taxa = ["taxid_berghia", "taxid1", "taxid2"]
    n_seqs = 5  # sequences per taxon

    # --- Transcriptomes (protein .aa files) ---
    for taxid in taxa:
        entries = [(f"{taxid}_prot_{i}", random_protein_seq(200 + i * 10))
                   for i in range(1, n_seqs + 1)]
        write_fasta(os.path.join(transcriptomes, f"{taxid}.aa"), entries)
        # Also write nucleotide versions
        nuc_entries = [(f"{taxid}_prot_{i}", random_dna_seq(600 + i * 30))
                       for i in range(1, n_seqs + 1)]
        write_fasta(os.path.join(transcriptomes, f"{taxid}.mrna"), nuc_entries)

    # Berghia-specific transcriptome (config.sh expects taxid_berghia_berghia.aa)
    berghia_entries = [(f"taxid_berghia_prot_{i}", random_protein_seq(250))
                       for i in range(1, n_seqs + 1)]
    write_fasta(os.path.join(transcriptomes, "taxid_berghia_berghia.aa"), berghia_entries)

    # --- Genomes (2 needed for synteny) ---
    for taxid in ["taxid_berghia", "taxid1"]:
        # Genome FASTA
        scaffolds = [(f"scaffold_{i}", random_dna_seq(3000)) for i in range(1, 4)]
        genome_name = f"{taxid}_berghia.fasta" if taxid == "taxid_berghia" else f"{taxid}.fasta"
        write_fasta(os.path.join(genomes, genome_name), scaffolds)
        # Predicted proteins
        prots = [(f"{taxid}_gene_{i}", random_protein_seq(150)) for i in range(1, 4)]
        write_fasta(os.path.join(genomes, f"{taxid}.proteins.fa"), prots)
        # GFF annotation
        gff_path = os.path.join(genomes, f"{taxid}.gff")
        os.makedirs(os.path.dirname(gff_path), exist_ok=True)
        with open(gff_path, "w") as f:
            f.write("##gff-version 3\n")
            for i in range(1, 4):
                f.write(f"scaffold_{i}\tpred\tgene\t{i * 1000}\t{i * 1000 + 500}\t.\t+\t.\t"
                        f"ID={taxid}_gene_{i}\n")

    # --- Reference sequences (Nath et al. structure) ---
    # Step 01 expects nath_et_al/{one_to_one_ortholog,lse}/{group}/Species.faa
    ref_dir = os.path.join(references, "nath_et_al")
    ref_entries = []
    for i in range(1, 8):
        seq_id = f"ref_species{i}_prot1"
        ref_entries.append((seq_id, random_protein_seq(300)))
    # Create subdirectory structure matching step 01's find commands
    conserved_dir = os.path.join(ref_dir, "one_to_one_ortholog", "TestGroup")
    lse_dir = os.path.join(ref_dir, "lse", "TestGroup")
    os.makedirs(conserved_dir, exist_ok=True)
    os.makedirs(lse_dir, exist_ok=True)
    write_fasta(os.path.join(conserved_dir, "10001_Conserved_Species1.faa"), ref_entries[:2])
    write_fasta(os.path.join(conserved_dir, "10002_Conserved_Species2.faa"), ref_entries[2:4])
    write_fasta(os.path.join(lse_dir, "10003_LSE_Species3.faa"), ref_entries[4:6])
    write_fasta(os.path.join(lse_dir, "10004_LSE_Species4.faa"), ref_entries[6:])

    # Also write CDS references
    ref_cds = [(eid, random_dna_seq(900)) for eid, _ in ref_entries]
    cds_conserved = os.path.join(ref_dir, "one_to_one_ortholog", "TestGroup")
    cds_lse = os.path.join(ref_dir, "lse", "TestGroup")
    write_fasta(os.path.join(cds_conserved, "10001_Conserved_Species1_cds.fna"), ref_cds[:2])
    write_fasta(os.path.join(cds_conserved, "10002_Conserved_Species2_cds.fna"), ref_cds[2:4])
    write_fasta(os.path.join(cds_lse, "10003_LSE_Species3_cds.fna"), ref_cds[4:6])
    write_fasta(os.path.join(cds_lse, "10004_LSE_Species4_cds.fna"), ref_cds[6:])

    # --- Pre-populated results that steps expect ---
    # Step 01 outputs (so step 02 can run)
    write_fasta(os.path.join(ref_results, "all_references.fa"), ref_entries)
    write_fasta(os.path.join(ref_results, "conserved_references.fa"), ref_entries[:4])
    write_fasta(os.path.join(ref_results, "lse_references.fa"), ref_entries[4:])
    os.makedirs(os.path.join(ref_results, "cds"), exist_ok=True)
    write_fasta(os.path.join(ref_results, "cds", "all_references_cds.fna"), ref_cds)

    # ID map CSV
    id_map_rows = []
    for i, (seq_id, seq) in enumerate(ref_entries, 1):
        id_map_rows.append([seq_id, f"ref_{i}", f"species{i}", "reference",
                            f"Species {i}", "", str(len(seq))])
    write_csv(os.path.join(ref_results, "id_map.csv"),
              "original_id,short_id,taxid,source_type,species_name,gene_name,seq_length",
              id_map_rows)

    # HMM directory (empty files - stubs will handle)
    hmm_dir = os.path.join(results, "hmms")
    os.makedirs(hmm_dir, exist_ok=True)
    for hmm in ["conserved.hmm", "lse.hmm"]:
        open(os.path.join(hmm_dir, hmm), "w").close()

    # Expression data CSV (simple format)
    expr_path = os.path.join(base_dir, "expression_data.csv")
    with open(expr_path, "w") as f:
        f.write("gene_id,rhinophore,oral_veil,foot,digestive\n")
        for i in range(1, n_seqs + 1):
            vals = [random.uniform(0, 100) for _ in range(4)]
            f.write(f"taxid_berghia_prot_{i},{','.join(f'{v:.1f}' for v in vals)}\n")

    # Salmon quant directory structure (for expression analysis)
    for tissue in ["rhinophore", "oral_veil", "foot"]:
        quant_dir = os.path.join(base_dir, "expression_data", tissue)
        os.makedirs(quant_dir, exist_ok=True)
        with open(os.path.join(quant_dir, "quant.sf"), "w") as f:
            f.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
            for i in range(1, n_seqs + 1):
                tpm = random.uniform(0, 200)
                reads = int(tpm * 10)
                f.write(f"taxid_berghia_prot_{i}\t600\t550\t{tpm:.2f}\t{reads}\n")

    # --- Create empty placeholder dirs ---
    for d in ["chemogpcrs", "orthogroups", "orthogroups/input",
              "phylogenies/protein", "phylogenies/visualizations",
              "selective_pressure", "selective_pressure/nucleotide",
              "asr", "synteny", "synteny/blast", "synteny/gff",
              "ranking", "structural_analysis", "busco",
              "busco/single_copy", "busco/alignments", "busco/gene_trees",
              "lse_classification", "clustering", "classification",
              "candidates", "mapping", "cafe", "notung",
              "report", "logs", "checkpoints", "provenance",
              "hhdb", "gproteins", "ecl_analysis"]:
        os.makedirs(os.path.join(results, d), exist_ok=True)

    # Create custom_hmms directory (referenced by config)
    custom_hmms = os.path.join(base_dir, "custom_hmms")
    os.makedirs(custom_hmms, exist_ok=True)

    print(f"Generated synthetic test data in {base_dir}", file=sys.stderr)
    print(f"  Taxa: {taxa}", file=sys.stderr)
    print(f"  Sequences per taxon: {n_seqs}", file=sys.stderr)
    print(f"  Genomes: 2 (for synteny)", file=sys.stderr)
    return True


if __name__ == "__main__":
    base = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    generate_all(base)
