import io
import csv
from Bio.Blast import NCBIWWW

# DNA to RNA transcription
def transcribe_dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

# RNA to protein translation
def translate_rna_to_protein(rna_sequence):
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    protein_sequence = ''
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) == 3:  # Ensure codon is of length 3
            protein_sequence += codon_table.get(codon, '*')
    return protein_sequence

# GC content calculation
def calculate_gc_content(dna_sequence):
    cleaned_sequence = ''.join([n for n in dna_sequence.upper() if n in 'ATGC'])
    if len(cleaned_sequence) == 0:
        return 0.0
    gc_count = cleaned_sequence.count('G') + cleaned_sequence.count('C')
    return (gc_count / len(cleaned_sequence)) * 100

# Nucleotide distribution
def get_nucleotide_distribution(dna_sequence):
    nucleotides = ['A', 'T', 'C', 'G']
    distribution = {nucleotide: dna_sequence.upper().count(nucleotide) for nucleotide in nucleotides}
    return distribution

# Generate FASTA format
def generate_fasta(dna_sequence, sequence_id="sequence"):
    return f">{sequence_id}\n{dna_sequence}"

# Generate CSV
def generate_csv(dna_sequence, rna_sequence, protein_sequence):
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["DNA Sequence", "RNA Sequence", "Protein Sequence"])
    writer.writerow([dna_sequence, rna_sequence, protein_sequence])
    return output.getvalue()

# Query NCBI BLAST
def query_ncbi_blast(sequence):
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        return result_handle.read()
    except Exception as e:
        return f"An error occurred while querying BLAST: {e}"
