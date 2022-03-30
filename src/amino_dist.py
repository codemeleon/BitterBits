from Bio import SeqIO, Seq
# https://www.biostars.org/p/9480935/


def amino_dist(fasta_file, codon_table):
    """
    Calculates the amino acid distribution of a given fasta file.
    :param fasta_file: The fasta file to be analyzed.
    :param codon_table: The codon table to be used.
    :return: A dictionary with the amino acid distribution.
    """
    # Initialize the dictionary
    amino_dist = {}
    # Iterate over the fasta file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Get the sequence
        seq = str(record.seq)
        # Iterate over the sequence
        for i in range(0, len(seq), 3):
            # Get the codon
            codon = seq[i:i+3]
            # Get the amino acid
            amino_acid = codon_table[codon]
            # Check if the amino acid is in the dictionary
            if amino_acid in amino_dist:
                # If it is, increment the count
                amino_dist[amino_acid] += 1
            else:
                # If it is not, add it to the dictionary
                amino_dist[amino_acid] = 1
    # Return the dictionary
    return amino_dist


codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

amino_acid
