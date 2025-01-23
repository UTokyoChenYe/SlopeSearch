from Bio import SeqIO
from typing import List

def load_sequences(file_path: str) -> List[str]:
    """
    Loading DNA sequences from a FASTA file

    Input: fasta file path
    File example:
        >seq1
        ACGTACGTACGT
        >seq2
        TGCATGCATGCA
        >seq3
        GGGAAAACCCGGG

    Output: a list of DNA sequences
    Output example: ['ACGTACGTACGT', 'TGCATGCATGCA', 'GGGAAAACCCGGG']
    """
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required")
    return sequences