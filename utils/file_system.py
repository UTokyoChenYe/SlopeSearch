from Bio import SeqIO
from typing import List, Tuple

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
    return 

def load_sequences_for_evaluation(file_path: str) -> Tuple[List[str], List[str]]:
    records = list(SeqIO.parse(file_path, "fasta"))
    names = [rec.id for rec in records]
    seqs = [str(rec.seq).upper() for rec in records]
    return names, seqs