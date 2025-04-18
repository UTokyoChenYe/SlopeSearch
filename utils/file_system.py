from Bio import SeqIO
from typing import List, Tuple
import os
import glob

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

def load_sequences_for_evaluation(file_path: str) -> Tuple[List[str], List[str]]:
    records = list(SeqIO.parse(file_path, "fasta"))
    names = [rec.id for rec in records]
    seqs = [str(rec.seq).upper() for rec in records]
    return names, seqs

def load_sequences_for_evaluation_from_multiple_files(directory: str) -> Tuple[List[str], List[str]]:
    """
    Load sequences from all FASTA files in a directory and its subdirectories.
    
    :param directory: The root directory to search for FASTA files.
    :return: A tuple containing a list of sequence names and a list of sequences.
    """
    names = []
    seqs = []

    # Use glob to find all FASTA files in the directory and subdirectories
    fasta_files = glob.glob(os.path.join(directory, '**', '*.fasta'), recursive=True)

    for file_path in fasta_files:
        # Parse each FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            names.append(record.id)
            seqs.append(str(record.seq).upper())

    return names, seqs