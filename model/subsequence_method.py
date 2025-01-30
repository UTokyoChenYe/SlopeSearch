import random
from collections import Counter
from typing import List

# 0. tools
def reverse_complement(sequence: str) -> str:
    """generate the reverse of seq"""
    complement = str.maketrans("ACGT", "TGCA") 
    return sequence.translate(complement)[::-1]  

def count_kmers(sequences: List[str], k: int) -> Counter:
    """calculate k-mer frequencies in a given sequence"""
    kmers = [sequence[i:i+k] for sequence in sequences for i in range(len(sequence) - k + 1)]
    # kmers = s[i:i+k] for s in sequences for i in range(len(s) - k + 1)
    return Counter(kmers)

def generate_pattern(length):
    """Generate a random binary pattern of given length."""
    return ''.join(random.choice(['0', '1']) for _ in range(length))

def extract_spaced_word(sequence: str, pattern: str) -> List[str]:
    words = []
    pattern_length = len(pattern)
    for i in range(len(sequence) - pattern_length + 1):
        spaced_word = ''.join(sequence[i+j] for j in range(pattern_length) if pattern[j] == '1')
        words.append(spaced_word)
    return words

# 1. basic k-mer matches
def basic_kmer_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse = reverse_complement(seq1)
        seq1 = [seq1, seq1_reverse]
        seq2 = [seq2] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = count_kmers(seq1, k)
    kmer_count2 = count_kmers(seq2, k)
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += kmer_count1[kmer] * kmer_count2[kmer]
    return matches

# 2. spaced-word matches
def spaced_word_matches(seq1: str, seq2: str, k: int) -> int:
    """calculate the number of spaced-word matches between two sequences"""

    pattern = generate_pattern(k)

    words1 = extract_spaced_word(seq1, pattern)
    words2 = extract_spaced_word(seq2, pattern)

    word_count1 = Counter(words1)
    word_count2 = Counter(words2)
    matches = 0
    for word in word_count1:
        if word in word_count2:
            matches += word_count1[word] * word_count2[word]
    return matches