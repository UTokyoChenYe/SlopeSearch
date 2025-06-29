import random
from collections import Counter
from typing import List

from utils.sequence_tool import *

# 0. tools
def count_kmers(sequences: List[str], k: int) -> Counter:
    """calculate k-mer frequencies in a given sequence"""
    kmers = [sequence[i:i+k] for sequence in sequences for i in range(len(sequence) - k + 1)]
    # kmers = s[i:i+k] for s in sequences for i in range(len(s) - k + 1)
    return Counter(kmers)

def count_kmers_start_ry(sequences: List[str], k: int) -> Counter:
    """calculate k-mer frequencies in a given sequence"""
    purines = {'A', 'G'} 
    pyrimidines = {'C', 'T'}
    kmers = []
    for sequence in sequences:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]  
            if kmer[0] in purines and kmer[1] in pyrimidines:
                kmers.append(kmer)
    return Counter(kmers)

def count_kmers_start_rr(sequences: List[str], k: int) -> Counter:
    """计算以两个嘌呤（A或G）开头的k-mer频率"""
    purines = {'A', 'G'}
    kmers = []
    for sequence in sequences:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            # 检查前两个碱基是否为嘌呤
            if kmer[0] in purines and kmer[1] in purines:
                kmers.append(kmer)
    return Counter(kmers)

# def count_kmers_start_ry_4_6(sequences: List[str], k: int) -> Counter:
#     """计算以特定 RY46 模式开头的 k-mer 频率"""
#     patterns = [
#         'RRRRRY', 'RRYRRY', 'RRYRYY', 'RYRRRR', 'RYRRRY',
#         'RYRYRY', 'RYYRRR', 'RYYRRY', 'RYYRYR', 'RYYRYY',
#         'RYYYRY', 'RYYYYY', 'YRYRRY', 'YYYRRR', 'YYYRRY',
#         'YYYYRY'
#     ]
#     purines = {'A', 'G'}
#     pyrimidines = {'C', 'T'}
#     kmers = []
#     for sequence in sequences:
#         for i in range(len(sequence) - k + 1):
#             kmer = sequence[i:i + k]
#             # Convert kmer to RY pattern
#             ry_pattern = ''.join(['R' if base in purines else 'Y' for base in kmer])
#             # Check if the k-mer starts with any of the patterns
#             if any(ry_pattern.startswith(pattern) for pattern in patterns):
#                 kmers.append(kmer)
#     return Counter(kmers)


def count_kmers_start_ry_4_6(sequences: List[str], k: int) -> Counter:
    """计算以特定 RY46 模式开头的 k-mer 频率"""
    patterns = [
        'RRRRRY', 'RRYRRY', 'RRYRYY', 'RYRRRR', 'RYRRRY',
        'RYRYRY', 'RYYRRR', 'RYYRRY', 'RYYRYR', 'RYYRYY',
        'RYYYRY', 'RYYYYY', 'YRYRRY', 'YYYRRR', 'YYYRRY',
        'YYYYRY'
    ]
    table = str.maketrans("ACGT", "RYRY")
    word_length = len(next(iter(patterns)))

    for sequence in sequences:
        ry_sequence = sequence.translate(table)
        for i in range(len(sequence) - k + 1):
            ry_pattern = ry_sequence[i:i + word_length]
            if ry_pattern in patterns:
                yield sequence[i:i + k]



# def count_kmers_start_ry_4_9(sequences: List[str], k: int) -> Counter:
#     """计算以特定 RY49 模式开头的 k-mer 频率"""
#     patterns = [
#         'rrrrrrrry', 'rrrrrryry', 'rrrrryrrr', 'rrrrryrry',
#         'rrrrryyrr', 'rrrrryyry', 'rrrrryyyr', 'rrrrryyyy',
#         'rrryryrrr', 'rrryryrry', 'rrryryyrr', 'rrryryyry',
#         'rrryryyyr', 'rrryryyyy', 'rryrrrrrr', 'rryrrrrry',
#         'rryrrryry', 'rryrryrrr', 'rryrryrry', 'rryrryryr',
#         'rryrryyrr', 'rryrryyry', 'rryrryyyy', 'rryryryry',
#         'rryyryrrr', 'rryyryrry', 'rryyryyry', 'rryyryyyy',
#         'rryyyyrrr', 'rryyyyrry', 'rryyyyryr', 'rryyyyyrr',
#         'ryrrrryrr', 'ryrrrryry', 'ryrrrryyr', 'ryrrrryyy',
#         'ryrrryrrr', 'ryrrryrry', 'ryrrryyrr', 'ryrrryyry',
#         'ryrrryyyr', 'ryrrryyyy', 'ryrryryyr', 'ryrryryyy',
#         'ryrryyyrr', 'ryrryyyry', 'ryryryrrr', 'ryryryrry',
#         'ryryryyrr', 'ryryryyry', 'ryryryyyr', 'ryryryyyy',
#         'ryyrrrrrr', 'ryyrrrrry', 'ryyrrrryr', 'ryyrrryrr',
#         'ryyrrryry', 'ryyrrryyr', 'ryyrrryyy', 'ryyrryrrr',
#         'ryyrryrry', 'ryyrryryr', 'ryyrryyrr', 'ryyrryyry',
#         'ryyrryyyr', 'ryyrryyyy', 'ryyryryrr', 'ryyryryry',
#         'ryyryryyr', 'ryyryryyy', 'ryyryyrrr', 'ryyryyrry',
#         'ryyryyyrr', 'ryyryyyry', 'ryyyryrrr', 'ryyyryrry',
#         'ryyyryyrr', 'ryyyryyry', 'ryyyryyyr', 'ryyyryyyy',
#         'ryyyyryyr', 'ryyyyryyy', 'ryyyyyryr', 'ryyyyyyrr',
#         'ryyyyyyyy', 'yryrrrrrr', 'yryrrrrry', 'yryrrryry',
#         'yryrryrrr', 'yryrryrry', 'yryrryryr', 'yryrryyrr',
#         'yryrryyry', 'yryrryyyy', 'yryyryrrr', 'yryyryrry',
#         'yryyryyry', 'yryyryyyy', 'yryyyyrrr', 'yryyyyrry',
#         'yryyyyryr', 'yryyyyyrr', 'yyrrrryyr', 'yyrrrryyy',
#         'yyrryryyr', 'yyrryryyy', 'yyyrrrrrr', 'yyyrrrrry',
#         'yyyrrrryr', 'yyyrrryrr', 'yyyrrryry', 'yyyrrryyr',
#         'yyyrrryyy', 'yyyrryrrr', 'yyyrryrry', 'yyyrryryr',
#         'yyyrryyrr', 'yyyrryyry', 'yyyrryyyr', 'yyyrryyyy',
#         'yyyryryrr', 'yyyryryry', 'yyyryryyr', 'yyyryryyy',
#         'yyyyyryyr', 'yyyyyryyy', 'yyyyyyryr', 'yyyyyyyrr'
#     ]
#     purines = {'A', 'G'}
#     pyrimidines = {'C', 'T'}
#     kmers = []
#     for sequence in sequences:
#         for i in range(len(sequence) - k + 1):
#             kmer = sequence[i:i + k]
#             # Convert kmer to RY pattern
#             ry_pattern = ''.join(['r' if base in purines else 'y' for base in kmer])
#             # Check if the k-mer starts with any of the patterns
#             if any(ry_pattern.startswith(pattern) for pattern in patterns):
#                 kmers.append(kmer)
#     return Counter(kmers)

def count_kmers_start_ry_4_9(sequences: List[str], k: int) -> Counter:
    """计算以特定 RY49 模式开头的 k-mer 频率"""
    patterns = [
        'rrrrrrrry', 'rrrrrryry', 'rrrrryrrr', 'rrrrryrry',
        'rrrrryyrr', 'rrrrryyry', 'rrrrryyyr', 'rrrrryyyy',
        'rrryryrrr', 'rrryryrry', 'rrryryyrr', 'rrryryyry',
        'rrryryyyr', 'rrryryyyy', 'rryrrrrrr', 'rryrrrrry',
        'rryrrryry', 'rryrryrrr', 'rryrryrry', 'rryrryryr',
        'rryrryyrr', 'rryrryyry', 'rryrryyyy', 'rryryryry',
        'rryyryrrr', 'rryyryrry', 'rryyryyry', 'rryyryyyy',
        'rryyyyrrr', 'rryyyyrry', 'rryyyyryr', 'rryyyyyrr',
        'ryrrrryrr', 'ryrrrryry', 'ryrrrryyr', 'ryrrrryyy',
        'ryrrryrrr', 'ryrrryrry', 'ryrrryyrr', 'ryrrryyry',
        'ryrrryyyr', 'ryrrryyyy', 'ryrryryyr', 'ryrryryyy',
        'ryrryyyrr', 'ryrryyyry', 'ryryryrrr', 'ryryryrry',
        'ryryryyrr', 'ryryryyry', 'ryryryyyr', 'ryryryyyy',
        'ryyrrrrrr', 'ryyrrrrry', 'ryyrrrryr', 'ryyrrryrr',
        'ryyrrryry', 'ryyrrryyr', 'ryyrrryyy', 'ryyrryrrr',
        'ryyrryrry', 'ryyrryryr', 'ryyrryyrr', 'ryyrryyry',
        'ryyrryyyr', 'ryyrryyyy', 'ryyryryrr', 'ryyryryry',
        'ryyryryyr', 'ryyryryyy', 'ryyryyrrr', 'ryyryyrry',
        'ryyryyyrr', 'ryyryyyry', 'ryyyryrrr', 'ryyyryrry',
        'ryyyryyrr', 'ryyyryyry', 'ryyyryyyr', 'ryyyryyyy',
        'ryyyyryyr', 'ryyyyryyy', 'ryyyyyryr', 'ryyyyyyrr',
        'ryyyyyyyy', 'yryrrrrrr', 'yryrrrrry', 'yryrrryry',
        'yryrryrrr', 'yryrryrry', 'yryrryryr', 'yryrryyrr',
        'yryrryyry', 'yryrryyyy', 'yryyryrrr', 'yryyryrry',
        'yryyryyry', 'yryyryyyy', 'yryyyyrrr', 'yryyyyrry',
        'yryyyyryr', 'yryyyyyrr', 'yyrrrryyr', 'yyrrrryyy',
        'yyrryryyr', 'yyrryryyy', 'yyyrrrrrr', 'yyyrrrrry',
        'yyyrrrryr', 'yyyrrryrr', 'yyyrrryry', 'yyyrrryyr',
        'yyyrrryyy', 'yyyrryrrr', 'yyyrryrry', 'yyyrryryr',
        'yyyrryyrr', 'yyyrryyry', 'yyyrryyyr', 'yyyrryyyy',
        'yyyryryrr', 'yyyryryry', 'yyyryryyr', 'yyyryryyy',
        'yyyyyryyr', 'yyyyyryyy', 'yyyyyyryr', 'yyyyyyyrr'
    ]
    table = str.maketrans("ACGT", "ryry")
    word_length = len(next(iter(patterns)))

    for sequence in sequences:
        ry_sequence = sequence.translate(table)
        for i in range(len(sequence) - k + 1):
            ry_pattern = ry_sequence[i:i + word_length]
            if ry_pattern in patterns:
                yield sequence[i:i + k]




# def count_kmers_start_ry_4_push(sequences: List[str], k: int) -> Counter:
#     """计算以特定 RY4_push 模式开头的 k-mer 频率"""
#     patterns = [
#         'rrrrryrrr', 'rrrrryrry', 'rrrrryryr', 'rrrrryryy',
#         'rrrrryyrr', 'rrrrryyry', 'rrrrryyyy', 'rrrryyyrr',
#         'rrrryyyry', 'rryrrrrry', 'rryrryrrr', 'rryrryrry',
#         'rryrryryr', 'rryryrrrr', 'rryryrrry', 'rryryrryr',
#         'rryryrryy', 'rryryryry', 'rryryyrrr', 'rryryyrry',
#         'rryryyryr', 'rryryyryy', 'rryryyyrr', 'rryryyyry',
#         'rryryyyyr', 'rryryyyyy', 'rryyyyrrr', 'rryyyyrry',
#         'rryyyyryr', 'rryyyyryy', 'ryrrrrrrr', 'ryrrrrrry',
#         'ryrrrryrr', 'ryrrrryry', 'ryrrrryyr', 'ryrrryrrr',
#         'ryrrryrry', 'ryrrryryr', 'ryrrryryy', 'ryrrryyrr',
#         'ryrrryyry', 'ryrrryyyy', 'ryrryyrrr', 'ryrryyrry',
#         'ryrryyryr', 'ryrryyryy', 'ryrryyyrr', 'ryrryyyry',
#         'ryrryyyyy', 'ryryryrrr', 'ryryryrry', 'ryryryyrr',
#         'ryryryyry', 'ryryryyyy', 'ryyrrrrrr', 'ryyrrrrry',
#         'ryyrrryry', 'ryyrryrrr', 'ryyrryrry', 'ryyrryryr',
#         'ryyrryryy', 'ryyrryyrr', 'ryyrryyry', 'ryyrryyyy',
#         'ryyryrrrr', 'ryyryrrry', 'ryyryrryy', 'ryyryryry',
#         'ryyryyrrr', 'ryyryyrry', 'ryyryyryr', 'ryyryyryy',
#         'ryyryyyrr', 'ryyryyyry', 'ryyryyyyy', 'ryyyyyrrr',
#         'ryyyyyrry', 'ryyyyyryr', 'ryyyyyryy', 'ryyyyyyrr',
#         'ryyyyyyry', 'yrrrryyyy', 'yrrryyyrr', 'yrrryyyry',
#         'yryryyyrr', 'yryryyyry', 'yyrrrryrr', 'yyrrrryry',
#         'yyrrrryyr', 'yyrrryrrr', 'yyrrryrry', 'yyrrryyrr',
#         'yyrrryyry', 'yyrrryyyy', 'yyrryyyrr', 'yyrryyyry',
#         'yyryrryrr', 'yyryrryry', 'yyryryrrr', 'yyryryrry',
#         'yyryryyrr', 'yyryryyry', 'yyryryyyy', 'yyryyyyrr',
#         'yyryyyyry', 'yyyrrrrrr', 'yyyrrrrry', 'yyyrrryry',
#         'yyyrryrrr', 'yyyrryrry', 'yyyrryryr', 'yyyrryryy',
#         'yyyrryyrr', 'yyyrryyry', 'yyyrryyyy', 'yyyryrrrr',
#         'yyyryrrry', 'yyyryrryy', 'yyyryryry', 'yyyryyrrr',
#         'yyyryyrry', 'yyyryyryr', 'yyyryyryy', 'yyyryyyrr',
#         'yyyryyyry', 'yyyryyyyy', 'yyyyyyyrr', 'yyyyyyyry'
#     ]
#     purines = {'A', 'G'}
#     pyrimidines = {'C', 'T'}
#     kmers = []
#     for sequence in sequences:
#         for i in range(len(sequence) - k + 1):
#             kmer = sequence[i:i + k]
#             # Convert kmer to RY pattern
#             ry_pattern = ''.join(['r' if base in purines else 'y' for base in kmer])
#             # Check if the k-mer starts with any of the patterns
#             if any(ry_pattern.startswith(pattern) for pattern in patterns):
#                 kmers.append(kmer)
#     return Counter(kmers)


def count_kmers_start_ry_4_push(sequences: List[str], k: int) -> Counter:
    """计算以特定 RY4_push 模式开头的 k-mer 频率"""
    patterns = [
        'rrrrryrrr', 'rrrrryrry', 'rrrrryryr', 'rrrrryryy',
        'rrrrryyrr', 'rrrrryyry', 'rrrrryyyy', 'rrrryyyrr',
        'rrrryyyry', 'rryrrrrry', 'rryrryrrr', 'rryrryrry',
        'rryrryryr', 'rryryrrrr', 'rryryrrry', 'rryryrryr',
        'rryryrryy', 'rryryryry', 'rryryyrrr', 'rryryyrry',
        'rryryyryr', 'rryryyryy', 'rryryyyrr', 'rryryyyry',
        'rryryyyyr', 'rryryyyyy', 'rryyyyrrr', 'rryyyyrry',
        'rryyyyryr', 'rryyyyryy', 'ryrrrrrrr', 'ryrrrrrry',
        'ryrrrryrr', 'ryrrrryry', 'ryrrrryyr', 'ryrrryrrr',
        'ryrrryrry', 'ryrrryryr', 'ryrrryryy', 'ryrrryyrr',
        'ryrrryyry', 'ryrrryyyy', 'ryrryyrrr', 'ryrryyrry',
        'ryrryyryr', 'ryrryyryy', 'ryrryyyrr', 'ryrryyyry',
        'ryrryyyyy', 'ryryryrrr', 'ryryryrry', 'ryryryyrr',
        'ryryryyry', 'ryryryyyy', 'ryyrrrrrr', 'ryyrrrrry',
        'ryyrrryry', 'ryyrryrrr', 'ryyrryrry', 'ryyrryryr',
        'ryyrryryy', 'ryyrryyrr', 'ryyrryyry', 'ryyrryyyy',
        'ryyryrrrr', 'ryyryrrry', 'ryyryrryy', 'ryyryryry',
        'ryyryyrrr', 'ryyryyrry', 'ryyryyryr', 'ryyryyryy',
        'ryyryyyrr', 'ryyryyyry', 'ryyryyyyy', 'ryyyyyrrr',
        'ryyyyyrry', 'ryyyyyryr', 'ryyyyyryy', 'ryyyyyyrr',
        'ryyyyyyry', 'yrrrryyyy', 'yrrryyyrr', 'yrrryyyry',
        'yryryyyrr', 'yryryyyry', 'yyrrrryrr', 'yyrrrryry',
        'yyrrrryyr', 'yyrrryrrr', 'yyrrryrry', 'yyrrryyrr',
        'yyrrryyry', 'yyrrryyyy', 'yyrryyyrr', 'yyrryyyry',
        'yyryrryrr', 'yyryrryry', 'yyryryrrr', 'yyryryrry',
        'yyryryyrr', 'yyryryyry', 'yyryryyyy', 'yyryyyyrr',
        'yyryyyyry', 'yyyrrrrrr', 'yyyrrrrry', 'yyyrrryry',
        'yyyrryrrr', 'yyyrryrry', 'yyyrryryr', 'yyyrryryy',
        'yyyrryyrr', 'yyyrryyry', 'yyyrryyyy', 'yyyryrrrr',
        'yyyryrrry', 'yyyryrryy', 'yyyryryry', 'yyyryyrrr',
        'yyyryyrry', 'yyyryyryr', 'yyyryyryy', 'yyyryyyrr',
        'yyyryyyry', 'yyyryyyyy', 'yyyyyyyrr', 'yyyyyyyry'
    ]
    table = str.maketrans("ACGT", "ryry")
    word_length = len(next(iter(patterns)))

    for sequence in sequences:
        ry_sequence = sequence.translate(table)
        for i in range(len(sequence) - k + 1):
            ry_pattern = ry_sequence[i:i + word_length]
            if ry_pattern in patterns:
                yield sequence[i:i + k]




# def count_kmers_start_ry_4_pull(sequences: List[str], k: int) -> Counter:
#     """计算以特定 RY4_pull 模式开头的 k-mer 频率"""
#     patterns = [
#         'rrrrrrrrr', 'rrrrrryrr', 'rrrrrryry', 'rrrrrryyr',
#         'rrrrrryyy', 'rrryryrrr', 'rrryryyrr', 'rrryryyry',
#         'rrryryyyr', 'rrryryyyy', 'rryrrrrrr', 'rryrrryrr',
#         'rryrrryry', 'rryrrryyr', 'rryrryrrr', 'rryrryrry',
#         'rryrryryr', 'rryrryyrr', 'rryrryyry', 'rryrryyyr',
#         'rryrryyyy', 'rryryryyr', 'rryryryyy', 'rryyryrrr',
#         'rryyryyrr', 'rryyryyry', 'rryyryyyr', 'rryyryyyy',
#         'rryyyryrr', 'rryyyyrrr', 'rryyyyyrr', 'rryyyyyyr',
#         'ryrrrryrr', 'ryrrrryry', 'ryrrrryyr', 'ryrrrryyy',
#         'ryrrryyyr', 'ryrrryyyy', 'ryrryryyr', 'ryrryryyy',
#         'ryryrryrr', 'ryryrryry', 'ryryrryyr', 'ryryrryyy',
#         'ryryryrrr', 'ryryryrry', 'ryryryryr', 'ryyrrrrrr',
#         'ryyrrrrry', 'ryyrrryrr', 'ryyrrryry', 'ryyrrryyr',
#         'ryyrrryyy', 'ryyrryrrr', 'ryyrryrry', 'ryyrryryr',
#         'ryyrryyrr', 'ryyrryyry', 'ryyrryyyr', 'ryyrryyyy',
#         'ryyryryrr', 'ryyryryry', 'ryyryryyr', 'ryyryryyy',
#         'ryyyrrrrr', 'ryyyrrrry', 'ryyyrryrr', 'ryyyrryry',
#         'ryyyrryyr', 'ryyyrryyy', 'ryyyryryr', 'ryyyryyrr',
#         'ryyyryyry', 'ryyyryyyr', 'ryyyryyyy', 'ryyyyryrr',
#         'ryyyyryry', 'ryyyyryyr', 'ryyyyryyy', 'ryyyyyyyr',
#         'yrrrrryrr', 'yrrrrryry', 'yrrrrryyr', 'yrrrrryyy',
#         'yrryryrrr', 'yryrrrrrr', 'yryrrryrr', 'yryrrryry',
#         'yryrrryyr', 'yryryryyr', 'yryryryyy', 'yryyryrrr',
#         'yryyryyrr', 'yryyryyry', 'yryyryyyr', 'yryyryyyy',
#         'yryyyryrr', 'yryyyyrrr', 'yryyyyyrr', 'yryyyyyyr',
#         'yyrrrryrr', 'yyrrrryry', 'yyrrrryyr', 'yyrrrryyy',
#         'yyrryryyr', 'yyrryryyy', 'yyryrryrr', 'yyryrryry',
#         'yyryrryyr', 'yyryrryyy', 'yyryryrrr', 'yyyrrryrr',
#         'yyyrrryry', 'yyyrrryyr', 'yyyrrryyy', 'yyyryryyr',
#         'yyyryryyy', 'yyyyrrrrr', 'yyyyrryrr', 'yyyyrryry',
#         'yyyyrryyr', 'yyyyrryyy', 'yyyyyryrr', 'yyyyyryry',
#         'yyyyyryyr', 'yyyyyryyy', 'yyyyyyyyr', 'yyyyyyyyy'
#     ]
#     purines = {'A', 'G'}
#     pyrimidines = {'C', 'T'}
#     kmers = []
#     for sequence in sequences:
#         for i in range(len(sequence) - k + 1):
#             kmer = sequence[i:i + k]
#             # Convert kmer to RY pattern
#             ry_pattern = ''.join(['r' if base in purines else 'y' for base in kmer])
#             # Check if the k-mer starts with any of the patterns
#             if any(ry_pattern.startswith(pattern) for pattern in patterns):
#                 kmers.append(kmer)
#     return Counter(kmers)


def count_kmers_start_ry_4_pull(sequences: List[str], k: int) -> Counter:
    """计算以特定 RY4_pull 模式开头的 k-mer 频率"""
    patterns = [
        'rrrrrrrrr', 'rrrrrryrr', 'rrrrrryry', 'rrrrrryyr',
        'rrrrrryyy', 'rrryryrrr', 'rrryryyrr', 'rrryryyry',
        'rrryryyyr', 'rrryryyyy', 'rryrrrrrr', 'rryrrryrr',
        'rryrrryry', 'rryrrryyr', 'rryrryrrr', 'rryrryrry',
        'rryrryryr', 'rryrryyrr', 'rryrryyry', 'rryrryyyr',
        'rryrryyyy', 'rryryryyr', 'rryryryyy', 'rryyryrrr',
        'rryyryyrr', 'rryyryyry', 'rryyryyyr', 'rryyryyyy',
        'rryyyryrr', 'rryyyyrrr', 'rryyyyyrr', 'rryyyyyyr',
        'ryrrrryrr', 'ryrrrryry', 'ryrrrryyr', 'ryrrrryyy',
        'ryrrryyyr', 'ryrrryyyy', 'ryrryryyr', 'ryrryryyy',
        'ryryrryrr', 'ryryrryry', 'ryryrryyr', 'ryryrryyy',
        'ryryryrrr', 'ryryryrry', 'ryryryryr', 'ryyrrrrrr',
        'ryyrrrrry', 'ryyrrryrr', 'ryyrrryry', 'ryyrrryyr',
        'ryyrrryyy', 'ryyrryrrr', 'ryyrryrry', 'ryyrryryr',
        'ryyrryyrr', 'ryyrryyry', 'ryyrryyyr', 'ryyrryyyy',
        'ryyryryrr', 'ryyryryry', 'ryyryryyr', 'ryyryryyy',
        'ryyyrrrrr', 'ryyyrrrry', 'ryyyrryrr', 'ryyyrryry',
        'ryyyrryyr', 'ryyyrryyy', 'ryyyryryr', 'ryyyryyrr',
        'ryyyryyry', 'ryyyryyyr', 'ryyyryyyy', 'ryyyyryrr',
        'ryyyyryry', 'ryyyyryyr', 'ryyyyryyy', 'ryyyyyyyr',
        'yrrrrryrr', 'yrrrrryry', 'yrrrrryyr', 'yrrrrryyy',
        'yrryryrrr', 'yryrrrrrr', 'yryrrryrr', 'yryrrryry',
        'yryrrryyr', 'yryryryyr', 'yryryryyy', 'yryyryrrr',
        'yryyryyrr', 'yryyryyry', 'yryyryyyr', 'yryyryyyy',
        'yryyyryrr', 'yryyyyrrr', 'yryyyyyrr', 'yryyyyyyr',
        'yyrrrryrr', 'yyrrrryry', 'yyrrrryyr', 'yyrrrryyy',
        'yyrryryyr', 'yyrryryyy', 'yyryrryrr', 'yyryrryry',
        'yyryrryyr', 'yyryrryyy', 'yyryryrrr', 'yyyrrryrr',
        'yyyrrryry', 'yyyrrryyr', 'yyyrrryyy', 'yyyryryyr',
        'yyyryryyy', 'yyyyrrrrr', 'yyyyrryrr', 'yyyyrryry',
        'yyyyrryyr', 'yyyyrryyy', 'yyyyyryrr', 'yyyyyryry',
        'yyyyyryyr', 'yyyyyryyy', 'yyyyyyyyr', 'yyyyyyyyy'
    ]
    table = str.maketrans("ACGT", "ryry")
    word_length = len(next(iter(patterns)))

    for sequence in sequences:
        ry_sequence = sequence.translate(table)
        for i in range(len(sequence) - k + 1):
            ry_pattern = ry_sequence[i:i + word_length]
            if ry_pattern in patterns:
                yield sequence[i:i + k]



# 1. basic k-mer matches
def basic_kmer_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = count_kmers(seq1, k)
    kmer_count2 = count_kmers(seq2, k)
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += kmer_count1[kmer] * kmer_count2[kmer]
    return matches

# # 2. spaced-word matches
# def spaced_word_matches(seq1: str, seq2: str, k: int) -> int:
#     """calculate the number of spaced-word matches between two sequences"""

#     pattern = generate_pattern(k)

#     words1 = extract_spaced_word(seq1, pattern)
#     words2 = extract_spaced_word(seq2, pattern)

#     word_count1 = Counter(words1)
#     word_count2 = Counter(words2)
#     matches = 0
#     for word in word_count1:
#         if word in word_count2:
#             matches += word_count1[word] * word_count2[word]
#     return matches

# 2. start_ry_matches
def start_ry_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq2_reverse_comple = reverse_complement(seq2)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2, seq2_reverse_comple] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = count_kmers_start_ry(seq1, k)
    kmer_count2 = count_kmers_start_ry(seq2, k)
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            # matches += kmer_count1[kmer] * kmer_count2[kmer]
            matches += 0.5 * min(kmer_count1[kmer], kmer_count2[kmer])
    return matches

# 3. start_rr_matches
def start_rr_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq2_reverse_comple = reverse_complement(seq2)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2, seq2_reverse_comple] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = count_kmers_start_rr(seq1, k)
    kmer_count2 = count_kmers_start_rr(seq2, k)
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += 0.5 * min(kmer_count1[kmer], kmer_count2[kmer])
    return matches 

# 4. start_ry_4_6_matches
def start_ry_4_6_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq2_reverse_comple = reverse_complement(seq2)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2, seq2_reverse_comple] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = Counter(count_kmers_start_ry_4_6(seq1, k))
    kmer_count2 = Counter(count_kmers_start_ry_4_6(seq2, k))
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += 0.5 * min(kmer_count1[kmer], kmer_count2[kmer])
    return matches

def start_ry_4_9_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq2_reverse_comple = reverse_complement(seq2)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2, seq2_reverse_comple] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = Counter(count_kmers_start_ry_4_9(seq1, k))
    kmer_count2 = Counter(count_kmers_start_ry_4_9(seq2, k))
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += 0.5 * min(kmer_count1[kmer], kmer_count2[kmer])
    return matches

def start_ry_4_push_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq2_reverse_comple = reverse_complement(seq2)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2, seq2_reverse_comple] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = Counter(count_kmers_start_ry_4_push(seq1, k))
    kmer_count2 = Counter(count_kmers_start_ry_4_push(seq2, k))
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += 0.5 * min(kmer_count1[kmer], kmer_count2[kmer])
    return matches

def start_ry_4_pull_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of k-mer matches between two sequences"""
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq2_reverse_comple = reverse_complement(seq2)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2, seq2_reverse_comple] # only one time reverse is okay, two reverse will make it slower
    kmer_count1 = Counter(count_kmers_start_ry_4_pull(seq1, k))
    kmer_count2 = Counter(count_kmers_start_ry_4_pull(seq2, k))
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += 0.5 * min(kmer_count1[kmer], kmer_count2[kmer])
    return matches


def spaced_word_matches(seq1: str, seq2: str, k: int, single_seq: bool) -> int:
    """calculate the number of spaced-word matches between two sequences"""
    # 根据 k 生成 pattern
    pattern = generate_pattern(k)
    
    if single_seq == 1:
        seq1 = [seq1]
        seq2 = [seq2]
    else:
        seq1_reverse_comple = reverse_complement(seq1)
        seq2_reverse_comple = reverse_complement(seq2)
        seq1 = [seq1, seq1_reverse_comple]
        seq2 = [seq2, seq2_reverse_comple] # only one time reverse is okay, two reverse will make it slower
    
    # 提取 spaced words
    spaced_words1 = []
    spaced_words2 = []
    
    for sequence in seq1:
        spaced_words1.extend(extract_spaced_word(sequence, pattern))
    
    for sequence in seq2:
        spaced_words2.extend(extract_spaced_word(sequence, pattern))
    
    # 计算匹配数
    word_count1 = Counter(spaced_words1)
    word_count2 = Counter(spaced_words2)
    matches = 0
    for word in word_count1:
        if word in word_count2:
            matches += 0.5 * min(word_count1[word], word_count2[word])
    return matches