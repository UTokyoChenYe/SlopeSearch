import random
from typing import List

def generate_random_sequence(length: int) -> str:
    """generate a random DNA sequence of a given length"""
    return ''.join(random.choice('ACGT') for _ in range(length))


def mutate_sequence(sequence: str, mutation_rate: float) -> str:
    """generate a mutated sequence with a given mutation rate"""
    mutated_sequence = list(sequence)
    for i in range(len(mutated_sequence)):
        if random.random() < mutation_rate:
            mutated_sequence[i] = random.choice('ACGT'.replace(mutated_sequence[i], ''))
    return ''.join(mutated_sequence)

def reverse_complement(sequence: str) -> str:
    """generate the reverse of seq"""
    complement = str.maketrans("ACGT", "TGCA") 
    return sequence.translate(complement)[::-1]  

# For spaced-word matches
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