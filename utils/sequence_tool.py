import random

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