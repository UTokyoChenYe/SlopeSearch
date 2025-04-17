import numpy as np
import matplotlib.pyplot as plt
import random
import math
from typing import List
from collections import Counter
import pandas as pd
from tqdm import tqdm

# 1. generate random DNA sequences and mutate them
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

# 2. calculate k-mer frequencies
def count_kmers(sequence: str, k: int) -> Counter:
    """calculate k-mer frequencies in a given sequence"""
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return Counter(kmers)

def kmer_matches(seq1: str, seq2: str, k: int) -> int:
    """calculate the number of k-mer matches between two sequences"""
    kmer_count1 = count_kmers(seq1, k)
    kmer_count2 = count_kmers(seq2, k)
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += kmer_count1[kmer] * kmer_count2[kmer]
    return matches

# 3. calculate spaced-word matches
def spaced_word_matches(seq1: str, seq2: str, pattern: str) -> int:
    """calculate the number of spaced-word matches between two sequences"""
    def extract_spaced_word(sequence: str, pattern: str) -> List[str]:
        words = []
        pattern_length = len(pattern)
        for i in range(len(sequence) - pattern_length + 1):
            spaced_word = ''.join(sequence[i+j] for j in range(pattern_length) if pattern[j] == '1')
            words.append(spaced_word)
        return words

    words1 = extract_spaced_word(seq1, pattern)
    words2 = extract_spaced_word(seq2, pattern)

    word_count1 = Counter(words1)
    word_count2 = Counter(words2)
    matches = 0
    for word in word_count1:
        if word in word_count2:
            matches += word_count1[word] * word_count2[word]
    return matches

# 4. calculate match probability at different k values
def generate_pattern(length):
    """Generate a random binary pattern of given length."""
    return ''.join(random.choice(['0', '1']) for _ in range(length))

def calculate_match_probability(seq1: str, seq2: str, method: str = "kmer") -> float:
    """
    Calculate the match probability for different k values using k-mer or spaced-word methods.
    :param seq1: str, first DNA sequence
    :param seq2: str, second DNA sequence
    :param method: str, method to use ('kmer' or 'spaced_word')
    :return: float, estimated match probability
    """
    L = len(seq1)  # Assume seq1 and seq2 have the same length
    
    # Calculate k_min and k_max based on the formulas from the paper
    k_min = math.ceil((math.log(L) + 0.69) / 0.875)
    k_max = math.floor(math.log(L) / 0.634)
    
    match_counts = []
    for k in range(k_min, k_max + 1):
        if method == "kmer":
            matches = kmer_matches(seq1, seq2, k)
        elif method == "spaced_word":
            # Generate a random pattern of length k
            pattern_k = generate_pattern(k)
            matches = spaced_word_matches(seq1, seq2, pattern_k)
        else:
            raise ValueError("Invalid method. Use 'kmer' or 'spaced_word'.")
        
        match_counts.append(matches)
    
    # Calculate match probability p, using the slope of F(k) = ln(N_k) between k_min and k_max to estimate the match probability
    x = np.array([k_min, k_max])
    y = np.array([match_counts[0], match_counts[-1]]) + 1e-10  # Prevent log(0) by adding a small value
    slope, _ = np.polyfit(x, np.log(y), 1)  # Calculate slope by fitting log of match counts between k_min and k_max
    
    # Estimate match probability p
    p_hat = np.exp(slope)
    return p_hat




# 5. estimate Jukes-Cantor distance from match probability
def estimate_jukes_cantor_distance(p_hat: float) -> float:
    """estimate Jukes-Cantor distance from match probability"""
    if p_hat >= 1:
        p_hat = 0.999  # prevent log(0) when p_hat = 1
    try:
        d = -3/4 * np.log(1 - 4/3 * (1 - p_hat))
        return max(0, d)  # prevent negative distances
    except ValueError:
        return 0  # return 0 if the formula fails

# 6. simulate sequence pairs and calculate estimated distances
def simulate_and_calculate():
    # num_pairs = 20000  # number of sequence pairs to simulate
    # length = 100000  # sequence length
    num_pairs = 200  # number of sequence pairs to simulate
    length = 1000  # sequence length
    # distances = np.linspace(0.05, 1.0, 10)  # true Jukes-Cantor distances to simulate
    distances = np.linspace(0.00, 1.0, 21)
    

    results_kmer = []
    results_spaced = []

    for dist in tqdm(distances, desc="Simulating sequence pairs"):
        for _ in tqdm(range(num_pairs), desc="pairs for 1 distance: "):
            # generate random sequence and mutate it
            seq1 = generate_random_sequence(length)
            seq2 = mutate_sequence(seq1, dist)

            # calculate k-mer mode match probability
            p_hat_kmer = calculate_match_probability(seq1, seq2, method="kmer")
            estimated_distance_kmer = estimate_jukes_cantor_distance(p_hat_kmer)
            results_kmer.append((dist, estimated_distance_kmer))

            # calculate spaced-word mode match probability
            p_hat_spaced = calculate_match_probability(seq1, seq2, method="spaced_word")
            estimated_distance_spaced = estimate_jukes_cantor_distance(p_hat_spaced)
            results_spaced.append((dist, estimated_distance_spaced))

    # transform results into pandas DataFrame
    df_kmer = pd.DataFrame(results_kmer, columns=['true_distance', 'estimated_distance'])
    df_spaced = pd.DataFrame(results_spaced, columns=['true_distance', 'estimated_distance'])

    # calculate average and standard deviation of estimated distances for each true distance
    average_estimates_kmer = df_kmer.groupby('true_distance')['estimated_distance'].mean()
    std_estimates_kmer = df_kmer.groupby('true_distance')['estimated_distance'].std()

    average_estimates_spaced = df_spaced.groupby('true_distance')['estimated_distance'].mean()
    std_estimates_spaced = df_spaced.groupby('true_distance')['estimated_distance'].std()

    # plot the results
    plt.figure(figsize=(12, 8))

    # k-mer mode plot
    plt.errorbar(average_estimates_kmer.index, average_estimates_kmer, yerr=std_estimates_kmer,
                 fmt='o-', color='b', ecolor='r', capsize=5, label='k-mer Mode')

    # spaced-word mode plot
    plt.errorbar(average_estimates_spaced.index, average_estimates_spaced, yerr=std_estimates_spaced,
                 fmt='o-', color='g', ecolor='orange', capsize=5, label='Spaced-word Mode')
    
    # ideal line plot
    plt.plot([0, 1], [0, 1], linestyle='--', color='black', label='Ideal Line (y = x)')

    # set plot properties
    plt.ylim(0, 1)

    plt.xlabel('True Jukes-Cantor Distance')
    plt.ylabel('Average Estimated Distance')
    plt.title('Average Estimated Distance vs True Jukes-Cantor Distance')
    plt.grid(True)
    plt.legend()
    plt.savefig('estimated_distance_vs_true_distance.png')
    plt.show()

if __name__ == "__main__":
    simulate_and_calculate()