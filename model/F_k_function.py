import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)

import math
import numpy as np
from model.subsequence_method import basic_kmer_matches, spaced_word_matches


def calculate_match_probability(seq1: str, seq2: str, method: str = "basic_kmer") -> float:
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
        if method == "basic_kmer":
            matches = basic_kmer_matches(seq1, seq2, k)
        elif method == "spaced_word":
            matches = spaced_word_matches(seq1, seq2, k)
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