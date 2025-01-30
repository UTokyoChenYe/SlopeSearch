import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)

import math
import numpy as np
from tqdm import tqdm
from model.subsequence_method import basic_kmer_matches, spaced_word_matches


def calculate_match_probability(seq1: str, seq2: str, show_all_F_k: bool, single_seq: bool, method: str = "basic_kmer") -> float:
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

    if show_all_F_k == True:
        F_k = []
        k_values = list(range(1, 25))
        for k in tqdm(k_values, desc="Calculating F(k) for different k values"):
            if method == "basic_kmer":
                matches = basic_kmer_matches(seq1, seq2, k, single_seq)
            elif method == "spaced_word":
                matches = spaced_word_matches(seq1, seq2, k, single_seq)
            else:
                raise ValueError("Invalid align-free method.")
            
            match_counts.append(matches)
            F_k_value = np.log(matches)
            F_k.append(F_k_value)
    else:
        for k in range(k_min, k_max + 1):
            if method == "basic_kmer":
                matches = basic_kmer_matches(seq1, seq2, k)
            elif method == "spaced_word":
                matches = spaced_word_matches(seq1, seq2, k)
            else:
                raise ValueError("Invalid align-free method.")
            
            match_counts.append(matches)
    
    # Calculate match probability p, using the slope of F(k) = ln(N_k) between k_min and k_max to estimate the match probability
    x = np.array([k_min, k_max])
    y = np.array([match_counts[0], match_counts[-1]]) + 1e-10  # Prevent log(0) by adding a small value
    slope, _ = np.polyfit(x, np.log(y), 1)  # Calculate slope by fitting log of match counts between k_min and k_max
    
    # Estimate match probability p
    p_hat = np.exp(slope)

    if show_all_F_k == True:
        return F_k, p_hat
    else:
        return p_hat