import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)

import math
import numpy as np
from tqdm import tqdm

from model.subsequence_method import *
from utils.sequence_tool import *

class F_k_function_paper:
    def __init__(self, seq1: str, seq2: str, method, args):
        self.method = method
        self.single_seq = args.get("single_seq", True)
        self.L_1 = len(seq1)
        self.L_2 = len(seq2)
        self.L = (self.L_1 + self.L_2) / 2  # Average length of the two sequences

        # Calculate k_min and k_max based on the formulas from the paper
        self.k_min = math.ceil((math.log(self.L) + 0.69) / 0.875)
        self.k_max = math.floor(math.log(self.L) / 0.634)

        self.F_k = []

        self.k_values = list(range(2, 25))

        for k in tqdm(self.k_values, desc="Calculating F(k) for different k values"):
            if self.method == "basic_kmer_matches":
                matches = basic_kmer_matches(seq1, seq2, k, self.single_seq)
            # elif self.method == "spaced_word":
            #     matches = spaced_word_matches(seq1, seq2, k, self.single_seq)
            elif self.method == "start_ry_matches":
                matches = start_ry_matches(seq1, seq2, k, self.single_seq)
            else:
                raise ValueError("Invalid align-free method.")
            
            F_k_value = np.log(matches - 2 * self.L_1 * self.L_2 * (1/4)**k)
            
            self.F_k.append(F_k_value)
    
    def get_F_k(self):
        return self.F_k


class F_k_function_v1:
    def __init__(self, seq1: str, seq2: str, method, args):
        self.method = method
        self.single_seq = args.get("single_seq", True)
        self.L_1 = len(seq1)
        self.L_2 = len(seq2)
        self.L = (self.L_1 + self.L_2) / 2  # Average length of the two sequences

        # Calculate k_min and k_max based on the formulas from the paper
        self.k_min = math.ceil((math.log(self.L) + 0.69) / 0.875)
        self.k_max = math.floor(math.log(self.L) / 0.634)

        self.F_k = []

        self.k_values = list(range(2, 25))

        for k in tqdm(self.k_values, desc="Calculating F(k) for different k values"):
            if self.method == "basic_kmer_matches":
                matches = basic_kmer_matches(seq1, seq2, k, self.single_seq)
            # elif self.method == "spaced_word":
            #     matches = spaced_word_matches(seq1, seq2, k, self.single_seq)
            elif self.method == "start_ry_matches":
                matches = start_ry_matches(seq1, seq2, k, self.single_seq)
            else:
                raise ValueError("Invalid align-free method.")
            
            F_k_value = np.log(matches)
            
            self.F_k.append(F_k_value)
    
    def get_F_k(self):
        return self.F_k


class F_k_function_v2:
    def __init__(self, seq1: str, seq2: str, method,  args):
        self.method = method
        self.single_seq = args.get("single_seq", True)
        self.L_1 = len(seq1)
        self.L_2 = len(seq2)
        self.L = (self.L_1 + self.L_2) / 2  # Average length of the two sequences

        # Calculate k_min and k_max based on the formulas from the paper
        self.k_min = math.ceil((math.log(self.L) + 0.69) / 0.875)
        self.k_max = math.floor(math.log(self.L) / 0.634)

        self.F_k = []

        self.k_values = list(range(2, 25))


        for k in tqdm(self.k_values, desc="Calculating F(k) for different k values"):
            if self.method == "basic_kmer_matches":
                matches = basic_kmer_matches(seq1, seq2, k, self.single_seq)
            # elif self.method == "spaced_word":
            #     matches = spaced_word_matches(seq1, seq2, k, self.single_seq)
            elif self.method == "start_ry_matches":
                matches = start_ry_matches(seq1, seq2, k, self.single_seq)
            else:
                raise ValueError("Invalid align-free method.")
            
            if self.method == "start_ry_matches":
                non_homologous_matches = basic_kmer_matches(seq1, reverse(seq2), k, self.single_seq) / 4
            else:
                non_homologous_matches = basic_kmer_matches(seq1, reverse(seq2), k, self.single_seq)
            
            diff = matches - non_homologous_matches
            if diff <= 0:
                self.F_k.append(np.nan)
            else:
                F_k_value = np.log(diff)
                self.F_k.append(F_k_value)
    
    def get_F_k(self):
        return self.F_k




# def calculate_match_probability(seq1: str, seq2: str, show_all_F_k: bool, single_seq: bool, method: str = "basic_kmer") -> float:
#     """
#     Calculate the match probability for different k values using k-mer or spaced-word methods.
#     :param seq1: str, first DNA sequence
#     :param seq2: str, second DNA sequence
#     :param method: str, method to use ('kmer' or 'spaced_word')
#     :return: float, estimated match probability
#     """
#     L_1 = len(seq1)
#     L_2 = len(seq2)
#     L = (L_1 + L_2) / 2  # Average length of the two sequences
    
#     # Calculate k_min and k_max based on the formulas from the paper
#     k_min = math.ceil((math.log(L) + 0.69) / 0.875)
#     k_max = math.floor(math.log(L) / 0.634)
    
#     match_counts = []

#     if show_all_F_k == True:
#         F_k = []
#         # k_values = list(range(1, 25))
#         k_values = list(range(2, 25))
#         for k in tqdm(k_values, desc="Calculating F(k) for different k values"):
#             if method == "basic_kmer_matches":
#                 matches = basic_kmer_matches(seq1, seq2, k, single_seq)
#             elif method == "spaced_word":
#                 matches = spaced_word_matches(seq1, seq2, k, single_seq)
#             elif method == "start_ry_matches":
#                 matches = start_ry_matches(seq1, seq2, k, single_seq)
#             else:
#                 raise ValueError("Invalid align-free method.")
            
#             match_counts.append(matches)
#             F_k_value = np.log(matches)
#             F_k.append(F_k_value)
#     else:
#         for k in range(k_min, k_max + 1):
#             if method == "basic_kmer":
#                 matches = basic_kmer_matches(seq1, seq2, k)
#             elif method == "spaced_word":
#                 matches = spaced_word_matches(seq1, seq2, k)
#             else:
#                 raise ValueError("Invalid align-free method.")
            
#             match_counts.append(matches)
    
#     # Calculate match probability p, using the slope of F(k) = ln(N_k) between k_min and k_max to estimate the match probability
#     x = np.array([k_min, k_max])
#     y = np.array([match_counts[0], match_counts[-1]]) + 1e-10  # Prevent log(0) by adding a small value
#     slope, _ = np.polyfit(x, np.log(y), 1)  # Calculate slope by fitting log of match counts between k_min and k_max
    
#     # Estimate match probability p
#     p_hat = np.exp(slope)

#     if show_all_F_k == True:
#         return F_k, p_hat
#     else:
#         return p_hat