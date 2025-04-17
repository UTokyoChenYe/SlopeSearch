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
from utils.logger import setup_logger

logger = setup_logger()

class F_k_funtion:
    def __init__(self, seq1:str, seq2:str, args:dict):
        self.k_mers_method = args.get("k_mers_method", "basic_kmer_matches") # method to calculate F_k
        self.bool_use_single_seq = args.get("bool_use_single_seq", True) # whether taking reverse complement into account
        self.bool_use_empirical_formula = args.get("bool_use_empirical_formula", True) # whether using empirical formula to calculate k_min and k_max
        self.background_matches_method = args.get("background_matches_method", "basic_kmer_matches") # method to calculate background matches
        self.k_show_values = args.get("k_show_values", list(range(2, 25))) # k values to show
        self.seq1 = seq1
        self.seq2 = seq2
        self.L_1 = len(seq1)
        self.L_2 = len(seq2)
        self.L = (self.L_1 + self.L_2) / 2  # Average length of the two sequences

        if self.bool_use_empirical_formula:
            logger.info("Using empirical formula to calculate k_min and k_max")
            # using Empirical formula to calculate k_min and k_max
            self.k_min = math.ceil((math.log(self.L) + 0.69) / 0.875)
            self.k_max = math.floor(math.log(self.L) / 0.634)
            logger.info(f"k_min: {self.k_min}, k_max: {self.k_max}")
        else:
            logger.info("Not using empirical formula to calculate k_min and k_max")
            # TODO: reproduct John's idea to calculate k_min and k_max
            pass

        # calculate F_k for different k values
        self.F_k_p_hat = []
        self.F_k_show = []

        
    def calculate_p_hat(self):
        logger.info("Calculating p_hat")
        for k in tqdm(range(self.k_min, self.k_max + 1), desc="Calculating F(k) for different k values"):
            if self.k_mers_method == "basic_kmer_matches":
                logger.info("Using basic kmer matches to calculate F(k)")
                matches = basic_kmer_matches(self.seq1, self.seq2, k, self.bool_use_single_seq)
            elif self.k_mers_method == "start_ry_matches":
                logger.info("Using start_ry_matches to calculate F(k)")
                matches = start_ry_matches(self.seq1, self.seq2, k, self.bool_use_single_seq)
            else:
                logger.error("Invalid align-free method.")
                raise ValueError("Invalid align-free method.")
            
            if self.background_matches_method == "basic_kmer_matches":
                logger.info("Using basic kmer matches to calculate background matches")
                background_matches = basic_kmer_matches(self.seq1, reverse(self.seq2), k, self.bool_use_single_seq)
            elif self.background_matches_method == "static_method":
                logger.info("Using static method to calculate background matches")
                background_matches = 2 * self.L_1 * self.L_2 * (1/4)**k
            elif self.background_matches_method == "no_background_matches":
                logger.info("No background matches")
                background_matches = 0
            else:
                logger.error("Invalid background matches method.")
                raise ValueError("Invalid background matches method.")
            
            self.F_k_p_hat.append(np.log(matches - background_matches))

        # calculate p_hat
        x = np.array([self.k_min, self.k_max])
        y = np.array([self.F_k_p_hat[0], self.F_k_p_hat[-1]]) + 1e-10  # Prevent log(0) by adding a small value
        slope, _ = np.polyfit(x, np.log(y), 1)  # Calculate slope by fitting log of match counts between k_min and k_max
        self.p_hat = np.exp(slope)
        logger.info(f"p_hat: {self.p_hat}")

        return self.p_hat


    def show_F_k_curve(self):
        logger.info("Showing F(k) curve")
        for k in tqdm(self.k_show_values, desc="Showing F(k) curve"):
            if self.k_mers_method == "basic_kmer_matches":
                logger.info("Using basic kmer matches to calculate F(k)")
                matches = basic_kmer_matches(self.seq1, self.seq2, k, self.bool_use_single_seq)
            elif self.k_mers_method == "start_ry_matches":
                logger.info("Using start_ry_matches to calculate F(k)")
                matches = start_ry_matches(self.seq1, self.seq2, k, self.bool_use_single_seq)
            else:
                logger.error("Invalid align-free method.")
                raise ValueError("Invalid align-free method.")

            if self.background_matches_method == "basic_kmer_matches":
                logger.info("Using basic kmer matches to calculate background matches")
                background_matches = basic_kmer_matches(self.seq1, reverse(self.seq2), k, self.bool_use_single_seq)
            elif self.background_matches_method == "static_method":
                logger.info("Using static method to calculate background matches")
                background_matches = 2 * self.L_1 * self.L_2 * (1/4)**k
            elif self.background_matches_method == "no_background_matches":
                logger.info("No background matches")
                background_matches = 0
            else:
                logger.error("Invalid background matches method.")
                raise ValueError("Invalid background matches method.")
            
            self.F_k_show.append(np.log(matches - background_matches))

        return self.F_k_show


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