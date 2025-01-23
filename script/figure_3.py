import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)

import yaml
from easydict import EasyDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

from utils.sequence_tool import generate_random_sequence, mutate_sequence
from model.F_k_function import calculate_match_probability
from model.jukes_cantor_model import estimate_jukes_cantor_distance

def main():
    config_path = "../config/run_fig_3.yaml"
    with open(config_path) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    args = EasyDict(config)

    num_pairs = args.get("num_pairs", 200)
    length = args.get("length", 1000)
    show_distance_min = args.get("show_distances_min", 0.00)
    show_distance_max = args.get("show_distances_max", 21)
    show_distance_step = args.get("show_distances_step", 1.0)
    figure_output_path = args.get("figure_output_path", "../result/figure/")
    show_all_F_k = args.get("show_all_F_k", False)

    distances = np.linspace(show_distance_min, show_distance_step, show_distance_max)

    # if you have new methods, you can add them here
    results_kmer = []
    results_spaced = []

    for dist in tqdm(distances, desc="Simulating sequence pairs in different distances: "):
        for _ in tqdm(range(num_pairs), desc="Calculating one distances value by different pairs: "):
            # generate random sequence and mutate it
            seq1 = generate_random_sequence(length)
            seq2 = mutate_sequence(seq1, dist)

            # calculate k-mer mode match probability
            p_hat_kmer = calculate_match_probability(seq1, seq2, show_all_F_k, method="basic_kmer")
            estimated_distance_kmer = estimate_jukes_cantor_distance(p_hat_kmer)
            results_kmer.append((dist, estimated_distance_kmer))

            # calculate spaced-word mode match probability
            p_hat_spaced = calculate_match_probability(seq1, seq2, show_all_F_k, method="spaced_word")
            estimated_distance_spaced = estimate_jukes_cantor_distance(p_hat_spaced)
            results_spaced.append((dist, estimated_distance_spaced))
    
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
    plt.savefig(figure_output_path + f'estimated_distance_vs_true_distance_L_{length}_pairs_{num_pairs}.png')
    plt.show()

if __name__ == "__main__":
    main()


