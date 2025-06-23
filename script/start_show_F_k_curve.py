import sys
import os
from datetime import datetime
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)

import yaml
from easydict import EasyDict
from utils.file_system import load_sequences_for_evaluation
from utils.logger import setup_logger

import matplotlib.pyplot as plt

from model.F_k_function import F_k_funtion

# Initialize logger
logger = setup_logger()

def main():
    try:
        logger.info("Starting show F(k) curve")
        
        config_path = "../config/show_F_k_curve.yaml"
        logger.info(f"Loading configuration from {config_path}")
        with open(config_path) as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        args = EasyDict(config)

        data_path = args.get("data_path", "../dataset/test.fasta")
        k_show_values = args.get("k_show_values", list(range(2, 25)))

        # test method parameters
        test_method = args.get("test method", {})
        test_k_mers_method = test_method.get("k_mers_method", "basic_kmer_matches")
        test_background_matches_method = test_method.get("background_matches_method", "static_method")
        test_bool_use_single_seq = test_method.get("bool_use_single_seq", True)
        test_bool_use_empirical_formula = test_method.get("bool_use_empirical_formula", True)

        # control method parameters
        control_method = args.get("control method", {})

        logger.info(f"Loading sequences from {data_path}")
        names, seqs = load_sequences_for_evaluation(data_path)
        logger.info(f"Loaded {len(seqs)} sequences")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        save_path = args.get("output_path", "../result") + f"/show_F_k_curve_{timestamp}_{test_k_mers_method}_{test_background_matches_method}_{test_bool_use_single_seq}_{test_bool_use_empirical_formula}.png"
        logger.info(f"Saving results to {save_path}")

        # 1. calculate F(k)
        F_k_function_test = F_k_funtion(seqs[0], seqs[1], test_method)
        F_k_function_control = F_k_funtion(seqs[0], seqs[1], control_method)

        F_k_show_test = F_k_function_test.show_F_k_curve()
        F_k_show_control = F_k_function_control.show_F_k_curve()

        # 2. plot F(k) curve
        plt.figure(figsize=(10, 6))
        plt.plot(k_show_values[:len(F_k_show_test)], F_k_show_test, marker='o', linestyle='-', color='b', label='F_k_ry')
        plt.plot(k_show_values[:len(F_k_show_control)], F_k_show_control, marker='s', linestyle='-', color='r', label='F_k_basic')

        # Add start_rr_matches plot
        F_k_function_rr = F_k_funtion(seqs[0], seqs[1], {"k_mers_method": "start_rr_matches", "bool_use_single_seq": test_bool_use_single_seq})
        F_k_show_rr = F_k_function_rr.show_F_k_curve()
        plt.plot(k_show_values[:len(F_k_show_rr)], F_k_show_rr, marker='^', linestyle='-', color='g', label='F_k_rr')

        plt.xlabel('Word Length (k)')
        plt.ylabel('F(k)')
        plt.title('Relationship between F(k) and Word Length (k)')
        plt.grid(True)
        plt.legend()
        plt.savefig(save_path)
    
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
        
