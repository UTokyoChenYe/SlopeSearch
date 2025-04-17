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

from model.upper_model import compute_distance

# Initialize logger
logger = setup_logger()

def main():
    try:
        logger.info("Starting AF project evaluation")
        
        config_path = "../config/af_project_evaluation.yaml"
        logger.info(f"Loading configuration from {config_path}")
        with open(config_path) as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        args = EasyDict(config)

        data_path = args.get("data_path", "../dataset/test.fasta")
        k_mers_method = args.get("k_mers_method", "basic_kmer_matches")
        background_matches_method = args.get("background_matches_method", "static_method")
        bool_use_single_seq = args.get("bool_use_single_seq", True)
        bool_use_empirical_formula = args.get("bool_use_empirical_formula", True)

        logger.info(f"Loading sequences from {data_path}")
        names, seqs = load_sequences_for_evaluation(data_path)
        logger.info(f"Loaded {len(seqs)} sequences")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        save_path = args.get("output_path", "../result")+f"/af_project_evaluation_{timestamp}_{k_mers_method}_{background_matches_method}_{bool_use_single_seq}_{bool_use_empirical_formula}.tsv"
        logger.info(f"Saving results to {save_path}")

        with open(save_path, "w") as f:
            for i in range(len(seqs)):
                for j in range(len(seqs)):
                    dist = compute_distance(seqs[i], seqs[j], args)
                    f.write(f"{names[i]}\t{names[j]}\t{dist:.6f}\n")
        
        logger.info("Evaluation completed successfully")
    
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()