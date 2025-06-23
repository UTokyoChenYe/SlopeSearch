import sys
import os
from datetime import datetime
import concurrent.futures
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..', '..'))
sys.path.append(project_root)

import yaml
from easydict import EasyDict
from utils.file_system import *
from utils.logger import setup_logger

from model.upper_model import compute_distance

# Initialize logger
logger = setup_logger()

# Function to compute distance between two sequences
# This will be executed in parallel

def compute_distance_parallel(args):
    seq1, seq2, config_args = args
    return compute_distance(seq1, seq2, config_args)

def main():
    try:
        logger.info("Starting AF project evaluation with start_ry_4_6_matches")
        
        config_path = "../../config/assembled-fish_mito/af_project_evaluation_ry_4_6.yaml"
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
        names, seqs = load_sequences_for_evaluation_from_multiple_files(data_path)
        logger.info(f"Loaded {len(seqs)} sequences")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        save_path_tsv = args.get("output_path", "../result")+f"/af_project_evaluation_{timestamp}_{k_mers_method}_{background_matches_method}_{bool_use_single_seq}_{bool_use_empirical_formula}.tsv"
        save_path_phy = args.get("output_path", "../result")+f"/af_project_evaluation_{timestamp}_{k_mers_method}_{background_matches_method}_{bool_use_single_seq}_{bool_use_empirical_formula}.phy"

        # # TSV file format:
        # # name1 name2 distance
        # logger.info(f"Saving results to {save_path_tsv}")
        # with open(save_path_tsv, "w") as f:
        #     for i in range(len(seqs)):
        #         for j in range(len(seqs)):
        #             dist = compute_distance(seqs[i], seqs[j], args)
        #             f.write(f"{names[i]}\t{names[j]}\t{dist:.6f}\n")
        

        # PHYLIP file format:
        logger.info(f"Saving results to {save_path_phy}")
        N = len(seqs)
        distance_matrix = [[0.0 for _ in range(N)] for _ in range(N)]

        # Use ProcessPoolExecutor for parallel computation
        with concurrent.futures.ProcessPoolExecutor() as executor:
            # Prepare arguments for parallel execution
            tasks = ((seqs[i], seqs[j], args) for i in range(N) for j in range(N) if i != j)
            results = executor.map(compute_distance_parallel, tasks)

            # Fill the distance matrix with results
            for i in range(N):
                for j in range(N):
                    if i != j:
                        distance_matrix[i][j] = next(results)

        # 写入 PHYLIP 格式
        with open(save_path_phy, "w") as f:
            f.write(f"{N}\n")
            for i in range(N):
                name = names[i][:10].ljust(10)  # PHYLIP 格式要求名字最多10字符
                row = " ".join(f"{distance_matrix[i][j]:.6f}" for j in range(N))
                f.write(f"{name}{row}\n")
        logger.info("PHYLIP evaluation completed successfully")

        logger.info("Evaluation completed successfully")
    
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()