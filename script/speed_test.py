import sys
import os
import time
from datetime import datetime
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)

import yaml
from easydict import EasyDict
from utils.file_system import load_sequences_for_evaluation_from_multiple_files
from utils.logger import setup_logger
import pandas as pd

from model.upper_model import compute_distance

# Initialize logger
logger = setup_logger()

def measure_matrix_build_time(seqs, names, method_config, method_name):
    """测量构建距离矩阵的时间"""
    logger.info(f"Testing method: {method_name}")
    
    start_time = time.time()
    
    try:
        N = len(seqs)
        distance_matrix = [[0.0 for _ in range(N)] for _ in range(N)]
        
        # 构建距离矩阵
        for i in range(N):
            for j in range(N):
                if i == j:
                    distance_matrix[i][j] = 0.0
                else:
                    distance_matrix[i][j] = compute_distance(seqs[i], seqs[j], method_config)
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        logger.info(f"{method_name} completed in {execution_time:.4f} seconds")
        logger.info(f"Matrix size: {N}x{N}, Total pairs: {N*(N-1)//2}")
        
        return execution_time, distance_matrix
        
    except Exception as e:
        end_time = time.time()
        execution_time = end_time - start_time
        logger.error(f"Error in {method_name}: {str(e)}")
        return execution_time, None

def main():
    try:
        logger.info("Starting speed test for distance matrix construction")
        
        # 加载配置文件
        config_path = "../config/speed_test/speed_test.yaml"
        logger.info(f"Loading configuration from {config_path}")
        with open(config_path) as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        args = EasyDict(config)
        
        # 获取数据路径
        data_path = args.get("data_path", "../dataset/af_project_dataset/AF-reference_datasets190511/genome-std-assembled-fish_mito/dataset/assembled-fish_mito")
        
        # 修复路径：如果路径是相对路径，需要相对于项目根目录解析
        if not os.path.isabs(data_path):
            # 相对于项目根目录解析路径
            data_path = os.path.abspath(os.path.join(project_root, data_path))
        
        logger.info(f"Resolved data path: {data_path}")
        
        # 检查路径是否存在
        if not os.path.exists(data_path):
            logger.error(f"Data path does not exist: {data_path}")
            raise FileNotFoundError(f"Data path not found: {data_path}")
        
        # 加载序列数据
        logger.info(f"Loading sequences from {data_path}")
        names, seqs = load_sequences_for_evaluation_from_multiple_files(data_path)
        logger.info(f"Loaded {len(seqs)} sequences")
        
        # 检查序列数量
        if len(seqs) < 2:
            logger.error(f"Need at least 2 sequences for distance matrix construction. Found {len(seqs)} sequences.")
            raise ValueError(f"Insufficient sequences: {len(seqs)} (need at least 2)")
        
        # 显示序列信息
        logger.info(f"Sequence names: {names[:5]}...")  # 显示前5个序列名
        logger.info(f"Sequence lengths: {[len(seq) for seq in seqs[:3]]}...")  # 显示前3个序列长度
        logger.info(f"Total pairs to calculate: {len(seqs) * (len(seqs) - 1) // 2}")
        
        # 从配置文件获取方法配置
        methods_config = args.get("methods", {})
        
        # 执行速度测试
        results = []
        for method_key, method_config in methods_config.items():
            method_name = method_config.get("name", method_key)
            logger.info(f"Testing method: {method_name}")
            
            # 提取方法参数
            method_params = {
                "k_mers_method": method_config.get("k_mers_method"),
                "bool_use_single_seq": method_config.get("bool_use_single_seq", True),
                "bool_use_empirical_formula": method_config.get("bool_use_empirical_formula", True),
                "background_matches_method": method_config.get("background_matches_method", "no_background_matches")
            }
            
            execution_time, distance_matrix = measure_matrix_build_time(seqs, names, method_params, method_name)
            
            # 计算平均距离（排除对角线）
            avg_distance = None
            if distance_matrix:
                total_distance = sum(sum(row) for row in distance_matrix)
                avg_distance = total_distance / (len(seqs) * len(seqs) - len(seqs))  # 排除对角线
            
            results.append({
                "Method": method_name,
                "Execution_Time_Seconds": execution_time,
                "Matrix_Size": f"{len(seqs)}x{len(seqs)}",
                "Total_Pairs": len(seqs) * (len(seqs) - 1) // 2,
                "Average_Distance": avg_distance
            })
        
        # 创建结果DataFrame
        df = pd.DataFrame(results)
        
        # 打印结果
        logger.info("\n" + "="*70)
        logger.info("DISTANCE MATRIX CONSTRUCTION SPEED TEST RESULTS")
        logger.info("="*70)
        for _, row in df.iterrows():
            logger.info(f"{row['Method']:<15}: {row['Execution_Time_Seconds']:.4f}s | {row['Matrix_Size']} | {row['Total_Pairs']} pairs")
        logger.info("="*70)
        
        # 保存结果到CSV文件
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = args.get("output_path", "../result/speed_test")
        
        # 修复输出路径：如果路径是相对路径，需要相对于项目根目录解析
        if not os.path.isabs(output_path):
            output_path = os.path.abspath(os.path.join(project_root, output_path))
        
        logger.info(f"Resolved output path: {output_path}")
        
        # 确保输出目录存在
        os.makedirs(output_path, exist_ok=True)
        
        csv_file = f"{output_path}/matrix_speed_test_results_{timestamp}.csv"
        df.to_csv(csv_file, index=False)
        logger.info(f"Results saved to {csv_file}")
        
        # 找出最快和最慢的方法
        fastest_method = df.loc[df['Execution_Time_Seconds'].idxmin()]
        slowest_method = df.loc[df['Execution_Time_Seconds'].idxmax()]
        
        logger.info(f"\nFastest method: {fastest_method['Method']} ({fastest_method['Execution_Time_Seconds']:.4f}s)")
        logger.info(f"Slowest method: {slowest_method['Method']} ({slowest_method['Execution_Time_Seconds']:.4f}s)")
        
        # 计算相对速度
        logger.info("\nRelative speed comparison:")
        for _, row in df.iterrows():
            relative_speed = fastest_method['Execution_Time_Seconds'] / row['Execution_Time_Seconds']
            logger.info(f"{row['Method']:<15}: {relative_speed:.2f}x slower than fastest")
        
        # 计算每对序列的平均时间
        logger.info("\nAverage time per sequence pair:")
        for _, row in df.iterrows():
            if row['Total_Pairs'] > 0:
                avg_time_per_pair = row['Execution_Time_Seconds'] / row['Total_Pairs']
                logger.info(f"{row['Method']:<15}: {avg_time_per_pair:.6f}s per pair")
            else:
                logger.info(f"{row['Method']:<15}: No pairs to calculate average time")
        
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
