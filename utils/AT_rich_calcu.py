import os
from Bio import SeqIO
from typing import Dict, List, Tuple

def calculate_at_richness(sequence: str) -> float:
    """
    Calculate AT richness (percentage of A and T nucleotides) in a sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        float: AT richness percentage
    """
    sequence = sequence.upper()
    total_length = len(sequence)
    if total_length == 0:
        return 0.0
    
    at_count = sequence.count('A') + sequence.count('T')
    return (at_count / total_length) * 100

def process_fasta_file(fasta_path: str) -> Dict[str, float]:
    """
    Process a FASTA file and calculate AT richness for each sequence.
    
    Args:
        fasta_path (str): Path to the FASTA file
        
    Returns:
        Dict[str, float]: Dictionary with sequence IDs as keys and AT richness as values
    """
    results = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            at_richness = calculate_at_richness(str(record.seq))
            results[record.id] = at_richness
    except Exception as e:
        print(f"Error processing {fasta_path}: {str(e)}")
    return results

def analyze_directory(directory_path: str) -> List[Tuple[str, Dict[str, float]]]:
    """
    Analyze all FASTA files in a directory.
    
    Args:
        directory_path (str): Path to the directory containing FASTA files
        
    Returns:
        List[Tuple[str, Dict[str, float]]]: List of tuples containing file names and their AT richness results
    """
    results = []
    for filename in os.listdir(directory_path):
        if filename.endswith(('.fasta', '.fa')):
            file_path = os.path.join(directory_path, filename)
            at_richness_results = process_fasta_file(file_path)
            results.append((filename, at_richness_results))
    return results

def main():
    # Path to the directory containing FASTA files
    # directory_path = "/home/chenye/project/slope-spam/slope-spam-python/dataset/af_project_dataset/AF-reference_datasets190511/genome-std-assembled-fish_mito/dataset/assembled-fish_mito"
    directory_path = "/home/chenye/project/slope-spam/slope-spam-python/dataset/af_project_dataset/AF-reference_datasets190511/genome-std-assembled-ecoli/dataset/assembled-ecoli"

    # Analyze all FASTA files in the directory
    results = analyze_directory(directory_path)
    
    # Print results
    for filename, at_richness_dict in results:
        print(f"\nResults for {filename}:")
        for seq_id, at_richness in at_richness_dict.items():
            print(f"Sequence {seq_id}: AT richness = {at_richness:.2f}%")

if __name__ == "__main__":
    main()
