# /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/codon_optimization.py

import pandas as pd
import numpy as np
import sys

def load_input(csv_file):
    """Load input sequences from a CSV file."""
    return pd.read_csv(csv_file)

def optimize_codon(sequence):
    """Perform codon optimization on the given sequence."""
    # Placeholder for actual optimization logic
    optimized_sequence = sequence  # Replace with actual optimization logic
    score = np.random.rand()  # Random score for demonstration
    return optimized_sequence, score

def save_results(task_id, iteration, aa_seq, dna_seq, score):
    """Save the results of the optimization to a log file."""
    log_file = f"optimization_log_{task_id}.txt"
    with open(log_file, 'a') as f:
        f.write(f"Iteration: {iteration}, AA Seq: {aa_seq}, DNA Seq: {dna_seq}, Score: {score}\n")

def main(task_id):
    # Load input sequences from CSV
    input_data = load_input('input_sequences.csv')
    
    # Iterate through the sequences
    for iteration, row in input_data.iterrows():
        aa_seq = row['AA_Seq']
        dna_seq, score = optimize_codon(aa_seq)
        save_results(task_id, iteration + 1, aa_seq, dna_seq, score)

if __name__ == "__main__":
    task_id = int(sys.argv[1])
    main(task_id)