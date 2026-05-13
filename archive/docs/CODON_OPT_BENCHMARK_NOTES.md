### 1. Python Script: `codon_optimization.py`

```python
import pandas as pd
import numpy as np
import sys

def load_input_data(csv_file):
    """Load input data from a CSV file."""
    return pd.read_csv(csv_file)

def optimize_codon(sequence):
    """Perform codon optimization on the given sequence."""
    # Placeholder for actual optimization logic
    optimized_sequence = sequence  # Replace with actual optimization logic
    score = np.random.rand()  # Replace with actual scoring logic
    return optimized_sequence, score

def save_results(task_id, iteration, aa_seq, dna_seq, score):
    """Save the results of the optimization to a log file."""
    log_file = f"/home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/log_{task_id}.txt"
    with open(log_file, 'a') as f:
        f.write(f"Iteration: {iteration}, AA Seq: {aa_seq}, DNA Seq: {dna_seq}, Score: {score}\n")

def main(task_id, csv_file):
    """Main function to run codon optimization."""
    data = load_input_data(csv_file)
    for iteration, row in data.iterrows():
        aa_seq = row['AA_seq']
        dna_seq, score = optimize_codon(aa_seq)
        save_results(task_id, iteration + 1, aa_seq, dna_seq, score)

if __name__ == "__main__":
    task_id = sys.argv[1]
    csv_file = sys.argv[2]
    main(task_id, csv_file)
```

### 2. Shell Script: `batch_optimization.sh`

```bash
#!/bin/bash
#SBATCH --job-name=codon_opt
#SBATCH --output=codon_opt_%A_%a.out
#SBATCH --error=codon_opt_%A_%a.err
#SBATCH --array=1-$(wc -l < batch_optimization_tasks.csv)
#SBATCH --time=01:00:00
#SBATCH --mem=4G

# Load necessary modules if required
# module load python/3.x

# Get the task ID
TASK_ID=$SLURM_ARRAY_TASK_ID

# Define the input CSV file
CSV_FILE="/home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/batch_optimization_tasks.csv"

# Run the Python script
python3 /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/codon_optimization.py $TASK_ID $CSV_FILE
```

### 3. Task File: `batch_optimization_tasks.csv`

```csv
AA_seq
MKTIIALSYIFCLVFADYKDDDDK
MKTIIALSYIFCLVFADYKDDDDK
MKTIIALSYIFCLVFADYKDDDDK
```

### Instructions to Create the Files

You can create these files in the specified directory using the following commands in your terminal:

```bash
mkdir -p /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark

# Create the Python script
cat << 'EOF' > /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/codon_optimization.py
import pandas as pd
import numpy as np
import sys

def load_input_data(csv_file):
    """Load input data from a CSV file."""
    return pd.read_csv(csv_file)

def optimize_codon(sequence):
    """Perform codon optimization on the given sequence."""
    # Placeholder for actual optimization logic
    optimized_sequence = sequence  # Replace with actual optimization logic
    score = np.random.rand()  # Replace with actual scoring logic
    return optimized_sequence, score

def save_results(task_id, iteration, aa_seq, dna_seq, score):
    """Save the results of the optimization to a log file."""
    log_file = f"/home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/log_{task_id}.txt"
    with open(log_file, 'a') as f:
        f.write(f"Iteration: {iteration}, AA Seq: {aa_seq}, DNA Seq: {dna_seq}, Score: {score}\n")

def main(task_id, csv_file):
    """Main function to run codon optimization."""
    data = load_input_data(csv_file)
    for iteration, row in data.iterrows():
        aa_seq = row['AA_seq']
        dna_seq, score = optimize_codon(aa_seq)
        save_results(task_id, iteration + 1, aa_seq, dna_seq, score)

if __name__ == "__main__":
    task_id = sys.argv[1]
    csv_file = sys.argv[2]
    main(task_id, csv_file)
EOF

# Create the shell script
cat << 'EOF' > /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/batch_optimization.sh
#!/bin/bash
#SBATCH --job-name=codon_opt
#SBATCH --output=codon_opt_%A_%a.out
#SBATCH --error=codon_opt_%A_%a.err
#SBATCH --array=1-$(wc -l < batch_optimization_tasks.csv)
#SBATCH --time=01:00:00
#SBATCH --mem=4G

# Load necessary modules if required
# module load python/3.x

# Get the task ID
TASK_ID=$SLURM_ARRAY_TASK_ID

# Define the input CSV file
CSV_FILE="/home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/batch_optimization_tasks.csv"

# Run the Python script
python3 /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/codon_optimization.py $TASK_ID $CSV_FILE
EOF

# Create the CSV task file
cat << 'EOF' > /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/batch_optimization_tasks.csv
AA_seq
MKTIIALSYIFCLVFADYKDDDDK
MKTIIALSYIFCLVFADYKDDDDK
MKTIIALSYIFCLVFADYKDDDDK
EOF

# Make the shell script executable
chmod +x /home/wenzhao/github_repo/clone_repeat_protein/codon_opt_benchmark/batch_optimization.sh
```

### Notes:
- Make sure to replace the placeholder logic in the `optimize_codon` function with your actual optimization logic.
- Adjust the parameters in the SLURM script as needed for your specific computational environment.
- The CSV file currently contains three identical sequences for demonstration purposes; you can modify it with your actual sequences.