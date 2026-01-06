#!/bin/bash
# Sequentially run M1..M15 codon optimization with increasing repeat counts
set -euo pipefail

# Activate conda environment
source /home/wenzhao/miniconda3/etc/profile.d/conda.sh
conda activate codon_opt

cd "$(dirname "$0")"

MOTIF="NEQIQAVIDAGALPALVQLLSSPNEQILQEALWALSNIASGG"
MUTATION_GRID="0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.0010,0.0011,0.0012,0.0013,0.0014,0.0015,0.0016,0.0017,0.0018,0.0019,0.0020,0.0021,0.0022,0.0023,0.0024,0.0025"

for i in $(seq 1 15); do
  name="M${i}"
  aa_seq="$(printf "%0.s${MOTIF}" $(seq 1 "$i"))"
  out_dir="./codon_opt_results_${name,,}"
  log_file="${out_dir}/${name,,}_run.log"

  mkdir -p "$out_dir"

  echo "[$(date +"%F %T")] Starting $name with ${i} motif repeats -> $out_dir"
  python codon_optimizer.py \
    --name "$name" \
    --aa-seq "$aa_seq" \
    --output-dir "$out_dir" \
    --mutation-grid "$MUTATION_GRID" \
    --verbose | tee "$log_file"
  echo "[$(date +"%F %T")] Finished $name" 
  echo "---------------------------------------------"
done
