#!/usr/bin/env bash
# Submit this file as one Digs/taskrunner command with 16 requested CPU cores.
# It rebuilds the historical notebook pattern population in memory, benchmarks
# serial versus 16-process execution, writes exactly two result Parquet files,
# validates every 7--60AA point against the committed notebook results in
# two-copy mode, or performs an exact two-copy/three-copy comparison in
# three-copy mode. Only the validated two-copy run removes superseded shards.

set -euo pipefail

repo=/home/wendai/projects/hurdler/clone_repeat_protein
python_bin=/home/wendai/.conda/envs/hurdler/bin/python
runner="$repo/scripts/run_success_landscape_single_files.py"
scan_copies="${1:-3}"
if [[ "$scan_copies" != 2 && "$scan_copies" != 3 ]]; then
    echo "Scan copies must be 2 or 3; received $scan_copies" >&2
    exit 2
fi

two_copy_root=/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run03_single_file_16core/raw
if [[ "$scan_copies" -eq 2 ]]; then
    output_root="$two_copy_root"
    figure_dir="$repo/studies/hurdler_validation/step02_success_landscape/figures"
else
    output_root=/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run06_three_copy_16core/raw
    figure_dir="$repo/studies/hurdler_validation/step02_success_landscape/figures/scan_3x"
fi
short_output="$output_root/short_motifs_1_5.parquet"
random_output="$output_root/random_modules_6_60.parquet"

requested_workers="${SLURM_CPUS_PER_TASK:-16}"
if [[ "$requested_workers" -ne 16 ]]; then
    echo "This production workflow requires exactly 16 allocated CPU cores; received $requested_workers" >&2
    exit 2
fi

"$python_bin" "$runner" \
    --repo-dir "$repo" \
    --benchmark \
    --scan-copies "$scan_copies" \
    --benchmark-workers 1,16

"$python_bin" "$runner" \
    --repo-dir "$repo" \
    --short-output "$short_output" \
    --random-output "$random_output" \
    --figure-dir "$figure_dir" \
    --workers 16 \
    --scan-copies "$scan_copies" \
    --short-max-length 5 \
    --random-min-length 6 \
    --random-max-length 60 \
    --tests 1000

# A second independent read proves that both completed files are intact before
# any legacy data is removed.
"$python_bin" "$runner" \
    --repo-dir "$repo" \
    --short-output "$short_output" \
    --random-output "$random_output" \
    --figure-dir "$figure_dir" \
    --workers 16 \
    --scan-copies "$scan_copies" \
    --short-max-length 5 \
    --random-min-length 6 \
    --random-max-length 60 \
    --tests 1000 \
    --validate-only

legacy_motif_shards=/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run01_production/raw/short_shards/motif_shards
legacy_hit_shards=/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/step02_success_landscape/runs/run01_production/raw/short_shards/short_motif_hits.parquet

if [[ "$scan_copies" -eq 2 ]]; then
    if [[ -d "$legacy_motif_shards" ]]; then
        rm -r -- "$legacy_motif_shards"
    fi
    if [[ -d "$legacy_hit_shards" ]]; then
        rm -r -- "$legacy_hit_shards"
    fi
    echo "Validated the two-copy Parquet files and removed any superseded legacy shard directories."
else
    "$python_bin" "$runner" \
        --compare-two-short "$two_copy_root/short_motifs_1_5.parquet" \
        --compare-two-random "$two_copy_root/random_modules_6_60.parquet" \
        --compare-three-short "$short_output" \
        --compare-three-random "$random_output" \
        --comparison-output "$repo/studies/hurdler_validation/step02_success_landscape/tables/scan_copy_comparison_2x_vs_3x.csv"
    echo "Validated the three-copy files and wrote the exact two-copy versus three-copy comparison."
fi
