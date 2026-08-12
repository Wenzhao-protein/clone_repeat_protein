#!/usr/bin/env bash
set -euo pipefail

input_catalog="$1"
output_catalog="$2"
output_directory="$(dirname "$output_catalog")"
workers="${3:-12}"

/home/wendai/.conda/envs/hurdler/bin/python \
  /home/wendai/projects/hurdler/clone_repeat_protein/studies/hurdler_validation/scripts/apply_natural_middle_units.py \
  --input "$input_catalog" \
  --output "$output_catalog" \
  --workers "$workers" \
  --target-count 100

/home/wendai/.conda/envs/hurdler/bin/hurdler infer-boundaries \
  --input "$output_catalog" \
  --output "$output_directory/natural100_full_sequence_scan_audit.parquet" \
  --candidates-output "$output_directory/natural100_full_sequence_period_candidates.parquet" \
  --units-output "$output_directory/natural100_full_sequence_scan_units.parquet" \
  --positions-output "$output_directory/natural100_full_sequence_scan_positions.parquet" \
  --workers "$workers"
