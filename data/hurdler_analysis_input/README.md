# data/hurdler_analysis_input

Inputs consumed directly by `src/hurdler/pipeline.py` when run with the
default `HURDLER_INPUT_DIR`. The contents mirror the subset of
`../reference_output/` that the pipeline reads:

- `methylation_check.csv`
- `neb_buffer_activity_cleaned.csv`
- `plasmid_digest_check.csv`
- `seamless_insert.csv`
- `slient_mutation.csv`
- `codon_usage.csv`

Refresh this folder if the pipeline changes the schema or column names
it expects.
