# legacy-optimized-v1 complete lookup

This is the complete, immutable lookup bundled for the interactive designer:

- 1,335,463 distinct protein-pattern keys;
- 6,549,904 normalized plasmid-specific enzyme-pair candidates;
- all eight maintained plasmids;
- every retained Site-I/Site-II pair and its complete Site-III alternatives.

`pattern_index.npz` is the compact query accelerator. The partitioned Parquet
catalog materializes every retained candidate so the notebook can display and
filter all routes rather than only the best pair. `metadata.json` freezes the
rule profile and upstream input hashes. `SHA256SUMS` verifies every binary.

Do not regenerate this directory during an interactive run. Build a separately
versioned artifact with `hurdler lookup build` when the scientific rule profile
or reference tables change.
