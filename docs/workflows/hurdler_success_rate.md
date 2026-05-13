# Workflow: HURDLER Success-Rate Analysis

> **Layout note (May 2026 reconstruction).**
> The pipeline now lives at [`../../src/hurdler/pipeline.py`](../../src/hurdler/pipeline.py).
> Inputs come from [`../../data/hurdler_analysis_input/`](../../data/hurdler_analysis_input/)
> (or set `HURDLER_INPUT_DIR`); outputs go to [`../../output/`](../../output/)
> (or set `HURDLER_OUTPUT_DIR`). Run from the repository root with
> `PYTHONPATH=src python -m hurdler.pipeline`.

# HURDLER Success Rate Analysis

Complete pipeline for analyzing HURDLER cloning strategy success rates across different plasmids and sequence lengths.

## Overview

This analysis calculates the probability of finding valid HURDLER enzyme combinations for random amino acid sequences of varying lengths (7-60 AA) across 8 common expression plasmids.

## Directory Structure

```
hurdler_analysis/
├── hurdler_analysis.py      # Main analysis script
├── README.md                 # This file
├── input/                    # Required input files
│   ├── methylation_check.csv
│   ├── neb_buffer_activity_cleaned.csv
│   ├── plasmid_digest_check.csv
│   ├── codon_usage.csv
│   ├── seamless_insert.csv
│   └── slient_mutation.csv
└── output/                   # Generated results (created automatically)
    ├── hurdler_three_site_combinations_df1.csv
    ├── hurdler_three_site_combinations_df2_optimized.csv
    ├── hurdler_lookup_optimized.pkl
    ├── hurdler_success_rate_results.csv
    ├── hurdler_success_rate_summary.csv
    ├── hurdler_success_rate_vs_sequence_length.pdf
    └── hurdler_success_rate_statistics.csv
```

## Requirements

- Python 3.7+
- Required packages:
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - biopython
  - tqdm

Install dependencies:
```bash
pip install pandas numpy matplotlib seaborn biopython tqdm
```

## Usage

Simply run the main script:

```bash
cd hurdler_analysis
python hurdler_analysis.py
```

The analysis will:
1. Generate three-site enzyme combinations (df1)
2. Expand with 3mer AA sequences and codon optimization (df2)
3. Create lookup dictionary with new structure
4. Calculate success rates for 7-60 AA sequences (1000 tests per length)
5. Generate visualizations and statistics

## Analysis Pipeline

### Step 1: Generate Three-Site Enzyme Combinations (df1)
- **Site I**: Seamless insert (regular enzymes, methylation compatible)
- **Site II**: Silent mutation (regular enzymes, methylation compatible)
- **Site III**: Type IIS enzymes (BsaI, BbsI, etc., no methylation requirement)

Quality filters applied:
- Commercial availability
- NEB availability
- Good ligation efficiency
- No star activity
- Plasmid compatibility
- Methylation compatibility (Sites I & II only)

### Step 2: Generate Optimized df2
Expands df1 by adding 3mer AA sequences:
- Filters out stop codons (*)
- Optimizes codon usage for E.coli
- Determines mutation direction (left/right search)
- Site I: Uses original sequence codon usage
- Site II: Uses mutated sequence codon usage

### Step 3: Create Lookup Dictionary
New structure: `{(3mer_i, 3mer_ii): [(site_i_info, site_ii_info, site_iii_enzyme, plasmid_dict), ...]}`

Where:
- **Key**: Tuple of two 3mer AA (no ordering)
- **Value**: List of tuples containing:
  - `site_i_info`: {'enzyme', '3mer_aa', '9mer_bp'}
  - `site_ii_info`: {'enzyme', '3mer_aa', '9mer_bp_original', '9mer_bp_mutated', 'direction'}
  - `site_iii_enzyme`: Enzyme name
  - `plasmid_dict`: {'plasmid_name': True/False, ...}

### Step 4: Calculate Success Rates
- Tests random AA sequences from 7-60 residues
- 1000 tests per length
- 8 plasmids tested
- Total: 432,000 tests
- Uses doubled sequence (ABCD → ABCDABCD) for 3mer extraction
- Distance constraint: 5 < distance < module_length

### Step 5: Visualization
Generates:
- Main plot: Success rate vs sequence length (PDF)
- Statistics table: Summary by plasmid

## Output Files

1. **hurdler_three_site_combinations_df1.csv**: Three-site enzyme combinations
2. **hurdler_three_site_combinations_df2_optimized.csv**: Expanded with 3mer AA (large file ~6GB)
3. **hurdler_lookup_optimized.pkl**: Fast lookup dictionary
4. **hurdler_success_rate_results.csv**: Raw test results (432,000 rows)
5. **hurdler_success_rate_summary.csv**: Summary statistics
6. **hurdler_success_rate_vs_sequence_length.pdf**: Visualization
7. **hurdler_success_rate_statistics.csv**: Plasmid comparison table

## Key Findings

- Success rates increase with sequence length
- Best performing plasmids at 60 AA:
  - pUC18: ~97-99%
  - pCold_I: ~97-99%
  - pQE-3: ~97-99%
- Lower performing plasmids:
  - pGEX-4T-1: ~40-50%
  - pMAL-c5X: ~40-50%

## Computation Time

- Total runtime: ~30-60 minutes (depending on system)
- Memory usage: ~8-10 GB peak
- Disk space: ~6-7 GB for df2 file

## Notes

- Lookup dictionary uses symmetric keys: both (A, B) and (B, A) point to the same entries
- Site orders in values maintain biological meaning: Site I → Site II → Site III
- All sequences tested are random 20-amino acid standard set (no stop codons)
- Random seed set to 42 for reproducibility
