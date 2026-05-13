# Workflow: HURDLER Three-Site Combinations

> **Layout note (May 2026 reconstruction).**
> Canonical code now lives in [`../../src/hurdler/`](../../src/hurdler/):
> `pipeline.py`, `query.py`, `validate.py`, `utils.py`. The reference
> inputs that used to live under `utils/output/` and `hurdler_analysis/input/`
> are now under [`../../data/reference_output/`](../../data/reference_output/)
> and [`../../data/hurdler_analysis_input/`](../../data/hurdler_analysis_input/).
> Generated artifacts go to [`../../output/`](../../output/) (gitignored).
> Earlier paths mentioned below have been updated where possible; if you
> see a stale reference, the canonical location is in the table above.

# HURDLER Three-Site Combination Analysis

This toolkit identifies and queries valid three-site restriction enzyme combinations for the HURDLER (Hidden Unremovable RE sites in Disposable Linker for Efficient Repeat protein cloning) method.

## Overview

The HURDLER method requires three restriction enzyme sites with specific properties:
- **Site I**: Seamless insertion, methylation insensitive
- **Site II**: Allows silent mutation, methylation insensitive, same overhang as Site III
- **Site III**: Type IIS enzyme (cuts outside recognition sequence), same overhang as Site II

## Files

### Main Analysis Notebook
- **hurdler_site_combination_analysis.ipynb**: Main analysis notebook that generates all data

### Scripts
- **hurdler_query.py**: Command-line tool for querying valid combinations
- **hurdler_validate.py**: Data validation and statistics generator

### Output Files
- **hurdler_three_site_combinations_df1.csv**: All valid site combinations with plasmid compatibility
- **hurdler_three_site_combinations_df2.csv**: Site combinations with 3mer AA sequences
- **hurdler_3mer_aa_lookup.csv**: Quick lookup table for 3mer AA queries
- **hurdler_plasmid_statistics.png**: Plasmid compatibility visualization
- **hurdler_overhang_distribution.png**: Overhang pattern distributions
- **hurdler_top_enzymes.png**: Most frequently used enzymes

## Requirements

```bash
pip install pandas numpy biopython matplotlib seaborn tqdm openpyxl
```

## Usage

### 1. Generate Data (First Time)

Run the Jupyter notebook to generate all combinations:

```bash
jupyter notebook hurdler_site_combination_analysis.ipynb
```

Execute all cells. This will:
1. Load and validate input data from `utils/output/`
2. Filter enzymes based on HURDLER requirements
3. Generate all valid three-site combinations (df1)
4. Map 3mer AA sequences to combinations (df2)
5. Create lookup tables and visualizations

**Estimated runtime**: 10-30 minutes depending on system

### 2. Query Valid Combinations

#### Single Query

Query specific 3mer AA sequences for a plasmid:

```bash
python hurdler_query.py --site-i-aa "ABC" --site-ii-aa "DEF" --plasmid "pET-28a(+)"
```

Example output:
```
Found 5 valid combination(s):
  Site I 3mer AA:  ABC
  Site II 3mer AA: DEF
  Plasmid:         pET-28a(+)

Combination 1:
  Site I:   EcoRI           (ovhg:  -4, frame: 2)
            DNA: GAATTC
  Site II:  BamHI           (ovhg:  -4, frame: 1, silent mutation)
            DNA: GGATCC
            Mutated: GGATCT
  Site III: BsaI            (ovhg:  -4, Type IIS)
```

#### Batch Query

Create an input CSV file (`queries.csv`):
```csv
site_i_3mer_aa,site_ii_3mer_aa,plasmid
ABC,DEF,pET-28a(+)
GHI,JKL,pGEX-4T-1
MNO,PQR,pMAL-c5X
```

Run batch query:
```bash
python hurdler_query.py --batch queries.csv --output results.csv
```

#### List Available Options

```bash
python hurdler_query.py --list
```

### 3. Validate Data

Run validation checks and generate statistics:

```bash
python hurdler_validate.py
```

Export summary to file:
```bash
python hurdler_validate.py --export validation_summary.txt
```

### 4. Use in Python Scripts

```python
import pandas as pd

# Load data
df2 = pd.read_csv('./output/hurdler_three_site_combinations_df2.csv')

# Query function
def find_hurdler_sites(site_i_3mer_aa, site_ii_3mer_aa, plasmid):
    plasmid_col = f'{plasmid}_compatible'
    results = df2[
        (df2['site_i_3mer_aa'] == site_i_3mer_aa) &
        (df2['site_ii_3mer_aa'] == site_ii_3mer_aa) &
        (df2[plasmid_col] == True)
    ]
    return results

# Example usage
combinations = find_hurdler_sites('ABC', 'DEF', 'pET-28a(+)')
print(f"Found {len(combinations)} valid combinations")

# Get enzyme details
for _, row in combinations.iterrows():
    print(f"Site I: {row['site_i']}, Site II: {row['site_ii']}, Site III: {row['site_iii']}")
```

## Data Structure

### df1 Structure (hurdler_three_site_combinations_df1.csv)

| Column | Description |
|--------|-------------|
| site_i | Site I enzyme name |
| site_ii | Site II enzyme name |
| site_iii | Site III enzyme name |
| ovhg_i | Site I overhang length (signed) |
| ovhg_ii | Site II overhang length (signed) |
| ovhg_iii | Site III overhang length (signed) |
| pGEX-4T-1_compatible | Boolean: compatible with pGEX-4T-1 |
| pMAL-c5X_compatible | Boolean: compatible with pMAL-c5X |
| ... | (one column per plasmid) |
| compatible_plasmids | List of compatible plasmids |
| num_compatible_plasmids | Number of compatible plasmids |

### df2 Structure (hurdler_three_site_combinations_df2.csv)

| Column | Description |
|--------|-------------|
| site_i | Site I enzyme name |
| site_i_3mer_aa | 3-amino acid sequence for Site I |
| site_i_dna | DNA sequence encoding the 3mer AA |
| site_i_frame | Frame shift value |
| site_i_codon_usage | Codon usage frequency |
| site_ii | Site II enzyme name |
| site_ii_3mer_aa | 3-amino acid sequence for Site II |
| site_ii_dna | DNA sequence encoding the 3mer AA |
| site_ii_dna_mutated | Mutated DNA (inactivates RE site) |
| site_ii_frame | Frame shift value |
| site_ii_codon_usage | Codon usage frequency |
| site_iii | Site III enzyme name |
| ovhg_i, ovhg_ii, ovhg_iii | Overhang lengths |
| (plasmid compatibility columns) | Boolean values |

## Supported Plasmids

1. pGEX-4T-1 (GST fusion, BamHI/EagI)
2. pMAL-c5X (MBP fusion, NdeI/HindIII)
3. pET-21a(+) (T7 expression, NdeI/XhoI)
4. pET-28a(+) (His-tag, BamHI/XhoI)
5. pET-28a(+)_start_codon (His-tag with start, NcoI/XhoI)
6. pCold_I (Cold shock, NdeI/BspMI)
7. pUC18 (Cloning vector, EcoRI/HindIII)
8. pQE-3 (His-tag, BamHI/HindIII)

## Selection Criteria Details

### All Sites (I, II, III)
- No ambiguous bases in recognition sequence
- Overhang length: 2-5 bp (Golden Gate compatible)
- No star activity
- Commercially available

### Site I Requirements
- Not sensitive to 6mA/5mC methylation (DH5α compatible)
- Allows seamless insertion
- Does not cut plasmid backbone

### Site II Requirements
- Not sensitive to 6mA/5mC methylation
- Allows silent mutation to inactivate site
- Same overhang pattern as Site III
- Does not cut plasmid backbone
- Different overhang from Site I OR orthogonal (non-complementary)

### Site III Requirements
- Type IIS enzyme (cuts outside recognition sequence)
- Same overhang pattern as Site II

## Troubleshooting

### No results found

1. Check if 3mer AA sequences exist in the database:
   ```bash
   python hurdler_query.py --list
   ```

2. Try different plasmids - not all combinations work with all plasmids

3. Check that input files are present in `utils/output/`:
   - methylation_check.csv
   - neb_buffer_activity_cleaned.csv
   - plasmid_digest_check.csv
   - restriction_enzyme_slient_mutation.csv
   - restriction_enzyme_seamless_insert.csv
   - orthogonality.csv

### Data validation errors

Run validation script:
```bash
python hurdler_validate.py
```

If files are missing, re-run the analysis notebook from the beginning.

## Citation

If you use this toolkit, please cite:
[Add your publication information here]

## Contact

For questions or issues, please contact:
[Add contact information]

## License

[Add license information]
