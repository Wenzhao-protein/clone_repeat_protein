# HURDLER Success Rate Analysis - Quick Start Guide

## 🚀 Quick Start

### Prerequisites
1. Ensure all data files are in place:
   ```
   ./utils/output/methylation_check.csv
   ./utils/output/neb_buffer_activity_cleaned.csv
   ./utils/output/plasmid_digest_check.csv
   ./utils/output/seamless_insert.csv
   ./utils/output/slient_mutation.csv
   ./utils/output/codon_usage.csv
   ./utils/input/golden_gate_fidelity/FileS03_T4_18h_25C.xlsx
   ```

2. Required Python packages:
   ```bash
   pip install pandas numpy matplotlib seaborn biopython tqdm openpyxl
   ```

### Running the Analysis

#### Option 1: Run All Cells (Recommended for first time)
```python
# In Jupyter/VS Code:
1. Open hurdler_success_rate_analysis.ipynb
2. Select "Run All" or execute cells sequentially
3. Total runtime: ~20-40 minutes
```

#### Option 2: Run Specific Steps
```python
# After Step 1-3 are completed once, you can re-run Step 4 with different parameters:

# In Step 4 cell, modify:
MATCHING_MODE = 'any'  # or 'greedy'
num_tests_per_length = 1000  # increase for more statistical power

# Then run Step 4-7 only
```

## 📊 Key Outputs

### Must-Check Files
1. **RE Pairing Matrices** (Step 1.5)
   - `./output/site_i_site_ii_pairing_matrix.csv`
   - `./output/site_i_site_ii_pairing_heatmap.pdf` ⭐
   - `./output/site_ii_site_iii_pairing_matrix.csv`
   - `./output/site_ii_site_iii_pairing_heatmap.pdf` ⭐

2. **Lookup Dictionary** (Step 3)
   - `./output/hurdler_lookup_pattern_based.pkl`

3. **Success Rate Results** (Step 4)
   - `./output/hurdler_success_rate_7_60aa_pattern.csv` ⭐
   - `./output/hurdler_success_rate_7_60aa_pattern.pdf` ⭐

## 🔍 Understanding the Results

### RE Pairing Matrix Interpretation
```python
# Load and inspect:
import pandas as pd
matrix = pd.read_csv('./output/site_i_site_ii_pairing_matrix.csv', index_col=0)

# True (1) = Can pair together
# False (0) = Cannot pair (not orthogonal enough)
```

### Success Rate Analysis
```python
# Load results:
df = pd.read_csv('./output/hurdler_success_rate_7_60aa_pattern.csv')

# Key columns:
# - module_length: Sequence length in amino acids
# - success_rate: Fraction of sequences with valid HURDLER sites (0-1)
# - num_successes: Number of successful tests
# - total_combinations: Total valid enzyme combinations found

# Quick stats:
print(f"Mean success rate: {df['success_rate'].mean():.1%}")
print(f"Success rate at 20 AA: {df[df['module_length']==20]['success_rate'].values[0]:.1%}")
print(f"Success rate at 60 AA: {df[df['module_length']==60]['success_rate'].values[0]:.1%}")
```

## ⚙️ Configuration Options

### Matching Mode (Step 4)
```python
# Greedy mode: Count all possible combinations (slower, more information)
MATCHING_MODE = 'greedy'

# Any mode: Stop at first match (faster, for existence checking)
MATCHING_MODE = 'any'
```

### Statistical Power (Step 4)
```python
# Fewer tests = faster but less reliable
num_tests_per_length = 100

# More tests = slower but more reliable
num_tests_per_length = 1000  # default

# High precision (warning: very slow)
num_tests_per_length = 10000
```

### Memory Management (Step 2 & 3)
```python
# Smaller chunks = less memory but slower
chunk_size = 100000

# Larger chunks = more memory but faster
chunk_size = 500000  # default

# For systems with 32GB+ RAM
chunk_size = 1000000
```

## 🐛 Troubleshooting

### Issue: Out of Memory
**Solution 1**: Reduce chunk size
```python
# In Step 2 and Step 3, change:
chunk_size = 100000  # or even 50000
```

**Solution 2**: Process in stages
```python
# Run Steps 1-2, restart kernel, then run Step 3-4
```

### Issue: Fidelity File Not Found
```python
# Error: FileNotFoundError: golden_gate_fidelity/FileS03_T4_18h_25C.xlsx

# Solution: Check file path
import os
assert os.path.exists('./utils/input/golden_gate_fidelity/FileS03_T4_18h_25C.xlsx')
```

### Issue: No Matches Found in Step 4
**Possible causes**:
1. Lookup dictionary is empty (check Step 3 output)
2. Distance constraints too strict
3. RE pairing too restrictive

**Debug**:
```python
# Check lookup size:
import pickle
with open('./output/hurdler_lookup_pattern_based.pkl', 'rb') as f:
    lookup = pickle.load(f)
print(f"Lookup patterns: {len(lookup)}")
print(f"Total entries: {sum(len(v) for v in lookup.values())}")
```

### Issue: Slow Performance
**Speed up tips**:
1. Use `MATCHING_MODE = 'any'` instead of 'greedy'
2. Reduce `num_tests_per_length` to 100 for quick testing
3. Use fewer module lengths:
   ```python
   module_lengths = list(range(10, 61, 5))  # Every 5 instead of every 1
   ```

## 📈 Interpreting Success Rates

### What the Numbers Mean
- **0%**: No sequences of this length have valid HURDLER sites
- **50%**: Half of random sequences have valid sites
- **100%**: All sequences have valid sites

### Expected Patterns
- Success rate increases with sequence length
- Longer sequences have more 3mer pairs to search
- Typical range: 20-90% for 10-60 AA sequences

### Red Flags
- Success rate decreases with length → Check lookup dictionary
- All success rates are 0% → RE pairing may be too restrictive
- All success rates are 100% → Constraints may be too loose

## 🎯 Next Steps

### After Successful Run
1. **Visualize Results**: Open the PDF in `./output/`
2. **Analyze Patterns**: Look at which sequence lengths have best success rates
3. **Check RE Pairing**: Inspect heatmaps to understand enzyme compatibility
4. **Optimize Parameters**: Try different matching modes and test counts

### Advanced Usage
1. **Custom Amino Acid Sets**:
   ```python
   # In Step 4, modify:
   amino_acids = ['A', 'C', 'D', 'E']  # Only test specific AAs
   ```

2. **Specific Module Lengths**:
   ```python
   # In Step 4, modify:
   module_lengths = [20, 30, 40]  # Test only these lengths
   ```

3. **Export for External Analysis**:
   ```python
   import pickle
   import pandas as pd
   
   # Export lookup to JSON for inspection
   with open('./output/hurdler_lookup_pattern_based.pkl', 'rb') as f:
       lookup = pickle.load(f)
   
   # Sample 10 patterns
   sample = dict(list(lookup.items())[:10])
   import json
   with open('./output/lookup_sample.json', 'w') as f:
       json.dump(sample, f, indent=2)
   ```

## 📚 Additional Resources

### Key Documentation Files
- `HURDLER_REWRITE_SUMMARY.md` - Detailed technical changes
- `FILE_INDEX.md` - Complete file listing
- `HURDLER_QUICKSTART.md` - Original quickstart (may be outdated)

### Related Notebooks
- `utils/re_pair_fidelity.ipynb` - Orthogonality calculation details
- `hurdler_query.py` - Query existing HURDLER database

### Support
For issues or questions:
1. Check the troubleshooting section above
2. Review the full summary in `HURDLER_REWRITE_SUMMARY.md`
3. Inspect intermediate outputs in `./output/` directory

---

**Last Updated**: 2026-01-11  
**Version**: 2.0 (Pattern-Based)
