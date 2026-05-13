#!/usr/bin/env python3
"""HURDLER data validation and statistics.

Sanity-checks the data products of :mod:`hurdler.pipeline` (or its
notebook equivalent) and prints a structured statistics report.

Verifies that:

- All reference inputs under ``data/reference_output/`` exist.
- ``df1`` contains the required columns and that Site II / Site III
  overhangs match the HURDLER constraint.
- ``df2`` contains the expected 3-mer AA / DNA columns and no
  unexpected NaNs.
- Per-plasmid compatibility counts are consistent.

Run from the repository root::

    PYTHONPATH=src python -m hurdler.validate
    PYTHONPATH=src python -m hurdler.validate --export output/validation.txt

See ``docs/workflows/hurdler_site_combinations.md`` for the broader
workflow.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys


def check_files_exist():
    """Verify the curated reference inputs required by the HURDLER pipeline.

    Returns ``True`` if every CSV under ``data/reference_output/`` that the
    pipeline depends on is present, ``False`` otherwise.
    """
    required_files = [
        "./data/reference_output/methylation_check.csv",
        "./data/reference_output/neb_buffer_activity_cleaned.csv",
        "./data/reference_output/plasmid_digest_check.csv",
        "./data/reference_output/restriction_enzyme_slient_mutation.csv",
        "./data/reference_output/restriction_enzyme_seamless_insert.csv",
        "./data/reference_output/orthogonality.csv",
    ]
    
    missing_files = []
    for file in required_files:
        if not Path(file).exists():
            missing_files.append(file)
    
    if missing_files:
        print("Error: Missing required input files:")
        for file in missing_files:
            print(f"  - {file}")
        return False
    
    return True


def validate_df1(df1_path):
    """Validate df1 structure and content"""
    print("\n=== Validating df1 ===")
    
    if not Path(df1_path).exists():
        print(f"Error: {df1_path} not found")
        return False
    
    df1 = pd.read_csv(df1_path)
    
    # Check required columns
    required_cols = ['site_i', 'site_ii', 'site_iii', 'ovhg_i', 'ovhg_ii', 'ovhg_iii']
    missing_cols = [col for col in required_cols if col not in df1.columns]
    
    if missing_cols:
        print(f"Error: Missing columns in df1: {', '.join(missing_cols)}")
        return False
    
    # Validate overhang patterns
    print(f"✓ Total combinations: {len(df1):,}")
    
    # Check Site II and III have same ovhg pattern
    same_ovhg = (df1['ovhg_ii'] == df1['ovhg_iii']).all()
    if same_ovhg:
        print("✓ All Site II and III have matching ovhg patterns")
    else:
        n_mismatch = (df1['ovhg_ii'] != df1['ovhg_iii']).sum()
        print(f"✗ Warning: {n_mismatch} combinations have mismatched Site II/III ovhg")
    
    # Check Site I and II have different ovhg OR are orthogonal
    same_ovhg_i_ii = (df1['ovhg_i'] == df1['ovhg_ii']).sum()
    diff_ovhg_i_ii = (df1['ovhg_i'] != df1['ovhg_ii']).sum()
    print(f"✓ Site I/II same ovhg (must be orthogonal): {same_ovhg_i_ii:,}")
    print(f"✓ Site I/II different ovhg: {diff_ovhg_i_ii:,}")
    
    # Check plasmid compatibility columns
    plasmid_cols = [col for col in df1.columns if col.endswith('_compatible')]
    if plasmid_cols:
        print(f"✓ Plasmid compatibility columns: {len(plasmid_cols)}")
        for col in plasmid_cols:
            n_compatible = df1[col].sum()
            print(f"  - {col.replace('_compatible', '')}: {n_compatible:,} compatible")
    
    return True


def validate_df2(df2_path):
    """Validate df2 structure and content"""
    print("\n=== Validating df2 ===")
    
    if not Path(df2_path).exists():
        print(f"Error: {df2_path} not found")
        return False
    
    df2 = pd.read_csv(df2_path)
    
    # Check required columns
    required_cols = [
        'site_i', 'site_i_3mer_aa', 'site_i_dna',
        'site_ii', 'site_ii_3mer_aa', 'site_ii_dna', 'site_ii_dna_mutated',
        'site_iii'
    ]
    missing_cols = [col for col in required_cols if col not in df2.columns]
    
    if missing_cols:
        print(f"Error: Missing columns in df2: {', '.join(missing_cols)}")
        return False
    
    print(f"✓ Total site-3mer AA combinations: {len(df2):,}")
    print(f"✓ Unique Site I 3mer AA: {df2['site_i_3mer_aa'].nunique():,}")
    print(f"✓ Unique Site II 3mer AA: {df2['site_ii_3mer_aa'].nunique():,}")
    print(f"✓ Unique (Site I AA, Site II AA) pairs: {df2[['site_i_3mer_aa', 'site_ii_3mer_aa']].drop_duplicates().shape[0]:,}")
    
    # Check for NaN values in critical columns
    nan_counts = df2[required_cols].isna().sum()
    if nan_counts.sum() > 0:
        print("\n✗ Warning: Found NaN values:")
        for col, count in nan_counts[nan_counts > 0].items():
            print(f"  - {col}: {count} NaN values")
    else:
        print("✓ No NaN values in critical columns")
    
    # Validate 3mer AA sequences (should be 3 characters)
    invalid_site_i = df2[df2['site_i_3mer_aa'].str.len() != 3]
    invalid_site_ii = df2[df2['site_ii_3mer_aa'].str.len() != 3]
    
    if len(invalid_site_i) > 0:
        print(f"✗ Warning: {len(invalid_site_i)} invalid Site I 3mer AA (not 3 characters)")
    else:
        print("✓ All Site I 3mer AA are 3 characters")
    
    if len(invalid_site_ii) > 0:
        print(f"✗ Warning: {len(invalid_site_ii)} invalid Site II 3mer AA (not 3 characters)")
    else:
        print("✓ All Site II 3mer AA are 3 characters")
    
    return True


def generate_statistics(df1_path, df2_path):
    """Generate detailed statistics"""
    print("\n" + "="*80)
    print("DETAILED STATISTICS")
    print("="*80)
    
    df1 = pd.read_csv(df1_path)
    df2 = pd.read_csv(df2_path)
    
    # Enzyme usage statistics
    print("\n1. ENZYME USAGE")
    print(f"   Site I unique enzymes: {df1['site_i'].nunique()}")
    print(f"   Site II unique enzymes: {df1['site_ii'].nunique()}")
    print(f"   Site III unique enzymes: {df1['site_iii'].nunique()}")
    
    print("\n   Top 5 most used Site I enzymes:")
    for enzyme, count in df1['site_i'].value_counts().head(5).items():
        print(f"     {enzyme}: {count:,} combinations ({count/len(df1)*100:.1f}%)")
    
    print("\n   Top 5 most used Site II enzymes:")
    for enzyme, count in df1['site_ii'].value_counts().head(5).items():
        print(f"     {enzyme}: {count:,} combinations ({count/len(df1)*100:.1f}%)")
    
    print("\n   Top 5 most used Site III enzymes:")
    for enzyme, count in df1['site_iii'].value_counts().head(5).items():
        print(f"     {enzyme}: {count:,} combinations ({count/len(df1)*100:.1f}%)")
    
    # Overhang pattern statistics
    print("\n2. OVERHANG PATTERNS")
    print("   Site I overhang distribution:")
    for ovhg, count in df1['ovhg_i'].value_counts().sort_index().items():
        print(f"     {ovhg:>3} bp: {count:,} ({count/len(df1)*100:.1f}%)")
    
    print("\n   Site II/III overhang distribution:")
    for ovhg, count in df1['ovhg_ii'].value_counts().sort_index().items():
        print(f"     {ovhg:>3} bp: {count:,} ({count/len(df1)*100:.1f}%)")
    
    # 3mer AA coverage
    print("\n3. 3MER AA COVERAGE")
    
    # Calculate theoretical maximum
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    theoretical_max = len(amino_acids) ** 3
    
    site_i_coverage = df2['site_i_3mer_aa'].nunique() / theoretical_max * 100
    site_ii_coverage = df2['site_ii_3mer_aa'].nunique() / theoretical_max * 100
    
    print(f"   Theoretical maximum 3mer AA combinations: {theoretical_max:,}")
    print(f"   Site I 3mer AA coverage: {site_i_coverage:.2f}%")
    print(f"   Site II 3mer AA coverage: {site_ii_coverage:.2f}%")
    
    print("\n   Top 10 most common Site I 3mer AA:")
    for aa, count in df2['site_i_3mer_aa'].value_counts().head(10).items():
        print(f"     {aa}: {count:,} combinations")
    
    print("\n   Top 10 most common Site II 3mer AA:")
    for aa, count in df2['site_ii_3mer_aa'].value_counts().head(10).items():
        print(f"     {aa}: {count:,} combinations")
    
    # Plasmid compatibility matrix
    print("\n4. PLASMID COMPATIBILITY MATRIX")
    plasmid_cols = [col for col in df1.columns if col.endswith('_compatible')]
    
    if plasmid_cols:
        print("\n   Combinations per plasmid:")
        for col in sorted(plasmid_cols):
            plasmid = col.replace('_compatible', '')
            n_compatible = df1[col].sum()
            pct = n_compatible / len(df1) * 100
            print(f"     {plasmid:<25}: {n_compatible:>8,} ({pct:>5.1f}%)")
        
        # Cross-plasmid compatibility
        print("\n   Multi-plasmid compatibility:")
        compatibility_counts = df1[[col for col in plasmid_cols]].sum(axis=1)
        for n in sorted(compatibility_counts.unique()):
            count = (compatibility_counts == n).sum()
            print(f"     Compatible with {n} plasmid(s): {count:,}")


def export_summary(df1_path, df2_path, output_path):
    """Export summary to text file"""
    import sys
    from io import StringIO
    
    # Redirect stdout to capture print output
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    
    # Generate all output
    generate_statistics(df1_path, df2_path)
    
    # Get the output
    output = mystdout.getvalue()
    sys.stdout = old_stdout
    
    # Write to file
    with open(output_path, 'w') as f:
        f.write(output)
    
    print(f"\nSummary exported to: {output_path}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Validate HURDLER data and generate statistics')
    parser.add_argument('--df1', default='./output/hurdler_three_site_combinations_df1.csv',
                        help='Path to df1 CSV file')
    parser.add_argument('--df2', default='./output/hurdler_three_site_combinations_df2.csv',
                        help='Path to df2 CSV file')
    parser.add_argument('--export', type=str, help='Export summary to text file')
    
    args = parser.parse_args()
    
    print("="*80)
    print("HURDLER DATA VALIDATION AND STATISTICS")
    print("="*80)
    
    # Check input files
    print("\n=== Checking Required Input Files ===")
    if not check_files_exist():
        print("\nPlease ensure all required input files are present.")
        sys.exit(1)
    print("✓ All required input files found")
    
    # Validate df1
    if not validate_df1(args.df1):
        sys.exit(1)
    
    # Validate df2
    if not validate_df2(args.df2):
        sys.exit(1)
    
    # Generate statistics
    generate_statistics(args.df1, args.df2)
    
    # Export if requested
    if args.export:
        export_summary(args.df1, args.df2, args.export)
    
    print("\n" + "="*80)
    print("Validation complete!")
    print("="*80)


if __name__ == '__main__':
    main()
