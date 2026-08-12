#!/usr/bin/env python3
"""
Quick test script to verify HURDLER analysis system is working correctly.
Run this after generating data to ensure everything is set up properly.
"""

import sys
from pathlib import Path
import pandas as pd

# Historical executable smoke script, not a pytest module.
__test__ = False


def test_file_exists(filepath, description):
    """Test if a file exists"""
    if Path(filepath).exists():
        print(f"✅ {description}: {filepath}")
        return True
    else:
        print(f"❌ {description} NOT FOUND: {filepath}")
        return False


def test_dataframe_structure(df, expected_cols, df_name):
    """Test if DataFrame has expected columns"""
    missing = [col for col in expected_cols if col not in df.columns]
    if missing:
        print(f"❌ {df_name} missing columns: {', '.join(missing)}")
        return False
    else:
        print(f"✅ {df_name} has all expected columns")
        return True


def main():
    print("="*80)
    print("HURDLER ANALYSIS SYSTEM - QUICK TEST")
    print("="*80)
    
    all_passed = True
    
    # Test 1: Check input files
    print("\n[TEST 1] Checking input files...")
    input_files = [
        ("./utils/output/methylation_check.csv", "Methylation check data"),
        ("./utils/output/neb_buffer_activity_cleaned.csv", "NEB activity data"),
        ("./utils/output/plasmid_digest_check.csv", "Plasmid compatibility"),
        ("./utils/output/restriction_enzyme_slient_mutation.csv", "Silent mutation data"),
        ("./utils/output/restriction_enzyme_seamless_insert.csv", "Seamless insert data"),
        ("./utils/output/orthogonality.csv", "Orthogonality matrix"),
    ]
    
    for filepath, desc in input_files:
        if not test_file_exists(filepath, desc):
            all_passed = False
    
    # Test 2: Check output files
    print("\n[TEST 2] Checking output files...")
    output_files = [
        ("./output/hurdler_three_site_combinations_df1.csv", "df1 (site combinations)"),
        ("./output/hurdler_three_site_combinations_df2.csv", "df2 (with 3mer AA)"),
        ("./output/hurdler_3mer_aa_lookup.csv", "Lookup table"),
    ]
    
    output_exists = True
    for filepath, desc in output_files:
        if not test_file_exists(filepath, desc):
            all_passed = False
            output_exists = False
    
    if not output_exists:
        print("\n⚠️  Output files not found. Please run hurdler_site_combination_analysis.ipynb first.")
        sys.exit(1)
    
    # Test 3: Validate df1 structure
    print("\n[TEST 3] Validating df1 structure...")
    try:
        df1 = pd.read_csv("./output/hurdler_three_site_combinations_df1.csv")
        print(f"✅ df1 loaded: {len(df1):,} rows, {len(df1.columns)} columns")
        
        expected_cols = ['site_i', 'site_ii', 'site_iii', 'ovhg_i', 'ovhg_ii', 'ovhg_iii']
        if not test_dataframe_structure(df1, expected_cols, "df1"):
            all_passed = False
        
        # Check ovhg patterns
        site_ii_iii_match = (df1['ovhg_ii'] == df1['ovhg_iii']).all()
        if site_ii_iii_match:
            print("✅ All Site II and III have matching ovhg patterns")
        else:
            print("❌ Some Site II/III ovhg patterns don't match")
            all_passed = False
        
        # Check plasmid columns
        plasmid_cols = [col for col in df1.columns if col.endswith('_compatible')]
        if len(plasmid_cols) > 0:
            print(f"✅ Found {len(plasmid_cols)} plasmid compatibility columns")
        else:
            print("❌ No plasmid compatibility columns found")
            all_passed = False
            
    except Exception as e:
        print(f"❌ Error loading df1: {e}")
        all_passed = False
    
    # Test 4: Validate df2 structure
    print("\n[TEST 4] Validating df2 structure...")
    try:
        df2 = pd.read_csv("./output/hurdler_three_site_combinations_df2.csv")
        print(f"✅ df2 loaded: {len(df2):,} rows, {len(df2.columns)} columns")
        
        expected_cols = [
            'site_i', 'site_i_3mer_aa', 'site_i_dna',
            'site_ii', 'site_ii_3mer_aa', 'site_ii_dna', 'site_ii_dna_mutated',
            'site_iii'
        ]
        if not test_dataframe_structure(df2, expected_cols, "df2"):
            all_passed = False
        
        # Check for NaN in critical columns
        critical_cols = ['site_i_3mer_aa', 'site_ii_3mer_aa']
        nan_count = df2[critical_cols].isna().sum().sum()
        if nan_count == 0:
            print("✅ No NaN values in 3mer AA columns")
        else:
            print(f"❌ Found {nan_count} NaN values in 3mer AA columns")
            all_passed = False
        
        # Check 3mer AA lengths
        site_i_valid = (df2['site_i_3mer_aa'].str.len() == 3).all()
        site_ii_valid = (df2['site_ii_3mer_aa'].str.len() == 3).all()
        
        if site_i_valid and site_ii_valid:
            print("✅ All 3mer AA sequences are 3 characters long")
        else:
            if not site_i_valid:
                print("❌ Some Site I 3mer AA are not 3 characters")
            if not site_ii_valid:
                print("❌ Some Site II 3mer AA are not 3 characters")
            all_passed = False
        
        # Report statistics
        print(f"\n📊 Statistics:")
        print(f"   Unique Site I 3mer AA: {df2['site_i_3mer_aa'].nunique():,}")
        print(f"   Unique Site II 3mer AA: {df2['site_ii_3mer_aa'].nunique():,}")
        print(f"   Unique (Site I, Site II) pairs: {df2[['site_i_3mer_aa', 'site_ii_3mer_aa']].drop_duplicates().shape[0]:,}")
        
    except Exception as e:
        print(f"❌ Error loading df2: {e}")
        all_passed = False
    
    # Test 5: Test query function
    print("\n[TEST 5] Testing query function...")
    try:
        # Get a sample 3mer AA pair
        sample_site_i_aa = df2['site_i_3mer_aa'].value_counts().index[0]
        sample_site_ii_aa = df2['site_ii_3mer_aa'].value_counts().index[0]
        
        # Find compatible plasmid
        plasmid_cols = [col for col in df2.columns if col.endswith('_compatible')]
        test_plasmid = None
        for col in plasmid_cols:
            test_plasmid = col.replace('_compatible', '')
            break
        
        if test_plasmid:
            # Query
            results = df2[
                (df2['site_i_3mer_aa'] == sample_site_i_aa) &
                (df2['site_ii_3mer_aa'] == sample_site_ii_aa) &
                (df2[f'{test_plasmid}_compatible'] == True)
            ]
            
            if len(results) > 0:
                print(f"✅ Query function works!")
                print(f"   Test query: Site I AA='{sample_site_i_aa}', Site II AA='{sample_site_ii_aa}', Plasmid='{test_plasmid}'")
                print(f"   Found {len(results)} valid combination(s)")
            else:
                print("⚠️  Query returned no results (this may be normal)")
        else:
            print("❌ No plasmid columns found for testing")
            all_passed = False
            
    except Exception as e:
        print(f"❌ Error testing query function: {e}")
        all_passed = False
    
    # Test 6: Check scripts are executable
    print("\n[TEST 6] Checking scripts...")
    scripts = [
        ("./hurdler_query.py", "Query tool"),
        ("./hurdler_validate.py", "Validation tool"),
    ]
    
    for filepath, desc in scripts:
        if test_file_exists(filepath, desc):
            if Path(filepath).stat().st_mode & 0o111:
                print(f"   ✅ {desc} is executable")
            else:
                print(f"   ⚠️  {desc} is not executable (run: chmod +x {filepath})")
    
    # Final summary
    print("\n" + "="*80)
    if all_passed:
        print("✅ ALL TESTS PASSED!")
        print("\nYour HURDLER analysis system is ready to use.")
        print("\nNext steps:")
        print("  1. Try a query: python hurdler_query.py --list")
        print("  2. Run validation: python hurdler_validate.py")
        print("  3. Read the docs: HURDLER_QUICKSTART.md")
    else:
        print("❌ SOME TESTS FAILED")
        print("\nPlease check the errors above and:")
        print("  1. Make sure you've run hurdler_site_combination_analysis.ipynb")
        print("  2. Check that all input files are in utils/output/")
        print("  3. Review the error messages for specific issues")
        sys.exit(1)
    print("="*80)


if __name__ == '__main__':
    main()
