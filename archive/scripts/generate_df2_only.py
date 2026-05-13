#!/usr/bin/env python3
"""
Generate df2 from existing df1 by adding 3mer AA sequences.
"""

import pandas as pd
from pathlib import Path
from tqdm import tqdm
import sys

def generate_df2_from_df1(df1_path, df_silent_mutation, df_seamless_insert, output_dir='./output'):
    """
    Generate df2: add 3mer AA sequences to df1
    """
    print("\n" + "="*80)
    print("GENERATING DF2: Adding 3mer AA Sequences")
    print("="*80)
    
    # Load df1
    print(f"\nLoading df1 from {df1_path}...")
    df1 = pd.read_csv(df1_path)
    print(f"Loaded {len(df1):,} three-site combinations")
    
    # Create lookup dictionaries
    print("\n1. Creating lookup dictionaries...")
    
    # Site I: seamless insert (enzyme -> 3mer AA mapping)
    site_i_lookup = {}
    for _, row in df_seamless_insert.iterrows():
        enzyme = row['name']
        three_mer_aa = row['re_site_shifted_tl']
        dna_seq = row['re_site_shifted']
        
        if enzyme not in site_i_lookup:
            site_i_lookup[enzyme] = []
        site_i_lookup[enzyme].append({
            '3mer_aa': three_mer_aa,
            'dna_sequence': dna_seq,
            'has_mutation': False  # Seamless insert doesn't have mutation
        })
    
    print(f"   Site I lookup: {len(site_i_lookup)} enzymes")
    
    # Site II/III: silent mutation (enzyme -> 3mer AA mapping)
    site_ii_lookup = {}
    for _, row in df_silent_mutation.iterrows():
        enzyme = row['name']
        three_mer_aa = row['re_site_shifted_tl']
        dna_seq = row['re_site_mutate_shifted']
        
        if enzyme not in site_ii_lookup:
            site_ii_lookup[enzyme] = []
        site_ii_lookup[enzyme].append({
            '3mer_aa': three_mer_aa,
            'dna_sequence': dna_seq,
            'has_mutation': True  # Silent mutation has mutation
        })
    
    print(f"   Site II/III lookup: {len(site_ii_lookup)} enzymes")
    
    # Expand df1 to df2 using vectorized operations
    print("\n2. Expanding combinations with 3mer AA sequences...")
    expanded_rows = []
    
    # Pre-compute for efficiency
    df1_values = df1.values
    df1_cols = df1.columns.tolist()
    site_i_idx = df1_cols.index('site_i_enzyme')
    site_ii_idx = df1_cols.index('site_ii_enzyme')
    site_iii_idx = df1_cols.index('site_iii_enzyme')
    
    batch_size = 1000
    for batch_start in tqdm(range(0, len(df1), batch_size), desc="Processing batches"):
        batch_end = min(batch_start + batch_size, len(df1))
        
        for idx in range(batch_start, batch_end):
            row_values = df1_values[idx]
            site_i = row_values[site_i_idx]
            site_ii = row_values[site_ii_idx]
            site_iii = row_values[site_iii_idx]
            
            # Get possible 3mer AA for each site
            site_i_options = site_i_lookup.get(site_i, [])
            site_ii_options = site_ii_lookup.get(site_ii, [])
            site_iii_options = site_ii_lookup.get(site_iii, [])  # Site III uses silent mutation
            
            if not site_i_options or not site_ii_options or not site_iii_options:
                continue
            
            # Generate all combinations for this row
            for i_option in site_i_options:
                for ii_option in site_ii_options:
                    for iii_option in site_iii_options:
                        # Build row as list for efficiency
                        expanded_row = list(row_values) + [
                            i_option['3mer_aa'],
                            i_option['dna_sequence'],
                            i_option['has_mutation'],
                            ii_option['3mer_aa'],
                            ii_option['dna_sequence'],
                            ii_option['has_mutation'],
                            iii_option['3mer_aa'],
                            iii_option['dna_sequence'],
                            iii_option['has_mutation']
                        ]
                        expanded_rows.append(expanded_row)
    
    # Create DataFrame with proper columns
    new_columns = df1_cols + [
        'site_i_3mer_aa', 'site_i_dna', 'site_i_has_mutation',
        'site_ii_3mer_aa', 'site_ii_dna', 'site_ii_has_mutation',
        'site_iii_3mer_aa', 'site_iii_dna', 'site_iii_has_mutation'
    ]
    df2 = pd.DataFrame(expanded_rows, columns=new_columns)
    print(f"\n3. Generated {len(df2):,} combinations with 3mer AA sequences")
    
    # Save df2
    output_path = Path(output_dir) / 'hurdler_three_site_combinations_df2.csv'
    df2.to_csv(output_path, index=False)
    print(f"\n✓ df2 saved to: {output_path}")
    
    # Create optimized lookup dictionary
    print("\n4. Creating optimized lookup dictionary...")
    lookup_dict = {}
    
    for _, row in tqdm(df2.iterrows(), total=len(df2), desc="Building lookup"):
        site_i_aa = row['site_i_3mer_aa']
        site_ii_aa = row['site_ii_3mer_aa']
        
        # Create key as frozenset (unordered pair)
        key = frozenset([site_i_aa, site_ii_aa])
        
        if key not in lookup_dict:
            lookup_dict[key] = []
        
        # Store complete information
        lookup_dict[key].append([
            row['site_i_enzyme'],
            row['site_ii_enzyme'],
            row['site_iii_enzyme'],
            row['site_i_3mer_aa'],
            row['site_ii_3mer_aa'],
            row['site_i_dna'],
            row['site_ii_dna'],
            # Store plasmid compatibility info
            {col: row[col] for col in row.index if '_compatible' in col}
        ])
    
    print(f"   Created lookup with {len(lookup_dict):,} unique 3mer AA pairs")
    
    # Save lookup dictionary as pickle for fast loading
    import pickle
    lookup_path = Path(output_dir) / 'hurdler_3mer_lookup.pkl'
    with open(lookup_path, 'wb') as f:
        pickle.dump(lookup_dict, f)
    print(f"\n✓ Lookup dictionary saved to: {lookup_path}")
    
    return df2, lookup_dict


def main():
    """Main execution function"""
    print("="*80)
    print("HURDLER DF2 GENERATION FROM DF1")
    print("="*80)
    
    # Check input files
    input_dir = Path('./utils/output')
    output_dir = Path('./output')
    df1_path = output_dir / 'hurdler_three_site_combinations_df1.csv'
    
    if not df1_path.exists():
        print(f"\nError: df1 not found at {df1_path}")
        print("Please run generate_hurdler_data.py first to generate df1.")
        sys.exit(1)
    
    print("\nChecking input files...")
    
    required_files = [
        'restriction_enzyme_slient_mutation.csv',
        'restriction_enzyme_seamless_insert.csv'
    ]
    
    for filename in required_files:
        filepath = input_dir / filename
        if not filepath.exists():
            print(f"✗ Missing: {filepath}")
            sys.exit(1)
        print(f"✓ Found: {filename}")
    
    # Load data
    print("\n" + "="*80)
    print("LOADING DATA")
    print("="*80)
    
    df_silent_mutation = pd.read_csv(input_dir / 'restriction_enzyme_slient_mutation.csv')
    print(f"✓ Silent mutation options: {df_silent_mutation.shape}")
    
    df_seamless_insert = pd.read_csv(input_dir / 'restriction_enzyme_seamless_insert.csv')
    print(f"✓ Seamless insert options: {df_seamless_insert.shape}")
    
    # Generate df2
    df2, lookup_dict = generate_df2_from_df1(df1_path, df_silent_mutation, df_seamless_insert, output_dir)
    
    print("\n" + "="*80)
    print("DF2 GENERATION COMPLETE")
    print("="*80)
    print(f"\nGenerated:")
    print(f"  df2: {len(df2):,} combinations with 3mer AA sequences")
    print(f"  Lookup: {len(lookup_dict):,} unique 3mer AA pairs")
    print(f"\nFiles saved to: {output_dir}/")


if __name__ == '__main__':
    main()
