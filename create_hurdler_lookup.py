#!/usr/bin/env python3
"""
Create optimized HURDLER lookup dictionary directly from df1.
This skips the df2 generation and creates the lookup dictionary directly.
"""

import pandas as pd
from pathlib import Path
from tqdm import tqdm
import pickle
import sys
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import numpy as np

def process_chunk(args):
    """
    Process a chunk of df1 rows in parallel.
    
    Parameters:
    -----------
    args : tuple
        (chunk_df, site_i_map, site_ii_map, plasmid_cols)
    
    Returns:
    --------
    dict : lookup dictionary for this chunk
    """
    chunk_df, site_i_map, site_ii_map, plasmid_cols = args
    
    chunk_lookup = defaultdict(list)
    
    for _, row in chunk_df.iterrows():
        site_i_enzyme = row['site_i_enzyme']
        site_ii_enzyme = row['site_ii_enzyme']
        site_iii_enzyme = row['site_iii_enzyme']
        
        # Get plasmid compatibility
        plasmid_compat = {col: row[col] for col in plasmid_cols}
        
        # Get all 3mer options for each site
        site_i_options = site_i_map.get(site_i_enzyme, [])
        site_ii_options = site_ii_map.get(site_ii_enzyme, [])
        
        # For each pair of Site I and Site II 3mers
        for i_3mer, i_dna in site_i_options:
            for ii_3mer, ii_dna in site_ii_options:
                # Create key
                key = frozenset([i_3mer, ii_3mer])
                
                # Add entry
                chunk_lookup[key].append([
                    site_i_enzyme,
                    site_ii_enzyme,
                    site_iii_enzyme,
                    i_3mer,
                    ii_3mer,
                    i_dna,
                    ii_dna,
                    plasmid_compat
                ])
    
    return dict(chunk_lookup)


def merge_lookups(lookup_list):
    """
    Merge multiple lookup dictionaries from parallel processing.
    
    Parameters:
    -----------
    lookup_list : list of dict
        List of lookup dictionaries from each process
    
    Returns:
    --------
    dict : merged lookup dictionary
    """
    merged = defaultdict(list)
    
    for lookup in lookup_list:
        for key, entries in lookup.items():
            merged[key].extend(entries)
    
    return dict(merged)


def create_lookup_from_df1(df1_path, df_silent_mutation, df_seamless_insert, output_dir='./output', n_processes=None):
    """
    Create optimized lookup dictionary directly from df1.
    
    Lookup structure:
    - Key: frozenset of two 3mer AA sequences (unordered pair)
    - Value: list of [site_i_enzyme, site_ii_enzyme, site_iii_enzyme, 
                      site_i_3mer_aa, site_ii_3mer_aa, 
                      site_i_dna (9mer), site_ii_dna (9mer),
                      plasmid_compatibility_dict]
    
    Parameters:
    -----------
    n_processes : int, optional
        Number of processes to use. If None, uses all available CPUs.
    """
    print("\n" + "="*80)
    print("CREATING OPTIMIZED LOOKUP DICTIONARY (MULTIPROCESSING)")
    print("="*80)
    
    # Determine number of processes
    if n_processes is None:
        n_processes = cpu_count()
    print(f"\nUsing {n_processes} CPU cores")
    
    # Load df1
    print(f"\nLoading df1 from {df1_path}...")
    df1 = pd.read_csv(df1_path)
    print(f"Loaded {len(df1):,} three-site combinations")
    
    # Create enzyme to 3mer mappings
    print("\n1. Creating enzyme-to-3mer mappings...")
    
    # Site I: seamless insert
    site_i_map = {}
    for _, row in df_seamless_insert.iterrows():
        enzyme = row['name']
        three_mer_aa = row['re_site_shifted_tl']
        dna_seq = row['re_site_shifted']
        
        if enzyme not in site_i_map:
            site_i_map[enzyme] = []
        site_i_map[enzyme].append((three_mer_aa, dna_seq))
    
    print(f"   Site I: {len(site_i_map)} enzymes")
    
    # Site II/III: silent mutation
    site_ii_map = {}
    for _, row in df_silent_mutation.iterrows():
        enzyme = row['name']
        three_mer_aa = row['re_site_shifted_tl']
        dna_seq = row['re_site_mutate_shifted']
        
        if enzyme not in site_ii_map:
            site_ii_map[enzyme] = []
        site_ii_map[enzyme].append((three_mer_aa, dna_seq))
    
    print(f"   Site II/III: {len(site_ii_map)} enzymes")
    
    # Build lookup dictionary using multiprocessing
    print(f"\n2. Building optimized lookup dictionary (multiprocessing with {n_processes} cores)...")
    
    plasmid_cols = [col for col in df1.columns if '_compatible' in col]
    
    # Split df1 into chunks for parallel processing
    chunk_size = len(df1) // n_processes + 1
    chunks = [df1.iloc[i:i+chunk_size].copy() for i in range(0, len(df1), chunk_size)]
    print(f"   Split into {len(chunks)} chunks of ~{chunk_size} rows each")
    
    # Prepare arguments for each process
    chunk_args = [(chunk, site_i_map, site_ii_map, plasmid_cols) for chunk in chunks]
    
    # Process chunks in parallel
    print("   Processing chunks in parallel...")
    with Pool(processes=n_processes) as pool:
        chunk_lookups = list(tqdm(
            pool.imap(process_chunk, chunk_args),
            total=len(chunk_args),
            desc="   Chunks processed"
        ))
    
    # Merge results
    print("   Merging results from all processes...")
    lookup_dict = merge_lookups(chunk_lookups)
    
    print(f"\n3. Created lookup with {len(lookup_dict):,} unique 3mer AA pairs")
    
    # Calculate statistics
    total_entries = sum(len(v) for v in lookup_dict.values())
    avg_entries = total_entries / len(lookup_dict) if lookup_dict else 0
    print(f"   Total entries: {total_entries:,}")
    print(f"   Average entries per pair: {avg_entries:.1f}")
    
    # Save full lookup dictionary
    lookup_path = Path(output_dir) / 'hurdler_3mer_lookup.pkl'
    with open(lookup_path, 'wb') as f:
        pickle.dump(lookup_dict, f)
    print(f"\n✓ Full lookup dictionary saved to: {lookup_path}")
    
    # Create lightweight lookup for success rate testing
    # Only store: key (3mer AA pair) -> set of compatible plasmids
    print("\n4. Creating lightweight lookup for success rate testing...")
    lightweight_lookup = {}
    
    for key, entries in tqdm(lookup_dict.items(), desc="   Processing pairs"):
        # For each 3mer AA pair, collect all compatible plasmids
        compatible_plasmids = {}
        
        for entry in entries:
            plasmid_compat = entry[7]  # 8th element is plasmid dict
            for plasmid_col, is_compatible in plasmid_compat.items():
                if is_compatible:
                    plasmid_name = plasmid_col.replace('_compatible', '')
                    if plasmid_name not in compatible_plasmids:
                        compatible_plasmids[plasmid_name] = True
        
        # Only store if at least one plasmid is compatible
        if compatible_plasmids:
            lightweight_lookup[key] = compatible_plasmids
    
    print(f"   Lightweight lookup: {len(lightweight_lookup):,} 3mer pairs with plasmid info")
    
    # Save lightweight lookup
    lightweight_path = Path(output_dir) / 'hurdler_3mer_lightweight_lookup.pkl'
    with open(lightweight_path, 'wb') as f:
        pickle.dump(lightweight_lookup, f)
    print(f"\n✓ Lightweight lookup saved to: {lightweight_path}")
    print(f"   (This file is much smaller and faster to load for testing)")
    
    # Also save a summary CSV for inspection
    summary_data = []
    for key, plasmid_dict in list(lightweight_lookup.items())[:100]:  # First 100 for inspection
        key_list = list(key)
        summary_data.append({
            '3mer_aa_1': key_list[0] if len(key_list) > 0 else '',
            '3mer_aa_2': key_list[1] if len(key_list) > 1 else key_list[0],
            'n_compatible_plasmids': len(plasmid_dict),
            'compatible_plasmids': ','.join(sorted(plasmid_dict.keys()))
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_path = Path(output_dir) / 'hurdler_lookup_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"✓ Summary saved to: {summary_path}")
    
    return lookup_dict, lightweight_lookup


def main():
    """Main execution function"""
    print("="*80)
    print("HURDLER OPTIMIZED LOOKUP CREATION")
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
    print(f"✓ Found: df1")
    
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
    
    # Create lookup
    lookup_dict, lightweight_lookup = create_lookup_from_df1(
        df1_path, df_silent_mutation, df_seamless_insert, output_dir, n_processes=None
    )
    
    print("\n" + "="*80)
    print("LOOKUP CREATION COMPLETE")
    print("="*80)
    print(f"\nFiles created:")
    print(f"  1. hurdler_3mer_lookup.pkl - Full lookup (for detailed queries)")
    print(f"  2. hurdler_3mer_lightweight_lookup.pkl - Lightweight (for success rate testing)")
    print(f"  3. hurdler_lookup_summary.csv - Summary for inspection")
    print(f"\nYou can now run hurdler_success_rate_test.py")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Create optimized HURDLER lookup dictionary with multiprocessing'
    )
    parser.add_argument('--n-processes', type=int, default=None,
                        help='Number of processes to use (default: all CPUs)')
    
    args = parser.parse_args()
    
    # Update the function call to pass n_processes
    original_main = main
    def main_with_args():
        """Main execution function with command line args"""
        print("="*80)
        print("HURDLER OPTIMIZED LOOKUP CREATION")
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
        print(f"✓ Found: df1")
        
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
        
        # Create lookup with specified number of processes
        lookup_dict, lightweight_lookup = create_lookup_from_df1(
            df1_path, df_silent_mutation, df_seamless_insert, output_dir, n_processes=args.n_processes
        )
        
        print("\n" + "="*80)
        print("LOOKUP CREATION COMPLETE")
        print("="*80)
        print(f"\nFiles created:")
        print(f"  1. hurdler_3mer_lookup.pkl - Full lookup (for detailed queries)")
        print(f"  2. hurdler_3mer_lightweight_lookup.pkl - Lightweight (for success rate testing)")
        print(f"  3. hurdler_lookup_summary.csv - Summary for inspection")
        print(f"\nYou can now run hurdler_success_rate_test.py")
    
    main_with_args()
