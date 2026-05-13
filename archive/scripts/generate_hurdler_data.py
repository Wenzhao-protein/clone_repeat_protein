#!/usr/bin/env python3
"""
Generate HURDLER three-site combination data directly.
This is a standalone version that doesn't require running the notebook.
"""

import pandas as pd
import numpy as np
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes
import itertools
from tqdm import tqdm
import warnings
import sys
from pathlib import Path

warnings.filterwarnings('ignore')


def check_methylation_compatible(enzyme_name, df_methylation):
    """Check if enzyme is not sensitive to 6mA/5mC methylation (DH5α compatible)"""
    dh5a_compatible = set(
        df_methylation.loc[~df_methylation['6mA_5mC_sensitive'], 'enzyme'].dropna().tolist() +
        df_methylation.loc[~df_methylation['6mA_5mC_sensitive'], 'prototype'].dropna().tolist()
    )
    return enzyme_name in dh5a_compatible


def check_neb_quality(enzyme_name, df_neb):
    """Check if enzyme has good NEB characteristics"""
    if enzyme_name not in df_neb['enzyme'].values:
        return False
    
    enzyme_data = df_neb[df_neb['enzyme'] == enzyme_name].iloc[0]
    
    # Check ligation efficiency is not low
    if enzyme_data['ligation_efficiencies'] == 'low':
        return False
    
    # Check no star activity
    if enzyme_data['star_activity'] == True:
        return False
    
    return True


def is_type_iis(enzyme_name):
    """Check if enzyme is Type IIS (cuts outside recognition sequence)
    
    In BioPython:
    - fst5: cut position on 5' strand (positive = from 5' end)
    - fst3: cut position on 3' strand (negative = from 3' end back)
    
    Type IIS: fst5 < 0 (left cut) OR fst5 > site_length (right cut)
    """
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        site_length = len(enzyme.site)
        # Type IIS: cut position is outside the recognition sequence
        return enzyme.fst5 < 0 or enzyme.fst5 > site_length
    except:
        return False


def is_regular_enzyme(enzyme_name):
    """Check if enzyme is regular (cuts inside recognition sequence)
    
    Regular enzymes: fst5 <= site_length (cuts within recognition sequence)
    """
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        site_length = len(enzyme.site)
        # Regular enzymes: 5' cut position is within or at the recognition sequence
        return enzyme.fst5 <= site_length
    except:
        return False


def get_ovhg_pattern(enzyme_name):
    """Get overhang pattern (length and direction)"""
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        return enzyme.ovhg  # Returns overhang length with sign
    except:
        return None


def has_degenerate_bases(enzyme_name):
    """Check if enzyme recognition site has degenerate bases"""
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        site = str(enzyme.site)
        # Standard DNA bases
        standard_bases = set('ATCG')
        # Check if all bases are standard
        return not all(base in standard_bases for base in site)
    except:
        return True  # If can't check, assume it has


def check_ovhg_length(enzyme_name):
    """Check if overhang length is 2-5 bp (Golden Gate compatible)"""
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        ovhg_len = abs(enzyme.ovhg) if enzyme.ovhg else 0
        return 2 <= ovhg_len <= 5
    except:
        return False


def check_plasmid_compatible(enzyme_name, df_plasmid_check):
    """Check if enzyme doesn't cut plasmid backbone"""
    if enzyme_name not in df_plasmid_check.index:
        return False
    # At least one plasmid should be compatible
    return df_plasmid_check.loc[enzyme_name].any()


def check_orthogonality(enzyme1, enzyme2, df_orthogonality, min_score=2):
    """Check if two enzymes are orthogonal (non-complementary sticky ends)
    
    For Type IIS enzymes not in orthogonality database:
    - If either enzyme is not in the database, assume orthogonal (return True)
    - Type IIS enzymes typically have non-complementary overhangs
    """
    # Look up orthogonality score
    result = df_orthogonality[
        ((df_orthogonality['re1'] == enzyme1) & (df_orthogonality['re2'] == enzyme2)) |
        ((df_orthogonality['re1'] == enzyme2) & (df_orthogonality['re2'] == enzyme1))
    ]
    
    if len(result) > 0:
        return result.iloc[0]['orthogonality'] >= min_score
    
    # If not in database (e.g., Type IIS enzymes), assume orthogonal
    # This is reasonable since Type IIS enzymes have flexible overhangs
    return True


def generate_df1(df_methylation, df_neb, df_plasmid_check, df_silent_mutation, 
                  df_seamless_insert, df_orthogonality, output_dir='./output'):
    """
    Generate df1: all valid three-site combinations with plasmid compatibility
    """
    print("\n" + "="*80)
    print("GENERATING DF1: Three-Site Combinations")
    print("="*80)
    
    # Get candidate enzymes for each site
    print("\n1. Filtering Site I candidates (seamless insert, regular enzyme)...")
    site_i_candidates = df_seamless_insert['name'].unique()
    site_i_candidates = [e for e in site_i_candidates if check_methylation_compatible(e, df_methylation)]
    site_i_candidates = [e for e in site_i_candidates if is_regular_enzyme(e)]
    print(f"   Found {len(site_i_candidates)} Site I candidates (regular enzyme)")
    
    print("\n2. Filtering Site II candidates (silent mutation, regular enzyme)...")
    site_ii_candidates = df_silent_mutation['name'].unique()
    site_ii_candidates = [e for e in site_ii_candidates if check_methylation_compatible(e, df_methylation)]
    site_ii_candidates = [e for e in site_ii_candidates if is_regular_enzyme(e)]
    print(f"   Found {len(site_ii_candidates)} Site II candidates (regular enzyme)")
    
    print("\n3. Filtering Site III candidates (Type IIS only, no methylation filter)...")
    # Site III: MUST be Type IIS enzymes (cuts outside recognition site)
    # Requirements: commercial availability, no degenerate bases, 2-5bp overhang,
    # plasmid compatible, NEB quality, NO methylation compatibility required
    from Bio import Restriction
    site_iii_candidates = [str(enz) for enz in Restriction.AllEnzymes]
    
    # Initial count
    print(f"   Starting with {len(site_iii_candidates)} enzymes from AllEnzymes")
    
    # Filter 1: Valid ovhang information
    site_iii_candidates = [e for e in site_iii_candidates if get_ovhg_pattern(e) is not None]
    print(f"   After ovhang filter: {len(site_iii_candidates)} enzymes")
    
    # Filter 2: No degenerate bases in recognition site
    site_iii_candidates = [e for e in site_iii_candidates if not has_degenerate_bases(e)]
    print(f"   After degenerate bases filter: {len(site_iii_candidates)} enzymes")
    
    # Filter 3: Overhang length 2-5 bp (Golden Gate compatible)
    site_iii_candidates = [e for e in site_iii_candidates if check_ovhg_length(e)]
    print(f"   After ovhang length (2-5bp) filter: {len(site_iii_candidates)} enzymes")
    
    # Filter 4: NEB quality (good ligation, no star activity)
    site_iii_candidates = [e for e in site_iii_candidates if check_neb_quality(e, df_neb)]
    print(f"   After NEB quality filter: {len(site_iii_candidates)} enzymes")
    
    # Filter 5: Plasmid compatibility (at least one plasmid)
    site_iii_candidates = [e for e in site_iii_candidates if check_plasmid_compatible(e, df_plasmid_check)]
    print(f"   After plasmid compatibility filter: {len(site_iii_candidates)} enzymes")
    
    # Filter 6: MUST be Type IIS (cuts outside recognition site)
    site_iii_candidates = [e for e in site_iii_candidates if is_type_iis(e)]
    print(f"   After Type IIS filter: {len(site_iii_candidates)} enzymes")
    print(f"   → Final Site III candidates (all Type IIS): {site_iii_candidates}")
    
    # Generate all combinations
    print("\n4. Generating all three-site combinations...")
    combinations = []
    
    for site_i in tqdm(site_i_candidates, desc="Site I"):
        ovhg_i = get_ovhg_pattern(site_i)
        if ovhg_i is None:
            continue
        
        for site_ii in site_ii_candidates:
            ovhg_ii = get_ovhg_pattern(site_ii)
            if ovhg_ii is None:
                continue
            
            # Check Site I and Site II orthogonality
            if not check_orthogonality(site_i, site_ii, df_orthogonality):
                continue
            
            for site_iii in site_iii_candidates:
                # Site III must be different from Site II to avoid overlapping cut sites
                if site_iii == site_ii:
                    continue
                
                # Site III must have same overhang as Site II
                ovhg_iii = get_ovhg_pattern(site_iii)
                if ovhg_iii != ovhg_ii:
                    continue
                
                # Check Site I and Site III orthogonality
                if not check_orthogonality(site_i, site_iii, df_orthogonality):
                    continue
                
                # Don't need to check Site II and Site III orthogonality since they have same overhang
                
                # Check plasmid compatibility for all three sites
                plasmid_compat = {}
                for plasmid in df_plasmid_check.columns:
                    all_compatible = (
                        df_plasmid_check.loc[site_i, plasmid] if site_i in df_plasmid_check.index else False
                    ) and (
                        df_plasmid_check.loc[site_ii, plasmid] if site_ii in df_plasmid_check.index else False
                    ) and (
                        df_plasmid_check.loc[site_iii, plasmid] if site_iii in df_plasmid_check.index else False
                    )
                    plasmid_compat[f'{plasmid}_compatible'] = all_compatible
                
                # Add combination
                combination = {
                    'site_i_enzyme': site_i,
                    'site_ii_enzyme': site_ii,
                    'site_iii_enzyme': site_iii,
                    'site_i_ovhg': ovhg_i,
                    'site_ii_ovhg': ovhg_ii,
                    'site_iii_ovhg': ovhg_iii,
                    **plasmid_compat
                }
                combinations.append(combination)
    
    df1 = pd.DataFrame(combinations)
    print(f"\n5. Generated {len(df1):,} valid three-site combinations")
    
    # Save df1
    output_path = Path(output_dir) / 'hurdler_three_site_combinations_df1.csv'
    df1.to_csv(output_path, index=False)
    print(f"\n✓ df1 saved to: {output_path}")
    
    return df1


def generate_df2(df1, df_silent_mutation, df_seamless_insert, output_dir='./output'):
    """
    Generate df2: add 3mer AA sequences to df1
    """
    print("\n" + "="*80)
    print("GENERATING DF2: Adding 3mer AA Sequences")
    print("="*80)
    
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
    
    # Site II: silent mutation (enzyme -> 3mer AA mapping)
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
    
    # Expand df1 to df2
    print("\n2. Expanding combinations with 3mer AA sequences...")
    print("   Note: Site III (Type IIS) doesn't require 3mer AA encoding")
    print("         Recognition site is outside the repeat unit")
    expanded_rows = []
    
    for _, row in tqdm(df1.iterrows(), total=len(df1), desc="Processing"):
        site_i = row['site_i_enzyme']
        site_ii = row['site_ii_enzyme']
        site_iii = row['site_iii_enzyme']
        
        # Get possible 3mer AA for Site I and Site II
        site_i_options = site_i_lookup.get(site_i, [])
        site_ii_options = site_ii_lookup.get(site_ii, [])
        
        # Site III is Type IIS - doesn't need 3mer AA
        # Its recognition site is outside the repeat unit
        
        # Generate all combinations
        for i_option in site_i_options:
            for ii_option in site_ii_options:
                expanded_row = row.to_dict()
                expanded_row.update({
                    'site_i_3mer_aa': i_option['3mer_aa'],
                    'site_i_dna': i_option['dna_sequence'],
                    'site_i_has_mutation': i_option['has_mutation'],
                    'site_ii_3mer_aa': ii_option['3mer_aa'],
                    'site_ii_dna': ii_option['dna_sequence'],
                    'site_ii_has_mutation': ii_option['has_mutation'],
                    # Site III info: only enzyme name and overhang (already in row)
                    'site_iii_3mer_aa': 'N/A',  # Type IIS doesn't encode AA
                    'site_iii_dna': 'N/A',  # Recognition site is outside repeat unit
                    'site_iii_has_mutation': False  # Not applicable
                })
                expanded_rows.append(expanded_row)
    
    df2 = pd.DataFrame(expanded_rows)
    print(f"\n3. Generated {len(df2):,} combinations with 3mer AA sequences")
    
    # Save df2
    output_path = Path(output_dir) / 'hurdler_three_site_combinations_df2.csv'
    df2.to_csv(output_path, index=False)
    print(f"\n✓ df2 saved to: {output_path}")
    
    return df2


def main():
    """Main execution function"""
    print("="*80)
    print("HURDLER THREE-SITE COMBINATION DATA GENERATION")
    print("="*80)
    
    # Check input files
    input_dir = Path('./utils/output')
    output_dir = Path('./output')
    output_dir.mkdir(exist_ok=True)
    
    required_files = [
        'methylation_check.csv',
        'neb_buffer_activity_cleaned.csv',
        'plasmid_digest_check.csv',
        'restriction_enzyme_slient_mutation.csv',
        'restriction_enzyme_seamless_insert.csv',
        'orthogonality.csv'
    ]
    
    print("\nChecking input files...")
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
    
    df_methylation = pd.read_csv(input_dir / 'methylation_check.csv')
    print(f"✓ Methylation check: {df_methylation.shape}")
    
    df_neb = pd.read_csv(input_dir / 'neb_buffer_activity_cleaned.csv')
    print(f"✓ NEB characteristics: {df_neb.shape}")
    
    df_plasmid_check = pd.read_csv(input_dir / 'plasmid_digest_check.csv', index_col=0)
    print(f"✓ Plasmid compatibility: {df_plasmid_check.shape}")
    
    df_silent_mutation = pd.read_csv(input_dir / 'restriction_enzyme_slient_mutation.csv')
    print(f"✓ Silent mutation options: {df_silent_mutation.shape}")
    
    df_seamless_insert = pd.read_csv(input_dir / 'restriction_enzyme_seamless_insert.csv')
    print(f"✓ Seamless insert options: {df_seamless_insert.shape}")
    
    df_orthogonality = pd.read_csv(input_dir / 'orthogonality.csv')
    print(f"✓ Orthogonality matrix: {df_orthogonality.shape}")
    
    # Generate df1
    df1 = generate_df1(
        df_methylation, df_neb, df_plasmid_check, df_silent_mutation,
        df_seamless_insert, df_orthogonality, output_dir
    )
    
    # Generate df2
    df2 = generate_df2(df1, df_silent_mutation, df_seamless_insert, output_dir)
    
    print("\n" + "="*80)
    print("DATA GENERATION COMPLETE")
    print("="*80)
    print(f"\nGenerated:")
    print(f"  df1: {len(df1):,} three-site combinations")
    print(f"  df2: {len(df2):,} combinations with 3mer AA sequences")
    print(f"\nFiles saved to: {output_dir}/")


if __name__ == '__main__':
    main()
