#!/usr/bin/env python
"""
Test new HURDLER data generation logic
"""
import pandas as pd
import numpy as np
import time
import sys
import Bio.Restriction as Restriction

print("="*80)
print("TESTING NEW HURDLER DATA GENERATION LOGIC")
print("="*80)

# ============================================================================
# Load base data
# ============================================================================
print("\n1. Loading base data...")
df_seamless_insert = pd.read_csv('./utils/output/restriction_enzyme_seamless_insert.csv')
df_silent_mutation = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')
df_plasmid_check = pd.read_csv('./utils/output/plasmid_digest_check.csv', index_col=0)
df_neb = pd.read_csv('./utils/output/neb_buffer_activity_cleaned.csv')

print(f"   ✓ Seamless insert: {len(df_seamless_insert)} rows, {df_seamless_insert['name'].nunique()} enzymes")
print(f"   ✓ Silent mutation: {len(df_silent_mutation)} rows, {df_silent_mutation['name'].nunique()} enzymes")
print(f"   ✓ Plasmid check: {df_plasmid_check.shape}")
print(f"   ✓ NEB data: {len(df_neb)} enzymes")

# ============================================================================
# Generate Site III data
# ============================================================================
print("\n2. Generating Site III enzyme data...")
site_iii_data = []

for enzyme_name in Restriction.AllEnzymes:
    enzyme_name = str(enzyme_name)
    try:
        enzyme = getattr(Restriction, enzyme_name)
        
        # Must be Type IIS
        site_length = len(enzyme.site)
        if not (enzyme.fst5 < 0 or enzyme.fst5 > site_length):
            continue
        
        # Must have ovhg in {-4, +2}
        if enzyme.ovhg not in [-4, 2]:
            continue
        
        # No degenerate bases
        site = str(enzyme.site)
        if not all(base in 'ATCG' for base in site):
            continue
        
        # Check NEB quality
        if enzyme_name in df_neb['enzyme'].values:
            enzyme_data = df_neb[df_neb['enzyme'] == enzyme_name].iloc[0]
            if enzyme_data['ligation_efficiencies'] == 'low':
                continue
            if enzyme_data['star_activity'] == True:
                continue
        else:
            continue
        
        # Check plasmid compatibility
        if enzyme_name in df_plasmid_check.index:
            if not df_plasmid_check.loc[enzyme_name].any():
                continue
        else:
            continue
        
        site_iii_data.append({
            'enzyme': enzyme_name,
            'ovhg': enzyme.ovhg,
            'recognition_site': str(enzyme.site)
        })
    except:
        continue

df_site_iii = pd.DataFrame(site_iii_data)
print(f"   ✓ Site III: {len(df_site_iii)} enzymes")
print(f"     Ovhg -4: {len(df_site_iii[df_site_iii['ovhg'] == -4])}")
print(f"     Ovhg +2: {len(df_site_iii[df_site_iii['ovhg'] == 2])}")

# ============================================================================
# Build Site I DataFrame
# ============================================================================
print("\n3. Building Site I DataFrame...")
start_time = time.time()
site_i_rows = []

for enzyme_name in df_seamless_insert['name'].unique():
    enzyme_data = df_seamless_insert[df_seamless_insert['name'] == enzyme_name]
    first_row = enzyme_data.iloc[0]
    enzyme_site = first_row.get('site', '')
    enzyme_ovhg = first_row['ovhg']
    
    # Get plasmid compatibility
    plasmid_compat = {}
    if enzyme_name in df_plasmid_check.index:
        for plasmid in df_plasmid_check.columns:
            plasmid_compat[f'site_i_{plasmid}_compatible'] = bool(df_plasmid_check.loc[enzyme_name, plasmid])
    else:
        for plasmid in df_plasmid_check.columns:
            plasmid_compat[f'site_i_{plasmid}_compatible'] = False
    
    # Add all 3mer AA combinations
    for _, row in enzyme_data.iterrows():
        site_i_rows.append({
            'site_i_enzyme': enzyme_name,
            'site_i_recognition_site': enzyme_site,
            'site_i_ovhg': enzyme_ovhg,
            'site_i_3mer_aa': row['re_site_shifted_tl'],
            'site_i_9mer_bp': row.get('re_site_shifted', ''),
            **plasmid_compat
        })

df_site_i = pd.DataFrame(site_i_rows)
elapsed = time.time() - start_time
print(f"   ✓ Site I: {len(df_site_i):,} rows in {elapsed:.2f}s")
print(f"     Unique enzymes: {df_site_i['site_i_enzyme'].nunique()}")
print(f"     Unique 3mer AA: {df_site_i['site_i_3mer_aa'].nunique()}")

# ============================================================================
# Build Site II DataFrame
# ============================================================================
print("\n4. Building Site II DataFrame...")
start_time = time.time()
site_ii_rows = []

for enzyme_name in df_silent_mutation['name'].unique():
    enzyme_data = df_silent_mutation[df_silent_mutation['name'] == enzyme_name]
    first_row = enzyme_data.iloc[0]
    enzyme_site = first_row.get('site', '')
    enzyme_ovhg = first_row['ovhg']
    
    # Get plasmid compatibility
    plasmid_compat = {}
    if enzyme_name in df_plasmid_check.index:
        for plasmid in df_plasmid_check.columns:
            plasmid_compat[f'site_ii_{plasmid}_compatible'] = bool(df_plasmid_check.loc[enzyme_name, plasmid])
    else:
        for plasmid in df_plasmid_check.columns:
            plasmid_compat[f'site_ii_{plasmid}_compatible'] = False
    
    # Add all 3mer AA combinations
    for _, row in enzyme_data.iterrows():
        site_ii_rows.append({
            'site_ii_enzyme': enzyme_name,
            'site_ii_recognition_site': enzyme_site,
            'site_ii_ovhg': enzyme_ovhg,
            'site_ii_3mer_aa': row['re_site_shifted_tl'],
            'site_ii_9mer_bp_original': row.get('re_site_shifted', ''),
            'site_ii_9mer_bp_mutated': row.get('re_site_mutate_shifted', ''),
            'site_ii_search_direction': 'right',
            **plasmid_compat
        })

df_site_ii = pd.DataFrame(site_ii_rows)
elapsed = time.time() - start_time
print(f"   ✓ Site II: {len(df_site_ii):,} rows in {elapsed:.2f}s")
print(f"     Unique enzymes: {df_site_ii['site_ii_enzyme'].nunique()}")
print(f"     Unique 3mer AA: {df_site_ii['site_ii_3mer_aa'].nunique()}")

# ============================================================================
# Combine with HURDLER pairing rules
# ============================================================================
print("\n5. Combining with HURDLER pairing rules...")
print("   Rules:")
print("     - Site I and Site II must have different overhangs")
print("     - Site III must have same overhang as Site II")
print("     - Both enzymes must be compatible with plasmid")
print("     - Pattern = Site_I_3mer.*?Site_II_3mer")

start_time = time.time()
combined_rows = []
skipped_same_ovhg = 0
skipped_same_pattern = 0
skipped_no_site_iii = 0
skipped_no_plasmid = 0

site_i_enzymes = df_site_i['site_i_enzyme'].unique()
site_ii_enzymes = df_site_ii['site_ii_enzyme'].unique()

print(f"\n   Processing {len(site_i_enzymes)} × {len(site_ii_enzymes)} = {len(site_i_enzymes) * len(site_ii_enzymes)} enzyme pairs...")

for i, site_i_enzyme in enumerate(site_i_enzymes):
    site_i_subset = df_site_i[df_site_i['site_i_enzyme'] == site_i_enzyme]
    site_i_ovhg = site_i_subset.iloc[0]['site_i_ovhg']
    
    for site_ii_enzyme in site_ii_enzymes:
        site_ii_subset = df_site_ii[df_site_ii['site_ii_enzyme'] == site_ii_enzyme]
        site_ii_ovhg = site_ii_subset.iloc[0]['site_ii_ovhg']
        
        # Rule 1: Different overhangs
        if site_i_ovhg == site_ii_ovhg:
            skipped_same_ovhg += 1
            continue
        
        # Rule 2: Find Site III with same ovhg as Site II
        compatible_site_iii = df_site_iii[df_site_iii['ovhg'] == site_ii_ovhg]
        if len(compatible_site_iii) == 0:
            skipped_no_site_iii += 1
            continue
        
        site_iii_list = compatible_site_iii['enzyme'].tolist()
        site_iii_str = ','.join(site_iii_list)
        
        # Combine all 3mer AA combinations
        for _, row_i in site_i_subset.iterrows():
            for _, row_ii in site_ii_subset.iterrows():
                # Skip same 3mer AA
                if row_i['site_i_3mer_aa'] == row_ii['site_ii_3mer_aa']:
                    skipped_same_pattern += 1
                    continue
                
                # Create pattern
                pattern = f"{row_i['site_i_3mer_aa']}.*?{row_ii['site_ii_3mer_aa']}"
                
                # Rule 3: Plasmid compatibility
                plasmid_compat = {}
                at_least_one_compatible = False
                
                for plasmid in df_plasmid_check.columns:
                    site_i_compat = row_i[f'site_i_{plasmid}_compatible']
                    site_ii_compat = row_ii[f'site_ii_{plasmid}_compatible']
                    
                    # Check Site III
                    site_iii_compat = False
                    for site_iii_enzyme in site_iii_list:
                        if site_iii_enzyme in df_plasmid_check.index:
                            if df_plasmid_check.loc[site_iii_enzyme, plasmid]:
                                site_iii_compat = True
                                break
                    
                    compat = site_i_compat and site_ii_compat and site_iii_compat
                    plasmid_compat[f'{plasmid}_compatible'] = compat
                    
                    if compat:
                        at_least_one_compatible = True
                
                if not at_least_one_compatible:
                    skipped_no_plasmid += 1
                    continue
                
                # Add row
                combined_rows.append({
                    'pattern': pattern,
                    'site_i_enzyme': row_i['site_i_enzyme'],
                    'site_i_recognition_site': row_i['site_i_recognition_site'],
                    'site_i_ovhg': row_i['site_i_ovhg'],
                    'site_i_3mer_aa': row_i['site_i_3mer_aa'],
                    'site_i_9mer_bp': row_i['site_i_9mer_bp'],
                    'site_ii_enzyme': row_ii['site_ii_enzyme'],
                    'site_ii_recognition_site': row_ii['site_ii_recognition_site'],
                    'site_ii_ovhg': row_ii['site_ii_ovhg'],
                    'site_ii_3mer_aa': row_ii['site_ii_3mer_aa'],
                    'site_ii_9mer_bp_original': row_ii['site_ii_9mer_bp_original'],
                    'site_ii_9mer_bp_mutated': row_ii['site_ii_9mer_bp_mutated'],
                    'site_ii_relative_to_i': 'right',
                    'compatible_site_iii_enzymes': site_iii_str,
                    'num_compatible_site_iii': len(site_iii_list),
                    **plasmid_compat
                })
    
    if (i + 1) % 10 == 0:
        elapsed = time.time() - start_time
        rate = (i + 1) / elapsed
        eta = (len(site_i_enzymes) - i - 1) / rate if rate > 0 else 0
        print(f"     {i+1}/{len(site_i_enzymes)} enzymes, {len(combined_rows):,} rows, ETA: {eta:.0f}s")

df_combined = pd.DataFrame(combined_rows)
elapsed = time.time() - start_time

print("\n" + "="*80)
print("✓ RESULTS")
print("="*80)
print(f"  Total rows: {len(df_combined):,}")
print(f"  Unique patterns: {df_combined['pattern'].nunique():,}")
print(f"  Unique Site I enzymes: {df_combined['site_i_enzyme'].nunique()}")
print(f"  Unique Site II enzymes: {df_combined['site_ii_enzyme'].nunique()}")
print(f"  Time: {elapsed:.1f}s")

print(f"\nSkipped:")
print(f"  Same overhang: {skipped_same_ovhg:,}")
print(f"  Same 3mer AA: {skipped_same_pattern:,}")
print(f"  No Site III: {skipped_no_site_iii:,}")
print(f"  No plasmid: {skipped_no_plasmid:,}")

# Show sample
print(f"\nSample (first 3 rows):")
print(df_combined[['pattern', 'site_i_enzyme', 'site_i_ovhg', 'site_ii_enzyme', 'site_ii_ovhg']].head(3))

print("\n" + "="*80)
print("TEST COMPLETE!")
print("="*80)
