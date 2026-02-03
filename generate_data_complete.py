#!/usr/bin/env python3
"""
Complete data generation script following the optimized notebook workflow
"""

import pandas as pd
import numpy as np
from Bio import Restriction
import gc
import sys
import time

print("="*80)
print("HURDLER DATA GENERATION - COMPLETE WORKFLOW")
print("="*80)

# ============================================================================
# STEP 1: Load Base Data
# ============================================================================
print("\nSTEP 1: Loading base data...")
df_seamless_insert = pd.read_csv('./utils/output/restriction_enzyme_seamless_insert.csv')
df_silent_mutation = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')
df_plasmid_check = pd.read_csv('./utils/output/plasmid_digest_check.csv', index_col=0)
df_neb = pd.read_csv('./utils/output/neb_buffer_activity_cleaned.csv')
df_fidelity = pd.read_csv('./utils/output/orthogonality.csv')

plasmid_names = [
    'pGEX-4T-1', 'pMAL-c5X', 'pET-21a(+)', 'pET-28a(+)',
    'pET-28a(+)_start_codon', 'pCold_I', 'pUC18', 'pQE-3'
]

print(f"  ✓ Seamless insert: {len(df_seamless_insert)} rows, {df_seamless_insert['name'].nunique()} enzymes")
print(f"  ✓ Silent mutation: {len(df_silent_mutation)} rows, {df_silent_mutation['name'].nunique()} enzymes")
print(f"  ✓ Fidelity/orthogonality: {len(df_fidelity)} pairs")

# ============================================================================
# STEP 2: Generate Site III Data
# ============================================================================
print("\n" + "="*80)
print("STEP 2: Generating Site III data...")
print("="*80)

site_iii_data = []
for enzyme in Restriction.AllEnzymes:
    try:
        enzyme_name = enzyme.__name__
        site = str(enzyme.site)
        ovhg = enzyme.ovhg
        fst5 = enzyme.fst5
        site_length = len(site)
        
        # Filter criteria
        if ovhg is None or fst5 is None:
            continue
        if not all(base in 'ATCG' for base in site):
            continue
        if ovhg not in [-4, 2]:
            continue
        
        # Type IIS check
        is_type_iis = fst5 < 0 or fst5 > site_length
        if not is_type_iis:
            continue
        
        # NEB quality check
        if enzyme_name not in df_neb['enzyme'].values:
            continue
        enzyme_neb = df_neb[df_neb['enzyme'] == enzyme_name].iloc[0]
        if enzyme_neb['ligation_efficiencies'] == 'low':
            continue
        if enzyme_neb['star_activity'] == True:
            continue
        
        site_iii_data.append({
            'enzyme': enzyme_name,
            'ovhg': ovhg,
            'site': site,
            'site_length': site_length,
            'is_type_iis': is_type_iis
        })
    except:
        continue

df_site_iii = pd.DataFrame(site_iii_data)
df_site_iii.to_csv('./output/hurdler_site_iii_data.csv', index=False)
print(f"✓ Generated {len(df_site_iii)} Site III enzymes")
print(f"  Overhang distribution: {dict(df_site_iii['ovhg'].value_counts().sort_index())}")

# ============================================================================
# STEP 3: Build RE Pairing Matrix
# ============================================================================
print("\n" + "="*80)
print("STEP 3: Building RE pairing matrix...")
print("="*80)

site_i_enzymes = sorted(df_seamless_insert['name'].unique())
site_ii_enzymes = sorted(df_silent_mutation['name'].unique())

print(f"Site I enzymes: {len(site_i_enzymes)}")
print(f"Site II enzymes: {len(site_ii_enzymes)}")

def get_enzyme_ovhg(enzyme_name):
    try:
        enzyme = getattr(Restriction, enzyme_name)
        return enzyme.ovhg
    except:
        return None

def check_enzyme_pairing(enzyme1, enzyme2, df_fidelity, threshold=1):
    ovhg1 = get_enzyme_ovhg(enzyme1)
    ovhg2 = get_enzyme_ovhg(enzyme2)
    
    if ovhg1 != ovhg2:
        return True
    
    result = df_fidelity[
        ((df_fidelity['re1'] == enzyme1) & (df_fidelity['re2'] == enzyme2)) |
        ((df_fidelity['re1'] == enzyme2) & (df_fidelity['re2'] == enzyme1))
    ]
    
    if len(result) > 0:
        return result.iloc[0]['orthogonality'] >= threshold
    
    return True

print("Building matrix...")
site_i_ii_matrix = pd.DataFrame(
    index=site_i_enzymes,
    columns=site_ii_enzymes,
    dtype=bool
)

for i, enzyme_i in enumerate(site_i_enzymes):
    for j, enzyme_ii in enumerate(site_ii_enzymes):
        site_i_ii_matrix.loc[enzyme_i, enzyme_ii] = check_enzyme_pairing(
            enzyme_i, enzyme_ii, df_fidelity, threshold=1
        )
    
    if (i + 1) % 10 == 0:
        print(f"  Processed {i+1}/{len(site_i_enzymes)} Site I enzymes")

total_pairs = len(site_i_enzymes) * len(site_ii_enzymes)
compatible_pairs = site_i_ii_matrix.sum().sum()
print(f"✓ Total pairs: {total_pairs:,}, Compatible: {compatible_pairs:,} ({compatible_pairs/total_pairs*100:.1f}%)")

site_i_ii_matrix.to_csv('./output/site_i_site_ii_pairing_matrix.csv')
print("✓ Matrix saved")

# ============================================================================
# STEP 4: Generate Site I and Site II Dataframes
# ============================================================================
print("\n" + "="*80)
print("STEP 4: Generating Site I and Site II dataframes...")
print("="*80)

# Site I
df_site_i = df_seamless_insert.sort_values('codon_usage', ascending=False).groupby(
    ['name', 're_site_shifted_tl']
).first().reset_index()

df_site_i = df_site_i[['name', 're_site_shifted_tl', 're_site_shifted', 'codon_usage']].rename(columns={
    'name': 'enzyme',
    're_site_shifted_tl': '3mer_aa',
    're_site_shifted': 'dna_seq',
    'codon_usage': 'codon_usage_freq'
})

print(f"✓ Site I: {len(df_site_i):,} unique (enzyme, 3mer_AA) pairs")

# Site II
def get_mutation_direction(orig, mutated):
    for i in range(len(orig)):
        if orig[i] != mutated[i]:
            return 'left' if i < 4 else 'right'
    return 'unknown'

df_silent_mutation['search_direction'] = df_silent_mutation.apply(
    lambda row: get_mutation_direction(row['re_site_shifted'], row['re_site_mutate_shifted']), 
    axis=1
)

df_site_ii = df_silent_mutation.sort_values('codon_usage_mutate', ascending=False).groupby(
    ['name', 're_site_shifted_tl']
).first().reset_index()

df_site_ii = df_site_ii[['name', 're_site_shifted_tl', 're_site_shifted', 're_site_mutate_shifted', 
                          'codon_usage_mutate', 'search_direction']].rename(columns={
    'name': 'enzyme',
    're_site_shifted_tl': '3mer_aa',
    're_site_shifted': 'dna_seq_original',
    're_site_mutate_shifted': 'dna_seq_mutated',
    'codon_usage_mutate': 'codon_usage_freq'
})

print(f"✓ Site II: {len(df_site_ii):,} unique (enzyme, 3mer_AA) pairs")

df_site_i.to_csv('./output/hurdler_site_i_data.csv', index=False)
df_site_ii.to_csv('./output/hurdler_site_ii_data.csv', index=False)
print("✓ Dataframes saved")

# ============================================================================
# STEP 5: Generate Combined Dataframe
# ============================================================================
print("\n" + "="*80)
print("STEP 5: Generating combined dataframe...")
print("="*80)

# Add overhang info
df_site_i['ovhg'] = df_site_i['enzyme'].apply(get_enzyme_ovhg)
df_site_ii['ovhg'] = df_site_ii['enzyme'].apply(get_enzyme_ovhg)

# Get compatible pairs
compatible_pairs = []
for enzyme_i in site_i_ii_matrix.index:
    for enzyme_ii in site_i_ii_matrix.columns:
        if site_i_ii_matrix.loc[enzyme_i, enzyme_ii]:
            compatible_pairs.append((enzyme_i, enzyme_ii))

print(f"Compatible enzyme pairs: {len(compatible_pairs):,}")

# Generate combinations in batches
output_path = './output/hurdler_combined_dataframe.csv'
batch_size = 100
total_rows = 0
first_batch = True

num_batches = (len(compatible_pairs) + batch_size - 1) // batch_size
print(f"Processing in {num_batches} batches...")

batch_start_time = time.time()

for batch_idx in range(0, len(compatible_pairs), batch_size):
    batch_end = min(batch_idx + batch_size, len(compatible_pairs))
    pairs_batch = compatible_pairs[batch_idx:batch_end]
    
    batch_num = batch_idx // batch_size + 1
    
    expanded_list = []
    
    for enzyme_i, enzyme_ii in pairs_batch:
        site_i_filtered = df_site_i[df_site_i['enzyme'] == enzyme_i]
        site_ii_filtered = df_site_ii[df_site_ii['enzyme'] == enzyme_ii]
        
        if len(site_i_filtered) == 0 or len(site_ii_filtered) == 0:
            continue
        
        # Plasmid compatibility
        plasmid_compat = {}
        for plasmid in plasmid_names:
            enzyme_i_ok = enzyme_i in df_plasmid_check.index and df_plasmid_check.loc[enzyme_i, plasmid]
            enzyme_ii_ok = enzyme_ii in df_plasmid_check.index and df_plasmid_check.loc[enzyme_ii, plasmid]
            plasmid_compat[f"{plasmid}_compatible"] = enzyme_i_ok and enzyme_ii_ok
        
        if not any(plasmid_compat.values()):
            continue
        
        # Site III compatibility
        ovhg_ii = site_ii_filtered.iloc[0]['ovhg']
        site_iii_compatible = ','.join(
            df_site_iii[df_site_iii['ovhg'] == ovhg_ii]['enzyme'].tolist()
        )
        
        if not site_iii_compatible:
            continue
        
        # Create all 3mer combinations
        for _, site_i_row in site_i_filtered.iterrows():
            for _, site_ii_row in site_ii_filtered.iterrows():
                if site_ii_row['search_direction'] == 'right':
                    pattern = f"{site_i_row['3mer_aa']}.*?{site_ii_row['3mer_aa']}"
                else:
                    pattern = f"{site_ii_row['3mer_aa']}.*?{site_i_row['3mer_aa']}"
                
                new_row = {
                    'pattern': pattern,
                    'enzyme_i': enzyme_i,
                    '3mer_aa_i': site_i_row['3mer_aa'],
                    'dna_seq_i': site_i_row['dna_seq'],
                    'codon_usage_freq_i': site_i_row['codon_usage_freq'],
                    'ovhg_i': site_i_row['ovhg'],
                    'enzyme_ii': enzyme_ii,
                    '3mer_aa_ii': site_ii_row['3mer_aa'],
                    'dna_seq_original': site_ii_row['dna_seq_original'],
                    'dna_seq_mutated': site_ii_row['dna_seq_mutated'],
                    'codon_usage_freq_ii': site_ii_row['codon_usage_freq'],
                    'search_direction': site_ii_row['search_direction'],
                    'ovhg_ii': site_ii_row['ovhg'],
                    'site_iii_compatible': site_iii_compatible
                }
                new_row.update(plasmid_compat)
                expanded_list.append(new_row)
    
    if len(expanded_list) > 0:
        df_batch = pd.DataFrame(expanded_list)
        
        if first_batch:
            df_batch.to_csv(output_path, index=False, mode='w')
            first_batch = False
        else:
            df_batch.to_csv(output_path, index=False, mode='a', header=False)
        
        total_rows += len(df_batch)
        del df_batch
    
    del expanded_list
    gc.collect()
    
    if batch_num % 10 == 0 or batch_num == num_batches:
        elapsed = time.time() - batch_start_time
        eta_seconds = (elapsed / batch_end) * (len(compatible_pairs) - batch_end) if batch_end > 0 else 0
        print(f"  Batch {batch_num}/{num_batches}: {total_rows:,} rows, {elapsed:.1f}s, ETA {eta_seconds:.1f}s")
        sys.stdout.flush()

print(f"\n✓ Generated {total_rows:,} combinations")
print(f"✓ Saved to: {output_path}")

# Display sample
if total_rows > 0:
    df_sample = pd.read_csv(output_path, nrows=min(100, total_rows))
    print(f"\nSample statistics:")
    print(f"  Unique patterns: {df_sample['pattern'].nunique()}")
    print(f"  Unique Site I enzymes: {df_sample['enzyme_i'].nunique()}")
    print(f"  Unique Site II enzymes: {df_sample['enzyme_ii'].nunique()}")
    print(f"\nFirst 3 rows:")
    print(df_sample[['pattern', 'enzyme_i', '3mer_aa_i', 'enzyme_ii', '3mer_aa_ii']].head(3))

print("\n" + "="*80)
print("✓ DATA GENERATION COMPLETE")
print("="*80)
