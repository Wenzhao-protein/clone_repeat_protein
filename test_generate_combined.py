#!/usr/bin/env python3
"""
Test script to generate combined dataframe following the optimized method
"""

import pandas as pd
import numpy as np
from Bio import Restriction
import sys

print("="*80)
print("TESTING COMBINED DATAFRAME GENERATION")
print("="*80)

# Step 1: Load base data
print("\n1. Loading base data...")
df_seamless = pd.read_csv('./utils/output/restriction_enzyme_seamless_insert.csv')
df_silent = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')
df_plasmid_check = pd.read_csv('./utils/output/plasmid_digest_check.csv', index_col=0)
df_neb = pd.read_csv('./utils/output/neb_buffer_activity_cleaned.csv')
codon_freq = pd.read_csv('./utils/output/codon_usage.csv')

print(f"  Seamless insert: {len(df_seamless)} rows, {df_seamless['name'].nunique()} enzymes")
print(f"  Silent mutation: {len(df_silent)} rows, {df_silent['name'].nunique()} enzymes")
print(f"  Plasmid check: {df_plasmid_check.shape}")

# Step 2: Process Site I and Site II
print("\n2. Processing Site I and Site II...")

codon_freq_dict = dict(zip(codon_freq['codon'], codon_freq['frequency']))

def calculate_codon_usage_freq(dna_seq, codon_freq_dict):
    if len(dna_seq) != 9:
        return 0
    codon1 = dna_seq[0:3]
    codon2 = dna_seq[3:6]
    codon3 = dna_seq[6:9]
    freq1 = codon_freq_dict.get(codon1, 0)
    freq2 = codon_freq_dict.get(codon2, 0)
    freq3 = codon_freq_dict.get(codon3, 0)
    return (freq1 + freq2 + freq3) / 3

# Site I
df_seamless['codon_usage_freq'] = df_seamless['codon_usage']  # Already calculated
df_site_i = df_seamless.sort_values('codon_usage_freq', ascending=False).groupby(
    ['name', 're_site_shifted_tl']
).first().reset_index()
df_site_i = df_site_i[['name', 're_site_shifted_tl', 're_site_shifted', 'codon_usage_freq']].rename(columns={
    'name': 'enzyme',
    're_site_shifted_tl': '3mer_aa',
    're_site_shifted': 'dna_seq'
})

print(f"  Site I: {len(df_site_i)} unique (enzyme, 3mer_AA) pairs")

# Site II
def get_mutation_direction(orig, mutated):
    for i in range(len(orig)):
        if orig[i] != mutated[i]:
            return 'left' if i < 4 else 'right'
    return 'unknown'

df_silent['search_direction'] = df_silent.apply(
    lambda row: get_mutation_direction(row['re_site_shifted'], row['re_site_mutate_shifted']), 
    axis=1
)
df_silent['codon_usage_freq'] = df_silent['codon_usage_mutate']  # Already calculated
df_site_ii = df_silent.sort_values('codon_usage_freq', ascending=False).groupby(
    ['name', 're_site_shifted_tl']
).first().reset_index()
df_site_ii = df_site_ii[['name', 're_site_shifted_tl', 're_site_shifted', 're_site_mutate_shifted', 
                          'codon_usage_freq', 'search_direction']].rename(columns={
    'name': 'enzyme',
    're_site_shifted_tl': '3mer_aa',
    're_site_shifted': 'dna_seq_original',
    're_site_mutate_shifted': 'dna_seq_mutated'
})

print(f"  Site II: {len(df_site_ii)} unique (enzyme, 3mer_AA) pairs")

# Step 3: Generate combined dataframe (simple cartesian product)
print("\n3. Generating combined dataframe (cartesian product)...")

# Add overhang info
def get_enzyme_ovhg(enzyme_name):
    try:
        enzyme = getattr(Restriction, enzyme_name)
        return enzyme.ovhg
    except:
        return None

df_site_i['ovhg'] = df_site_i['enzyme'].apply(get_enzyme_ovhg)
df_site_ii['ovhg'] = df_site_ii['enzyme'].apply(get_enzyme_ovhg)

# Merge
df_combined = df_site_i.merge(df_site_ii, how='cross', suffixes=('_i', '_ii'))

print(f"  Combined: {len(df_combined):,} combinations")
print(f"  Memory: {df_combined.memory_usage(deep=True).sum()/1e6:.1f} MB")

# Add pattern
def generate_pattern(row):
    if row['search_direction'] == 'right':
        return f"{row['3mer_aa_i']}.*?{row['3mer_aa_ii']}"
    else:
        return f"{row['3mer_aa_ii']}.*?{row['3mer_aa_i']}"

df_combined['pattern'] = df_combined.apply(generate_pattern, axis=1)

print(f"  Unique patterns: {df_combined['pattern'].nunique():,}")

# Save
output_path = './output/hurdler_combined_test.csv'
df_combined.to_csv(output_path, index=False)
print(f"\n✓ Saved to: {output_path}")
print(f"✓ Shape: {df_combined.shape}")

# Show sample
print("\n4. Sample data:")
print(df_combined[['pattern', 'enzyme_i', '3mer_aa_i', 'enzyme_ii', '3mer_aa_ii', 'search_direction']].head(5))

print("\n" + "="*80)
print("TEST COMPLETE")
print("="*80)
