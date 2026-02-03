"""
HURDLER Success Rate Analysis
Complete pipeline for analyzing HURDLER cloning strategy success rates
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import pickle
import random
import gc
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
from Bio.Restriction import AllEnzymes

# Set plotting style
plt.style.use('default')
sns.set_palette("husl")

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Create output directory
output_dir = './output'
os.makedirs(output_dir, exist_ok=True)

print("="*80)
print("HURDLER SUCCESS RATE ANALYSIS")
print("="*80)

# ============================================================================
# STEP 1: GENERATE THREE-SITE ENZYME COMBINATIONS (DF1)
# ============================================================================

print("\n" + "="*80)
print("STEP 1: GENERATING THREE-SITE ENZYME COMBINATIONS (DF1)")
print("="*80)

# Load methylation sensitivity data
print("\n1. Loading methylation sensitivity data...")
df_methylation = pd.read_csv('./input/methylation_check.csv')
dh5a_compatible = set(
    df_methylation.loc[~df_methylation['6mA_5mC_sensitive'], 'enzyme'].dropna().tolist() +
    df_methylation.loc[~df_methylation['6mA_5mC_sensitive'], 'prototype'].dropna().tolist()
)
print(f"   DH5α compatible enzymes: {len(dh5a_compatible)}")

# Load NEB quality data
print("\n2. Loading NEB quality data...")
df_neb = pd.read_csv('./input/neb_buffer_activity_cleaned.csv')
print(f"   NEB enzymes in database: {len(df_neb)}")

# Load plasmid compatibility data
print("\n3. Loading plasmid compatibility data...")
df_plasmid = pd.read_csv('./input/plasmid_digest_check.csv', index_col=0)
plasmid_cols = ['pGEX-4T-1', 'pMAL-c5X', 'pET-21a(+)', 'pET-28a(+)', 
                'pET-28a(+)_start_codon', 'pCold_I', 'pUC18', 'pQE-3']

# Build enzyme database from Bio.Restriction
print("\n4. Building enzyme database from Bio.Restriction...")
enzyme_data = []

for enzyme in AllEnzymes:
    try:
        name = enzyme.__name__
        site = str(enzyme.site)
        ovhg = enzyme.ovhg
        fst5 = enzyme.fst5
        fst3 = enzyme.fst3
        suppl = enzyme.suppl
        is_comm = enzyme.is_comm()
        
        if ovhg is None or fst5 is None:
            continue
        
        if not all(base in 'ATCG' for base in site):
            continue
        
        if not is_comm:
            continue
        
        if abs(ovhg) not in [2, 3, 4, 5]:
            continue
        
        has_neb = 'N' in str(suppl)
            
        enzyme_data.append({
            'enzyme': name,
            'site': site,
            'ovhg': ovhg,
            'fst5': fst5,
            'fst3': fst3,
            'site_length': len(site),
            'suppl': suppl,
            'has_neb': has_neb,
            'methylation_compatible': name in dh5a_compatible
        })
    except:
        continue

df_enzymes = pd.DataFrame(enzyme_data)
print(f"   Total enzymes extracted: {len(df_enzymes)}")

# Define Type IIS detection function
def is_type_iis_enzyme(row):
    """Check if enzyme is Type IIS (cuts outside recognition site)"""
    return row['fst5'] < 0 or row['fst5'] > row['site_length']

df_enzymes['is_type_iis'] = df_enzymes.apply(is_type_iis_enzyme, axis=1)

# Add NEB quality checks
def check_neb_quality(enzyme_name, df_neb):
    """Check if enzyme has good NEB characteristics"""
    if enzyme_name not in df_neb['enzyme'].values:
        return False, False
    
    enzyme_data = df_neb[df_neb['enzyme'] == enzyme_name].iloc[0]
    ligation_ok = enzyme_data['ligation_efficiencies'] != 'low'
    no_star = enzyme_data['star_activity'] == False
    
    return ligation_ok, no_star

neb_quality = df_enzymes['enzyme'].apply(lambda x: check_neb_quality(x, df_neb))
df_enzymes['ligation_ok'] = neb_quality.apply(lambda x: x[0])
df_enzymes['no_star_activity'] = neb_quality.apply(lambda x: x[1])

# Filter Site I enzymes (seamless insert)
print("\n5. Filtering Site I enzymes (seamless insert)...")
site_i = df_enzymes[
    (df_enzymes['methylation_compatible'] == True) &
    (df_enzymes['is_type_iis'] == False) &
    (df_enzymes['has_neb'] == True) &
    (df_enzymes['ligation_ok'] == True) &
    (df_enzymes['no_star_activity'] == True)
].copy()

site_i = site_i.merge(df_plasmid, left_on='enzyme', right_index=True, how='left')
site_i = site_i[site_i[plasmid_cols].any(axis=1)]
print(f"   Site I enzymes: {len(site_i)}")

# Filter Site II enzymes (silent mutation)
print("\n6. Filtering Site II enzymes (silent mutation)...")
site_ii = df_enzymes[
    (df_enzymes['methylation_compatible'] == True) &
    (df_enzymes['is_type_iis'] == False) &
    (df_enzymes['has_neb'] == True) &
    (df_enzymes['ligation_ok'] == True) &
    (df_enzymes['no_star_activity'] == True)
].copy()

site_ii = site_ii.merge(df_plasmid, left_on='enzyme', right_index=True, how='left')
site_ii = site_ii[site_ii[plasmid_cols].any(axis=1)]
print(f"   Site II enzymes: {len(site_ii)}")

# Filter Site III enzymes (Type IIS)
print("\n7. Filtering Site III enzymes (Type IIS)...")
site_iii = df_enzymes[
    (df_enzymes['is_type_iis'] == True) &
    (df_enzymes['has_neb'] == True) &
    (df_enzymes['ligation_ok'] == True) &
    (df_enzymes['no_star_activity'] == True)
].copy()

site_iii = site_iii.merge(df_plasmid, left_on='enzyme', right_index=True, how='left')
site_iii = site_iii[site_iii[plasmid_cols].any(axis=1)]
print(f"   Site III enzymes (Type IIS): {len(site_iii)}")

# Filter Site II for overhang compatibility with Site III
if len(site_iii) > 0:
    print("\n8. Filtering Site II for overhang compatibility with Site III...")
    site_iii_ovhgs = set(site_iii['ovhg'].unique())
    site_ii_before = len(site_ii)
    site_ii = site_ii[site_ii['ovhg'].isin(site_iii_ovhgs)]
    print(f"   Site II after overhang filter: {len(site_ii)} (removed {site_ii_before - len(site_ii)})")

# Generate all three-site combinations
print("\n9. Generating three-site combinations...")
combinations = []

if len(site_i) > 0 and len(site_ii) > 0 and len(site_iii) > 0:
    for _, i_row in site_i.iterrows():
        for _, ii_row in site_ii.iterrows():
            for _, iii_row in site_iii.iterrows():
                for plasmid in plasmid_cols:
                    if i_row[plasmid] and ii_row[plasmid] and iii_row[plasmid]:
                        combinations.append({
                            'site_i_enzyme': i_row['enzyme'],
                            'site_ii_enzyme': ii_row['enzyme'],
                            'site_iii_enzyme': iii_row['enzyme'],
                            'site_i_ovhg': i_row['ovhg'],
                            'site_ii_ovhg': ii_row['ovhg'],
                            'site_iii_ovhg': iii_row['ovhg'],
                            f"{plasmid}_compatible": True
                        })

if len(combinations) > 0:
    df1 = pd.DataFrame(combinations)
    
    df1 = df1.groupby(['site_i_enzyme', 'site_ii_enzyme', 'site_iii_enzyme',
                       'site_i_ovhg', 'site_ii_ovhg', 'site_iii_ovhg']).agg({
        f"{plasmid}_compatible": 'any' for plasmid in plasmid_cols
    }).reset_index()
    
    for plasmid in plasmid_cols:
        col_name = f"{plasmid}_compatible"
        df1[col_name] = df1[col_name].fillna(False)
    
    print(f"   Generated {len(df1):,} three-site combinations")
    
    df1_path = os.path.join(output_dir, 'hurdler_three_site_combinations_df1.csv')
    df1.to_csv(df1_path, index=False)
    print(f"   ✓ df1 saved to: {df1_path}")
else:
    print("   ⚠ No valid three-site combinations found!")
    sys.exit(1)

# ============================================================================
# STEP 2: GENERATE OPTIMIZED DF2 WITH CODON USAGE FILTERING
# ============================================================================

print("\n" + "="*80)
print("STEP 2: GENERATING OPTIMIZED DF2")
print("="*80)

# Load E.coli codon usage
print("\n1. Loading E.coli codon usage...")
codon_freq = pd.read_csv('./input/codon_usage.csv')
codon_freq_dict = dict(zip(codon_freq['codon'], codon_freq['frequency']))
print(f"   Loaded {len(codon_freq_dict)} codon frequencies")

def calculate_codon_usage_freq(dna_seq, codon_dict):
    """Calculate 9mer DNA sequence codon usage frequency product"""
    if len(dna_seq) != 9:
        return 0.0
    freq = 1.0
    for i in range(0, 9, 3):
        codon = dna_seq[i:i+3]
        freq *= codon_dict.get(codon, 0.0)
    return freq

# Load restriction enzyme data
print("\n2. Loading restriction enzyme data...")
df_seamless = pd.read_csv('./input/seamless_insert.csv')
df_silent = pd.read_csv('./input/slient_mutation.csv')

# Process Site I (seamless insert)
print("\n3. Pre-processing Site I (seamless insert)...")
df_seamless = df_seamless[~df_seamless['re_3aa_site'].str.contains('\\*', na=False)].copy()
print(f"   After stop codon filter: {len(df_seamless)} rows")

df_seamless['codon_usage_freq'] = df_seamless['re_9bp_site'].apply(
    lambda x: calculate_codon_usage_freq(x, codon_freq_dict)
)

df_seamless = df_seamless.sort_values('codon_usage_freq', ascending=False).groupby(
    ['name', 're_3aa_site']
).first().reset_index()
print(f"   After deduplication: {len(df_seamless)} unique pairs")

site_i_data = df_seamless[['name', 're_3aa_site', 're_9bp_site', 'codon_usage_freq']].rename(columns={
    'name': 'enzyme',
    're_3aa_site': '3mer_aa',
    're_9bp_site': 'dna_seq'
})

# Process Site II (silent mutation)
print("\n4. Pre-processing Site II (silent mutation)...")
df_silent = df_silent[~df_silent['re_3aa_site'].str.contains('\\*', na=False)].copy()
print(f"   After stop codon filter: {len(df_silent)} rows")

df_silent['mutated_codon_usage_freq'] = df_silent['re_9bp_site_mutated'].apply(
    lambda x: calculate_codon_usage_freq(x, codon_freq_dict)
)

def get_mutation_direction(orig, mutated):
    """Determine if mutation is in first 4bp (search left) or last 4bp (search right)"""
    for i in range(len(orig)):
        if orig[i] != mutated[i]:
            if i < 4:
                return 'left'
            else:
                return 'right'
    return 'unknown'

df_silent['search_direction'] = df_silent.apply(
    lambda row: get_mutation_direction(row['re_9bp_site'], row['re_9bp_site_mutated']), axis=1
)

df_silent = df_silent.sort_values('mutated_codon_usage_freq', ascending=False).groupby(
    ['name', 're_3aa_site']
).first().reset_index()
print(f"   After deduplication: {len(df_silent)} unique pairs")

site_ii_data = df_silent[['name', 're_3aa_site', 're_9bp_site', 
                          're_9bp_site_mutated', 'mutated_codon_usage_freq', 'search_direction']].rename(columns={
    'name': 'enzyme',
    're_3aa_site': '3mer_aa',
    're_9bp_site': 'dna_seq_original',
    're_9bp_site_mutated': 'dna_seq_mutated',
    'mutated_codon_usage_freq': 'codon_usage_freq'
})

print(f"   Search direction distribution:")
print(f"     Left: {(site_ii_data['search_direction'] == 'left').sum()}")
print(f"     Right: {(site_ii_data['search_direction'] == 'right').sum()}")

# Release memory
del df_seamless, df_silent
gc.collect()

# Streaming expansion of df1 to df2
print("\n5. Streaming df1 to df2...")
output_path = './output/hurdler_three_site_combinations_df2_optimized.csv'
batch_size = 100
total_rows = 0
first_batch = True

num_batches = (len(df1) + batch_size - 1) // batch_size
print(f"   Processing {len(df1)} df1 rows in {num_batches} batches...")

for batch_idx in range(0, len(df1), batch_size):
    batch_end = min(batch_idx + batch_size, len(df1))
    df1_batch = df1.iloc[batch_idx:batch_end]
    
    batch_num = batch_idx // batch_size + 1
    
    expanded_list = []
    
    for _, row in df1_batch.iterrows():
        site_i_filtered = site_i_data[site_i_data['enzyme'] == row['site_i_enzyme']]
        site_ii_filtered = site_ii_data[site_ii_data['enzyme'] == row['site_ii_enzyme']]
        
        if len(site_i_filtered) == 0 or len(site_ii_filtered) == 0:
            continue
        
        temp_df = pd.DataFrame([row] * len(site_i_filtered))
        temp_df = temp_df.reset_index(drop=True)
        site_i_copy = site_i_filtered.reset_index(drop=True)
        
        temp_df['site_i_3mer_aa'] = site_i_copy['3mer_aa'].values
        temp_df['site_i_dna'] = site_i_copy['dna_seq'].values
        temp_df['site_i_codon_usage_freq'] = site_i_copy['codon_usage_freq'].values
        
        temp_expanded = []
        for _, temp_row in temp_df.iterrows():
            for _, site_ii_row in site_ii_filtered.iterrows():
                new_row = temp_row.to_dict()
                new_row['site_ii_3mer_aa'] = site_ii_row['3mer_aa']
                new_row['site_ii_dna_original'] = site_ii_row['dna_seq_original']
                new_row['site_ii_dna_mutated'] = site_ii_row['dna_seq_mutated']
                new_row['site_ii_codon_usage_freq'] = site_ii_row['codon_usage_freq']
                new_row['site_ii_search_direction'] = site_ii_row['search_direction']
                new_row['site_iii_3mer_aa'] = 'N/A'
                new_row['site_iii_dna'] = 'N/A'
                temp_expanded.append(new_row)
        
        expanded_list.extend(temp_expanded)
    
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
    
    if batch_num % 100 == 0:
        print(f"   Batch {batch_num}/{num_batches}: {total_rows:,} rows")
        sys.stdout.flush()

print(f"\n✓ Generated {total_rows:,} total combinations")
print(f"✓ df2 saved to: {output_path}")

# ============================================================================
# STEP 3: CREATE LOOKUP DICTIONARY (OPTIMIZED)
# ============================================================================

print("\n" + "="*80)
print("STEP 3: CREATING LOOKUP DICTIONARY (OPTIMIZED)")
print("="*80)

plasmid_cols_lookup = [
    'pGEX-4T-1_compatible',
    'pMAL-c5X_compatible',
    'pET-21a(+)_compatible',
    'pET-28a(+)_compatible',
    'pET-28a(+)_start_codon_compatible',
    'pCold_I_compatible',
    'pUC18_compatible',
    'pQE-3_compatible'
]

# Step 3.1: Build Site I lookup
print("\n1. Building Site I lookup...")
print("   Structure: {(enzyme, 3mer_aa): {'9mer_bp': ..., 'plasmids': {...}}}")
site_i_lookup = {}

chunk_size = 1000000
chunk_num = 0
total_rows = 0

for chunk in pd.read_csv('./output/hurdler_three_site_combinations_df2_optimized.csv',
                         chunksize=chunk_size,
                         usecols=['site_i_enzyme', 'site_i_3mer_aa', 'site_i_dna'] + plasmid_cols_lookup):
    chunk_num += 1
    total_rows += len(chunk)
    
    for _, row in chunk.iterrows():
        key = (row['site_i_enzyme'], row['site_i_3mer_aa'])
        
        if key not in site_i_lookup:
            plasmid_dict = {plasmid.replace('_compatible', ''): bool(row[plasmid]) 
                           for plasmid in plasmid_cols_lookup}
            site_i_lookup[key] = {
                '9mer_bp': row['site_i_dna'],
                'plasmids': plasmid_dict
            }
        else:
            # Update plasmid compatibility (OR operation)
            for plasmid in plasmid_cols_lookup:
                plasmid_name = plasmid.replace('_compatible', '')
                if row[plasmid]:
                    site_i_lookup[key]['plasmids'][plasmid_name] = True
    
    if chunk_num % 10 == 0:
        print(f"   Processed {total_rows:,} rows")
        sys.stdout.flush()
    
    del chunk
    gc.collect()

print(f"   Site I lookup size: {len(site_i_lookup):,} unique (enzyme, 3mer_aa) pairs")

# Step 3.2: Build Site II & III lookup
print("\n2. Building Site II & III lookup...")
print("   Structure: {(site_ii_enzyme, site_ii_3mer_aa, site_iii_enzyme): [info_dict, ...]}")
site_ii_iii_lookup = defaultdict(list)

chunk_num = 0
total_rows = 0

for chunk in pd.read_csv('./output/hurdler_three_site_combinations_df2_optimized.csv',
                         chunksize=chunk_size,
                         usecols=['site_ii_enzyme', 'site_ii_3mer_aa', 'site_ii_dna_original',
                                 'site_ii_dna_mutated', 'site_ii_search_direction',
                                 'site_iii_enzyme'] + plasmid_cols_lookup):
    chunk_num += 1
    total_rows += len(chunk)
    
    for _, row in chunk.iterrows():
        key = (row['site_ii_enzyme'], row['site_ii_3mer_aa'], row['site_iii_enzyme'])
        
        plasmid_dict = {plasmid.replace('_compatible', ''): bool(row[plasmid]) 
                       for plasmid in plasmid_cols_lookup}
        
        # Skip if no compatible plasmids
        if not any(plasmid_dict.values()):
            continue
        
        info = {
            '9mer_bp_original': row['site_ii_dna_original'],
            '9mer_bp_mutated': row['site_ii_dna_mutated'],
            'direction': row['site_ii_search_direction'],
            'plasmids': plasmid_dict
        }
        
        # Check if this exact info already exists
        if info not in site_ii_iii_lookup[key]:
            site_ii_iii_lookup[key].append(info)
    
    if chunk_num % 10 == 0:
        print(f"   Processed {total_rows:,} rows")
        sys.stdout.flush()
    
    del chunk
    gc.collect()

print(f"   Site II & III lookup size: {len(site_ii_iii_lookup):,} unique (enzyme_ii, 3mer_aa_ii, enzyme_iii) tuples")

# Step 3.3: Combine lookups
print("\n3. Combining lookups to create final dictionary...")
print("   Final structure: {(3mer_i, 3mer_ii): [(site_i_info, site_ii_info, site_iii_enzyme, plasmid_dict), ...]}")

final_lookup = defaultdict(list)

# Read df2 again to get site_i to site_ii/iii mappings
chunk_num = 0
total_rows = 0
entries_added = 0

for chunk in pd.read_csv('./output/hurdler_three_site_combinations_df2_optimized.csv',
                         chunksize=chunk_size,
                         usecols=['site_i_enzyme', 'site_i_3mer_aa', 
                                 'site_ii_enzyme', 'site_ii_3mer_aa',
                                 'site_iii_enzyme'] + plasmid_cols_lookup):
    chunk_num += 1
    total_rows += len(chunk)
    
    for _, row in chunk.iterrows():
        site_i_key = (row['site_i_enzyme'], row['site_i_3mer_aa'])
        site_ii_iii_key = (row['site_ii_enzyme'], row['site_ii_3mer_aa'], row['site_iii_enzyme'])
        
        # Check if both keys exist in their respective lookups
        if site_i_key not in site_i_lookup or site_ii_iii_key not in site_ii_iii_lookup:
            continue
        
        site_i_data = site_i_lookup[site_i_key]
        
        # For each site_ii_iii combination
        for site_ii_iii_data in site_ii_iii_lookup[site_ii_iii_key]:
            # Find compatible plasmids (intersection)
            plasmid_dict = {}
            has_compatible = False
            for plasmid_name in site_i_data['plasmids'].keys():
                is_compatible = (site_i_data['plasmids'][plasmid_name] and 
                               site_ii_iii_data['plasmids'][plasmid_name])
                plasmid_dict[plasmid_name] = is_compatible
                if is_compatible:
                    has_compatible = True
            
            # Skip if no compatible plasmids
            if not has_compatible:
                continue
            
            # Create site_i_info
            site_i_info = {
                'enzyme': row['site_i_enzyme'],
                '3mer_aa': row['site_i_3mer_aa'],
                '9mer_bp': site_i_data['9mer_bp']
            }
            
            # Create site_ii_info
            site_ii_info = {
                'enzyme': row['site_ii_enzyme'],
                '3mer_aa': row['site_ii_3mer_aa'],
                '9mer_bp_original': site_ii_iii_data['9mer_bp_original'],
                '9mer_bp_mutated': site_ii_iii_data['9mer_bp_mutated'],
                'direction': site_ii_iii_data['direction']
            }
            
            # Create entry
            entry = (site_i_info, site_ii_info, row['site_iii_enzyme'], plasmid_dict)
            
            # Store with both orderings (since key has no order)
            pair1 = (row['site_i_3mer_aa'], row['site_ii_3mer_aa'])
            pair2 = (row['site_ii_3mer_aa'], row['site_i_3mer_aa'])
            
            # Check if this exact entry already exists
            if entry not in final_lookup[pair1]:
                final_lookup[pair1].append(entry)
                entries_added += 1
            if entry not in final_lookup[pair2]:
                final_lookup[pair2].append(entry)
                entries_added += 1
    
    if chunk_num % 10 == 0:
        print(f"   Processed {total_rows:,} rows, added {entries_added:,} entries")
        sys.stdout.flush()
    
    del chunk
    gc.collect()

print(f"\n4. Total rows processed: {total_rows:,}")
print(f"   Total entries added: {entries_added:,}")

# Display statistics
print("\n5. Final lookup statistics:")
print("="*80)
print(f"  Total unique 3mer pairs: {len(final_lookup):,}")
if len(final_lookup) > 0:
    print(f"  Average entries per pair: {sum(len(v) for v in final_lookup.values()) / len(final_lookup):.1f}")
    
    # Sample entry
    sample_key = list(final_lookup.keys())[0]
    sample_entry = final_lookup[sample_key][0]
    print(f"\n6. Sample entry structure:")
    print(f"  Key: {sample_key}")
    print(f"  Value[0]:")
    print(f"    Site I: {sample_entry[0]}")
    print(f"    Site II: {sample_entry[1]}")
    print(f"    Site III: {sample_entry[2]}")
    print(f"    Plasmids: {sample_entry[3]}")

# Save lookup dictionary
lookup_path = './output/hurdler_lookup_optimized.pkl'
print(f"\n7. Saving lookup dictionary to {lookup_path}...")
with open(lookup_path, 'wb') as f:
    pickle.dump(dict(final_lookup), f)

# Clean up intermediate lookups to free memory
del site_i_lookup, site_ii_iii_lookup, final_lookup
gc.collect()

print("\n✓ Lookup dictionary created successfully!")

# ============================================================================
# STEP 4: CALCULATE SUCCESS RATES
# ============================================================================

print("\n" + "="*80)
print("STEP 4: TESTING HURDLER SUCCESS RATES")
print("="*80)

def generate_random_aa_sequence(length):
    """Generate random amino acid sequence"""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(random.choice(amino_acids) for _ in range(length))

def extract_3mer_pairs_with_positions(sequence):
    """Extract all 3mer AA from doubled sequence with their positions"""
    doubled_seq = sequence + sequence
    
    three_mer_positions = defaultdict(list)
    for i in range(len(doubled_seq) - 2):
        three_mer = doubled_seq[i:i+3]
        three_mer_positions[three_mer].append(i)
    
    return three_mer_positions

def check_hurdler_success(sequence, hurdler_lookup, plasmid):
    """Check if sequence has valid HURDLER solution"""
    module_length = len(sequence)
    three_mer_positions = extract_3mer_pairs_with_positions(sequence)
    
    three_mers = list(three_mer_positions.keys())
    
    for i in range(len(three_mers)):
        for j in range(i+1, len(three_mers)):
            site_i_3mer = three_mers[i]
            site_ii_3mer = three_mers[j]
            
            pair = (site_i_3mer, site_ii_3mer)
            if pair not in hurdler_lookup:
                pair = (site_ii_3mer, site_i_3mer)
                if pair not in hurdler_lookup:
                    continue
            
            entries = hurdler_lookup[pair]
            
            for entry in entries:
                site_i_info, site_ii_info, site_iii_enzyme, plasmid_dict = entry
                
                if not plasmid_dict.get(plasmid, False):
                    continue
                
                if site_i_info['3mer_aa'] == site_i_3mer:
                    site_i_positions = three_mer_positions[site_i_3mer]
                    site_ii_positions = three_mer_positions[site_ii_3mer]
                else:
                    site_i_positions = three_mer_positions[site_ii_3mer]
                    site_ii_positions = three_mer_positions[site_i_3mer]
                
                search_direction = site_ii_info['direction']
                
                for site_ii_pos in site_ii_positions:
                    if search_direction == 'left':
                        for site_i_pos in site_i_positions:
                            if site_i_pos < site_ii_pos:
                                distance = site_ii_pos - site_i_pos
                                if 5 < distance < module_length:
                                    return True, len(three_mer_positions)
                    else:
                        for site_i_pos in site_i_positions:
                            if site_i_pos > site_ii_pos:
                                distance = site_i_pos - site_ii_pos
                                if 5 < distance < module_length:
                                    return True, len(three_mer_positions)
    
    return False, len(three_mer_positions)

# Load lookup dictionary
print("\n1. Loading lookup dictionary...")
with open('./output/hurdler_lookup_optimized.pkl', 'rb') as f:
    hurdler_lookup = pickle.load(f)

plasmids = ['pGEX-4T-1', 'pMAL-c5X', 'pET-21a(+)', 'pET-28a(+)', 
            'pET-28a(+)_start_codon', 'pCold_I', 'pUC18', 'pQE-3']
print(f"   Loaded lookup with {len(hurdler_lookup):,} 3mer pairs")

# Test parameters
module_lengths = list(range(7, 61))
n_tests_per_length = 1000

print(f"\n2. Testing sequences...")
print(f"   Lengths: {min(module_lengths)} to {max(module_lengths)} AA")
print(f"   Tests per length: {n_tests_per_length:,}")
print(f"   Total tests: {len(module_lengths) * n_tests_per_length * len(plasmids):,}")

results = []

for length in tqdm(module_lengths, desc="Testing lengths"):
    for test_idx in range(n_tests_per_length):
        sequence = generate_random_aa_sequence(length)
        
        for plasmid in plasmids:
            success, n_pairs = check_hurdler_success(sequence, hurdler_lookup, plasmid)
            
            results.append({
                'module_length': length,
                'test_idx': test_idx,
                'plasmid': plasmid,
                'sequence': sequence,
                'success': success,
                'n_3mer_pairs': n_pairs
            })

df_results = pd.DataFrame(results)

# Save raw results
results_path = './output/hurdler_success_rate_results.csv'
df_results.to_csv(results_path, index=False)
print(f"\n✓ Raw results saved to: {results_path}")

# Calculate summary statistics
df_summary = df_results.groupby(['module_length', 'plasmid']).agg(
    n_success=('success', 'sum'),
    n_total=('success', 'count')
).reset_index()

df_summary['success_rate'] = df_summary['n_success'] / df_summary['n_total']
df_summary['success_rate_pct'] = df_summary['success_rate'] * 100

# Save summary
summary_path = './output/hurdler_success_rate_summary.csv'
df_summary.to_csv(summary_path, index=False)
print(f"✓ Summary saved to: {summary_path}")

# ============================================================================
# STEP 5: VISUALIZATION
# ============================================================================

print("\n" + "="*80)
print("STEP 5: GENERATING VISUALIZATIONS")
print("="*80)

# Process data for visualization
df_summary['Sequence_Length'] = df_summary['module_length']
df_summary['Probability'] = df_summary['success_rate_pct']

plasmid_dfs = {}
for plasmid in plasmids:
    plasmid_dfs[plasmid] = df_summary[df_summary['plasmid'] == plasmid][['Sequence_Length', 'Probability']].copy()
    plasmid_dfs[plasmid] = plasmid_dfs[plasmid].sort_values('Sequence_Length').reset_index(drop=True)

# Main plot
plt.figure(figsize=(8, 4))

for plasmid in plasmids:
    data = plasmid_dfs[plasmid]
    plt.plot(data['Sequence_Length'], data['Probability'],
             label=plasmid, linewidth=2)

plt.xlabel('Sequence Length (AA)')
plt.ylabel('Success Rate (%)')
plt.title('HURDLER Success Rate vs Sequence Length')
plt.legend(title="Plasmid")
plt.grid(True)
plt.tight_layout()

output_path = os.path.join(output_dir, "hurdler_success_rate_vs_sequence_length.pdf")
plt.savefig(output_path, format='pdf', dpi=600)
print(f"✓ Main plot saved to {output_path}")
plt.close()

# Statistics table
stats_list = []
for plasmid in plasmids:
    df = plasmid_dfs[plasmid]
    stats = {
        'Plasmid': plasmid,
        'Min (%)': df['Probability'].min(),
        'Max (%)': df['Probability'].max(),
        'Mean (%)': df['Probability'].mean(),
        '@20AA (%)': df[df['Sequence_Length'] == 20]['Probability'].values[0] if 20 in df['Sequence_Length'].values else np.nan,
        '@40AA (%)': df[df['Sequence_Length'] == 40]['Probability'].values[0] if 40 in df['Sequence_Length'].values else np.nan,
        '@60AA (%)': df[df['Sequence_Length'] == 60]['Probability'].values[0] if 60 in df['Sequence_Length'].values else np.nan
    }
    stats_list.append(stats)

df_stats = pd.DataFrame(stats_list)
df_stats = df_stats.sort_values('@60AA (%)', ascending=False).reset_index(drop=True)

stats_path = os.path.join(output_dir, "hurdler_success_rate_statistics.csv")
df_stats.to_csv(stats_path, index=False)
print(f"✓ Statistics saved to {stats_path}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nFiles generated:")
print(f"  1. {output_dir}/hurdler_three_site_combinations_df1.csv")
print(f"  2. {output_dir}/hurdler_three_site_combinations_df2_optimized.csv")
print(f"  3. {output_dir}/hurdler_lookup_optimized.pkl")
print(f"  4. {output_dir}/hurdler_success_rate_results.csv")
print(f"  5. {output_dir}/hurdler_success_rate_summary.csv")
print(f"  6. {output_dir}/hurdler_success_rate_vs_sequence_length.pdf")
print(f"  7. {output_dir}/hurdler_success_rate_statistics.csv")

print("\nStatistics Summary:")
print(df_stats.to_string(index=False))
