"""
Test the optimized Step 3 lookup creation
Only processes Step 3 to test memory efficiency
"""

import pandas as pd
import pickle
import gc
import sys
from collections import defaultdict

print("="*80)
print("TESTING OPTIMIZED STEP 3")
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
    
    if chunk_num % 5 == 0:
        print(f"   Processed {total_rows:,} rows")
        sys.stdout.flush()
    
    del chunk
    gc.collect()

print(f"   ✓ Site I lookup size: {len(site_i_lookup):,} unique (enzyme, 3mer_aa) pairs")

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
    
    if chunk_num % 5 == 0:
        print(f"   Processed {total_rows:,} rows")
        sys.stdout.flush()
    
    del chunk
    gc.collect()

print(f"   ✓ Site II & III lookup size: {len(site_ii_iii_lookup):,} unique (enzyme_ii, 3mer_aa_ii, enzyme_iii) tuples")

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
    
    if chunk_num % 5 == 0:
        print(f"   Processed {total_rows:,} rows, added {entries_added:,} entries")
        sys.stdout.flush()
    
    del chunk
    gc.collect()

print(f"\n4. ✓ Total rows processed: {total_rows:,}")
print(f"   ✓ Total entries added: {entries_added:,}")

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

print("\n" + "="*80)
print("✓ LOOKUP DICTIONARY CREATED SUCCESSFULLY!")
print("="*80)
