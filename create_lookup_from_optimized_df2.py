#!/usr/bin/env python
"""
从优化的df2创建lookup字典
只存储3mer AA对，用于快速检查成功率
"""
import pandas as pd
import pickle
import sys
from collections import defaultdict

print("="*80)
print("CREATING LOOKUP DICTIONARY FROM OPTIMIZED DF2")
print("="*80)

# 读取df2（分块读取以节省内存）
print("\n1. Reading df2 in chunks...")
lookup = defaultdict(set)  # {plasmid: {(3mer_aa_1, 3mer_aa_2)}}

chunk_size = 1000000
total_rows = 0
chunk_num = 0

plasmid_cols = [
    'pGEX-4T-1_compatible',
    'pMAL-c5X_compatible',
    'pET-21a(+)_compatible',
    'pET-28a(+)_compatible',
    'pET-28a(+)_start_codon_compatible',
    'pCold_I_compatible',
    'pUC18_compatible',
    'pQE-3_compatible'
]

for chunk in pd.read_csv('./output/hurdler_three_site_combinations_df2_optimized.csv', 
                         chunksize=chunk_size,
                         usecols=['site_i_3mer_aa', 'site_ii_3mer_aa'] + plasmid_cols):
    chunk_num += 1
    total_rows += len(chunk)
    
    # 对每个质粒，添加3mer AA对
    for plasmid in plasmid_cols:
        # 只保留兼容该质粒的组合
        compatible = chunk[chunk[plasmid] == True]
        
        if len(compatible) > 0:
            # 添加所有3mer AA对到set中
            for _, row in compatible.iterrows():
                pair = (row['site_i_3mer_aa'], row['site_ii_3mer_aa'])
                lookup[plasmid].add(pair)
    
    if chunk_num % 10 == 0:
        print(f"   Processed {total_rows:,} rows ({chunk_num} chunks)")
        sys.stdout.flush()

print(f"\n2. Total rows processed: {total_rows:,}")

# 统计信息
print("\n3. Lookup statistics:")
print("="*80)
for plasmid in plasmid_cols:
    plasmid_name = plasmid.replace('_compatible', '').replace('_start_codon', '')
    num_pairs = len(lookup[plasmid])
    print(f"  {plasmid_name:30s}: {num_pairs:,} unique 3mer AA pairs")

# 保存lookup字典
output_path = './output/hurdler_lookup_optimized.pkl'
print(f"\n4. Saving lookup dictionary to {output_path}...")
with open(output_path, 'wb') as f:
    pickle.dump(dict(lookup), f)

print("\n" + "="*80)
print("LOOKUP DICTIONARY CREATED SUCCESSFULLY")
print("="*80)
