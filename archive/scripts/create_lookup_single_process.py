#!/usr/bin/env python3
"""
单进程版本 - 逐行处理df2，创建lookup dictionary
"""

import pandas as pd
import pickle
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
import gc

print("="*80)
print("CREATING HURDLER LOOKUP DICTIONARY (Single Process)")
print("="*80)

# 1. 读取df2（分块）
print("\n1. Analyzing df2...")
df2_path = './output/hurdler_three_site_combinations_df2.csv'

# 先读取一小部分获取列名
df_sample = pd.read_csv(df2_path, nrows=100)
plasmid_cols = [col for col in df_sample.columns if '_compatible' in col]

print(f"   Found {len(plasmid_cols)} plasmids")

# 2. 逐块读取并构建lookup
print("\n2. Building lookup dictionary...")
print("   Processing in chunks...")

hurdler_lookup = defaultdict(lambda: defaultdict(set))
chunk_size = 100000
processed = 0

reader = pd.read_csv(df2_path, chunksize=chunk_size)

for chunk_idx, df_chunk in enumerate(reader):
    for _, row in df_chunk.iterrows():
        # 提取3mer AA
        site_i_aa = row['site_i_3mer_aa']
        site_ii_aa = row['site_ii_3mer_aa']
        
        # 跳过NaN
        if pd.isna(site_i_aa) or pd.isna(site_ii_aa):
            continue
        
        # 创建3mer AA对
        aa_pair = frozenset([site_i_aa, site_ii_aa])
        
        # 为每个兼容的plasmid添加
        for plasmid in plasmid_cols:
            if row[plasmid]:
                # 只存储是否存在（set），不存储详细enzyme信息
                hurdler_lookup[plasmid][aa_pair].add(True)
    
    processed += len(df_chunk)
    if (chunk_idx + 1) % 10 == 0:
        print(f"   Processed: {processed:,} rows ({processed/40980196*100:.1f}%)")
        gc.collect()

print(f"\n   Total processed: {processed:,} rows")

# 3. 转换为简化格式（只保留keys）
print("\n3. Converting to lite format (keys only)...")
hurdler_lookup_lite = {
    plasmid: set(aa_dict.keys())
    for plasmid, aa_dict in hurdler_lookup.items()
}

# 4. 统计
print("\n4. Statistics:")
print("="*80)
for plasmid in plasmid_cols:
    plasmid_name = plasmid.replace('_compatible', '')
    if plasmid in hurdler_lookup_lite:
        n_pairs = len(hurdler_lookup_lite[plasmid])
        print(f"{plasmid_name}: {n_pairs:,} unique 3mer AA pairs")

# 5. 保存
print("\n5. Saving...")
lite_path = Path('./output') / 'hurdler_lookup_lite.pkl'

with open(lite_path, 'wb') as f:
    pickle.dump(hurdler_lookup_lite, f, protocol=pickle.HIGHEST_PROTOCOL)

size_mb = lite_path.stat().st_size / 1024 / 1024
print(f"✓ Saved to: {lite_path} ({size_mb:.1f} MB)")

print("\n" + "="*80)
print("LOOKUP DICTIONARY CREATION COMPLETE")
print("="*80)
