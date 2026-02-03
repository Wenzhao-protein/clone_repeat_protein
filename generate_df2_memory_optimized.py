#!/usr/bin/env python3
"""
内存优化的df2生成 - 分批处理并直接写入文件
"""

import pandas as pd
from tqdm import tqdm
from pathlib import Path
import gc

print("="*80)
print("MEMORY-OPTIMIZED DF2 GENERATION")
print("="*80)

# 1. 读取df1
print("\n1. Loading df1...")
df1 = pd.read_csv('./output/hurdler_three_site_combinations_df1.csv')
print(f"   df1 shape: {df1.shape}")

# 2. 创建lookup dictionaries
print("\n2. Creating lookup dictionaries...")
df_seamless = pd.read_csv('./utils/output/restriction_enzyme_seamless_insert.csv')
df_silent = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')

# Site I lookup
site_i_lookup = df_seamless.groupby('name').apply(
    lambda x: x[['re_site_shifted_tl', 're_site_shifted']].values.tolist()
).to_dict()

# Site II lookup
site_ii_lookup = df_silent.groupby('name').apply(
    lambda x: x[['re_site_shifted_tl', 're_site_mutate_shifted']].values.tolist()
).to_dict()

print(f"   Site I lookup: {len(site_i_lookup)} enzymes")
print(f"   Site II lookup: {len(site_ii_lookup)} enzymes")

# 释放不需要的数据
del df_seamless, df_silent
gc.collect()

# 3. 分批扩展并写入
output_path = './output/hurdler_three_site_combinations_df2.csv'
batch_size = 500
total_rows = 0
first_batch = True

print(f"\n3. Expanding in batches of {batch_size}...")
print("   Note: Site III (Type IIS) doesn't encode 3mer AA")

for batch_start in range(0, len(df1), batch_size):
    batch_end = min(batch_start + batch_size, len(df1))
    df1_batch = df1.iloc[batch_start:batch_end]
    
    expanded_rows = []
    
    for _, row in df1_batch.iterrows():
        site_i = row['site_i_enzyme']
        site_ii = row['site_ii_enzyme']
        
        i_options = site_i_lookup.get(site_i, [])
        ii_options = site_ii_lookup.get(site_ii, [])
        
        if not i_options or not ii_options:
            continue
        
        for i_opt in i_options:
            for ii_opt in ii_options:
                new_row = row.to_dict()
                new_row.update({
                    'site_i_3mer_aa': i_opt[0],
                    'site_i_dna': i_opt[1],
                    'site_ii_3mer_aa': ii_opt[0],
                    'site_ii_dna': ii_opt[1],
                    'site_iii_3mer_aa': 'N/A',
                    'site_iii_dna': 'N/A'
                })
                expanded_rows.append(new_row)
    
    # 写入文件
    df_batch = pd.DataFrame(expanded_rows)
    
    if first_batch:
        df_batch.to_csv(output_path, index=False, mode='w')
        first_batch = False
    else:
        df_batch.to_csv(output_path, index=False, mode='a', header=False)
    
    total_rows += len(df_batch)
    print(f"   Batch {batch_start//batch_size + 1}: {len(df_batch):,} rows (total: {total_rows:,})")
    
    # 释放内存
    del expanded_rows, df_batch
    gc.collect()

print(f"\n4. Generated {total_rows:,} total combinations")
print(f"✓ df2 saved to: {output_path}")

# 5. 验证并统计
print("\n5. Verification...")
df2_sample = pd.read_csv(output_path, nrows=10)
print(f"   Columns: {df2_sample.columns.tolist()}")
print(f"   Sample rows: {len(df2_sample)}")
