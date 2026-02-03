#!/usr/bin/env python3
"""
优化的df2生成脚本 - 从df1扩展添加3mer AA序列
"""

import pandas as pd
from tqdm import tqdm
from pathlib import Path

print("="*80)
print("GENERATING DF2 from DF1")
print("="*80)

# 1. 读取df1
print("\n1. Loading df1...")
df1 = pd.read_csv('./output/hurdler_three_site_combinations_df1.csv')
print(f"   df1 shape: {df1.shape}")

# 2. 创建lookup dictionaries
print("\n2. Creating lookup dictionaries...")
df_seamless = pd.read_csv('./utils/output/restriction_enzyme_seamless_insert.csv')
df_silent = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')

# Site I lookup (seamless insert)
site_i_lookup = {}
for _, row in df_seamless.iterrows():
    enzyme = row['name']
    if enzyme not in site_i_lookup:
        site_i_lookup[enzyme] = []
    site_i_lookup[enzyme].append({
        '3mer_aa': row['re_site_shifted_tl'],
        'dna_seq': row['re_site_shifted']
    })

print(f"   Site I lookup: {len(site_i_lookup)} enzymes")

# Site II lookup (silent mutation)
site_ii_lookup = {}
for _, row in df_silent.iterrows():
    enzyme = row['name']
    if enzyme not in site_ii_lookup:
        site_ii_lookup[enzyme] = []
    site_ii_lookup[enzyme].append({
        '3mer_aa': row['re_site_shifted_tl'],
        'dna_seq': row['re_site_mutate_shifted']
    })

print(f"   Site II lookup: {len(site_ii_lookup)} enzymes")

# 3. 扩展df1到df2
print("\n3. Expanding df1 to df2...")
print("   Note: Site III (Type IIS) doesn't encode 3mer AA")

expanded_rows = []
for idx, row in tqdm(df1.iterrows(), total=len(df1), desc="Processing"):
    site_i = row['site_i_enzyme']
    site_ii = row['site_ii_enzyme']
    
    # 获取Site I和Site II的3mer AA选项
    i_options = site_i_lookup.get(site_i, [])
    ii_options = site_ii_lookup.get(site_ii, [])
    
    if not i_options or not ii_options:
        continue  # 跳过没有3mer AA数据的组合
    
    # 生成所有组合
    for i_opt in i_options:
        for ii_opt in ii_options:
            new_row = row.to_dict()
            new_row.update({
                'site_i_3mer_aa': i_opt['3mer_aa'],
                'site_i_dna': i_opt['dna_seq'],
                'site_ii_3mer_aa': ii_opt['3mer_aa'],
                'site_ii_dna': ii_opt['dna_seq'],
                # Site III (Type IIS) 不需要3mer AA
                'site_iii_3mer_aa': 'N/A',
                'site_iii_dna': 'N/A'
            })
            expanded_rows.append(new_row)

df2 = pd.DataFrame(expanded_rows)
print(f"\n4. Generated {len(df2):,} combinations with 3mer AA sequences")

# 4. 保存df2
output_path = Path('./output') / 'hurdler_three_site_combinations_df2.csv'
df2.to_csv(output_path, index=False)
print(f"\n✓ df2 saved to: {output_path}")

# 5. 统计信息
print("\n" + "="*80)
print("STATISTICS")
print("="*80)
print(f"df1 combinations: {len(df1):,}")
print(f"df2 combinations: {len(df2):,}")
print(f"Expansion factor: {len(df2)/len(df1):.1f}x")
print(f"\nUnique 3mer AA pairs:")
print(f"  Site I: {df2['site_i_3mer_aa'].nunique()}")
print(f"  Site II: {df2['site_ii_3mer_aa'].nunique()}")
print(f"  Combined: {df2.groupby(['site_i_3mer_aa', 'site_ii_3mer_aa']).ngroups}")
