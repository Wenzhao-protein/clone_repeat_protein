#!/usr/bin/env python3
"""
优化的df2生成 - 添加终止密码子过滤和codon usage优化
"""

import pandas as pd
from tqdm import tqdm
from pathlib import Path
import gc

def load_codon_usage():
    """加载E.coli codon usage frequency"""
    try:
        df_codon = pd.read_csv('./utils/output/codon_usage.csv')
        # 创建codon到frequency的映射
        codon_freq = {}
        for _, row in df_codon.iterrows():
            codon_freq[row['codon'].upper()] = row['frequency']
        return codon_freq
    except Exception as e:
        print(f"Warning: Could not load codon usage data: {e}")
        return {}

def calculate_codon_usage_freq(dna_seq, codon_freq):
    """计算9mer DNA序列的codon usage frequency（3个密码子的乘积）"""
    if len(dna_seq) != 9:
        return 0.0
    
    freq = 1.0
    for i in range(0, 9, 3):
        codon = dna_seq[i:i+3].upper()
        if codon in codon_freq:
            freq *= codon_freq[codon]
        else:
            return 0.0  # 未知密码子
    
    return freq

def has_stop_codon(aa_seq):
    """检查氨基酸序列是否包含终止密码子"""
    return '*' in str(aa_seq)

print("="*80)
print("OPTIMIZED DF2 GENERATION")
print("="*80)

# 1. 加载codon usage
print("\n1. Loading E.coli codon usage...")
codon_freq = load_codon_usage()
print(f"   Loaded {len(codon_freq)} codon frequencies")

# 2. 读取df1
print("\n2. Loading df1...")
df1 = pd.read_csv('./output/hurdler_three_site_combinations_df1.csv')
print(f"   df1 shape: {df1.shape}")

# 3. 读取原始数据
print("\n3. Loading restriction enzyme data...")
df_seamless = pd.read_csv('./utils/output/restriction_enzyme_seamless_insert.csv')
df_silent = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')

# 4. 预处理Site I数据
print("\n4. Pre-processing Site I (seamless insert)...")
print("   - Filtering stop codons (*)")
print("   - Selecting best codon usage for duplicate 3mer AA")

# 过滤终止密码子
df_seamless = df_seamless[~df_seamless['re_site_shifted_tl'].str.contains('\\*', na=False)]
print(f"   After stop codon filter: {len(df_seamless)} rows")

# 计算codon usage frequency
df_seamless['codon_usage_freq'] = df_seamless['re_site_shifted'].apply(
    lambda x: calculate_codon_usage_freq(x, codon_freq)
)

# 按enzyme和3mer AA分组，保留最高codon usage的
site_i_best = df_seamless.sort_values('codon_usage_freq', ascending=False).groupby(
    ['name', 're_site_shifted_tl']
).first().reset_index()

print(f"   After deduplication: {len(site_i_best)} unique (enzyme, 3mer_AA) pairs")

# 创建Site I lookup
site_i_lookup = {}
for _, row in site_i_best.iterrows():
    enzyme = row['name']
    if enzyme not in site_i_lookup:
        site_i_lookup[enzyme] = []
    site_i_lookup[enzyme].append({
        '3mer_aa': row['re_site_shifted_tl'],
        'dna_seq': row['re_site_shifted'],
        'codon_usage_freq': row['codon_usage_freq']
    })

print(f"   Site I lookup: {len(site_i_lookup)} enzymes")

# 5. 预处理Site II数据
print("\n5. Pre-processing Site II (silent mutation)...")
print("   - Filtering stop codons (*)")
print("   - Selecting best codon usage for duplicate 3mer AA")

# 过滤终止密码子
df_silent = df_silent[~df_silent['re_site_shifted_tl'].str.contains('\\*', na=False)]
print(f"   After stop codon filter: {len(df_silent)} rows")

# 计算突变后序列的codon usage frequency
df_silent['mutated_codon_usage_freq'] = df_silent['re_site_mutate_shifted'].apply(
    lambda x: calculate_codon_usage_freq(x, codon_freq)
)

# 按enzyme和3mer AA分组，保留最高mutated codon usage的
# 保留突变前后的序列
site_ii_best = df_silent.sort_values('mutated_codon_usage_freq', ascending=False).groupby(
    ['name', 're_site_shifted_tl']
).first().reset_index()

print(f"   After deduplication: {len(site_ii_best)} unique (enzyme, 3mer_AA) pairs")

# 创建Site II lookup
site_ii_lookup = {}
for _, row in site_ii_best.iterrows():
    enzyme = row['name']
    if enzyme not in site_ii_lookup:
        site_ii_lookup[enzyme] = []
    site_ii_lookup[enzyme].append({
        '3mer_aa': row['re_site_shifted_tl'],
        'dna_seq_original': row['re_site_shifted'],  # 突变前
        'dna_seq_mutated': row['re_site_mutate_shifted'],  # 突变后
        'codon_usage_freq': row['mutated_codon_usage_freq']
    })

print(f"   Site II lookup: {len(site_ii_lookup)} enzymes")

# 释放内存
del df_seamless, df_silent, site_i_best, site_ii_best
gc.collect()

# 6. 分批扩展df1到df2
print("\n6. Expanding df1 to df2...")
print("   Note: Site III (Type IIS) doesn't encode 3mer AA")

output_path = './output/hurdler_three_site_combinations_df2_optimized.csv'
batch_size = 500
total_rows = 0
first_batch = True

num_batches = (len(df1) + batch_size - 1) // batch_size
print(f"   Processing {len(df1)} df1 rows in {num_batches} batches...")

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
                    'site_i_3mer_aa': i_opt['3mer_aa'],
                    'site_i_dna': i_opt['dna_seq'],
                    'site_i_codon_usage_freq': i_opt['codon_usage_freq'],
                    'site_ii_3mer_aa': ii_opt['3mer_aa'],
                    'site_ii_dna_original': ii_opt['dna_seq_original'],
                    'site_ii_dna_mutated': ii_opt['dna_seq_mutated'],
                    'site_ii_codon_usage_freq': ii_opt['codon_usage_freq'],
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
    
    batch_num = batch_start // batch_size + 1
    total_rows += len(df_batch)
    
    # 更频繁的进度更新
    if batch_num % 1 == 0:
        avg_per_df1 = total_rows / batch_end if batch_end > 0 else 0
        print(f"   Batch {batch_num}/{num_batches}: {total_rows:,} rows total (avg {avg_per_df1:.1f} per df1 row)")
        sys.stdout.flush()
    
    del expanded_rows, df_batch
    gc.collect()

print(f"\n7. Generated {total_rows:,} total combinations")
print(f"✓ df2 saved to: {output_path}")

# 8. 统计信息
print("\n8. Statistics:")
print("="*80)

if total_rows > 0:
    df2_sample = pd.read_csv(output_path, nrows=1000)
    
    print(f"Unique 3mer AA:")
    print(f"  Site I: {df2_sample['site_i_3mer_aa'].nunique()}")
    print(f"  Site II: {df2_sample['site_ii_3mer_aa'].nunique()}")
    
    print(f"\nCodon usage frequency range:")
    print(f"  Site I: {df2_sample['site_i_codon_usage_freq'].min():.6f} - {df2_sample['site_i_codon_usage_freq'].max():.6f}")
    print(f"  Site II: {df2_sample['site_ii_codon_usage_freq'].min():.6f} - {df2_sample['site_ii_codon_usage_freq'].max():.6f}")

print("\n" + "="*80)
print("OPTIMIZED DF2 GENERATION COMPLETE")
print("="*80)
