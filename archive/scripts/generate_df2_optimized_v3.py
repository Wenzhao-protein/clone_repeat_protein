#!/usr/bin/env python
"""
优化版本3：使用更节省内存的方式生成df2
直接从原始数据流式处理，避免存储大型lookup字典
"""
import pandas as pd
import numpy as np
import sys
import gc

print("="*80)
print("OPTIMIZED DF2 GENERATION (V3 - MEMORY EFFICIENT)")
print("="*80)

# 1. 加载codon usage
print("\n1. Loading E.coli codon usage...")
codon_freq = pd.read_csv('./utils/output/codon_usage.csv')
codon_freq_dict = dict(zip(codon_freq['codon'], codon_freq['frequency']))
print(f"   Loaded {len(codon_freq_dict)} codon frequencies")

def calculate_codon_usage_freq(dna_seq, codon_dict):
    """计算9mer DNA序列的codon usage频率乘积"""
    if len(dna_seq) != 9:
        return 0.0
    freq = 1.0
    for i in range(0, 9, 3):
        codon = dna_seq[i:i+3]
        freq *= codon_dict.get(codon, 0.0)
    return freq

def has_stop_codon(aa_seq):
    """检查氨基酸序列是否包含终止密码子"""
    return '*' in aa_seq

# 2. 加载df1
print("\n2. Loading df1...")
df1 = pd.read_csv('./output/hurdler_three_site_combinations_df1.csv')
print(f"   df1 shape: {df1.shape}")

# 3. 加载RE数据
print("\n3. Loading restriction enzyme data...")
df_seamless = pd.read_csv('./utils/output/seamless_insert.csv')
df_silent = pd.read_csv('./utils/output/slient_mutation.csv')  # 注意：文件名拼写是slient

# 4. 预处理Site I
print("\n4. Pre-processing Site I (seamless insert)...")
print("   - Filtering stop codons (*)")
df_seamless = df_seamless[~df_seamless['re_3aa_site'].str.contains('\\*', na=False)].copy()
print(f"   After stop codon filter: {len(df_seamless)} rows")

print("   - Calculating codon usage frequency")
df_seamless['codon_usage_freq'] = df_seamless['re_9bp_site'].apply(
    lambda x: calculate_codon_usage_freq(x, codon_freq_dict)
)

print("   - Keeping best codon usage for duplicate (enzyme, 3mer_AA)")
df_seamless = df_seamless.sort_values('codon_usage_freq', ascending=False).groupby(
    ['name', 're_3aa_site']
).first().reset_index()
print(f"   After deduplication: {len(df_seamless)} unique pairs")

# 选择需要的列
site_i_data = df_seamless[['name', 're_3aa_site', 're_9bp_site', 'codon_usage_freq']].rename(columns={
    'name': 'enzyme',
    're_3aa_site': '3mer_aa',
    're_9bp_site': 'dna_seq'
})

# 5. 预处理Site II
print("\n5. Pre-processing Site II (silent mutation)...")
print("   - Filtering stop codons (*)")
df_silent = df_silent[~df_silent['re_3aa_site'].str.contains('\\*', na=False)].copy()
print(f"   After stop codon filter: {len(df_silent)} rows")

print("   - Calculating mutated codon usage frequency")
df_silent['mutated_codon_usage_freq'] = df_silent['re_9bp_site_mutated'].apply(
    lambda x: calculate_codon_usage_freq(x, codon_freq_dict)
)

print("   - Keeping best mutated codon usage for duplicate (enzyme, 3mer_AA)")
df_silent = df_silent.sort_values('mutated_codon_usage_freq', ascending=False).groupby(
    ['name', 're_3aa_site']
).first().reset_index()
print(f"   After deduplication: {len(df_silent)} unique pairs")

# 选择需要的列
site_ii_data = df_silent[['name', 're_3aa_site', 're_9bp_site', 
                          're_9bp_site_mutated', 'mutated_codon_usage_freq']].rename(columns={
    'name': 'enzyme',
    're_3aa_site': '3mer_aa',
    're_9bp_site': 'dna_seq_original',
    're_9bp_site_mutated': 'dna_seq_mutated',
    'mutated_codon_usage_freq': 'codon_usage_freq'
})

# 释放内存
del df_seamless, df_silent
gc.collect()

# 6. 流式处理生成df2
print("\n6. Streaming df1 to df2...")
print("   Note: Site III (Type IIS) doesn't encode 3mer AA")

output_path = './output/hurdler_three_site_combinations_df2_optimized.csv'
batch_size = 100  # 减小批次大小
total_rows = 0
first_batch = True

num_batches = (len(df1) + batch_size - 1) // batch_size
print(f"   Processing {len(df1)} df1 rows in {num_batches} batches...")
sys.stdout.flush()

for batch_idx in range(0, len(df1), batch_size):
    batch_end = min(batch_idx + batch_size, len(df1))
    df1_batch = df1.iloc[batch_idx:batch_end]
    
    batch_num = batch_idx // batch_size + 1
    
    # 对每个df1行，分别merge
    expanded_list = []
    
    for _, row in df1_batch.iterrows():
        # 获取该行对应的Site I和Site II数据
        site_i_filtered = site_i_data[site_i_data['enzyme'] == row['site_i_enzyme']]
        site_ii_filtered = site_ii_data[site_ii_data['enzyme'] == row['site_ii_enzyme']]
        
        if len(site_i_filtered) == 0 or len(site_ii_filtered) == 0:
            continue
        
        # 创建临时行数据
        temp_df = pd.DataFrame([row] * len(site_i_filtered))
        temp_df = temp_df.reset_index(drop=True)
        site_i_copy = site_i_filtered.reset_index(drop=True)
        
        # 添加Site I数据
        temp_df['site_i_3mer_aa'] = site_i_copy['3mer_aa'].values
        temp_df['site_i_dna'] = site_i_copy['dna_seq'].values
        temp_df['site_i_codon_usage_freq'] = site_i_copy['codon_usage_freq'].values
        
        # 扩展到Site II
        temp_expanded = []
        for _, temp_row in temp_df.iterrows():
            for _, site_ii_row in site_ii_filtered.iterrows():
                new_row = temp_row.to_dict()
                new_row['site_ii_3mer_aa'] = site_ii_row['3mer_aa']
                new_row['site_ii_dna_original'] = site_ii_row['dna_seq_original']
                new_row['site_ii_dna_mutated'] = site_ii_row['dna_seq_mutated']
                new_row['site_ii_codon_usage_freq'] = site_ii_row['codon_usage_freq']
                new_row['site_iii_3mer_aa'] = 'N/A'
                new_row['site_iii_dna'] = 'N/A'
                temp_expanded.append(new_row)
        
        expanded_list.extend(temp_expanded)
    
    # 写入批次
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
    
    # 打印进度
    avg_per_df1 = total_rows / batch_end if batch_end > 0 else 0
    print(f"   Batch {batch_num}/{num_batches}: {total_rows:,} rows total (avg {avg_per_df1:.1f} per df1 row)")
    sys.stdout.flush()

print(f"\n7. Generated {total_rows:,} total combinations")
print(f"✓ df2 saved to: {output_path}")

# 8. 统计信息
print("\n8. Statistics:")
print("="*80)

if total_rows > 0:
    df2_sample = pd.read_csv(output_path, nrows=min(1000, total_rows))
    
    print(f"Total rows: {total_rows:,}")
    print(f"\nSample statistics (first {len(df2_sample)} rows):")
    print(f"Unique 3mer AA:")
    print(f"  Site I: {df2_sample['site_i_3mer_aa'].nunique()}")
    print(f"  Site II: {df2_sample['site_ii_3mer_aa'].nunique()}")
    
    print(f"\nCodon usage frequency range:")
    print(f"  Site I: {df2_sample['site_i_codon_usage_freq'].min():.6f} - {df2_sample['site_i_codon_usage_freq'].max():.6f}")
    print(f"  Site II: {df2_sample['site_ii_codon_usage_freq'].min():.6f} - {df2_sample['site_ii_codon_usage_freq'].max():.6f}")

print("\n" + "="*80)
print("OPTIMIZED DF2 GENERATION COMPLETE")
print("="*80)
