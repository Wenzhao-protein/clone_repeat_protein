#!/usr/bin/env python3
"""
从df2创建HURDLER lookup dictionary用于成功率测试
"""

import pandas as pd
import pickle
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
import multiprocessing as mp

def create_lookup_batch(args):
    """处理一批df2数据，创建lookup entries"""
    df_batch, plasmid_cols = args
    
    batch_lookup = defaultdict(lambda: defaultdict(list))
    
    for _, row in df_batch.iterrows():
        # 提取site I和II的3mer AA
        site_i_aa = row['site_i_3mer_aa']
        site_ii_aa = row['site_ii_3mer_aa']
        
        if site_i_aa == 'N/A' or site_ii_aa == 'N/A':
            continue
        
        # 创建3mer AA对的frozenset (无序)
        aa_pair = frozenset([site_i_aa, site_ii_aa])
        
        # 为每个plasmid记录兼容的组合
        for plasmid in plasmid_cols:
            if row[plasmid]:  # 如果该plasmid兼容
                batch_lookup[plasmid][aa_pair].append({
                    'site_i_enzyme': row['site_i_enzyme'],
                    'site_ii_enzyme': row['site_ii_enzyme'],
                    'site_iii_enzyme': row['site_iii_enzyme'],
                    'site_i_ovhg': row['site_i_ovhg'],
                    'site_ii_ovhg': row['site_ii_ovhg'],
                    'site_iii_ovhg': row['site_iii_ovhg']
                })
    
    return batch_lookup

def main():
    print("="*80)
    print("CREATING HURDLER LOOKUP DICTIONARY")
    print("="*80)
    
    # 1. 读取df2
    print("\n1. Loading df2...")
    df2_path = './output/hurdler_three_site_combinations_df2.csv'
    
    if not Path(df2_path).exists():
        print(f"   ERROR: {df2_path} not found!")
        print("   Please run generate_df2_memory_optimized.py first.")
        return
    
    df2 = pd.read_csv(df2_path)
    print(f"   df2 shape: {df2.shape}")
    
    # 2. 识别plasmid列
    plasmid_cols = [col for col in df2.columns if '_compatible' in col]
    print(f"\n2. Found {len(plasmid_cols)} plasmids:")
    for pc in plasmid_cols:
        plasmid_name = pc.replace('_compatible', '')
        compatible_count = df2[pc].sum()
        print(f"   - {plasmid_name}: {compatible_count:,} compatible combinations")
    
    # 3. 创建lookup dictionary
    print("\n3. Creating lookup dictionary...")
    print("   Structure: {plasmid: {frozenset(3mer_AA_pair): [enzyme_combinations]}}")
    
    # 使用多进程处理
    n_workers = min(mp.cpu_count() - 1, 8)
    batch_size = len(df2) // (n_workers * 4)
    
    print(f"   Using {n_workers} workers, batch size: {batch_size}")
    
    # 分批
    batches = []
    for i in range(0, len(df2), batch_size):
        df_batch = df2.iloc[i:i+batch_size]
        batches.append((df_batch, plasmid_cols))
    
    # 并行处理
    with mp.Pool(n_workers) as pool:
        batch_results = list(tqdm(
            pool.imap(create_lookup_batch, batches),
            total=len(batches),
            desc="Processing batches"
        ))
    
    # 4. 合并结果
    print("\n4. Merging results...")
    hurdler_lookup = defaultdict(lambda: defaultdict(list))
    
    for batch_dict in tqdm(batch_results, desc="Merging"):
        for plasmid, aa_dict in batch_dict.items():
            for aa_pair, combinations in aa_dict.items():
                hurdler_lookup[plasmid][aa_pair].extend(combinations)
    
    # 转换为普通dict
    hurdler_lookup = {
        plasmid: dict(aa_dict)
        for plasmid, aa_dict in hurdler_lookup.items()
    }
    
    # 5. 统计信息
    print("\n5. Statistics:")
    print("="*80)
    for plasmid in plasmid_cols:
        plasmid_name = plasmid.replace('_compatible', '')
        if plasmid in hurdler_lookup:
            n_aa_pairs = len(hurdler_lookup[plasmid])
            total_combinations = sum(len(v) for v in hurdler_lookup[plasmid].values())
            print(f"\n{plasmid_name}:")
            print(f"  - Unique 3mer AA pairs: {n_aa_pairs}")
            print(f"  - Total enzyme combinations: {total_combinations:,}")
            
            # 样例
            sample_pair = list(hurdler_lookup[plasmid].keys())[0]
            sample_combs = hurdler_lookup[plasmid][sample_pair]
            print(f"  - Example AA pair: {set(sample_pair)}")
            print(f"    → {len(sample_combs)} enzyme combinations")
    
    # 6. 保存
    output_path = Path('./output') / 'hurdler_lookup_dict.pkl'
    print(f"\n6. Saving to {output_path}...")
    
    with open(output_path, 'wb') as f:
        pickle.dump(hurdler_lookup, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    print(f"✓ Lookup dictionary saved ({output_path.stat().st_size / 1024 / 1024:.1f} MB)")
    
    # 7. 创建简化版本(只有keys)
    print("\n7. Creating lightweight version (keys only)...")
    hurdler_lookup_lite = {
        plasmid: set(aa_dict.keys())
        for plasmid, aa_dict in hurdler_lookup.items()
    }
    
    lite_path = Path('./output') / 'hurdler_lookup_lite.pkl'
    with open(lite_path, 'wb') as f:
        pickle.dump(hurdler_lookup_lite, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    print(f"✓ Lite version saved ({lite_path.stat().st_size / 1024 / 1024:.1f} MB)")
    print("\n" + "="*80)
    print("LOOKUP DICTIONARY CREATION COMPLETE")
    print("="*80)

if __name__ == '__main__':
    main()
