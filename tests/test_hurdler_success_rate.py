#!/usr/bin/env python3
"""
计算HURDLER成功率：针对不同plasmid和module length
"""

import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from tqdm import tqdm
import random
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

# Historical executable analysis script, not a pytest module.
__test__ = False

def generate_random_aa_sequence(length):
    """生成随机氨基酸序列"""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(random.choice(amino_acids) for _ in range(length))

def extract_3mer_pairs(sequence):
    """提取序列中所有不重复的3mer AA"""
    # 处理circular sequence边界
    extended_seq = sequence[-2:] + sequence + sequence[:2]
    
    # 提取所有3mer
    three_mers = set()
    for i in range(len(sequence)):
        three_mer = extended_seq[i:i+3]
        if len(three_mer) == 3:
            three_mers.add(three_mer)
    
    return three_mers

def check_hurdler_success(sequence, hurdler_lookup, plasmid):
    """检查一个序列是否能找到HURDLER方案
    
    Args:
        sequence: 氨基酸序列
        hurdler_lookup: lookup字典
        plasmid: 质粒名称（已包含_compatible后缀）
    """
    three_mers = extract_3mer_pairs(sequence)
    
    if plasmid not in hurdler_lookup:
        return False, len(three_mers)
    
    plasmid_lookup = hurdler_lookup[plasmid]
    
    # 检查是否有任何3mer对在lookup中
    # lookup中存储的是(site_i_3mer, site_ii_3mer)对
    # 我们需要检查序列中的任意两个3mer是否在lookup中
    three_mer_list = list(three_mers)
    for i in range(len(three_mer_list)):
        for j in range(i+1, len(three_mer_list)):
            # 尝试两种顺序
            pair1 = (three_mer_list[i], three_mer_list[j])
            pair2 = (three_mer_list[j], three_mer_list[i])
            
            if pair1 in plasmid_lookup or pair2 in plasmid_lookup:
                return True, len(three_mers)
    
    return False, len(three_mers)

def test_success_rate(module_lengths, n_tests_per_length, hurdler_lookup, plasmids, output_dir='./output'):
    """测试成功率"""
    
    results = []
    
    print("\n" + "="*80)
    print("HURDLER SUCCESS RATE TESTING")
    print("="*80)
    print(f"\nModule lengths: {module_lengths}")
    print(f"Tests per length: {n_tests_per_length}")
    print(f"Plasmids: {plasmids}")
    
    for length in tqdm(module_lengths, desc="Testing lengths"):
        for test_idx in range(n_tests_per_length):
            # 生成随机序列
            sequence = generate_random_aa_sequence(length)
            
            # 测试每个plasmid
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
    
    # 保存原始结果
    results_path = Path(output_dir) / 'hurdler_success_rate_results.csv'
    df_results.to_csv(results_path, index=False)
    print(f"\n✓ Raw results saved to: {results_path}")
    
    return df_results

def analyze_results(df_results, output_dir='./output'):
    """分析结果并生成报告"""
    
    print("\n" + "="*80)
    print("RESULTS ANALYSIS")
    print("="*80)
    
    # 按module length和plasmid分组计算成功率
    success_rate = df_results.groupby(['module_length', 'plasmid'])['success'].agg([
        ('success_rate', 'mean'),
        ('n_success', 'sum'),
        ('n_total', 'count')
    ]).reset_index()
    
    success_rate['success_rate_pct'] = success_rate['success_rate'] * 100
    
    # 保存汇总结果
    summary_path = Path(output_dir) / 'hurdler_success_rate_summary.csv'
    success_rate.to_csv(summary_path, index=False)
    print(f"\n✓ Summary saved to: {summary_path}")
    
    # 打印结果
    print("\nSuccess Rate by Module Length and Plasmid:")
    print("="*80)
    for plasmid in success_rate['plasmid'].unique():
        print(f"\n{plasmid}:")
        plasmid_data = success_rate[success_rate['plasmid'] == plasmid]
        for _, row in plasmid_data.iterrows():
            print(f"  Length {row['module_length']:2d}: {row['success_rate_pct']:6.2f}% "
                  f"({int(row['n_success'])}/{int(row['n_total'])})")
    
    # 生成可视化
    print("\n" + "="*80)
    print("GENERATING PLOTS")
    print("="*80)
    
    # 1. 成功率 vs Module Length
    plt.figure(figsize=(12, 8))
    
    for plasmid in success_rate['plasmid'].unique():
        plasmid_data = success_rate[success_rate['plasmid'] == plasmid].sort_values('module_length')
        plt.plot(plasmid_data['module_length'], plasmid_data['success_rate_pct'], 
                marker='o', label=plasmid, linewidth=2)
    
    plt.xlabel('Module Length (amino acids)', fontsize=12)
    plt.ylabel('Success Rate (%)', fontsize=12)
    plt.title('HURDLER Success Rate vs Module Length', fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plot_path = Path(output_dir) / 'hurdler_success_rate_plot.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved to: {plot_path}")
    plt.close()
    
    # 2. Heatmap
    pivot_data = success_rate.pivot(index='module_length', columns='plasmid', values='success_rate_pct')
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(pivot_data, annot=True, fmt='.1f', cmap='YlGnBu', cbar_kws={'label': 'Success Rate (%)'})
    plt.title('HURDLER Success Rate Heatmap', fontsize=14, fontweight='bold')
    plt.xlabel('Plasmid', fontsize=12)
    plt.ylabel('Module Length (amino acids)', fontsize=12)
    plt.tight_layout()
    
    heatmap_path = Path(output_dir) / 'hurdler_success_rate_heatmap.png'
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    print(f"✓ Heatmap saved to: {heatmap_path}")
    plt.close()
    
    return success_rate

def main():
    # 参数设置
    module_lengths = list(range(4, 61))  # 4-60
    n_tests_per_length = 1000
    output_dir = './output'
    
    # 加载lookup dictionary
    print("\n" + "="*80)
    print("LOADING HURDLER LOOKUP DICTIONARY")
    print("="*80)
    
    lookup_path = Path(output_dir) / 'hurdler_lookup_optimized.pkl'
    
    if not lookup_path.exists():
        print(f"ERROR: {lookup_path} not found!")
        print("Please run create_lookup_from_optimized_df2.py first.")
        return
    
    with open(lookup_path, 'rb') as f:
        hurdler_lookup = pickle.load(f)
    
    plasmids = list(hurdler_lookup.keys())
    print(f"\nLoaded lookup for {len(plasmids)} plasmids:")
    for plasmid in plasmids:
        n_pairs = len(hurdler_lookup[plasmid])
        print(f"  - {plasmid}: {n_pairs:,} unique 3mer AA pairs")
    
    # 测试成功率
    df_results = test_success_rate(
        module_lengths=module_lengths,
        n_tests_per_length=n_tests_per_length,
        hurdler_lookup=hurdler_lookup,
        plasmids=plasmids,
        output_dir=output_dir
    )
    
    # 分析结果
    success_rate = analyze_results(df_results, output_dir=output_dir)
    
    print("\n" + "="*80)
    print("SUCCESS RATE TESTING COMPLETE")
    print("="*80)

if __name__ == '__main__':
    # 设置随机种子以确保可重复性
    random.seed(42)
    np.random.seed(42)
    
    main()
