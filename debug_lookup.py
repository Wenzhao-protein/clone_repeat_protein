#!/usr/bin/env python
"""
调试HURDLER lookup匹配问题
"""
import pickle

# 加载lookup
with open('./output/hurdler_lookup_optimized.pkl', 'rb') as f:
    lookup = pickle.load(f)

# 检查lookup数据结构
plasmid = 'pUC18_compatible'
pairs_sample = list(lookup[plasmid])[:10]

print("Lookup structure check:")
print(f"Plasmid: {plasmid}")
print(f"Total pairs: {len(lookup[plasmid])}")
print(f"\nFirst 10 pairs:")
for i, pair in enumerate(pairs_sample):
    print(f"  {i+1}. Type: {type(pair)}, Value: {pair}")

# 测试一个简单的序列
test_seq = "EAKI"
extended = test_seq[-2:] + test_seq + test_seq[:2]
print(f"\nTest sequence: {test_seq}")
print(f"Extended: {extended}")

# 提取3mer
three_mers = set()
for i in range(len(test_seq)):
    three_mer = extended[i:i+3]
    if len(three_mer) == 3:
        three_mers.add(three_mer)
        print(f"  3mer at {i}: {three_mer}")

print(f"\nUnique 3mers: {three_mers}")

# 创建3mer对
from itertools import combinations
three_mer_pairs = list(combinations(sorted(three_mers), 2))
print(f"\n3mer pairs:")
for pair in three_mer_pairs:
    print(f"  {pair} (type: {type(pair)})")
    
    # 检查是否在lookup中
    if pair in lookup[plasmid]:
        print(f"    -> FOUND in lookup!")
    else:
        print(f"    -> NOT found")

# 检查df2中实际存储的格式
print("\n" + "="*80)
print("Checking df2 for actual 3mer AA examples:")
import pandas as pd
df2_sample = pd.read_csv('./output/hurdler_three_site_combinations_df2_optimized.csv', nrows=100)
df2_uc18 = df2_sample[df2_sample['pUC18_compatible'] == True]

if len(df2_uc18) > 0:
    print(f"\nFound {len(df2_uc18)} pUC18-compatible rows in first 100")
    for idx, row in df2_uc18.head(3).iterrows():
        print(f"\nRow {idx}:")
        print(f"  Site I 3mer AA: {row['site_i_3mer_aa']}")
        print(f"  Site II 3mer AA: {row['site_ii_3mer_aa']}")
        pair = (row['site_i_3mer_aa'], row['site_ii_3mer_aa'])
        print(f"  Pair: {pair}")
        print(f"  In lookup? {pair in lookup[plasmid]}")
else:
    print("No pUC18-compatible rows in first 100")
