#!/usr/bin/env python
"""
深入调试匹配问题
"""
import pickle

# 加载lookup
with open('./output/hurdler_lookup_optimized.pkl', 'rb') as f:
    lookup = pickle.load(f)

plasmid = 'pUC18_compatible'

# 获取lookup中的一个实际3mer对
sample_pair = list(lookup[plasmid])[0]
print(f"Sample pair from lookup: {sample_pair}")
print(f"Type: {type(sample_pair)}")
print(f"site_i: {sample_pair[0]}, site_ii: {sample_pair[1]}")

# 创建一个包含这两个3mer的模拟序列
site_i_3mer = sample_pair[0]
site_ii_3mer = sample_pair[1]

# 简单拼接创建一个序列
test_seq = site_i_3mer + site_ii_3mer
print(f"\nTest sequence: {test_seq}")
print(f"Length: {len(test_seq)}")

# 提取3mer
extended = test_seq[-2:] + test_seq + test_seq[:2]
print(f"Extended: {extended}")

three_mers = set()
for i in range(len(test_seq)):
    three_mer = extended[i:i+3]
    if len(three_mer) == 3:
        three_mers.add(three_mer)
        print(f"  3mer at {i}: {three_mer}")

print(f"\nExtracted 3mers: {three_mers}")
print(f"Contains site_i ({site_i_3mer})? {site_i_3mer in three_mers}")
print(f"Contains site_ii ({site_ii_3mer})? {site_ii_3mer in three_mers}")

# 检查匹配
found = False
three_mer_list = list(three_mers)
for i in range(len(three_mer_list)):
    for j in range(i+1, len(three_mer_list)):
        pair1 = (three_mer_list[i], three_mer_list[j])
        pair2 = (three_mer_list[j], three_mer_list[i])
        
        if pair1 in lookup[plasmid]:
            print(f"\n✓ FOUND: {pair1}")
            found = True
            break
        if pair2 in lookup[plasmid]:
            print(f"\n✓ FOUND: {pair2}")
            found = True
            break
    if found:
        break

if not found:
    print(f"\n✗ NOT FOUND - checking all combinations:")
    for i in range(len(three_mer_list)):
        for j in range(i+1, len(three_mer_list)):
            pair1 = (three_mer_list[i], three_mer_list[j])
            pair2 = (three_mer_list[j], three_mer_list[i])
            print(f"  {pair1}: {pair1 in lookup[plasmid]}")
            print(f"  {pair2}: {pair2 in lookup[plasmid]}")

# 直接检查sample_pair
print(f"\nDirect check of original pair:")
print(f"  {sample_pair} in lookup? {sample_pair in lookup[plasmid]}")

# 检查df2中这个对的实际情况
import pandas as pd
print("\n" + "="*80)
print("Checking df2 for this specific pair:")
for chunk in pd.read_csv('./output/hurdler_three_site_combinations_df2_optimized.csv', 
                         chunksize=100000,
                         usecols=['site_i_3mer_aa', 'site_ii_3mer_aa', 'pUC18_compatible']):
    matches = chunk[
        (chunk['site_i_3mer_aa'] == site_i_3mer) &
        (chunk['site_ii_3mer_aa'] == site_ii_3mer) &
        (chunk['pUC18_compatible'] == True)
    ]
    if len(matches) > 0:
        print(f"\n✓ Found {len(matches)} matching rows in df2")
        print(matches.iloc[0])
        break
else:
    print("NOT found in df2")
