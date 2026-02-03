#!/usr/bin/env python
"""
简单测试：随机生成序列并检查匹配率
"""
import pickle
import random

# 加载lookup
with open('./output/hurdler_lookup_optimized.pkl', 'rb') as f:
    lookup = pickle.load(f)

plasmid = 'pUC18_compatible'
print(f"Plasmid: {plasmid}")
print(f"Total pairs in lookup: {len(lookup[plasmid]):,}")

# 生成10个随机序列并测试
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
n_tests = 100
n_success = 0

for test_idx in range(n_tests):
    # 生成长度10的随机序列
    length = 10
    sequence = ''.join(random.choice(amino_acids) for _ in range(length))
    
    # 提取3mer
    extended = sequence[-2:] + sequence + sequence[:2]
    three_mers = set()
    for i in range(len(sequence)):
        three_mer = extended[i:i+3]
        if len(three_mer) == 3:
            three_mers.add(three_mer)
    
    # 检查所有3mer对
    three_mer_list = list(three_mers)
    found = False
    for i in range(len(three_mer_list)):
        for j in range(i+1, len(three_mer_list)):
            pair1 = (three_mer_list[i], three_mer_list[j])
            pair2 = (three_mer_list[j], three_mer_list[i])
            
            if pair1 in lookup[plasmid] or pair2 in lookup[plasmid]:
                found = True
                n_success += 1
                print(f"Test {test_idx+1}: ✓ {sequence} - found pair!")
                break
        if found:
            break
    
    if not found and test_idx < 10:  # 只打印前10个失败的
        print(f"Test {test_idx+1}: ✗ {sequence} - no match (3mers: {three_mers})")

print(f"\nSuccess rate: {n_success}/{n_tests} = {n_success/n_tests*100:.1f}%")
