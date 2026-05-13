#!/usr/bin/env python3
"""
查找所有Type IIS酶并检查methylation compatibility
"""

import pandas as pd
from Bio import Restriction
import warnings
warnings.filterwarnings('ignore')

# 获取所有酶
all_enzymes = Restriction.AllEnzymes

print("="*100)
print("扫描所有限制酶，查找Type IIS")
print("="*100)

type_iis_all = []
for enz in all_enzymes:
    try:
        site_len = len(str(enz.site))
        
        # Type IIS: fst5 < 0 或 fst5 > site_length
        if enz.fst5 < 0 or enz.fst5 > site_len:
            type_iis_all.append({
                'enzyme': str(enz),
                'site': str(enz.site),
                'site_length': site_len,
                'fst5': enz.fst5,
                'fst3': enz.fst3,
                'ovhg': enz.ovhg
            })
    except:
        pass

df_all_iis = pd.DataFrame(type_iis_all)
print(f"\n在所有限制酶中找到 {len(df_all_iis)} 个Type IIS酶")
print(df_all_iis.to_string(index=False))

# 保存
df_all_iis.to_csv('output/all_type_iis_enzymes.csv', index=False)
print(f"\n已保存至: output/all_type_iis_enzymes.csv")

# 检查methylation compatibility
print("\n" + "="*100)
print("检查methylation compatibility")
print("="*100)

# 读取methylation compatible酶
df_seamless = pd.read_csv('./utils/output/restriction_enzyme_seamless_insert.csv')
df_silent = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')

seamless_names = set(df_seamless['name'].unique())
silent_names = set(df_silent['name'].unique())
all_methyl_compat = seamless_names | silent_names

print(f"\nMethylation compatible酶总数: {len(all_methyl_compat)}")
print(f"  - Seamless insert: {len(seamless_names)}")
print(f"  - Silent mutation: {len(silent_names)}")

# 检查Type IIS中有多少是methylation compatible
type_iis_names = set(df_all_iis['enzyme'])
type_iis_methyl_compat = type_iis_names & all_methyl_compat

print(f"\nType IIS酶中methylation compatible的数量: {len(type_iis_methyl_compat)}")
if type_iis_methyl_compat:
    print("酶名称:")
    for name in sorted(type_iis_methyl_compat):
        row = df_all_iis[df_all_iis['enzyme'] == name].iloc[0]
        print(f"  {name}: site={row['site']}, fst5={row['fst5']}, ovhg={row['ovhg']}")
        
        # 检查在哪个数据集中
        in_seamless = name in seamless_names
        in_silent = name in silent_names
        print(f"    → seamless: {in_seamless}, silent: {in_silent}")

# 结论
print("\n" + "="*100)
print("结论")
print("="*100)
print(f"1. 全库中共有 {len(df_all_iis)} 个Type IIS酶")
print(f"2. 其中只有 {len(type_iis_methyl_compat)} 个是methylation compatible")
print(f"3. 这解释了为什么在filtered数据中找不到Type IIS酶")
