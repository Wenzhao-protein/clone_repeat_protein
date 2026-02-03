#!/usr/bin/env python3
"""
诊断脚本：正确识别Type IIS酶
"""

import pandas as pd
from Bio import Restriction
from Bio.Restriction import Analysis, RestrictionBatch
import warnings
warnings.filterwarnings('ignore')

def check_type_iis_methods(enzyme_name):
    """使用多种方法检查是否为Type IIS"""
    try:
        enz = getattr(Restriction, enzyme_name)
        
        results = {
            'enzyme': enzyme_name,
            'site': str(enz.site),
            'site_length': len(str(enz.site)),
            'fst5': enz.fst5,
            'fst3': enz.fst3,
            'ovhg': enz.ovhg,
        }
        
        # 方法1: 检查类名或类型
        enzyme_class = enz.__class__.__name__
        results['class_name'] = enzyme_class
        
        # 方法2: 检查是否有特殊的cut特征
        # Type IIS通常：fst5 < 0 (左侧切割) 或 fst5 > site_length (右侧切割)
        results['method2_left_cut'] = enz.fst5 < 0
        results['method2_right_cut'] = enz.fst5 > len(str(enz.site))
        results['method2_is_typeIIS'] = results['method2_left_cut'] or results['method2_right_cut']
        
        # 方法3: 检查fst3
        results['method3_fst3_outside'] = abs(enz.fst3) > len(str(enz.site))
        
        # 方法4: 检查是否切割位点完全在识别位点外
        site_len = len(str(enz.site))
        cut_min = min(enz.fst5, enz.fst5 + enz.ovhg) if enz.ovhg else enz.fst5
        cut_max = max(enz.fst5, enz.fst5 + enz.ovhg) if enz.ovhg else enz.fst5
        
        # 切割区间：[cut_min, cut_max]
        # 识别区间：[0, site_len]
        results['cut_range'] = f"[{cut_min}, {cut_max}]"
        results['site_range'] = f"[0, {site_len}]"
        results['method4_cuts_outside'] = (cut_max <= 0) or (cut_min >= site_len)
        
        return results
    except:
        return None

# 已知的Type IIS酶
known_type_iis = ['BsaI', 'BsmBI', 'BbsI', 'SapI', 'BtgZI', 'FokI', 'AlwI']

print("="*100)
print("已知Type IIS酶的特征分析")
print("="*100)

results = []
for enz_name in known_type_iis:
    r = check_type_iis_methods(enz_name)
    if r:
        results.append(r)

df_known = pd.DataFrame(results)
print(df_known.to_string(index=False))

print("\n" + "="*100)
print("分析silent_mutation数据中的酶")
print("="*100)

# 读取数据
df_silent = pd.read_csv('./utils/output/restriction_enzyme_slient_mutation.csv')
unique_enzymes = df_silent['name'].unique()

print(f"\n总共 {len(unique_enzymes)} 个独特的酶")

# 检查所有酶
all_results = []
for enz_name in unique_enzymes:
    r = check_type_iis_methods(enz_name)
    if r:
        all_results.append(r)

df_all = pd.DataFrame(all_results)

# 按不同方法统计
print("\n方法2 (fst5 < 0 或 fst5 > site_length):")
type_iis_m2 = df_all[df_all['method2_is_typeIIS']]
print(f"  找到 {len(type_iis_m2)} 个Type IIS酶")
if len(type_iis_m2) > 0:
    print(f"  示例: {list(type_iis_m2['enzyme'].values[:10])}")

print("\n方法4 (切割范围完全在识别位点外):")
type_iis_m4 = df_all[df_all['method4_cuts_outside']]
print(f"  找到 {len(type_iis_m4)} 个Type IIS酶")
if len(type_iis_m4) > 0:
    print(f"  示例: {list(type_iis_m4['enzyme'].values[:10])}")

# 保存详细结果
df_all.to_csv('output/enzyme_type_analysis.csv', index=False)
print(f"\n详细结果已保存至: output/enzyme_type_analysis.csv")

# 显示一些边界案例
print("\n" + "="*100)
print("边界案例分析 (fst5 = 0, 1, site_length, site_length+1)")
print("="*100)
boundary_cases = df_all[
    (df_all['fst5'] == 0) | 
    (df_all['fst5'] == 1) | 
    (df_all['fst5'] == df_all['site_length']) |
    (df_all['fst5'] == df_all['site_length'] + 1)
]
if len(boundary_cases) > 0:
    print(boundary_cases[['enzyme', 'site', 'site_length', 'fst5', 'fst3', 'method2_is_typeIIS', 'method4_cuts_outside']].to_string(index=False))
else:
    print("没有找到边界案例")

print("\n" + "="*100)
print("负数fst5的案例（切割在识别位点左侧）")
print("="*100)
negative_fst5 = df_all[df_all['fst5'] < 0]
if len(negative_fst5) > 0:
    print(f"找到 {len(negative_fst5)} 个酶有负数fst5")
    print(negative_fst5[['enzyme', 'site', 'fst5', 'fst3', 'ovhg']].to_string(index=False))
else:
    print("没有找到负数fst5的酶")

print("\n" + "="*100)
print("fst5 > site_length的案例（切割在识别位点右侧）")
print("="*100)
beyond_site = df_all[df_all['fst5'] > df_all['site_length']]
if len(beyond_site) > 0:
    print(f"找到 {len(beyond_site)} 个酶的fst5 > site_length")
    print(beyond_site[['enzyme', 'site', 'site_length', 'fst5', 'fst3', 'ovhg']].to_string(index=False))
else:
    print("没有找到fst5 > site_length的酶")
