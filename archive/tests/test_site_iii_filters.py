#!/usr/bin/env python3
"""
测试Site III筛选条件
"""

import pandas as pd
import Bio
from Bio import Restriction
import warnings
warnings.filterwarnings('ignore')

# 加载数据
print("Loading data...")
df_neb = pd.read_csv('./utils/output/neb_buffer_activity_cleaned.csv')
df_plasmid_check = pd.read_csv('./utils/output/plasmid_digest_check.csv', index_col=0)
df_methylation = pd.read_csv('./utils/output/methylation_check.csv')

# 辅助函数
def has_degenerate_bases(enzyme_name):
    """Check if enzyme recognition site has degenerate bases"""
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        site = str(enzyme.site)
        standard_bases = set('ATCG')
        return not all(base in standard_bases for base in site)
    except:
        return True

def check_ovhg_length(enzyme_name):
    """Check if overhang length is 2-5 bp (Golden Gate compatible)"""
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        ovhg_len = abs(enzyme.ovhg) if enzyme.ovhg else 0
        return 2 <= ovhg_len <= 5
    except:
        return False

def check_neb_quality(enzyme_name, df_neb):
    """Check if enzyme has good NEB characteristics"""
    if enzyme_name not in df_neb['enzyme'].values:
        return False
    
    enzyme_data = df_neb[df_neb['enzyme'] == enzyme_name].iloc[0]
    
    if enzyme_data['ligation_efficiencies'] == 'low':
        return False
    
    if enzyme_data['star_activity'] == True:
        return False
    
    return True

def check_plasmid_compatible(enzyme_name, df_plasmid_check):
    """Check if enzyme doesn't cut plasmid backbone"""
    if enzyme_name not in df_plasmid_check.index:
        return False
    return df_plasmid_check.loc[enzyme_name].any()

def get_ovhg_pattern(enzyme_name):
    """Get overhang pattern (length and direction)"""
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        return enzyme.ovhg
    except:
        return None

def is_type_iis(enzyme_name):
    """Check if enzyme is Type IIS"""
    try:
        enzyme = getattr(Bio.Restriction, enzyme_name)
        site_length = len(enzyme.site)
        return enzyme.fst5 < 0 or enzyme.fst5 > site_length
    except:
        return False

def check_methylation_compatible(enzyme_name, df_methylation):
    """Check if enzyme is not sensitive to 6mA/5mC methylation"""
    dh5a_compatible = set(
        df_methylation.loc[~df_methylation['6mA_5mC_sensitive'], 'enzyme'].dropna().tolist() +
        df_methylation.loc[~df_methylation['6mA_5mC_sensitive'], 'prototype'].dropna().tolist()
    )
    return enzyme_name in dh5a_compatible

# 测试Site III筛选
print("\n" + "="*80)
print("Site III Filtering Process")
print("="*80)

all_enzymes = [str(enz) for enz in Restriction.AllEnzymes]
print(f"\n1. Total enzymes in AllEnzymes: {len(all_enzymes)}")

# Filter 1: Valid overhang
candidates = [e for e in all_enzymes if get_ovhg_pattern(e) is not None]
print(f"2. After valid overhang filter: {len(candidates)}")

# Filter 2: No degenerate bases
candidates = [e for e in candidates if not has_degenerate_bases(e)]
print(f"3. After no degenerate bases filter: {len(candidates)}")

# Filter 3: Overhang length 2-5 bp
candidates = [e for e in candidates if check_ovhg_length(e)]
print(f"4. After overhang length (2-5bp) filter: {len(candidates)}")

# Filter 4: NEB quality
candidates = [e for e in candidates if check_neb_quality(e, df_neb)]
print(f"5. After NEB quality filter: {len(candidates)}")

# Filter 5: Plasmid compatibility
candidates = [e for e in candidates if check_plasmid_compatible(e, df_plasmid_check)]
print(f"6. After plasmid compatibility filter: {len(candidates)}")

# 统计Type IIS
type_iis_enzymes = [e for e in candidates if is_type_iis(e)]
regular_enzymes = [e for e in candidates if not is_type_iis(e)]

print(f"\n7. Final Site III candidates: {len(candidates)}")
print(f"   - Type IIS: {len(type_iis_enzymes)}")
print(f"   - Regular: {len(regular_enzymes)}")

# 显示一些Type IIS示例
if type_iis_enzymes:
    print("\n" + "="*80)
    print("Type IIS enzymes found in Site III candidates:")
    print("="*80)
    
    for enz_name in type_iis_enzymes[:20]:  # 显示前20个
        enz = getattr(Bio.Restriction, enz_name)
        print(f"  {enz_name:15s} site={str(enz.site):15s} fst5={enz.fst5:3d} ovhg={enz.ovhg:3d}")

# 检查methylation compatibility
methyl_compat = [e for e in candidates if check_methylation_compatible(e, df_methylation)]
print(f"\n8. Methylation compatible (for comparison): {len(methyl_compat)}")
print(f"   - This means {len(candidates) - len(methyl_compat)} Site III candidates are NOT methylation compatible")

# 分析overhang分布
print("\n" + "="*80)
print("Overhang distribution in Site III candidates:")
print("="*80)
ovhg_dist = {}
for e in candidates:
    ovhg = get_ovhg_pattern(e)
    if ovhg not in ovhg_dist:
        ovhg_dist[ovhg] = []
    ovhg_dist[ovhg].append(e)

for ovhg in sorted(ovhg_dist.keys()):
    enzymes = ovhg_dist[ovhg]
    type_iis_count = sum(1 for e in enzymes if is_type_iis(e))
    print(f"  ovhg={ovhg:3d}: {len(enzymes):4d} enzymes ({type_iis_count} Type IIS)")

print("\n" + "="*80)
print("Summary")
print("="*80)
print(f"Site III can now use {len(candidates)} enzymes (vs 57 before)")
print(f"Including {len(type_iis_enzymes)} Type IIS enzymes that were previously excluded")
