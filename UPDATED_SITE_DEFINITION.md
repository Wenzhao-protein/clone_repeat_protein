# HURDLER Site Definition 更新分析

## Type IIS 酶分析结果

### 发现
- **全库中共有 173 个 Type IIS 酶**
- **methylation compatible的Type IIS酶数量: 0**
- **结论**: Type IIS酶都不是methylation compatible，无法用于HURDLER系统

### Type IIS 识别方法（已修正）
```python
def is_type_iis(enzyme):
    """正确的Type IIS识别方法"""
    site_len = len(str(enzyme.site))
    # Type IIS: 切割位点在识别位点外部
    # fst5 < 0 (左侧切割) 或 fst5 > site_length (右侧切割)
    return enzyme.fst5 < 0 or enzyme.fst5 > site_len
```

### 示例 Type IIS 酶
| Enzyme | Site | site_len | fst5 | fst3 | ovhg | 位置 |
|--------|------|----------|------|------|------|------|
| BsaI   | GGTCTC | 6 | 7 | 5 | -4 | 右侧切割 |
| BsmBI  | CGTCTC | 6 | 7 | 5 | -4 | 右侧切割 |
| BbsI   | GAAGAC | 6 | 8 | 6 | -4 | 右侧切割 |
| FokI   | GGATG  | 5 | 14 | 13 | -4 | 右侧切割 |
| PpiI   | GAACNNNNNCTC | 12 | -7 | -24 | 5 | 左侧切割 |

---

## 重新定义 Site III

### 原假设（错误）
- Site III = Type IIS酶
- **问题**: 所有Type IIS酶都不是methylation compatible

### 新理解

根据用户说明："识别位点和切割位点不重叠"可能指：

#### 选项1: Site II 和 Site III 的识别位点不重叠
```
含义: 
- Site II 和 Site III 是两个不同的酶
- 它们的识别序列不同（避免在同一位置识别）
- 但它们产生相同的粘性末端（overhang）用于连接重复单元
```

**实现**:
```python
# Site II
site_ii_enzyme = 'BamHI'  # GGATCC, ovhg=-4

# Site III  
site_iii_enzyme = 'EcoRI'  # GAATTC, ovhg=-4 (相同)

# 条件:
# 1. site_ii_enzyme != site_iii_enzyme ✓
# 2. ovhg(site_ii) == ovhg(site_iii) ✓
```

**可用组合数**:
- ovhg=-4: 36 enzymes → 36×35/2 = 630 pairs
- ovhg=-2: 14 enzymes → 14×13/2 = 91 pairs
- ovhg=2: 6 enzymes → 6×5/2 = 15 pairs
- **总计**: 736 pairs

---

## Type IIT 酶的考虑

### Type IIT 定义
- 识别**回文序列**的reverse complement
- 例如: EcoRI识别 GAATTC，其reverse complement也是 GAATTC
- **关键**: 对于非回文序列，正向和反向是两个独立的识别位点

### 示例分析

#### 回文序列（Palindromic）
```
EcoRI: GAATTC
Reverse complement: GAATTC (相同)
→ 只有1个独特的识别位点
```

#### 非回文序列（Non-palindromic）
```
BsaI: GGTCTC (如果是Type IIT)
Reverse complement: GAGACC
→ 有2个独立的识别位点

在双链DNA上:
正向: 5'-GGTCTC-3'
      3'-CCAGAG-5'

反向: 5'-GAGACC-3'
      3'-CTCTGG-5'
```

### 对HURDLER的影响

如果某些酶是Type IIT（识别双向），需要考虑：

1. **在序列搜索时**:
   - 需要同时检查序列和其reverse complement
   - 如果都存在，则有2个潜在识别位点

2. **在orthogonality检查时**:
   - Site I, II, III的识别位点都不应overlap
   - 需要检查正向和反向

3. **实现建议**:
```python
def check_site_in_sequence(enzyme_site, sequence):
    """检查酶识别位点是否在序列中"""
    from Bio.Seq import Seq
    
    # 正向
    if enzyme_site in sequence:
        return True
    
    # 反向complement（Type IIT考虑）
    site_rc = str(Seq(enzyme_site).reverse_complement())
    if site_rc in sequence:
        return True
    
    return False
```

---

## 建议的新Site定义

### Site I (Seamless Insert)
- **Pool**: Seamless insert enzymes (70)
- **Filter**: 
  - Methylation compatible ✓
  - Regular enzyme (fst5 ≤ site_length) ✓
- **功能**: 在重复单元前插入，无框移

### Site II (Silent Mutation)
- **Pool**: Silent mutation enzymes (57)
- **Filter**:
  - Methylation compatible ✓
  - Regular enzyme (fst5 ≤ site_length) ✓
- **功能**: 在重复单元内切割，产生粘性末端

### Site III (Silent Mutation)
- **Pool**: Silent mutation enzymes (57)
- **Filter**:
  - Methylation compatible ✓
  - Regular enzyme (fst5 ≤ site_length) ✓
  - **enzyme != Site II enzyme** ✓ (不同酶，避免识别位点重叠)
  - **ovhang == Site II ovhang** ✓ (相同粘性末端，用于连接)
- **功能**: 在重复单元内另一位置切割，产生相同粘性末端

### Orthogonality 检查
```python
def check_orthogonality(site_i, site_ii, site_iii, plasmid_seq):
    """
    检查三个site在plasmid上的正交性
    
    考虑Type IIT：需要检查正向和反向complement
    """
    from Bio.Seq import Seq
    
    # 获取所有识别位点（包括reverse complement）
    def get_all_sites(enzyme):
        site = str(enzyme.site)
        site_rc = str(Seq(site).reverse_complement())
        if site == site_rc:  # Palindromic
            return [site]
        else:  # Non-palindromic, Type IIT
            return [site, site_rc]
    
    sites_i = get_all_sites(site_i)
    sites_ii = get_all_sites(site_ii)
    sites_iii = get_all_sites(site_iii)
    
    # 查找所有位置
    positions_i = find_all_positions(sites_i, plasmid_seq)
    positions_ii = find_all_positions(sites_ii, plasmid_seq)
    positions_iii = find_all_positions(sites_iii, plasmid_seq)
    
    # 检查是否有overlap
    if has_overlap(positions_i, positions_ii, positions_iii):
        return False
    
    return True
```

---

## 需要确认的问题

1. ✅ Type IIS酶不可用（都不是methylation compatible）

2. ❓ **"识别位点和切割位点不重叠" 的正确含义**:
   - [ ] A. Site II和III是不同的酶（不同的识别序列）
   - [ ] B. Site II和III的识别位点在序列上不overlap
   - [ ] C. 其他含义

3. ❓ **Site II和III是否必须产生相同overhang**:
   - [ ] 是 - 用于连接重复单元
   - [ ] 否 - 其他原因

4. ❓ **Type IIT的处理**:
   - [ ] 需要考虑reverse complement（双向识别）
   - [ ] 不需要特殊处理
   - [ ] 只使用palindromic enzymes

5. ❓ **当前数据的适用性**:
   - [ ] 使用现有的regular enzymes (Site I: 70, Site II/III: 57)
   - [ ] 需要重新筛选包含Type IIS的数据（但它们不methylation compatible）

---

## 下一步行动

请明确以下信息，我将据此更新代码：

1. Site III的正确定义
2. 是否需要考虑Type IIT的双向识别
3. Site II和III是否必须产生相同overhang
