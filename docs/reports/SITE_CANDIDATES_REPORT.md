# HURDLER Site Candidates 分析报告

## 执行时间
2026-01-10

## 数据源
- Seamless insert: 14,180 rows, 70 unique enzymes
- Silent mutation: 7,727 rows, 57 unique enzymes

## Site I 候选酶（Seamless Insert）

### 统计
- **总数**: 70 enzymes
- **Type IIS**: 0 (所有酶的fst5 ≤ site_length)
- **Regular**: 70 (所有酶都是常规切割)

### 代表性酶
| Enzyme | Site | Length | fst5 | fst3 | Ovhang | Type |
|--------|------|--------|------|------|--------|------|
| AatII  | GACGTC | 6 | 5 | -5 | 4 | Regular |
| EcoRI  | GAATTC | 6 | 1 | -1 | -4 | Regular |
| BamHI  | GGATCC | 6 | 1 | -1 | -4 | Regular |

### 结论
✅ Site I：使用常规酶，seamless insert，methylation compatible

---

## Site II 候选酶（Silent Mutation）

### 统计
- **总数**: 57 enzymes
- **Type IIS**: 0 (所有酶的fst5 ≤ site_length)
- **Regular**: 57 (所有酶都是常规切割)

### fst5分布
- fst5=1: 38 enzymes (66.7%)
- fst5=2: 13 enzymes (22.8%)
- fst5=3-5: 6 enzymes (10.5%)

### Overhang分布
- ovhg=-4: 36 enzymes (63.2%) - 可形成630对
- ovhg=-2: 14 enzymes (24.6%) - 可形成91对
- ovhg=2: 6 enzymes (10.5%) - 可形成15对
- ovhg=-3: 1 enzyme (1.7%)

### 代表性酶
| Enzyme | Site | Length | fst5 | fst3 | Ovhang |
|--------|------|--------|------|------|--------|
| BamHI  | GGATCC | 6 | 1 | -1 | -4 |
| EcoRI  | GAATTC | 6 | 1 | -1 | -4 |
| XhoI   | CTCGAG | 6 | 1 | -1 | -4 |

### 结论
✅ Site II：使用常规酶，silent mutation，methylation compatible

---

## Site III 候选酶 - **需要明确定义**

### 当前发现
❌ **问题**: 数据中没有Type IIS酶（fst5 > site_length）

### 可能的方案

#### 方案A：Type IIS酶（原假设）
- **条件**: fst5 > site_length（切割在识别位点外）
- **数据**: silent mutation中找到0个
- **结论**: ❌ 不可行

#### 方案B：与Site II配对（推荐）
- **条件**: 
  1. 来自silent mutation pool
  2. Methylation compatible
  3. 与Site II enzyme不同（避免识别位点重叠）
  4. 产生相同的overhang（用于连接重复单元）
- **数据**: 
  - ovhg=-4: 36 enzymes → 可形成630个不同组合
  - ovhg=-2: 14 enzymes → 可形成91个不同组合
  - ovhg=2: 6 enzymes → 可形成15个不同组合
- **结论**: ✅ 可行，可产生大量组合

#### 方案C：其他条件？
- 需要您明确指定

---

## 关键问题

### 1. "识别位点和切割位点不重叠"的含义

**可能理解1**: Type IIS酶
- 切割位点在识别位点外部
- 但数据中不存在这样的酶

**可能理解2**: Site II和Site III识别位点不重叠
- 使用不同的酶（不同的识别序列）
- 避免在同一位置切割
- 但产生相同的粘性末端用于连接

### 2. HURDLER工作原理推测

```
重复单元: [A1-A2-A3]

Site I (seamless): 在单元前插入，无框移
Site II: 在单元内切割，产生粘性末端A
Site III: 在单元内另一位置切割，产生粘性末端A（相同）

连接: 
单元1-SiteII端 + 单元2-SiteIII端 → 形成重复
因为粘性末端相同，可以连接
```

---

## 建议

### 如果采用方案B（推荐）：

**Site I**: 
- Pool: Seamless insert enzymes
- Filter: Methylation compatible, regular enzyme
- Count: 70 enzymes

**Site II**: 
- Pool: Silent mutation enzymes  
- Filter: Methylation compatible, regular enzyme
- Count: 57 enzymes

**Site III**: 
- Pool: Silent mutation enzymes
- Filter: 
  - Methylation compatible
  - Regular enzyme
  - **enzyme != Site II enzyme**
  - **overhang == Site II overhang**
- Count: 56 enzymes per Site II choice (57-1)

**预期df1大小**:
- 约70 × (sum of valid II-III pairs) combinations
- 如果考虑orthogonality和plasmid compatibility，会进一步减少

---

## 需要确认的问题

**请明确以下信息：**

1. ✅ Site I使用常规酶（seamless insert）- 正确吗？

2. ✅ Site II使用常规酶（silent mutation）- 正确吗？

3. ❓ **Site III应该是：**
   - [ ] A. Type IIS酶（但数据中不存在）
   - [ ] B. 与Site II不同但相同overhang的酶
   - [ ] C. 其他条件：_________________

4. ❓ "识别位点和切割位点不重叠"具体指：
   - [ ] A. Type IIS特性
   - [ ] B. Site II和III使用不同酶（不同识别位点）
   - [ ] C. 其他：_________________

5. ❓ Site II和Site III必须产生相同overhang吗？
   - [ ] 是，用于连接重复单元
   - [ ] 否，其他原因：_________________

---

## 数据文件

详细信息已保存至：
- `output/site_i_candidates_detail.csv` (70 enzymes)
- `output/site_ii_candidates_detail.csv` (57 enzymes)
- `inspect_site_candidates_executed.ipynb` (完整分析notebook)

请查看这些文件并提供反馈。
