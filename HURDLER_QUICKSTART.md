# HURDLER Analysis - Quick Start Guide

## 🚀 快速开始

### 第一步：生成数据

在Jupyter中运行主分析notebook：

```bash
jupyter notebook hurdler_site_combination_analysis.ipynb
```

**注意**：首次运行会生成所有数据，可能需要10-30分钟。运行完成后会在 `./output/` 目录生成以下文件：
- hurdler_three_site_combinations_df1.csv
- hurdler_three_site_combinations_df2.csv
- hurdler_3mer_aa_lookup.csv
- 可视化图表（PNG文件）

### 第二步：查询组合

#### 方法1：单个查询（命令行）

```bash
python hurdler_query.py --site-i-aa "AAA" --site-ii-aa "DEF" --plasmid "pET-28a(+)"
```

#### 方法2：批量查询

1. 准备CSV文件（例如 `my_queries.csv`）：
```csv
site_i_3mer_aa,site_ii_3mer_aa,plasmid
AAA,DEF,pET-28a(+)
ABC,GHI,pGEX-4T-1
```

2. 运行批量查询：
```bash
python hurdler_query.py --batch my_queries.csv --output my_results.csv
```

#### 方法3：在Python中使用

```python
import pandas as pd

# 加载数据
df2 = pd.read_csv('./output/hurdler_three_site_combinations_df2.csv')

# 查询
results = df2[
    (df2['site_i_3mer_aa'] == 'AAA') &
    (df2['site_ii_3mer_aa'] == 'DEF') &
    (df2['pET-28a(+)_compatible'] == True)
]

print(f"找到 {len(results)} 个有效组合")
```

### 第三步：查看可用的3mer AA序列

```bash
python hurdler_query.py --list
```

## 📊 理解结果

### 查询结果包含的信息

对于每个有效组合，你会得到：

1. **Site I** (无缝插入位点)
   - 酶名称
   - 3mer AA序列
   - DNA序列
   - 粘性末端长度

2. **Site II** (允许沉默突变的位点)
   - 酶名称
   - 3mer AA序列
   - 原始DNA序列
   - 突变后DNA序列（失活RE位点）
   - 粘性末端长度（与Site III相同）

3. **Site III** (Type IIS酶)
   - 酶名称
   - 粘性末端长度（与Site II相同）

### 示例输出解读

```
Combination 1:
  Site I:   EcoRI           (ovhg:  -4, frame: 2)
            DNA: GAATTC
  Site II:  BamHI           (ovhg:  -4, frame: 1, silent mutation)
            DNA: GGATCC
            Mutated: GGATCT
  Site III: BsaI            (ovhg:  -4, Type IIS)
```

**解读**：
- Site I用EcoRI，产生-4bp粘性末端（5'突出）
- Site II用BamHI，也是-4bp粘性末端，可以通过沉默突变（GGATCC→GGATCT）失活
- Site III用BsaI（Type IIS酶），粘性末端与Site II相同
- 三个位点都不会切割指定质粒的骨架

## 🎯 应用场景

### 场景1：设计重复蛋白克隆策略

1. 确定你的重复单元包含的3mer AA序列
2. 使用查询工具找到合适的RE位点组合
3. 在DNA序列设计中编码这些RE位点

### 场景2：检查已有序列

```python
# 假设你的重复单元序列是 "NEQIQAVIDAGAL"
sequence = "NEQIQAVIDAGAL"

# 提取所有3mer
three_mers = [sequence[i:i+3] for i in range(len(sequence)-2)]
print(three_mers)  # ['NEQ', 'EQI', 'QIQ', 'IQA', ...]

# 对每个3mer查询是否有可用组合
for aa in three_mers:
    # 查询该3mer是否在数据库中
    ...
```

### 场景3：优化多质粒兼容性

```bash
# 查找同时兼容多个质粒的组合
python hurdler_validate.py --export summary.txt
```

查看输出中的"Multi-plasmid compatibility"部分。

## ⚙️ 高级选项

### 紧凑输出格式

```bash
python hurdler_query.py --site-i-aa "AAA" --site-ii-aa "DEF" --plasmid "pET-28a(+)" --compact
```

### 自定义数据路径

```bash
python hurdler_validate.py --df1 /path/to/df1.csv --df2 /path/to/df2.csv
```

### 导出验证报告

```bash
python hurdler_validate.py --export validation_report.txt
```

## 📁 文件结构

```
clone_repeat_protein/
├── hurdler_site_combination_analysis.ipynb  # 主分析notebook
├── hurdler_query.py                         # 查询工具
├── hurdler_validate.py                      # 验证工具
├── HURDLER_ANALYSIS_README.md              # 完整文档
├── HURDLER_QUICKSTART.md                   # 本文档
├── example_batch_query.csv                 # 批量查询示例
├── output/
│   ├── hurdler_three_site_combinations_df1.csv
│   ├── hurdler_three_site_combinations_df2.csv
│   ├── hurdler_3mer_aa_lookup.csv
│   └── *.png (可视化图表)
└── utils/output/
    ├── methylation_check.csv
    ├── neb_buffer_activity_cleaned.csv
    ├── plasmid_digest_check.csv
    └── ...
```

## 🐛 常见问题

### Q: 查询没有返回结果

**A**: 可能的原因：
1. 该3mer AA组合不存在：运行 `python hurdler_query.py --list` 查看可用序列
2. 指定的质粒不兼容：尝试其他质粒
3. 数据未生成：确保已运行主notebook生成数据

### Q: 运行notebook太慢

**A**: 
- 首次运行需要生成大量组合，这是正常的
- 可以在notebook中调整候选酶池的大小来加速
- 生成后的数据可以重复使用，不需要重新运行

### Q: 需要特定的粘性末端长度

**A**: 
```python
# 在df2中筛选特定ovhang长度
df2_filtered = df2[
    (df2['ovhg_i'] == -4) &  # Site I 4bp 5'突出
    (df2['ovhg_ii'] == -4)    # Site II/III 4bp 5'突出
]
```

### Q: 如何找到最常用的组合

**A**: 查看统计报告：
```bash
python hurdler_validate.py
```

看"Top enzymes"部分了解最常用的酶。

## 💡 提示

1. **优先选择常用酶**：EcoRI, BamHI, XhoI等更容易获得且性能更可靠
2. **考虑密码子优化**：df2中包含codon_usage信息，优先选择高频密码子
3. **检查多质粒兼容性**：如果需要在多个表达系统间转换，选择兼容多个质粒的组合
4. **保存查询结果**：批量查询结果保存为CSV，方便后续分析

## 📞 获取帮助

```bash
# 查看帮助信息
python hurdler_query.py --help
python hurdler_validate.py --help
```

查看完整文档：`HURDLER_ANALYSIS_README.md`

---

**祝你的HURDLER克隆项目顺利！** 🧬🔬
