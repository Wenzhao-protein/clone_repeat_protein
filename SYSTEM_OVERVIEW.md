# HURDLER三位点组合分析系统 - 创建说明

## 📦 已创建的文件

### 1. 主分析Notebook
**文件**: `hurdler_site_combination_analysis.ipynb`

**功能**：
- 完整的数据分析流程
- 筛选满足HURDLER要求的三个RE位点组合
- 生成df1和df2数据集
- 创建查询函数
- 生成可视化图表和统计分析

**输出**：
- `output/hurdler_three_site_combinations_df1.csv` - 基础组合数据
- `output/hurdler_three_site_combinations_df2.csv` - 包含3mer AA序列的完整数据
- `output/hurdler_3mer_aa_lookup.csv` - 快速查询表
- 多个PNG可视化图表

**关键特性**：
- ✅ Site I和II检查甲基化敏感性（DH5α兼容）
- ✅ Site III必须是Type IIS酶（切割位点在识别序列外）
- ✅ Site II和III有相同的ovhg pattern
- ✅ Site II允许沉默突变
- ✅ Site I和II要么ovhg不同，要么是正交的（不互补）
- ✅ 检查每个质粒的兼容性（不切割骨架）

---

### 2. 查询工具
**文件**: `hurdler_query.py`

**功能**：命令行查询工具，根据3mer AA序列和质粒查找有效组合

**使用方式**：

```bash
# 单个查询
python hurdler_query.py --site-i-aa "ABC" --site-ii-aa "DEF" --plasmid "pET-28a(+)"

# 批量查询
python hurdler_query.py --batch queries.csv --output results.csv

# 列出可用选项
python hurdler_query.py --list

# 紧凑输出
python hurdler_query.py --site-i-aa "ABC" --site-ii-aa "DEF" --plasmid "pET-28a(+)" --compact
```

**输出信息**：
- 每个site的酶名称
- 3mer AA序列
- DNA序列（原始和突变后）
- Ovhang信息
- Frame shift信息

---

### 3. 验证工具
**文件**: `hurdler_validate.py`

**功能**：数据验证和详细统计

**使用方式**：

```bash
# 运行验证
python hurdler_validate.py

# 导出报告
python hurdler_validate.py --export validation_summary.txt

# 自定义数据路径
python hurdler_validate.py --df1 /path/to/df1.csv --df2 /path/to/df2.csv
```

**验证内容**：
- ✅ 检查所有必需输入文件是否存在
- ✅ 验证df1和df2的结构完整性
- ✅ 检查Site II/III的ovhg是否匹配
- ✅ 验证3mer AA序列格式
- ✅ 统计酶使用频率
- ✅ 分析质粒兼容性
- ✅ 计算3mer AA覆盖率

---

### 4. 文档
**文件**: 
- `HURDLER_ANALYSIS_README.md` - 完整技术文档
- `HURDLER_QUICKSTART.md` - 快速入门指南（中文）

**内容包括**：
- 详细的使用说明
- 数据结构文档
- 示例代码
- 故障排除指南
- 应用场景示例

---

### 5. 示例文件
**文件**: `example_batch_query.csv`

批量查询的示例输入文件，展示CSV格式。

---

## 🎯 核心功能实现

### DataFrame 1 (df1)
包含所有有效的三位点组合及质粒兼容性：

**列**：
- `site_i`, `site_ii`, `site_iii` - 酶名称
- `ovhg_i`, `ovhg_ii`, `ovhg_iii` - 粘性末端长度
- `[plasmid]_compatible` - 每个质粒的兼容性布尔值
- `compatible_plasmids` - 兼容质粒列表
- `num_compatible_plasmids` - 兼容质粒数量

### DataFrame 2 (df2)
基于df1，添加3mer AA序列信息：

**新增列**：
- `site_i_3mer_aa` - Site I对应的3氨基酸序列
- `site_i_dna` - 编码DNA序列
- `site_i_frame` - Frame shift值
- `site_i_codon_usage` - 密码子使用频率
- `site_ii_3mer_aa` - Site II对应的3氨基酸序列
- `site_ii_dna` - 原始DNA序列
- `site_ii_dna_mutated` - 沉默突变后的DNA序列
- `site_ii_frame` - Frame shift值
- `site_ii_codon_usage` - 密码子使用频率

### 查询函数
`find_hurdler_sites(site_i_3mer_aa, site_ii_3mer_aa, plasmid, df2)`

**输入**：
- 两个3mer AA序列（Site I和Site II）
- 质粒名称

**返回**：
- DataFrame包含所有满足条件的三位点组合
- 如果没有找到，返回None

**特点**：
- ✅ 自动检查质粒兼容性
- ✅ 验证Site II允许沉默突变
- ✅ 确保Site II和III有相同ovhg
- ✅ 提供详细的酶和序列信息

---

## 📊 数据处理流程

```
输入数据源
├── methylation_check.csv (甲基化敏感性)
├── neb_buffer_activity_cleaned.csv (NEB酶特性)
├── plasmid_digest_check.csv (质粒兼容性)
├── restriction_enzyme_slient_mutation.csv (Site II候选)
├── restriction_enzyme_seamless_insert.csv (Site I候选)
└── orthogonality.csv (酶对正交性)
    ↓
筛选酶池
├── Site I: 甲基化不敏感 + 无星号活性
├── Site II: Site I条件 + 允许沉默突变
└── Site III: Type IIS酶
    ↓
生成组合 (df1)
├── 检查Site I和II的ovhg差异或正交性
├── 检查Site II和III的ovhg相同
└── 验证质粒兼容性
    ↓
添加3mer AA (df2)
├── 从seamless_insert获取Site I的3mer AA
├── 从silent_mutation获取Site II的3mer AA
└── 合并生成完整数据
    ↓
输出
├── df1.csv (基础组合)
├── df2.csv (完整数据)
├── lookup.csv (快速查询)
└── 可视化图表
```

---

## ✅ 满足的需求检查

### 1. 三个RE位点 ✅
- [x] Site I: 无缝插入，甲基化不敏感
- [x] Site II: 沉默突变，甲基化不敏感，与Site III相同ovhg
- [x] Site III: Type IIS，与Site II相同ovhg

### 2. 基本要求 ✅
- [x] 无简并碱基
- [x] 粘性末端2-5bp
- [x] 无星号活性
- [x] 商业可获得

### 3. 甲基化兼容性 ✅
- [x] Site I和II对6mA/5mC不敏感

### 4. Type IIS检查 ✅
- [x] Site III切割位点在识别序列外

### 5. Ovhang匹配 ✅
- [x] Site II和III相同ovhg pattern
- [x] Site I和II不同ovhg或正交

### 6. 沉默突变 ✅
- [x] Site II允许沉默突变失活

### 7. 质粒兼容性 ✅
- [x] Site I和II不切割质粒骨架
- [x] 针对8个质粒分别检查

### 8. 数据输出 ✅
- [x] df1: 基础组合 + 质粒兼容性
- [x] df2: 完整数据 + 3mer AA序列
- [x] 查询函数: 输入3mer AA和质粒，返回有效组合

---

## 🚀 使用工作流

### 首次使用
1. 运行notebook生成数据（10-30分钟）
2. 验证数据：`python hurdler_validate.py`
3. 查看可用选项：`python hurdler_query.py --list`

### 日常查询
```bash
# 命令行快速查询
python hurdler_query.py --site-i-aa "XXX" --site-ii-aa "YYY" --plasmid "pET-28a(+)"

# 或在Python中
import pandas as pd
df2 = pd.read_csv('./output/hurdler_three_site_combinations_df2.csv')
results = df2[(df2['site_i_3mer_aa']=='XXX') & ...]
```

### 批量分析
1. 准备CSV文件
2. `python hurdler_query.py --batch input.csv --output output.csv`
3. 分析结果CSV

---

## 📈 预期输出规模

基于repo中的数据估计：

- **Site I候选**: ~50-100个酶
- **Site II候选**: ~30-60个酶
- **Site III候选**: ~10-30个Type IIS酶
- **总组合数**: 15,000-180,000个（取决于筛选严格度）
- **df2行数**: 可能数百万行（每个酶-3mer AA对应一行）

**优化提示**：
- 数据预计算，避免重复查询
- 使用lookup表加速
- CSV格式便于筛选和分析

---

## 🔧 可扩展性

### 添加新质粒
在notebook的质粒列表中添加，重新运行即可。

### 修改筛选标准
在"Build Enzyme Pools"部分调整筛选条件。

### 自定义正交性阈值
在`check_orthogonality`函数中修改`min_score`参数。

### 添加新的输出格式
在df2基础上添加新列或创建新的衍生DataFrame。

---

## 📝 注意事项

1. **首次运行耗时**：生成所有组合需要较长时间，这是一次性操作
2. **文件大小**：df2可能较大（几百MB），正常现象
3. **内存使用**：处理大数据集时可能需要8GB+内存
4. **数据更新**：如果更新了input数据，需要重新运行notebook

---

## 🎉 总结

你现在拥有了一套完整的HURDLER三位点分析系统：

✅ **自动化分析** - 一键生成所有有效组合
✅ **灵活查询** - 命令行或Python查询
✅ **质量验证** - 内置数据验证工具
✅ **完整文档** - 中英文使用指南
✅ **可扩展** - 易于添加新质粒或修改标准

所有代码都有详细注释，便于理解和修改。祝你的HURDLER项目成功！🧬
