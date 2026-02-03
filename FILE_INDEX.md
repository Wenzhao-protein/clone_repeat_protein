# HURDLER分析工具包 - 文件索引

## 🎯 快速导航

### 我想...

#### 📖 了解系统
→ 阅读 [SYSTEM_OVERVIEW.md](SYSTEM_OVERVIEW.md) - 系统完整说明
→ 阅读 [HURDLER_QUICKSTART.md](HURDLER_QUICKSTART.md) - 5分钟快速入门

#### 🚀 开始使用
1. 打开 `hurdler_site_combination_analysis.ipynb` 并运行所有单元格
2. 运行 `python test_hurdler_system.py` 验证系统
3. 运行 `python hurdler_query.py --list` 查看可用选项

#### 🔍 查询RE组合
→ 使用 `hurdler_query.py --site-i-aa "XXX" --site-ii-aa "YYY" --plasmid "质粒名"`
→ 或查看 [查询示例](#查询示例)

#### 📊 查看统计数据
→ 运行 `python hurdler_validate.py`
→ 或在notebook中查看可视化图表

#### ❓ 解决问题
→ 查看 [HURDLER_ANALYSIS_README.md](HURDLER_ANALYSIS_README.md) 的"Troubleshooting"部分
→ 运行 `python test_hurdler_system.py` 诊断问题

---

## 📁 文件清单

### 核心文件（必需）

| 文件 | 类型 | 用途 | 首次使用？ |
|------|------|------|-----------|
| **hurdler_site_combination_analysis.ipynb** | Notebook | 生成所有数据 | ✅ 必须先运行 |
| **hurdler_query.py** | 脚本 | 查询工具 | 在notebook运行后使用 |
| **hurdler_validate.py** | 脚本 | 验证工具 | 用于检查数据质量 |
| **test_hurdler_system.py** | 脚本 | 系统测试 | 运行以验证安装 |

### 文档文件

| 文件 | 语言 | 内容 |
|------|------|------|
| **SYSTEM_OVERVIEW.md** | 中文 | 完整系统说明和创建文档 |
| **HURDLER_QUICKSTART.md** | 中文 | 5分钟快速入门指南 |
| **HURDLER_ANALYSIS_README.md** | 英文 | 详细技术文档 |
| **FILE_INDEX.md** | 中文 | 本文件 - 文件索引 |

### 示例文件

| 文件 | 用途 |
|------|------|
| **example_batch_query.csv** | 批量查询输入示例 |

---

## 🔄 工作流程

```
┌─────────────────────────────────────────────────────────┐
│ 1. 首次设置                                              │
│    运行: hurdler_site_combination_analysis.ipynb        │
│    耗时: 10-30分钟                                       │
│    生成: df1.csv, df2.csv, lookup.csv, 可视化图表       │
└─────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────┐
│ 2. 验证系统                                              │
│    运行: python test_hurdler_system.py                  │
│    确认所有测试通过                                      │
└─────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────┐
│ 3. 日常使用（选择一种方式）                              │
│                                                          │
│  A. 命令行查询                                           │
│     python hurdler_query.py --site-i-aa "XXX" ...       │
│                                                          │
│  B. 批量查询                                             │
│     python hurdler_query.py --batch input.csv ...       │
│                                                          │
│  C. Python脚本                                           │
│     import pandas as pd                                 │
│     df2 = pd.read_csv('output/...df2.csv')             │
│     # 自定义查询和分析                                   │
└─────────────────────────────────────────────────────────┘
```

---

## 📊 输出文件位置

### 生成的数据文件（在 `./output/` 目录）

| 文件 | 大小估计 | 内容 |
|------|----------|------|
| `hurdler_three_site_combinations_df1.csv` | 几MB | 所有有效的三位点组合 |
| `hurdler_three_site_combinations_df2.csv` | 几十MB-几百MB | 包含3mer AA序列的完整数据 |
| `hurdler_3mer_aa_lookup.csv` | <1MB | 快速查询索引 |
| `hurdler_plasmid_statistics.png` | <100KB | 质粒兼容性统计图 |
| `hurdler_overhang_distribution.png` | <100KB | Overhang分布图 |
| `hurdler_top_enzymes.png` | <100KB | 高频酶统计图 |

---

## 💡 常用命令速查

### 查询相关

```bash
# 查看可用的3mer AA序列
python hurdler_query.py --list

# 单个查询（详细输出）
python hurdler_query.py --site-i-aa "NEQ" --site-ii-aa "IQA" --plasmid "pET-28a(+)"

# 单个查询（紧凑输出）
python hurdler_query.py --site-i-aa "NEQ" --site-ii-aa "IQA" --plasmid "pET-28a(+)" --compact

# 批量查询
python hurdler_query.py --batch my_queries.csv --output results.csv
```

### 验证和统计

```bash
# 验证数据完整性
python hurdler_validate.py

# 导出验证报告
python hurdler_validate.py --export report.txt

# 测试系统
python test_hurdler_system.py
```

### 在Python中使用

```python
# 基础查询
import pandas as pd
df2 = pd.read_csv('./output/hurdler_three_site_combinations_df2.csv')

# 查询特定组合
results = df2[
    (df2['site_i_3mer_aa'] == 'NEQ') &
    (df2['site_ii_3mer_aa'] == 'IQA') &
    (df2['pET-28a(+)_compatible'] == True)
]

# 查看结果
print(f"找到 {len(results)} 个组合")
for _, row in results.iterrows():
    print(f"{row['site_i']} + {row['site_ii']} + {row['site_iii']}")
```

---

## 🎓 学习路径

### 初级用户（只需要查询）
1. ✅ 阅读 [HURDLER_QUICKSTART.md](HURDLER_QUICKSTART.md)
2. ✅ 运行 notebook 生成数据
3. ✅ 使用 `hurdler_query.py` 查询

### 中级用户（需要定制分析）
1. ✅ 完成初级步骤
2. ✅ 阅读 [HURDLER_ANALYSIS_README.md](HURDLER_ANALYSIS_README.md)
3. ✅ 学习在Python中使用df2进行自定义分析
4. ✅ 根据需要修改查询脚本

### 高级用户（需要修改系统）
1. ✅ 完成中级步骤
2. ✅ 阅读 [SYSTEM_OVERVIEW.md](SYSTEM_OVERVIEW.md)
3. ✅ 理解notebook中的数据处理流程
4. ✅ 根据需要修改筛选标准或添加新功能

---

## 📋 查询示例

### 示例1：查找特定蛋白序列的RE位点

```bash
# 假设你的重复单元包含 "NEQ" 和 "IQA"
python hurdler_query.py \
  --site-i-aa "NEQ" \
  --site-ii-aa "IQA" \
  --plasmid "pET-28a(+)"
```

### 示例2：批量查询多个组合

创建 `my_queries.csv`:
```csv
site_i_3mer_aa,site_ii_3mer_aa,plasmid
NEQ,IQA,pET-28a(+)
GAL,PAL,pGEX-4T-1
LPA,VQL,pMAL-c5X
```

运行:
```bash
python hurdler_query.py --batch my_queries.csv --output my_results.csv
```

### 示例3：在Python中筛选特定条件

```python
import pandas as pd

# 加载数据
df2 = pd.read_csv('./output/hurdler_three_site_combinations_df2.csv')

# 筛选：4bp overhang，兼容pET-28a(+)
filtered = df2[
    (df2['ovhg_i'] == -4) &
    (df2['ovhg_ii'] == -4) &
    (df2['pET-28a(+)_compatible'] == True)
]

print(f"找到 {len(filtered)} 个满足条件的组合")

# 查看使用了哪些酶
print("\n使用的Site I酶:")
print(filtered['site_i'].value_counts().head())
```

---

## ⚡ 性能提示

1. **首次生成数据很慢是正常的** - 一次性操作，10-30分钟
2. **df2文件可能很大** - 几百MB是正常的，包含了所有组合
3. **查询很快** - 使用pandas筛选，毫秒级响应
4. **内存使用** - 加载df2可能需要几GB内存

---

## 🔧 定制化

### 添加新质粒
在notebook的第5部分修改 `plasmids` 列表

### 修改筛选标准
在notebook的第3部分修改候选酶的筛选条件

### 更改正交性阈值
在 `check_orthogonality` 函数中修改 `min_score` 参数

### 添加新的分析
在notebook最后添加新的分析单元格

---

## 📞 获取帮助

```bash
# 脚本帮助信息
python hurdler_query.py --help
python hurdler_validate.py --help

# 测试系统状态
python test_hurdler_system.py

# 查看详细文档
cat HURDLER_ANALYSIS_README.md
```

---

## ✅ 检查清单

使用前确认：
- [ ] 已安装所需Python包（pandas, biopython, matplotlib, seaborn, tqdm）
- [ ] `utils/output/` 中有所有必需的输入文件
- [ ] 已运行 `hurdler_site_combination_analysis.ipynb`
- [ ] 运行 `test_hurdler_system.py` 所有测试通过
- [ ] 已创建 `output/` 目录并包含生成的CSV文件

---

**最后更新**: 2026年1月10日
**版本**: 1.0
**作者**: HURDLER Analysis System
