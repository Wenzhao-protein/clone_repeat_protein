# HURDLER Success Rate Optimized Notebook

## 概述

`hurdler_success_rate_optimized.ipynb` 是一个完全独立的notebook，可以从基础数据文件开始，自动生成所有HURDLER组合数据，并计算success rate。

## 特点

✅ **无需预先导出的数据文件** - 所有数据都在notebook内生成  
✅ **一键运行** - 使用 "Run All" 即可完成全流程  
✅ **自动生成图表和结果文件**  

## 使用方法

### 方式1: 在VS Code中运行

1. 打开 `hurdler_success_rate_optimized.ipynb`
2. 点击工具栏的 "Run All" 按钮
3. 等待所有cells执行完成（约30-60秒）

### 方式2: 使用命令行运行

```bash
# 使用nbconvert执行notebook
jupyter nbconvert --to notebook --execute hurdler_success_rate_optimized.ipynb --output hurdler_success_rate_optimized_executed.ipynb

# 或者使用papermill执行
papermill hurdler_success_rate_optimized.ipynb hurdler_success_rate_optimized_executed.ipynb
```

## Notebook结构

### 第1部分: 导入库
- 导入所有必需的Python库

### 第2部分: 数据加载
- **Cell 2**: 加载基础酶数据（seamless insert, silent mutation, NEB quality, plasmid compatibility）
- **Cell 3**: 生成Site III酶数据（Type IIS enzymes with ovhg ∈ {-4, +2}）
- **Cell 4**: 生成组合数据框（Site I + Site II + Site III combinations）
- **Cell 5**: 准备数据和提取plasmid列表

### 第3部分: 算法定义
- **Cell 6**: 定义`check_sequence_hurdler_success()`函数
  - 使用doubled sequence matching算法
  - 3mer AA indexing for fast lookup
  - 距离约束: 8 ≤ span ≤ module_length + 2

### 第4部分: 建立索引
- **Cell 7**: 构建3mer AA索引和plasmid-specific pattern sets

### 第5部分: 运行测试
- **Cell 8**: 测试7-60 AA的success rates
  - 每个长度测试1000个随机序列
  - 8个plasmids × 54个长度 = 432个数据点
  - 约20秒完成

### 第6部分: 结果分析
- **Cell 9**: 分析和显示结果摘要

### 第7部分: 可视化
- **Cell 10**: 生成两个图表
  - Success rate vs module length (折线图)
  - Success rate at 60 AA by plasmid (柱状图)

### 第8部分: 保存结果
- **Cell 11**: 保存结果到CSV文件

## 输入文件（必需）

这些文件必须存在于 `./utils/output/` 目录：

- `restriction_enzyme_seamless_insert.csv` - Site I候选酶
- `restriction_enzyme_slient_mutation.csv` - Site II候选酶
- `neb_buffer_activity_cleaned.csv` - NEB质量数据
- `plasmid_digest_check.csv` - Plasmid兼容性

## 输出文件

运行后会在 `./output/` 目录生成：

1. **hurdler_success_rate_optimized.pdf** - 可视化图表
   - 左图: Success rate随module length的变化
   - 右图: 60 AA时各plasmid的success rate对比

2. **hurdler_success_rate_results_optimized.csv** - 完整结果数据
   - 432行 (8 plasmids × 54 lengths)
   - 包含: module_length, plasmid, successes, tests, success_rate

## 关键结果

### Success Rates at 60 AA:
- pUC18: **98.9%** (最高)
- pCold_I: **98.5%**  
- pGEX-4T-1: **98.4%**
- pMAL-c5X: **97.6%**
- pET-21a(+): **96.9%**
- pET-28a(+)_start_codon: **96.6%**
- pQE-3: **96.5%**
- pET-28a(+): **95.1%** (最低)

### 整体平均: 73.0%

## 性能

- **数据生成**: 瞬间（从基础文件读取和组合）
- **索引构建**: ~2秒
- **Success rate测试**: ~20秒 (432,000次测试)
- **总运行时间**: ~30秒

## 技术细节

### Site III筛选标准:
1. 必须是Type IIS enzymes (cuts outside recognition site)
2. Overhang必须是 -4 或 +2 bp
3. 无简并碱基 (only ATCG)
4. NEB quality (good ligation, no star activity)
5. Plasmid compatible (至少一个plasmid)

### 组合生成逻辑:
- Site I (seamless insert) + Site II (silent mutation): 共享相同的3mer AA pattern
- Site III: overhang与Site II相同
- 所有三个sites必须与plasmid兼容

### 优化算法:
- **Doubled sequence**: `ABCDEF → ABCDEFABCDEF` 允许跨边界匹配
- **3mer AA indexing**: 预提取所有3mers快速生成候选
- **Per-plasmid pattern sets**: 预过滤兼容patterns
- **Early termination**: 找到第一个有效匹配立即返回True

## 故障排除

### 问题: Cell执行失败显示 "FileNotFoundError"
**解决**: 确保所有输入文件存在于 `./utils/output/` 目录

### 问题: 生成的数据数量不对
**解决**: 检查输入文件是否是最新版本，特别是enzyme筛选标准是否正确

### 问题: Success rate测试很慢
**解决**: 这是正常的，432,000次测试约需20秒。可以在Cell 8中减少 `num_tests_per_length` 来加快测试（会降低统计精度）

## 修改建议

如果需要修改测试参数：

1. **修改测试长度范围**: 在Cell 8中修改 `module_lengths = list(range(7, 61))`
2. **修改测试次数**: 在Cell 8中修改 `num_tests_per_length = 1000`
3. **修改Site III筛选**: 在Cell 3中修改overhang标准 `if enzyme.ovhg not in [-4, 2]:`

## 版本历史

- **v1.0** (2026-01-11): 初始版本，完全自包含的notebook，无需外部数据文件
