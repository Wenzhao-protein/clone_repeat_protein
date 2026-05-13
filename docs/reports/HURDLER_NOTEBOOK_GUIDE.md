# HURDLER Success Rate Optimized Notebook - 使用指南

## 📋 概述

这个notebook完整实现了HURDLER系统的数据生成、分析和导出流程。**完全独立运行，不依赖预先计算的文件**。

## ✨ 核心特性

1. **完整数据生成**：从原始酶数据生成所有组合
2. **RE配对矩阵**：基于orthogonality的智能酶配对过滤
3. **内存优化**：批处理避免109M行的内存溢出
4. **快速导出**：生成lookup文件用于实际HURDLER查询

## 🔄 Notebook结构

### Step 1: 加载基础数据
- 加载seamless insert, silent mutation, plasmid check, NEB数据
- **Cell 2-3**: 数据加载和Site III生成

### Step 1.5: 构建RE配对矩阵 ⭐
- **Cell 5**: 生成Site I-II酶配对矩阵
- 基于orthogonality score ≥ 1的规则
- 将3,990对过滤到3,749对兼容酶对
- **关键**：避免生成109M行，只生成有效组合

### Step 2: 处理Site I和Site II数据
- **Cell 6**: 提取和清理Site I/II数据
- 计算codon usage
- 确定mutation direction
- 生成14,180个Site I + 7,727个Site II条目

### Step 3: 生成Combined Dataframe ⭐⭐
- **Cell 7**: 组合Site I + Site II + Site III
- **只处理**配对矩阵允许的3,749个酶对
- 批处理生成（100对/批）
- 添加plasmid兼容性和Site III信息
- **输出**：4,016,521行组合（~750 MB）
- **时间**：约2-3分钟

### Step 4-6: 构建Lookup结构
- **Cell 8**: 加载combined data
- **Cell 9**: 构建pattern lookup字典
- **Cell 10**: 构建plasmid pattern sets和3mer index

### Step 7-9: 成功率分析
- 随机序列测试
- 不同module length的成功率
- 可视化分析

### Step 10: 最终导出 ⭐
- **Cell 最后**: 导出所有lookup文件
  - hurdler_combined_dataframe.csv
  - hurdler_pattern_lookup.pkl
  - hurdler_plasmid_patterns.pkl
  - hurdler_3mer_aa_index.pkl
  - hurdler_fast_match_package.pkl（整合包）
  - hurdler_lookup_summary.csv

## 🚀 使用方法

### 方法1：运行整个Notebook
```python
# 在Jupyter中：Run All Cells
# 总时间：约5-10分钟
```

### 方法2：使用生成的数据文件
如果已经运行过一次，后续可以直接从Step 4开始：
```python
# 跳过Cells 1-7（数据已生成）
# 从Cell 8开始运行（加载已有数据）
```

### 方法3：使用独立脚本
```bash
# 只生成数据，不运行分析
python generate_data_complete.py
# 时间：约2-3分钟
# 输出：所有CSV和PKL文件
```

## 📊 生成的数据

### 中间文件
- `hurdler_site_iii_data.csv`: 16个Type IIS酶
- `site_i_site_ii_pairing_matrix.csv`: 70×57酶配对矩阵
- `hurdler_site_i_data.csv`: 14,180个Site I条目
- `hurdler_site_ii_data.csv`: 7,727个Site II条目

### 主要输出
- `hurdler_combined_dataframe.csv`: **4,016,521行**（750 MB）
  - 包含所有有效组合
  - Pattern, 酶信息, DNA序列, plasmid兼容性

### Lookup文件（用于查询）
- `hurdler_pattern_lookup.pkl`: Pattern → 组合列表
- `hurdler_plasmid_patterns.pkl`: Plasmid → Pattern集合
- `hurdler_3mer_aa_index.pkl`: 所有3mer AA集合
- `hurdler_fast_match_package.pkl`: 整合查询包
- `hurdler_lookup_summary.csv`: 统计摘要

## 🔧 关键优化

### 1. RE配对矩阵过滤
**问题**：直接组合会生成14,180 × 7,727 = 109,516,860行
**解决**：使用orthogonality过滤到3,749对酶
**结果**：最终4M行（可管理的大小）

### 2. 批处理生成
- 每批100个酶对
- 及时释放内存
- 流式写入CSV
- 避免内存溢出

### 3. 数据结构优化
- Pattern lookup: O(1)查询
- Plasmid patterns: 快速plasmid过滤
- 3mer index: 快速序列预检查

## 📈 性能指标

- **数据生成时间**: 2-3分钟
- **Pattern lookup构建**: 2-3分钟
- **内存使用**: ~2 GB峰值
- **最终数据大小**: ~750 MB CSV + ~100 MB PKL

## 🎯 与原notebook的差异

### 改进
1. ✅ 完全独立，不依赖预计算文件
2. ✅ 添加RE配对矩阵生成逻辑
3. ✅ 批处理避免内存问题
4. ✅ 正确处理列名（re_site_shifted vs re_9bp_site）
5. ✅ 完整的导出功能

### 保持一致
1. ✅ 相同的酶过滤规则
2. ✅ 相同的配对逻辑（orthogonality ≥ 1）
3. ✅ 相同的pattern生成方法
4. ✅ 相同的plasmid兼容性检查

## 💡 使用建议

1. **首次运行**：Run All（生成所有数据）
2. **后续分析**：从Step 4开始（加载已有数据）
3. **只需查询文件**：运行到Step 10导出即可
4. **调试**：使用`generate_data_complete.py`独立生成数据

## 🐛 故障排查

### 内存不足
- 减小batch_size（Cell 7，默认100）
- 分段运行：先生成数据，再分析

### 时间太长
- 数据生成约2-3分钟是正常的
- 可以使用后台脚本：`nohup python generate_data_complete.py &`

### 缺少文件
- 确保`./utils/output/`下有所有输入文件
- 检查orthogonality.csv是否存在

## 📞 快速参考

```python
# 快速加载lookup用于查询
import pickle
with open('./output/hurdler_fast_match_package.pkl', 'rb') as f:
    package = pickle.load(f)

pattern_lookup = package['pattern_lookup']
plasmid_patterns = package['plasmid_patterns']
all_3mers = package['3mer_aa_index']

# 查询示例
pattern = "ADV.*?ALK"
if pattern in pattern_lookup:
    matches = pattern_lookup[pattern]
    print(f"Found {len(matches)} enzyme combinations")
```

## ✅ 验证成功

运行完成后应该看到：
- ✅ 5个CSV文件生成
- ✅ 4个PKL文件生成
- ✅ hurdler_combined_dataframe.csv = 4,016,521行
- ✅ 成功率分析完成
- ✅ 所有图表显示正常

---

**版本**: Optimized v2.0  
**日期**: 2026-01-11  
**状态**: ✅ 完整可运行
