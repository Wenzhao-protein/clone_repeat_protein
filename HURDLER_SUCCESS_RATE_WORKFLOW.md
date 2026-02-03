# HURDLER Success Rate Analysis Workflow

## 概述

本工作流程用于计算HURDLER三位点组合系统对不同质粒和内部模块长度的成功率。

## 系统设计

### 三位点定义

1. **Site I** (Seamless Insert)
   - 70个候选酶（常规酶）
   - 需要methylation compatible
   - 识别位点编码3mer AA，无框移

2. **Site II** (Silent Mutation)
   - 57个候选酶（常规酶）
   - 需要methylation compatible  
   - 识别位点通过沉默突变编码3mer AA

3. **Site III** (Type IIS)
   - 7个候选酶：BsaI, BbsI, BspMI, SapI, BsgI, BseRI, BsmI
   - **不需要**methylation compatible
   - 识别位点在重复单元外（Type IIS特性）
   - 必须满足：商业可获得、无简并碱基、2-5bp overhang、NEB质量、无星号活性

### 关键要求

- Site II和III必须产生**相同的overhang**（用于重复单元连接）
- Site II和III必须是**不同的酶**（避免识别位点重叠）
- 所有三个位点必须与质粒骨架**正交**（orthogonal）

## 数据生成流程

### 步骤1：生成df1（三位点组合）

```bash
python generate_hurdler_data.py --skip-df2
```

**输出：**
- `output/hurdler_three_site_combinations_df1.csv`
- 4,899个有效的三位点组合

**包含列：**
- site_i_enzyme, site_ii_enzyme, site_iii_enzyme
- site_i_ovhg, site_ii_ovhg, site_iii_ovhg
- [plasmid]_compatible (每个质粒一列)

### 步骤2：生成df2（添加3mer AA序列）

```bash
python generate_df2_memory_optimized.py
```

**特性：**
- 分批处理，内存优化
- Site I和II：添加re_site_shifted_tl (3mer AA)
- Site III：标记为'N/A'（Type IIS不编码AA）

**输出：**
- `output/hurdler_three_site_combinations_df2.csv`
- 预计数百万行（df1 × Site I选项 × Site II选项）

### 步骤3：创建lookup dictionary

```bash
python create_hurdler_lookup_from_df2.py
```

**功能：**
- 从df2创建快速查找字典
- 结构：`{plasmid: {frozenset(3mer_AA_pair): [enzyme_combinations]}}`
- 支持多进程加速

**输出：**
- `output/hurdler_lookup_dict.pkl` - 完整版本
- `output/hurdler_lookup_lite.pkl` - 轻量版（仅keys，用于成功率测试）

### 步骤4：计算成功率

```bash
python test_hurdler_success_rate.py
```

**测试参数：**
- Module lengths: 4-60氨基酸
- 每个长度：1000次随机测试
- 所有可用质粒

**输出：**
- `output/hurdler_success_rate_results.csv` - 原始结果
- `output/hurdler_success_rate_summary.csv` - 汇总统计
- `output/hurdler_success_rate_plot.png` - 成功率曲线图
- `output/hurdler_success_rate_heatmap.png` - 热力图

## 当前进度

✅ **已完成：**
1. df1生成：4,899个三位点组合
2. Site III筛选：7个Type IIS酶（不需methylation compatible）
3. Type IIS识别方法修正：`fst5 < 0` 或 `fst5 > site_length`

🔄 **进行中：**
- df2生成（分批处理，内存优化）

⏸️ **待执行：**
- 创建lookup dictionary
- 成功率测试
- 结果分析和可视化

## 关键发现

### Type IIS酶分析

在BioPython全库中找到173个Type IIS酶，但：
- **0个**是methylation compatible
- 经过质量筛选后（NEB、plasmid compatible等），只有**7个**满足Site III要求

### Site III候选酶

| Enzyme | Site | fst5 | Overhang | Site II配对数 |
|--------|------|------|----------|--------------|
| BsaI   | GGTCTC | 7 | -4 | 36 |
| BbsI   | GAAGAC | 8 | -4 | 36 |
| BspMI  | ACCTGC | 10 | -4 | 36 |
| SapI   | GCTCTTC | 8 | -3 | 1 |
| BsgI   | GTGCAG | 22 | 2 | 6 |
| BseRI  | GAGGAG | 16 | 2 | 6 |
| BsmI   | GAATGC | 7 | 2 | 6 |

### 配对统计

- ovhg=-4: 36 Site II × 3 Site III = **108 pairs**
- ovhg=-3: 1 Site II × 1 Site III = **1 pair**
- ovhg=2: 6 Site II × 3 Site III = **18 pairs**
- **总计**: 127个Site II-III配对

## 文件说明

### 生成脚本
- `generate_hurdler_data.py` - 主生成脚本（df1+df2）
- `generate_df2_memory_optimized.py` - 内存优化的df2生成
- `create_hurdler_lookup_from_df2.py` - 创建lookup dict
- `test_hurdler_success_rate.py` - 成功率测试

### 辅助脚本
- `find_all_type_iis.py` - 查找所有Type IIS酶
- `test_site_iii_filters.py` - 测试Site III筛选逻辑
- `diagnose_type_iis.py` - Type IIS酶诊断

### 输入数据
- `utils/output/restriction_enzyme_seamless_insert.csv`
- `utils/output/restriction_enzyme_slient_mutation.csv`
- `utils/output/methylation_check.csv`
- `utils/output/neb_buffer_activity_cleaned.csv`
- `utils/output/plasmid_digest_check.csv`
- `utils/output/orthogonality.csv`

## 注意事项

1. **内存要求**: df2生成需要大量内存（~6GB+），使用分批处理避免OOM
2. **Type IIS理解**: Type IIS酶的识别位点在重复单元外，因此不需要编码3mer AA
3. **Orthogonality**: Type IIS酶不在orthogonality数据库中，代码中假设为正交
4. **随机种子**: 成功率测试使用固定种子（42）以确保可重复性

## 性能优化

- 多进程处理（lookup dict创建）
- 分批读写（df2生成）
- 轻量级lookup（仅keys，用于测试）
- Progress bar（tqdm）

## 参考

- Type IIS酶定义：fst5 < 0 或 fst5 > site_length
- Golden Gate assembly: 2-5bp overhang
- Methylation compatibility: DH5α (6mA/5mC不敏感)
