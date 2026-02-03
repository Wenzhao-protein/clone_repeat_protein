# HURDLER成功率测试系统使用指南

## 系统架构

该系统用于测试HURDLER方法在不同长度重复蛋白序列上的成功率。

### 数据流程

```
df1 (酶组合) → Lookup Dictionary → 轻量级Lookup → 成功率测试
```

### 关键文件

1. **数据生成**（一次性运行）：
   - `generate_hurdler_data.py` - 生成df1（三位点酶组合）
   - `create_hurdler_lookup.py` - 从df1创建lookup字典
   
2. **测试脚本**：
   - `hurdler_success_rate_test.py` - 完整成功率测试（4-60长度）
   - `quick_test_hurdler.py` - 快速验证测试
   
3. **输出文件**：
   - `hurdler_three_site_combinations_df1.csv` - 57,415行三位点组合
   - `hurdler_3mer_lookup.pkl` - 完整lookup（用于详细查询）
   - `hurdler_3mer_lightweight_lookup.pkl` - **轻量级lookup（仅用于成功率测试）**
   - `hurdler_lookup_summary.csv` - 摘要表格

## 快速开始

### 第一步：生成数据（已完成）

```bash
# df1已生成
ls -lh output/hurdler_three_site_combinations_df1.csv
# 输出: 2.2M, 57,415 combinations
```

### 第二步：创建Lookup字典（运行中）

```bash
python create_hurdler_lookup.py
# 预计运行时间: 20-30分钟
# 进度会显示在终端
```

生成的文件：
- `hurdler_3mer_lookup.pkl` - 完整版本
- **`hurdler_3mer_lightweight_lookup.pkl`** - 轻量版本（测试用）

### 第三步：运行成功率测试

```bash
# 快速测试（验证系统）
python quick_test_hurdler.py

# 完整测试（4-60长度，每个长度1000次）
python hurdler_success_rate_test.py

# 自定义参数
python hurdler_success_rate_test.py --min-length 10 --max-length 50 --n-trials 500
```

## 优化说明

### 重复序列的环状特性

对于重复蛋白序列，考虑首尾相连的边界：

```python
# 标准3mers
for i in range(len(sequence) - 2):
    three_mers.add(sequence[i:i+3])

# 环状边界3mers
three_mers.add(sequence[-2:] + sequence[0])   # 末尾2个 + 开头1个
three_mers.add(sequence[-1] + sequence[:2])    # 末尾1个 + 开头2个
```

### Lookup字典结构

**轻量级Lookup**（用于成功率测试）：
```python
{
    frozenset(['ACE', 'DFG']): {
        'pET-28a(+)': True,
        'pGEX-4T-1': True,
        ...
    },
    ...
}
```

- **Key**: frozenset of 3mer AA对（无序）
- **Value**: 兼容plasmid的字典
- **优势**: 
  - 只需检查key是否存在
  - 内存占用小
  - 加载速度快

**完整Lookup**（用于详细查询）：
```python
{
    frozenset(['ACE', 'DFG']): [
        [site_i_enzyme, site_ii_enzyme, site_iii_enzyme,
         site_i_3mer_aa, site_ii_3mer_aa,
         site_i_dna_9mer, site_ii_dna_9mer,
         plasmid_compatibility_dict],
        ...
    ],
    ...
}
```

## 成功率测试算法

```python
for each sequence_length from 4 to 60:
    for trial in range(1000):
        # 1. 生成随机序列
        sequence = generate_random_aa(length)
        
        # 2. 提取3mers（含环状边界）
        three_mers = extract_3mers_circular(sequence)
        
        # 3. 检查每个plasmid
        for plasmid in plasmids:
            # 检查任意3mer对是否可行
            for pair in all_pairs(three_mers):
                key = frozenset(pair)
                if key in lookup and plasmid in lookup[key]:
                    success = True
                    break
```

## 输出结果

### 1. 原始结果
`hurdler_success_rate_raw_results.csv`:
```csv
length,trial,plasmid,n_unique_3mers,feasible,n_valid_pairs
4,0,pET-28a(+),4,True,2
4,1,pET-28a(+),4,False,0
...
```

### 2. 汇总结果
`hurdler_success_rate_summary.csv`:
```csv
length,plasmid,success_rate,avg_n_3mers,avg_n_pairs
4,pET-28a(+),0.45,4.0,0.9
5,pET-28a(+),0.52,5.0,1.2
...
```

### 3. 可视化图表

- **主图**: `hurdler_success_rate_by_length.png`
  - 横坐标: 序列长度（4-60）
  - 纵坐标: 成功率（0-100%）
  - 每个plasmid一条曲线

- **详细图**: `hurdler_success_rate_detailed.png`
  - 左：选定长度的柱状图对比
  - 右：热力图（plasmid × length）

## 性能指标

### 当前状态
- df1行数: 57,415
- 3mer对数量: ~10,000-15,000（预估）
- Lookup加载时间: <1秒
- 单个序列测试: <0.001秒

### 完整测试预估
- 总测试数: 57 × 1000 × 8 = 456,000次检查
- 预计时间: 5-10分钟

## 故障排查

### 问题: lookup文件不存在
```bash
# 检查create_hurdler_lookup.py是否完成
ps aux | grep create_hurdler_lookup.py

# 查看进度
tail -f create_lookup.log
```

### 问题: 测试速度慢
```bash
# 减少trials
python hurdler_success_rate_test.py --n-trials 100

# 缩小长度范围
python hurdler_success_rate_test.py --min-length 10 --max-length 30
```

### 问题: 内存不足
轻量级lookup已经优化，如仍有问题：
```python
# 分批测试
for length_range in [(4, 20), (21, 40), (41, 60)]:
    test_success_rate(min_length=length_range[0], 
                      max_length=length_range[1])
```

## 下一步

1. ✅ 等待`create_hurdler_lookup.py`完成
2. ✅ 运行`quick_test_hurdler.py`验证系统
3. ✅ 运行`hurdler_success_rate_test.py`完整测试
4. ✅ 分析结果图表

## 参数说明

### hurdler_success_rate_test.py参数

```bash
--min-length INT      # 最小序列长度 (默认: 4)
--max-length INT      # 最大序列长度 (默认: 60)
--n-trials INT        # 每个长度的测试次数 (默认: 1000)
--output-dir PATH     # 输出目录 (默认: ./output)
--skip-test           # 跳过测试，只重新绘图
```

### 示例

```bash
# 快速测试（降低trials）
python hurdler_success_rate_test.py --n-trials 100

# 聚焦特定长度范围
python hurdler_success_rate_test.py --min-length 20 --max-length 40

# 只重新绘制图表
python hurdler_success_rate_test.py --skip-test
```

## 文件大小参考

```
hurdler_three_site_combinations_df1.csv      ~2.2 MB
hurdler_3mer_lookup.pkl                      ~XX MB (待生成)
hurdler_3mer_lightweight_lookup.pkl          ~X MB (待生成)
hurdler_success_rate_raw_results.csv         ~XX MB (测试后)
hurdler_success_rate_summary.csv             ~XX KB (测试后)
```
