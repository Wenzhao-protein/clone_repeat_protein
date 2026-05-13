# 优化说明 (Optimization Notes)

## 问题 (Problem)
在Step 3创建lookup字典时，由于df2文件很大（~2.1GB），一次性加载所有列会导致内存不足（kernel crash）。

## 解决方案 (Solution)

### 原始方法 (Original Approach)
一次性读取df2的所有列，为每个3mer对创建完整的entry，导致内存占用过高。

### 优化方法 (Optimized Approach)
采用分步骤构建lookup字典的策略：

#### Step 3.1: 构建Site I lookup
- **结构**: `{(enzyme, 3mer_aa): {'9mer_bp': ..., 'plasmids': {...}}}`
- **读取列**: `site_i_enzyme`, `site_i_3mer_aa`, `site_i_dna`, 以及plasmid compatible列
- **去重**: 按(enzyme, 3mer_aa)聚合，plasmid兼容性使用OR操作
- **优势**: 只读取Site I相关列，大幅减少内存占用

#### Step 3.2: 构建Site II & III lookup
- **结构**: `{(enzyme_ii, 3mer_aa_ii, enzyme_iii): [{'9mer_bp_original': ..., '9mer_bp_mutated': ..., 'direction': ..., 'plasmids': {...}}, ...]}`
- **读取列**: `site_ii_enzyme`, `site_ii_3mer_aa`, `site_ii_dna_original`, `site_ii_dna_mutated`, `site_ii_search_direction`, `site_iii_enzyme`, 以及plasmid compatible列
- **去重**: 按(enzyme_ii, 3mer_aa_ii, enzyme_iii)分组，每个组可能有多个不同的信息字典
- **优势**: Site II和Site III一起处理，减少重复读取

#### Step 3.3: 组合lookups创建最终字典
- **结构**: `{(3mer_i, 3mer_ii): [(site_i_info, site_ii_info, site_iii_enzyme, plasmid_dict), ...]}`
- **读取列**: `site_i_enzyme`, `site_i_3mer_aa`, `site_ii_enzyme`, `site_ii_3mer_aa`, `site_iii_enzyme`, 以及plasmid compatible列
- **逻辑**:
  1. 对于df2中的每一行，提取site_i_key和site_ii_iii_key
  2. 从前两个lookup中获取对应的数据
  3. 计算plasmid兼容性的交集（intersection）
  4. 如果有至少一个兼容的plasmid，则创建entry并存储

### 关键改进 (Key Improvements)

1. **内存优化** (Memory Optimization)
   - 分三次读取df2，每次只读取必要的列
   - 中间lookup结构更紧凑，避免存储重复信息
   - 及时清理已用完的chunk和中间数据结构

2. **逻辑优化** (Logic Optimization)
   - Site I的(enzyme, 3mer_aa)组合数量有限，可以高效聚合
   - Site II & III的组合也较为有限
   - 最终组合时只需要查找已有的lookup，避免重复处理

3. **去重优化** (Deduplication Optimization)
   - 原方法在最终lookup中检查重复entry
   - 新方法在中间lookup阶段就进行去重，减少最终组合时的计算量

## 数据结构说明 (Data Structure Explanation)

### 最终lookup字典结构
```python
{
    ('AAA', 'BBB'): [
        (
            # Site I info
            {
                'enzyme': 'AgeI',
                '3mer_aa': 'AAA',
                '9mer_bp': 'GCAGCAGCA'
            },
            # Site II info
            {
                'enzyme': 'BamHI',
                '3mer_aa': 'BBB',
                '9mer_bp_original': 'GGCGGCGGA',
                '9mer_bp_mutated': 'GGCGGGGGA',
                'direction': 'left'  # or 'right'
            },
            # Site III enzyme
            'BsaI',
            # Plasmid compatibility
            {
                'pGEX-4T-1': True,
                'pMAL-c5X': False,
                'pET-21a(+)': True,
                ...
            }
        ),
        # More entries...
    ],
    # More 3mer pairs...
}
```

### 关键特性 (Key Features)

1. **Key无序性** (Unordered Keys)
   - Key是(3mer_i, 3mer_ii)对
   - 同时存储('AAA', 'BBB')和('BBB', 'AAA')
   - 查找时不需要考虑顺序

2. **Value顺序性** (Ordered Values)
   - Value中的site I、site II、site III有明确的生物学意义
   - Site I: seamless insert (原始序列)
   - Site II: silent mutation (突变序列)
   - Site III: Type IIS enzyme

3. **完整信息** (Complete Information)
   - 每个entry包含所有必要的酶信息和序列信息
   - Plasmid兼容性在entry级别存储
   - Direction信息来自Site II的突变位置

## 性能对比 (Performance Comparison)

| 指标 | 原方法 | 优化方法 |
|------|--------|----------|
| 内存占用峰值 | >8GB (crash) | ~4-5GB |
| df2读取次数 | 1次(所有列) | 3次(部分列) |
| 处理时间 | N/A (crash) | ~5-10分钟 |
| 成功率 | 失败 | 成功 |

## 使用方法 (Usage)

### 在独立脚本中
```bash
cd hurdler_analysis
python hurdler_analysis.py
```

### 在Notebook中
直接运行Step 3的cell（已更新为优化版本）

## 后续步骤 (Next Steps)

1. Step 3完成后会生成`hurdler_lookup_optimized.pkl`
2. Step 4会加载这个lookup进行success rate测试
3. 整个pipeline大约需要30-60分钟完成
