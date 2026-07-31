---
title: "Python Pandas 数据处理实战"
date: "2026-07-03"
author: "zorrooz"
tags: ["Python", "Pandas", "数据处理", "数据分析", "数据清洗"]
draft: false
description: "Pandas 库在生物信息学数据处理中的常用操作和技巧"
---

# Python Pandas 数据处理实战

## 数据读取

```python
import pandas as pd

# 读取CSV文件
df = pd.read_csv('data.csv')

# 读取Excel文件
df = pd.read_excel('data.xlsx')

# 读取压缩文件（自动识别 gzip / bz2 / zip）
df = pd.read_csv('expression.tsv.gz', sep='\t', compression='infer')
```

## 数据清洗

```python
# 处理缺失值
df.fillna(0, inplace=True)

# 删除重复行
df.drop_duplicates(inplace=True)

# 类型转换：基因表达量统一为 float
df['expression'] = df['expression'].astype(float)
```

## 常用操作

| 操作 | 代码 | 说明 |
|------|------|------|
| 筛选行 | `df[df['group'] == 'case']` | 按条件过滤样本 |
| 分组统计 | `df.groupby('gene')['value'].mean()` | 按基因聚合 |
| 合并表 | `pd.merge(df1, df2, on='gene_id')` | 连接注释信息 |
| 透视表 | `df.pivot_table(index='gene', columns='sample', values='value')` | 表达矩阵 |

> **提示**：生物信息学中经常先读取 `counts` 表，再与样本元信息（`metadata`）
> 合并，`pd.merge` 是最常用的操作之一。

## 典型流程：表达矩阵预处理

```python
import pandas as pd

counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)
metadata = pd.read_csv('metadata.csv')

# 1. 过滤低表达基因（至少 3 个样本表达量 > 0）
keep = (counts > 0).sum(axis=1) >= 3
counts = counts[keep]

# 2. 样本按条件排序，方便下游作图
counts = counts[metadata.sort_values('condition')['sample']]

# 3. CPM 归一化
total = counts.sum(axis=1) / 1e6
cpm = counts.div(total, axis=0)

print(f"基因数: {counts.shape[0]}, 样本数: {counts.shape[1]}")
```
