---
title: "Python 数据处理实战：文件 IO、正则与 NumPy/Pandas"
date: "2026-08-04"
author: "zorrooz"
tags: ["Python", "NumPy", "Pandas", "数据处理", "教程"]
draft: false
description: "文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景"
---

# Python 数据处理实战：文件 IO、正则与 NumPy/Pandas

数据处理是 Python 最核心的应用场景。本文讲解文件读写、正则表达式，以及科学计算的两大支柱 NumPy 与 Pandas，示例贴合生物信息学常见任务。

## 1. 文件读写

### 1.1 文本文件

推荐使用 `with` 语句，自动管理文件关闭：

```python
# 写入
with open("output.txt", "w", encoding="utf-8") as f:
    f.write("第一行\n")
    f.write("第二行\n")

# 读取全部
with open("output.txt", "r", encoding="utf-8") as f:
    content = f.read()
print(content)

# 逐行读取（大文件推荐）
with open("output.txt", "r", encoding="utf-8") as f:
    for line in f:
        print(line.strip())   # strip 去掉换行符
```

模式说明：`r` 读、`w` 写（覆盖）、`a` 追加、`rb`/`wb` 二进制。

### 1.2 解析 FASTA 文件

```python
def read_fasta(path):
    """读取 FASTA 文件，返回 {id: 序列} 字典。"""
    sequences = {}
    current_id = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:].split()[0]   # 取第一个字段作为 id
                sequences[current_id] = ""
            elif current_id is not None:
                sequences[current_id] += line
    return sequences

seqs = read_fasta("example.fasta")
for seq_id, seq in seqs.items():
    print(seq_id, len(seq))
```

### 1.3 CSV 文件

标准库 `csv` 模块：

```python
import csv

# 读取
with open("data.csv", "r", newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)          # 按表头读取为字典
    rows = list(reader)
    print(rows[0]["gene"], rows[0]["expression"])

# 写入
with open("out.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["gene", "value"])
    writer.writerows([["TP53", 12.3], ["BRCA1", 8.9]])
```

> 大数据场景推荐 Pandas 的 `read_csv`（见第 4 节）。

## 2. 正则表达式

正则用于文本模式匹配。`re` 模块常用函数：

```python
import re

text = "The protein P12345 has 500 amino acids, code: A1B2C3"

# match：从开头匹配
m = re.match(r"The", text)
print(m.group())          # The

# search：搜索第一个匹配
m = re.search(r"\d+", text)
print(m.group())          # 12345（第一个数字串）

# findall：返回所有匹配
print(re.findall(r"[A-Z]\d", text))   # ['P1', 'A1', 'B2', 'C3']

# sub：替换
print(re.sub(r"\d+", "#", text))
```

### 2.1 常用元字符

| 模式 | 含义 | 示例 |
|------|------|------|
| `.` | 任意字符（除换行） | `a.c` 匹配 `abc` |
| `*` | 前一个字符 0 次或多次 | `ab*c` 匹配 `ac`、`abc` |
| `+` | 前一个字符 1 次或多次 | `ab+c` 匹配 `abc` |
| `?` | 前一个字符 0 次或 1 次 | `colou?r` 匹配 color/colour |
| `^` / `$` | 行首 / 行尾 | `^ATG` 以 ATG 开头 |
| `[abc]` | 字符集 | `[ATGC]` 任意碱基 |
| `\d` / `\w` / `\s` | 数字 / 单词字符 / 空白 | |
| `{n,m}` | 重复 n 到 m 次 | `\d{3}` 三个数字 |
| `(…)` | 分组捕获 | `(ATG){2}` |

### 2.2 提取信息实战

从基因注释文本提取坐标与名称：

```python
line = 'gene=TP53;gene_id=ENSG00000141510;chromosome=17;start=7668402;end=7687550'

m = re.search(r"gene=(\w+);.*?start=(\d+);end=(\d+)", line)
if m:
    gene, start, end = m.groups()
    print(f"{gene}: {start}-{end} ({int(end) - int(start)} bp)")
# TP53: 7668402-7687550 (19148 bp)
```

### 2.3 编译正则提升性能

```python
motif = re.compile(r"^ATG.*(TAA|TAG|TGA)$")   # 完整 ORF 粗略匹配
for seq in ["ATGCCCTAA", "ATGGGG"]:
    print(seq, bool(motif.match(seq)))
```

## 3. NumPy：数值计算基础

### 3.1 数组创建

```python
import numpy as np

a = np.array([1, 2, 3, 4])
b = np.zeros((2, 3))          # 2x3 全零矩阵
c = np.ones(5)
d = np.arange(0, 1, 0.2)      # [0.  0.2 0.4 0.6 0.8]
e = np.linspace(0, 1, 5)      # 0 到 1 等分 5 份
f = np.random.rand(3, 3)      # 均匀分布随机数
```

### 3.2 向量化运算（核心优势）

```python
x = np.array([1.0, 2.0, 3.0])
y = np.array([4.0, 5.0, 6.0])

print(x + y)          # [5. 7. 9.]，无需循环
print(x * y)          # [4. 10. 18.]
print(np.sqrt(x))     # 逐元素开方
print(x.mean(), x.std())   # 2.0 0.816...

# 广播：不同形状自动扩展
matrix = np.array([[1, 2], [3, 4]])
print(matrix * 10)            # 每个元素乘 10
print(matrix + np.array([100, 200]))  # 行广播
```

### 3.3 索引与切片

```python
arr = np.arange(10)
print(arr[2:8:2])        # [2 4 6]，起始:结束:步长
print(arr[arr > 5])      # 布尔掩码 [6 7 8 9]

m = np.random.rand(4, 4)
print(m[0, :])           # 第一行
print(m[:, 1])           # 第二列
print(m[1:3, 1:3])       # 子矩阵
```

### 3.4 生物信息学示例：归一化表达矩阵

```python
# 模拟基因表达矩阵：4 基因 × 5 样本
expr = np.array([
    [12.3, 15.1, 9.8, 20.4, 11.2],
    [3.2, 4.1, 2.9, 5.0, 3.6],
    [80.5, 92.3, 75.1, 100.2, 88.7],
    [45.6, 41.2, 48.9, 39.8, 50.1],
])

# 按行（基因）做 Z-score 标准化：(x - mean) / std
means = expr.mean(axis=1, keepdims=True)
stds = expr.std(axis=1, keepdims=True)
z = (expr - means) / stds

print(np.round(z, 2))
```

## 4. Pandas：表格数据处理

### 4.1 Series 与 DataFrame

```python
import pandas as pd

# Series：一维带标签数组
s = pd.Series([1, 2, 3], index=["a", "b", "c"])

# DataFrame：二维表格
df = pd.DataFrame({
    "gene": ["TP53", "BRCA1", "EGFR", "MYC"],
    "expression": [12.3, 8.9, 45.1, 33.2],
    "chromosome": [17, 17, 7, 8],
})
print(df)
```

### 4.2 读取与写入

```python
df = pd.read_csv("expression.csv")
df = pd.read_csv("data.tsv", sep="\t")
df = pd.read_excel("data.xlsx", sheet_name="Sheet1")
df = pd.read_json("data.json")

df.to_csv("out.csv", index=False)
```

### 4.3 查看与筛选

```python
print(df.head())          # 前 5 行
print(df.info())          # 列类型与缺失值
print(df.describe())      # 数值列统计摘要

# 列选择
print(df["gene"])
print(df[["gene", "expression"]])

# 行筛选（布尔索引）
high = df[df["expression"] > 10]
print(high)

# 多条件
sel = df[(df["expression"] > 10) & (df["chromosome"] == 17)]
print(sel)

# isin 筛选
sel = df[df["gene"].isin(["TP53", "MYC"])]
```

### 4.4 分组与聚合

```python
# 按染色体分组，统计表达均值与数量
grouped = df.groupby("chromosome")["expression"].agg(["mean", "count", "max"])
print(grouped)

# 按分组转换（如组内标准化）
df["zscore"] = df.groupby("chromosome")["expression"].transform(
    lambda x: (x - x.mean()) / x.std()
)
```

### 4.5 缺失值处理

```python
df = pd.DataFrame({"a": [1, None, 3], "b": [4, 5, None]})

print(df.isna().sum())        # 每列缺失数量
print(df.dropna())            # 删除含缺失的行
print(df.fillna(0))           # 填充为 0
print(df.fillna(df.mean()))   # 用均值填充
```

### 4.6 实战：差异表达基因筛选

```python
# 模拟：基因 + 处理组/对照组表达 + p 值
results = pd.DataFrame({
    "gene": [f"GENE{i}" for i in range(100)],
    "log2fc": np.random.normal(0, 1.5, 100),
    "pvalue": np.random.rand(100),
})

# 筛选显著差异基因：|log2FC| > 1 且 p < 0.05
sig = results[(results["log2fc"].abs() > 1) & (results["pvalue"] < 0.05)]

# 按 log2FC 排序
sig = sig.sort_values("log2fc", ascending=False)
print(f"显著基因数: {len(sig)}")
print(sig.head())
```

## 5. 小结

- 文件：`with open` + 逐行读取处理大文件；CSV 用 `csv` 模块或 Pandas
- 正则：`re.search` / `findall` / `sub`，先写小模式再组合
- NumPy：向量化运算 + 广播，避免 Python 循环
- Pandas：`read_csv` → 筛选/分组/缺失值处理 → 导出，是表格数据的标准工作流

至此 Python 教程三部曲完成：入门 → 进阶 → 数据处理。
