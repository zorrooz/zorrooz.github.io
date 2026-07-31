---
title: "Python Pandas Data Processing in Practice"
date: "2026-07-03"
author: "zorrooz"
tags: ["Python", "Pandas", "Data Processing", "Analysis", "Data Cleaning"]
draft: false
description: "Common operations and techniques of the Pandas library in bioinformatics data processing"
---

# Python Pandas Data Processing in Practice

## Data Reading

```python
import pandas as pd

# Read CSV file
df = pd.read_csv('data.csv')

# Read Excel file
df = pd.read_excel('data.xlsx')

# Read compressed files (gzip / bz2 / zip auto-detected)
df = pd.read_csv('expression.tsv.gz', sep='\t', compression='infer')
```

## Data Cleaning

```python
# Handle missing values
df.fillna(0, inplace=True)

# Remove duplicate rows
df.drop_duplicates(inplace=True)

# Type conversion: unify gene expression values as float
df['expression'] = df['expression'].astype(float)
```

## Common Operations

| Operation | Code | Description |
|-----------|------|-------------|
| Filter rows | `df[df['group'] == 'case']` | Filter samples by condition |
| Group stats | `df.groupby('gene')['value'].mean()` | Aggregate by gene |
| Merge tables | `pd.merge(df1, df2, on='gene_id')` | Join annotation data |
| Pivot table | `df.pivot_table(index='gene', columns='sample', values='value')` | Expression matrix |

> **Note**: In bioinformatics you often load a `counts` table first and then merge it
> with sample metadata; `pd.merge` is one of the most frequently used operations.

## Typical Workflow: Expression Matrix Preprocessing

```python
import pandas as pd

counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)
metadata = pd.read_csv('metadata.csv')

# 1. Filter low-expression genes (expressed in at least 3 samples)
keep = (counts > 0).sum(axis=1) >= 3
counts = counts[keep]

# 2. Order samples by condition for downstream plotting
counts = counts[metadata.sort_values('condition')['sample']]

# 3. CPM normalization
total = counts.sum(axis=1) / 1e6
cpm = counts.div(total, axis=0)

print(f"Genes: {counts.shape[0]}, Samples: {counts.shape[1]}")
```
