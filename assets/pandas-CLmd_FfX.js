const n=`---
title: "Python Pandas 数据处理实战"
date: "2025-09-19"
author: "zorrooz"
tags: ["Python", "Pandas", "数据处理", "数据分析", "数据清洗"]
draft: false
description: "Pandas 库在生物信息学数据处理中的常用操作和技巧"
---

# Python Pandas 数据处理实战

## 数据读取

\`\`\`python
import pandas as pd

# 读取CSV文件
df = pd.read_csv('data.csv')

# 读取Excel文件
df = pd.read_excel('data.xlsx')
\`\`\`

## 数据清洗

\`\`\`python
# 处理缺失值
df.fillna(0, inplace=True)

# 删除重复行
df.drop_duplicates(inplace=True)`;export{n as default};
