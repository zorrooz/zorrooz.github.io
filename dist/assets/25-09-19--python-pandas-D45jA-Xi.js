const n=`---\r
title: "Python Pandas 数据处理实战"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Python", "Pandas", "数据处理", "数据分析", "数据清洗"]\r
draft: false\r
description: "Pandas 库在生物信息学数据处理中的常用操作和技巧"\r
---\r
\r
# Python Pandas 数据处理实战\r
\r
## 数据读取\r
\r
\`\`\`python\r
import pandas as pd\r
\r
# 读取CSV文件\r
df = pd.read_csv('data.csv')\r
\r
# 读取Excel文件\r
df = pd.read_excel('data.xlsx')\r
\`\`\`\r
\r
## 数据清洗\r
\r
\`\`\`python\r
# 处理缺失值\r
df.fillna(0, inplace=True)\r
\r
# 删除重复行\r
df.drop_duplicates(inplace=True)`;export{n as default};
