const n=`---\r
title: "用 Pandas 处理 FASTQ 质量统计"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Python", "Pandas", "FASTQ", "数据清洗", "生物信息", "入门"]\r
draft: false\r
description: "演示如何用 Python Pandas 快速分析测序数据质量报告。"\r
---\r
\r
# 用 Pandas 处理 FASTQ 质量统计\r
\r
简要示例：读取 \`fastqc_data.txt\` 并提取基本统计：\r
\r
\`\`\`python\r
import pandas as pd\r
df = pd.read_csv("fastqc_data.txt", sep="\\t", comment="#")\r
print(df.head())\r
\`\`\`\r
`;export{n as default};
