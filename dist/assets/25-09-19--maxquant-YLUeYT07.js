const n=`---\r
title: "MaxQuant 蛋白质组学分析工作流"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["MaxQuant", "蛋白质组学", "质谱分析", "生物信息学", "定量分析"]\r
draft: false\r
description: "详细记录 MaxQuant 在蛋白质组学数据分析中的标准流程与参数设置"\r
---\r
\r
# MaxQuant 蛋白质组学分析工作流\r
\r
## 基本配置\r
\r
1. **输入文件设置**\r
   - RAW 文件：质谱原始数据\r
   - FASTA 文件：蛋白质数据库\r
\r
2. **主要参数**\r
   - 酶切方式：Trypsin/P\r
   - 最大漏切位点：2\r
   - 前体质量容差：20 ppm（一级）\r
   - 碎片离子容差：0.5 Da（二级）\r
\r
## 分析流程\r
\r
\`\`\`bash\r
# 运行 MaxQuant\r
mono MaxQuantCmd.exe mqpar.xml\r
\`\`\`\r
\r
典型输出文件：\r
- evidence.txt：肽段鉴定结果\r
- proteinGroups.txt：蛋白质定量结果\r
- summary.txt：分析统计摘要`;export{n as default};
