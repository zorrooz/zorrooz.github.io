---
title: "MaxQuant 蛋白质组学分析工作流"
date: "2025-09-19"
author: "zorrooz"
tags: ["MaxQuant", "蛋白质组学", "质谱分析", "生物信息学", "定量分析"]
draft: false
description: "详细记录 MaxQuant 在蛋白质组学数据分析中的标准流程与参数设置"
---

# MaxQuant 蛋白质组学分析工作流

## 基本配置

1. **输入文件设置**
   - RAW 文件：质谱原始数据
   - FASTA 文件：蛋白质数据库

2. **主要参数**
   - 酶切方式：Trypsin/P
   - 最大漏切位点：2
   - 前体质量容差：20 ppm（一级）
   - 碎片离子容差：0.5 Da（二级）

## 分析流程

```bash
# 运行 MaxQuant
mono MaxQuantCmd.exe mqpar.xml
```

典型输出文件：
- evidence.txt：肽段鉴定结果
- proteinGroups.txt：蛋白质定量结果
- summary.txt：分析统计摘要