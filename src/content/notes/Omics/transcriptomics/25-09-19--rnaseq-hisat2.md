---
title: "RNA-Seq 分析：HISAT2 比对流程"
date: "2025-09-19"
tags: ["RNA-Seq", "HISAT2", "转录组", "比对", "生物信息学"]
draft: false
description: "使用 HISAT2 进行 RNA-Seq 数据比对的标准工作流程"
---

# RNA-Seq 分析：HISAT2 比对流程

## HISAT2 特点

- 支持拼接感知的比对
- 高效的索引结构
- 适用于各种长度的RNA-Seq reads

## 基本命令

```bash
# 构建索引
hisat2-build ref.fa ref_index

# 比对RNA-Seq数据
hisat2 -x ref_index -1 read1.fq -2 read2.fq -S output.sam