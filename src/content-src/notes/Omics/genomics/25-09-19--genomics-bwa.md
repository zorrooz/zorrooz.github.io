---
title: "BWA MEM 比对人类基因组实战"
date: "2025-09-19"
author: "zorrooz"
tags: ["BWA", "基因组比对", "NGS", "命令行", "人类基因组", "SAM/BAM"]
draft: false
description: "记录 BWA MEM 在人类 WGS 数据中的标准使用流程与参数调优"
---

# BWA MEM 比对人类基因组实战

## 基本命令

```bash
# 索引参考基因组
bwa index ref.fa

# 比对双端测序数据
bwa mem -t 8 ref.fa read1.fq read2.fq > output.sam

# 转换为BAM格式并排序
samtools view -bS output.sam | samtools sort -o output.bam
```

## 参数优化

- `-t`: 线程数，根据服务器配置调整
- `-M`: 标记较短的split hits为secondary
- `-R`: 添加read group信息