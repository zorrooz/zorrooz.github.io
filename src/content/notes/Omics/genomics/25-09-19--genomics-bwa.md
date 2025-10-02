---
title: "BWA MEM 比对人类基因组实战"
date: "2025-09-16"
author: "zorrooz"
tags: ["BWA", "基因组比对", "NGS", "命令行", "人类基因组", "SAM/BAM"]
draft: false
description: "记录 BWA MEM 在人类 WGS 数据中的标准使用流程与参数调优。"
---

# BWA MEM 比对人类基因组实战

常用命令：

```bash
bwa mem -t 8 ref.fa read1.fq read2.fq | samtools sort -o output.bam
```
