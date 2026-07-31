---
title: "BWA MEM 比对人类基因组实战"
date: "2026-06-19"
author: "zorrooz"
tags: ["BWA", "比对", "NGS", "命令行", "人类基因组", "SAM/BAM"]
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

# 标记重复（下游变异检测前推荐）
samtools markdup output.bam output.markdup.bam
```

> **提示**：`bwa index` 只需要执行一次；如果参考基因组不变，可以复用索引文件，
> 不必在每次比对前重复构建。

## 常用参数

| 参数 | 说明 | 建议值 |
|------|------|--------|
| `-t` | 线程数，根据服务器配置调整 | 8–32 |
| `-M` | 将较短的 split hits 标记为 secondary | 变异检测时建议开启 |
| `-R` | 添加 read group 信息（`@RG\tID:xxx\tSM:xxx`） | 多样本合并必需 |
| `-k` | 最小种子长度 | 默认 19，读长较短时可调低 |
| `-T` | 低于该比对分数的不输出 | 默认 30 |

## 参数优化

- `-t`: 线程数，根据服务器配置调整，I/O 密集时可配合 `-o` 直接输出 BAM
- `-M`: 标记较短的split hits为secondary，方便 GATK 等工具识别
- `-R`: 添加read group信息，多样本标记重复（`MarkDuplicates`）时需要

## 比对结果评估

比对完成后建议检查以下指标：

- 总 reads 数与比对率（通常 > 95% 视为合格）
- 重复率：过高（> 20%）需排查文库质量
- 插入片段分布：单峰且主峰符合文库设计

```bash
samtools flagstat output.bam
samtools stats output.bam | grep "insert size average"
```
