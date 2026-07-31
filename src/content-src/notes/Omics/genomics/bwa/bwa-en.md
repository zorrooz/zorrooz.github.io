---
title: "BWA MEM Alignment to Human Genome in Practice"
date: "2026-06-19"
author: "zorrooz"
tags: ["BWA", "Alignment", "NGS", "Command Line", "Human Genome", "SAM/BAM"]
draft: false
description: "Documenting the standard workflow and parameter tuning for BWA MEM in human WGS data"
---

# BWA MEM Alignment to Human Genome in Practice

## Basic Commands

```bash
# Index the reference genome
bwa index ref.fa

# Align paired-end sequencing data
bwa mem -t 8 ref.fa read1.fq read2.fq > output.sam

# Convert to BAM format and sort
samtools view -bS output.sam | samtools sort -o output.bam

# Mark duplicates (recommended before downstream variant calling)
samtools markdup output.bam output.markdup.bam
```

> **Note**: `bwa index` only needs to run once; if the reference genome is unchanged,
> you can reuse the index instead of rebuilding it before every alignment.

## Common Parameters

| Parameter | Description | Suggested Value |
|-----------|-------------|-----------------|
| `-t` | Number of threads, adjust based on server configuration | 8–32 |
| `-M` | Mark shorter split hits as secondary | Recommended for variant calling |
| `-R` | Add read group info (`@RG\tID:xxx\tSM:xxx`) | Required for multi-sample merging |
| `-k` | Minimum seed length | Default 19, lower for short reads |
| `-T` | Do not output alignments scoring below this | Default 30 |

## Parameter Optimization

- `-t`: Number of threads, adjust based on server configuration; pair with `-o` to write BAM directly on I/O-bound jobs
- `-M`: Mark shorter split hits as secondary so tools like GATK can recognize them
- `-R`: Add read group information, needed for multi-sample duplicate marking (`MarkDuplicates`)

## Evaluating Alignment Results

After alignment, it is recommended to check the following metrics:

- Total reads and mapping rate (usually > 95% is considered good)
- Duplicate rate: excessively high (> 20%) warrants library quality investigation
- Insert size distribution: unimodal with the main peak matching the library design

```bash
samtools flagstat output.bam
samtools stats output.bam | grep "insert size average"
```
