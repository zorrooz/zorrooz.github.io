---
title: "BWA MEM Human Genome Alignment Practice"
date: "2025-09-19"
tags: ["BWA", "Genome Alignment", "NGS", "Command Line", "Human Genome", "SAM/BAM"]
draft: false
description: "Recording standard usage workflow and parameter optimization of BWA MEM in human WGS data"
---

# BWA MEM Human Genome Alignment Practice

## Basic Commands

```bash
# Index reference genome
bwa index ref.fa

# Align paired-end sequencing data
bwa mem -t 8 ref.fa read1.fq read2.fq > output.sam

# Convert to BAM format and sort
samtools view -bS output.sam | samtools sort -o output.bam
```

## Parameter Optimization

- `-t`: Number of threads, adjust according to server configuration
- `-M`: Mark shorter split hits as secondary
- `-R`: Add read group information