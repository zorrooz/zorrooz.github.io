---
title: "BWA MEM Alignment to Human Genome in Practice"
date: "2025-09-19"
tags: ["BWA", "Genome Alignment", "NGS", "Command Line", "Human Genome", "SAM/BAM"]
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
```

## Parameter Optimization

- `-t`: Number of threads, adjust based on server configuration
- `-M`: Mark shorter split hits as secondary
- `-R`: Add read group information