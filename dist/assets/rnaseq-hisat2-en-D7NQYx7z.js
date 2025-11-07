const n=`---
title: "RNA-Seq Analysis: HISAT2 Alignment Workflow"
date: "2025-09-19"
author: "zorrooz"
tags: ["RNA-Seq", "HISAT2", "Transcriptome", "Alignment", "Bioinformatics"]
draft: false
description: "Standard workflow for aligning RNA-Seq data using HISAT2"
---

# RNA-Seq Analysis: HISAT2 Alignment Workflow

## HISAT2 Features

- Supports splice-aware alignment
- Efficient index structure
- Suitable for RNA-Seq reads of various lengths

## Basic Commands

\`\`\`bash
# Build index
hisat2-build ref.fa ref_index

# Align RNA-Seq data
hisat2 -x ref_index -1 read1.fq -2 read2.fq -S output.sam`;export{n as default};
