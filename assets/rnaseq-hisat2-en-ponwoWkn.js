const n=`---\r
title: "RNA-Seq Analysis: HISAT2 Alignment Workflow"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["RNA-Seq", "HISAT2", "Transcriptomics", "Alignment", "Bioinformatics"]\r
draft: false\r
description: "Standard workflow for RNA-Seq data alignment using HISAT2"\r
---\r
\r
# RNA-Seq Analysis: HISAT2 Alignment Workflow\r
\r
## HISAT2 Features\r
\r
- Supports splice-aware alignment\r
- Efficient index structure\r
- Suitable for RNA-Seq reads of various lengths\r
\r
## Basic Commands\r
\r
\`\`\`bash\r
# Build index\r
hisat2-build ref.fa ref_index\r
\r
# Align RNA-Seq data\r
hisat2 -x ref_index -1 read1.fq -2 read2.fq -S output.sam`;export{n as default};
