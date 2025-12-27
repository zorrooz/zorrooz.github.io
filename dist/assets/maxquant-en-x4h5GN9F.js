const n=`---
title: "MaxQuant Proteomics Analysis Workflow"
date: "2025-09-19"
author: "zorrooz"
tags: ["MaxQuant", "Proteomics", "Mass Spectrometry", "Bioinformatics", "Quantitative Analysis"]
draft: false
description: "Detailed documentation of standard procedures and parameter settings for MaxQuant in proteomics data analysis"
---

# MaxQuant Proteomics Analysis Workflow

## Basic Configuration

1. **Input File Settings**
   - RAW files: Mass spectrometry raw data
   - FASTA files: Protein database

2. **Main Parameters**
   - Enzyme: Trypsin/P
   - Maximum missed cleavages: 2
   - Precursor mass tolerance: 20 ppm (MS1)
   - Fragment mass tolerance: 0.5 Da (MS2)

## Analysis Workflow

\`\`\`bash
# Run MaxQuant
mono MaxQuantCmd.exe mqpar.xml
\`\`\`

Typical output files:
- evidence.txt: Peptide identification results
- proteinGroups.txt: Protein quantification results
- summary.txt: Analysis statistics summary`;export{n as default};
