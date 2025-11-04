const n=`---\r
title: "MaxQuant Proteomics Analysis Workflow"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["MaxQuant", "Proteomics", "Mass Spectrometry Analysis", "Bioinformatics", "Quantitative Analysis"]\r
draft: false\r
description: "Detailed recording of standard workflow and parameter settings for MaxQuant in proteomics data analysis"\r
---\r
\r
# MaxQuant Proteomics Analysis Workflow\r
\r
## Basic Configuration\r
\r
1. **Input File Settings**\r
   - RAW files: Mass spectrometry raw data\r
   - FASTA files: Protein database\r
\r
2. **Main Parameters**\r
   - Digestion method: Trypsin/P\r
   - Maximum missed cleavages: 2\r
   - Precursor mass tolerance: 20 ppm (MS1)\r
   - Fragment ion tolerance: 0.5 Da (MS2)\r
\r
## Analysis Workflow\r
\r
\`\`\`bash\r
# Run MaxQuant\r
mono MaxQuantCmd.exe mqpar.xml\r
\`\`\`\r
\r
Typical output files:\r
- evidence.txt: Peptide identification results\r
- proteinGroups.txt: Protein quantification results\r
- summary.txt: Analysis statistics summary`;export{n as default};
