const n=`---\r
title: "BWA MEM Human Genome Alignment Practice"\r
date: "2025-09-19"\r
tags: ["BWA", "Genome Alignment", "NGS", "Command Line", "Human Genome", "SAM/BAM"]\r
draft: false\r
description: "Recording standard usage workflow and parameter optimization of BWA MEM in human WGS data"\r
---\r
\r
# BWA MEM Human Genome Alignment Practice\r
\r
## Basic Commands\r
\r
\`\`\`bash\r
# Index reference genome\r
bwa index ref.fa\r
\r
# Align paired-end sequencing data\r
bwa mem -t 8 ref.fa read1.fq read2.fq > output.sam\r
\r
# Convert to BAM format and sort\r
samtools view -bS output.sam | samtools sort -o output.bam\r
\`\`\`\r
\r
## Parameter Optimization\r
\r
- \`-t\`: Number of threads, adjust according to server configuration\r
- \`-M\`: Mark shorter split hits as secondary\r
- \`-R\`: Add read group information`;export{n as default};
