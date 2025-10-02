const r=`---\r
title: "BWA MEM 比对人类基因组实战"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["BWA", "基因组比对", "NGS", "命令行", "人类基因组", "SAM/BAM"]\r
draft: false\r
description: "记录 BWA MEM 在人类 WGS 数据中的标准使用流程与参数调优"\r
---\r
\r
# BWA MEM 比对人类基因组实战\r
\r
## 基本命令\r
\r
\`\`\`bash\r
# 索引参考基因组\r
bwa index ref.fa\r
\r
# 比对双端测序数据\r
bwa mem -t 8 ref.fa read1.fq read2.fq > output.sam\r
\r
# 转换为BAM格式并排序\r
samtools view -bS output.sam | samtools sort -o output.bam\r
\`\`\`\r
\r
## 参数优化\r
\r
- \`-t\`: 线程数，根据服务器配置调整\r
- \`-M\`: 标记较短的split hits为secondary\r
- \`-R\`: 添加read group信息`;export{r as default};
