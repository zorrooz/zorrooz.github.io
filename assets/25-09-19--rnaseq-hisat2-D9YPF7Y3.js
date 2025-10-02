const r=`---\r
title: "RNA-Seq 分析：HISAT2 比对流程"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["RNA-Seq", "HISAT2", "转录组", "比对", "生物信息学"]\r
draft: false\r
description: "使用 HISAT2 进行 RNA-Seq 数据比对的标准工作流程"\r
---\r
\r
# RNA-Seq 分析：HISAT2 比对流程\r
\r
## HISAT2 特点\r
\r
- 支持拼接感知的比对\r
- 高效的索引结构\r
- 适用于各种长度的RNA-Seq reads\r
\r
## 基本命令\r
\r
\`\`\`bash\r
# 构建索引\r
hisat2-build ref.fa ref_index\r
\r
# 比对RNA-Seq数据\r
hisat2 -x ref_index -1 read1.fq -2 read2.fq -S output.sam`;export{r as default};
