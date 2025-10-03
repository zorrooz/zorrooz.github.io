const r=`---\r
title: "Biopython 生物信息学应用"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Python", "Biopython", "生物信息学", "序列分析", "基因组学"]\r
draft: false\r
description: "Biopython 库在序列分析和基因组数据处理中的应用"\r
---\r
\r
# Biopython 生物信息学应用\r
\r
## 序列操作\r
\r
\`\`\`python\r
from Bio import SeqIO\r
from Bio.Seq import Seq\r
\r
# 读取FASTA文件\r
for record in SeqIO.parse("sequences.fasta", "fasta"):\r
    print(record.id, len(record.seq))\r
\r
# 序列翻译\r
dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")\r
protein_seq = dna_seq.translate()`;export{r as default};
