const n=`---
title: "Biopython 生物信息学应用"
date: "2025-09-19"
author: "zorrooz"
tags: ["Python", "Biopython", "生物信息学", "序列分析", "基因组学"]
draft: false
description: "Biopython 库在序列分析和基因组数据处理中的应用"
---

# Biopython 生物信息学应用

## 序列操作

\`\`\`python
from Bio import SeqIO
from Bio.Seq import Seq

# 读取FASTA文件
for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(record.id, len(record.seq))

# 序列翻译
dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
protein_seq = dna_seq.translate()`;export{n as default};
