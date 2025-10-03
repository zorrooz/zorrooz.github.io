const r=`---\r
title: "BLAST 序列比对基础"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["BLAST", "序列比对", "数据库搜索", "同源性", "生物信息学"]\r
draft: false\r
description: "BLAST 工具在序列比对和同源性分析中的基本应用"\r
---\r
\r
# BLAST 序列比对基础\r
\r
## BLAST 类型\r
\r
- **blastn**: 核酸序列比对核酸数据库\r
- **blastp**: 蛋白序列比对蛋白数据库  \r
- **blastx**: 核酸序列翻译后比对蛋白数据库\r
- **tblastn**: 蛋白序列比对翻译的核酸数据库\r
\r
## 基本命令\r
\r
\`\`\`bash\r
# 本地BLAST搜索\r
blastp -query protein.fa -db nr -out results.txt`;export{r as default};
