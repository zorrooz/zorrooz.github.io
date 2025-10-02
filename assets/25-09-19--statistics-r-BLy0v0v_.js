const r=`---\r
title: "R语言统计分析基础"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["R语言", "统计分析", "假设检验", "回归分析", "数据科学"]\r
draft: false\r
description: "R语言在生物统计学和数据分析中的基本统计方法"\r
---\r
\r
# R语言统计分析基础\r
\r
## 描述性统计\r
\r
\`\`\`r\r
# 基本统计量\r
summary(data)\r
\r
# 相关性分析\r
cor(data$var1, data$var2)\r
\r
# t检验\r
t.test(group1, group2)\r
\`\`\`\r
\r
## 假设检验\r
\r
- 正态性检验：\`shapiro.test()\`\r
- 方差齐性检验：\`bartlett.test()\`\r
- 卡方检验：\`chisq.test()\``;export{r as default};
