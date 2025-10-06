---
title: "R语言统计分析基础"
date: "2025-09-19"
tags: ["R语言", "统计分析", "假设检验", "回归分析", "数据科学"]
draft: false
description: "R语言在生物统计学和数据分析中的基本统计方法"
---

# R语言统计分析基础

## 描述性统计

```r
# 基本统计量
summary(data)

# 相关性分析
cor(data$var1, data$var2)

# t检验
t.test(group1, group2)
```

## 假设检验

- 正态性检验：`shapiro.test()`
- 方差齐性检验：`bartlett.test()`
- 卡方检验：`chisq.test()`