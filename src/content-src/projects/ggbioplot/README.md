---
title: ggBioPlot 项目介绍
date: 2025-09-22
tags: [R语言, 可视化, ggplot2]
category: [项目, ggBioPlot]
---

# ggBioPlot 项目介绍

ggBioPlot 是一套面向生物数据可视化的 ggplot2 扩展主题与组件集合，强调简洁与发表级图形。

## 特性
- 统一主题与字体
- 常用生物学图型组件
- 语义化配色方案

## 示例
```r
library(ggplot2)
library(ggBioPlot)

ggplot(df, aes(x, y)) +
  geom_point() +
  theme_bio()