---
title: "R语言 ggplot2 数据可视化"
date: "2026-07-05"
author: "zorrooz"
tags: ["R语言", "ggplot2", "数据可视化", "生物信息学"]
draft: false
description: "ggplot2 在生物信息学数据可视化中的应用"
---

# R语言 ggplot2 数据可视化

## 基本语法

ggplot2 的核心思想是**图层叠加**：数据 → 映射 → 几何对象 → 统计变换 → 主题。

```r
library(ggplot2)

# 创建散点图
ggplot(data, aes(x = expression, y = condition)) +
  geom_point() +
  theme_minimal() +
  labs(title = "基因表达散点图")
```

## 常用图表类型

| 图表 | 几何对象 | 适用场景 |
|------|----------|----------|
| 箱线图 | `geom_boxplot()` | 分组比较分布 |
| 柱状图 | `geom_bar()` | 类别计数 |
| 线图 | `geom_line()` | 时间序列 / 曲线 |
| 热图 | `geom_tile()` | 表达矩阵可视化 |
| 小提琴图 | `geom_violin()` | 分布 + 密度 |

## 分面与主题

```r
# 按分组分面，观察不同条件下的趋势
ggplot(data, aes(x = time, y = expression)) +
  geom_line(aes(color = group)) +
  facet_wrap(~ gene) +
  theme_bw(base_size = 14)
```

> **提示**：生物信息学出图前先把数据整理为**长格式**
> （`tidyr::pivot_longer`），ggplot2 的分组、分面都依赖长格式数据。

## 保存图片

```r
# 支持多种格式，推荐 PDF 用于论文，PNG 用于预览
ggsave("plot.pdf", width = 6, height = 4, dpi = 300)
ggsave("plot.png", width = 6, height = 4, dpi = 300)
```
