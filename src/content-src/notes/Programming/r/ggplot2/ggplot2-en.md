---
title: "R Language ggplot2 Data Visualization"
date: "2026-07-05"
author: "zorrooz"
tags: ["R Language", "ggplot2", "Data Visualization", "Bioinformatics"]
draft: false
description: "Application of ggplot2 in Bioinformatics Data Visualization"
---

# R Language ggplot2 Data Visualization

## Basic Syntax

The core idea of ggplot2 is **layered grammar**: data → mapping → geometric objects → statistical transformation → theme.

```r
library(ggplot2)

# Create scatter plot
ggplot(data, aes(x = expression, y = condition)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Gene Expression Scatter Plot")
```

## Common Chart Types

| Chart | Geometric Object | Use Case |
|-------|------------------|----------|
| Box plot | `geom_boxplot()` | Compare distributions across groups |
| Bar chart | `geom_bar()` | Count categories |
| Line chart | `geom_line()` | Time series / curves |
| Heatmap | `geom_tile()` | Expression matrix visualization |
| Violin plot | `geom_violin()` | Distribution + density |

## Faceting & Themes

```r
# Facet by group to inspect trends across conditions
ggplot(data, aes(x = time, y = expression)) +
  geom_line(aes(color = group)) +
  facet_wrap(~ gene) +
  theme_bw(base_size = 14)
```

> **Note**: Reshape your data into **long format** first
> (`tidyr::pivot_longer`) — grouping and faceting in ggplot2 rely on long data.

## Saving Figures

```r
# Multiple formats: PDF for papers, PNG for preview
ggsave("plot.pdf", width = 6, height = 4, dpi = 300)
ggsave("plot.png", width = 6, height = 4, dpi = 300)
```
