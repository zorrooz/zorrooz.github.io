---
title: "R Language ggplot2 Data Visualization"
date: "2025-09-19"
author: "zorrooz"
tags: ["R Language", "ggplot2", "Data Visualization", "Statistical Charts", "Bioinformatics"]
draft: false
description: "Application of ggplot2 in bioinformatics data visualization"
---

# R Language ggplot2 Data Visualization

## Basic Syntax

```r
library(ggplot2)

# Create scatter plot
ggplot(data, aes(x = expression, y = condition)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Gene Expression Scatter Plot")
```

## Common Chart Types

- Box plot: `geom_boxplot()`
- Bar chart: `geom_bar()`
- Line chart: `geom_line()`