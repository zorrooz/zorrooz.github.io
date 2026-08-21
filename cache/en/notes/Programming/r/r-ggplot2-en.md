---
title: "R ggplot2 Data Visualization: Grammar of Graphics and Common Charts"
date: "2026-08-04"
author: "zorrooz"
tags: ["R Language", "ggplot2", "Data Visualization", "Tutorial"]
draft: false
description: "Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures"
---

# R ggplot2 Data Visualization: Grammar of Graphics and Common Charts

ggplot2 is R's most powerful plotting package, built on the **Grammar of Graphics**: data, mapping, geometric objects, statistical transformations, coordinates, facets, and themes, stacked layer by layer into a figure.

## 1. Basic Principles

```r
library(ggplot2)

# Skeleton of a plot
ggplot(data = dataset, aes(x = variable, y = variable)) +
  geom_xxx() +          # Geometric layer (point/line/bar/boxplot...)
  labs(...) +           # Labels
  theme_xxx()           # Theme
```

- `data`: data frame
- `aes()`: aesthetic mapping (x, y, color, fill, size, shape)
- `geom_*`: geometric objects, one per layer

## 2. Quick Start with Built-in Datasets

```r
# mpg: car fuel economy data; iris: iris flower data
head(mpg)
head(iris)
```

### 2.1 Scatter Plot

```r
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point()

# Color by category + smooth trend line
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6)
```

### 2.2 Boxplot

```r
ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5)   # Overlay points to show distribution
```

### 2.3 Histogram

```r
ggplot(iris, aes(x = Petal.Length)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white")
```

### 2.4 Bar Chart

```r
# Count bar chart
ggplot(mpg, aes(x = class)) +
  geom_bar()

# Value bar chart (requires aggregation first)
library(dplyr)
avg_hwy <- mpg %>%
  group_by(class) %>%
  summarise(mean_hwy = mean(hwy), .groups = "drop")

ggplot(avg_hwy, aes(x = reorder(class, mean_hwy), y = mean_hwy)) +
  geom_col(fill = "#3b82f6") +
  coord_flip()          # Horizontal, easier to read long category names
```

### 2.5 Line Chart

```r
# Time series example
df <- data.frame(
  time = seq(0, 20, by = 2),
  value = sin(seq(0, 20, by = 2)) + rnorm(11, 0, 0.1)
)

ggplot(df, aes(x = time, y = value)) +
  geom_line(linewidth = 1, color = "#3b82f6") +
  geom_point(size = 2)
```

## 3. Faceting (facet)

Split into multiple subplots by variables:

```r
# Facet by drv: one row
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_wrap(~drv)

# Two-variable faceting: grid
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_grid(cyl ~ drv)
```

## 4. Color Themes

### 4.1 Scale Control

```r
# Continuous color scale
ggplot(mpg, aes(x = displ, y = hwy, color = year)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red")

# Discrete color palette
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Set1")     # Built-in RColorBrewer palette

# Manually specified
scale_color_manual(values = c("#3b82f6", "#10b981", "#f59e0b"))
```

### 4.2 Theme Customization

```r
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 2.5)

# Classic white background theme
p + theme_bw()

# Minimal theme
p + theme_minimal()

# Custom details
p +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",                 # Legend position
    panel.grid.minor = element_blank(),      # Remove minor grid lines
    plot.title = element_text(face = "bold", size = 16)
  )
```

## 5. Labels and Export

```r
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 2.5) +
  labs(
    title = "Iris Measurement Data",
    subtitle = "Sepal size distribution",
    x = "Sepal Length (cm)",
    y = "Sepal Width (cm)",
    color = "Species"
  ) +
  theme_minimal()

# Save: ggsave automatically recognizes the extension
ggsave("iris_scatter.png", p, width = 8, height = 6, dpi = 300)
ggsave("iris_scatter.pdf", p, width = 8, height = 6)
```

> For journal submission, export **PDF vector graphics** (300 dpi PNG for web).

## 6. Bioinformatics Practice: Volcano Plot

Standard visualization for differential expression analysis:

```r
# Simulate differential expression results
set.seed(42)
de <- data.frame(
  gene = paste0("GENE", 1:500),
  log2fc = rnorm(500, 0, 1.8),
  pvalue = runif(500, 0, 1)
)

# Mark significant genes
de <- de %>%
  mutate(significance = case_when(
    abs(log2fc) > 1 & pvalue < 0.05 ~ "up/down",
    pvalue < 0.05 ~ "significant",
    TRUE ~ "ns"
  ))

ggplot(de, aes(x = log2fc, y = -log10(pvalue), color = significance)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("up/down" = "#ef4444", "significant" = "#3b82f6", "ns" = "#9ca3af")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  labs(title = "Differential Expression Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()
```

## 7. Summary

- Grammar of Graphics: `ggplot(data, aes()) + geom_xxx() + labs() + theme_xxx()`
- Common geometries: `geom_point` / `geom_boxplot` / `geom_histogram` / `geom_col` / `geom_line`
- Group by `aes(color/fill)`, split plots by `facet_wrap` / `facet_grid`
- `theme_minimal` + `labs` + `ggsave(dpi=300)` delivers publication-ready output

With this, the R tutorial trilogy is complete: Introduction → tidyverse → ggplot2.