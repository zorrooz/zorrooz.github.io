---
title: "R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning"
date: "2026-08-04"
author: "zorrooz"
tags: ["R Language", "tidyverse", "dplyr", "Tutorial"]
draft: false
description: "Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow."
---

# R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning

tidyverse is a collection of R packages with a consistent design philosophy, where `dplyr` handles data manipulation and `tidyr` handles data reshaping. The core idea is **chaining verbs with pipes**, making every step of data processing clear and readable.

## 1. Installation and Loading

```r
install.packages("tidyverse")   # install the whole suite at once

library(tidyverse)
# After loading, automatically includes: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats, tibble
```

## 2. The Pipe Operator %>%

The pipe passes the result on the left as the first argument of the function on the right, turning "nested calls" into a "pipeline":

```r
# Nested style (hard to read)
round(mean(sqrt(c(1, 2, 3, 4))), 2)

# Pipe style (read from left to right)
c(1, 2, 3, 4) |> sqrt() |> mean() |> round(2)
```

- Both `%>%` (magrittr) and `|>` (native, R 4.1+) are available; the tidyverse ecosystem conventionally uses `%>%`
- Shortcut: In RStudio, press `Ctrl+Shift+M` to insert `%>%`

## 3. Data Preparation

```r
# Create example data: gene expression experiment
expr_data <- tibble(
  gene      = rep(c("TP53", "BRCA1", "EGFR", "MYC"), each = 3),
  condition = rep(c("control", "treatment"), length.out = 12),
  replicate = rep(1:3, 4),
  expression = c(10.2, 11.1, 9.8,    # TP53 control
                 12.5, 13.0, 12.1,   # TP53 treatment
                 8.1, 7.9, 8.4,      # BRCA1
                 45.0, 44.2, 45.8,
                 33.5, 32.9, 33.1)
)
print(expr_data)
```

`tibble` is tidyverse's enhanced data frame: prints more friendlily, and does not automatically convert strings to factors.

## 4. Core Verbs

### 4.1 filter: Select Rows

```r
# Samples with expression greater than 20
expr_data %>%
  filter(expression > 20)

# Multiple conditions: & for AND, | for OR
expr_data %>%
  filter(condition == "treatment" & expression > 12)

# %in% for matching a set
expr_data %>%
  filter(gene %in% c("TP53", "MYC"))
```

### 4.2 select: Choose Columns

```r
expr_data %>%
  select(gene, expression)

# Exclude columns
expr_data %>%
  select(-replicate)

# Helper selectors
expr_data %>%
  select(starts_with("expr"))      # column names starting with "expr"
```

### 4.3 mutate: Add/Modify Columns

```r
expr_data %>%
  mutate(
    log2expr = log2(expression),
    group = paste(gene, condition, sep = "_")
  )

# Conditional column generation: case_when
expr_data %>%
  mutate(level = case_when(
    expression > 30 ~ "high",
    expression > 10 ~ "medium",
    TRUE            ~ "low"
  ))
```

### 4.4 arrange: Sort Rows

```r
expr_data %>%
  arrange(expression)              # ascending

expr_data %>%
  arrange(desc(expression))        # descending

expr_data %>%
  arrange(gene, desc(expression))  # multi-level sorting
```

### 4.5 summarise: Summary Statistics

```r
expr_data %>%
  summarise(
    mean_expr = mean(expression),
    sd_expr = sd(expression),
    n = n()
  )
```

### 4.6 group_by + summarise: Grouped Summaries (Most Common Combination)

```r
expr_data %>%
  group_by(gene, condition) %>%
  summarise(
    mean_expr = mean(expression),
    sd_expr = sd(expression),
    n = n(),
    .groups = "drop"
  )
```

### 4.7 A Full Pipe Chain

```r
# Average expression per gene under treatment, top 3 in descending order
top_genes <- expr_data %>%
  filter(condition == "treatment") %>%
  group_by(gene) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  arrange(desc(mean_expr)) %>%
  slice_max(mean_expr, n = 3)

print(top_genes)
```

## 5. Joining Multiple Tables

```r
gene_info <- tibble(
  gene = c("TP53", "BRCA1", "EGFR", "MYC"),
  chromosome = c(17, 17, 7, 8),
  function_desc = c("tumor suppressor", "DNA repair",
                    "receptor kinase", "transcription factor")
)

# Left join: keeps all rows from the left table
expr_data %>%
  left_join(gene_info, by = "gene")

# Inner join: keeps only rows found in both tables
expr_data %>%
  inner_join(gene_info, by = "gene")

# Right join / Full join
expr_data %>% right_join(gene_info, by = "gene")
expr_data %>% full_join(gene_info, by = "gene")
```

## 6. tidyr: Data Reshaping

### 6.1 pivot_longer: Wide to Long

```r
wide <- tibble(
  gene = c("TP53", "BRCA1"),
  sample1 = c(10.2, 8.1),
  sample2 = c(12.5, 7.9),
  sample3 = c(9.8, 8.4)
)

long <- wide %>%
  pivot_longer(cols = starts_with("sample"),
               names_to = "sample",
               values_to = "expression")
print(long)
```

### 6.2 pivot_wider: Long to Wide

```r
long %>%
  pivot_wider(names_from = sample, values_from = expression)
```

### 6.3 Handling Missing Values

```r
df <- tibble(a = c(1, NA, 3), b = c(NA, 5, 6))

df %>% drop_na()              # remove rows with missing values
df %>% replace_na(list(a = 0))  # fill specified columns
df %>% fill(a)                # fill downward
```

## 7. Practical Example: A Complete Differential Analysis Workflow

```r
# 1. Read data
# expr_data <- read_csv("expression.csv")

# 2. Data cleaning
cleaned <- expr_data %>%
  filter(!is.na(expression)) %>%              # remove missing values
  mutate(log2expr = log2(expression)) %>%     # transformation
  left_join(gene_info, by = "gene")           # annotation

# 3. Summary statistics
summary_table <- cleaned %>%
  group_by(gene, condition) %>%
  summarise(
    mean_log2 = mean(log2expr),
    n = n(),
    .groups = "drop"
  )

# 4. Calculate differences between treatment and control
diff_table <- summary_table %>%
  pivot_wider(names_from = condition, values_from = mean_log2) %>%
  mutate(log2fc = treatment - control) %>%
  arrange(desc(log2fc))

print(diff_table)
```

## 8. Summary

- **Pipes**: `%>%` lets data flow from left to right, one verb per step
- **The six key dplyr verbs**: `filter` / `select` / `mutate` / `arrange` / `summarise` / `group_by`
- **Join family**: `left_join` / `inner_join`, etc., merge tables by keys
- **tidyr**: `pivot_longer` / `pivot_wider` convert between wide and long tables

The next article will introduce ggplot2: creating publication-quality graphics using the grammar of layers.