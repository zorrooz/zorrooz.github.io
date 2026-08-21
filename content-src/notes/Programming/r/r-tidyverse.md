---
title: "R tidyverse 数据操作：dplyr 管道与数据清洗"
date: "2026-08-04"
author: "zorrooz"
tags: ["R语言", "tidyverse", "dplyr", "教程"]
draft: false
description: "用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程"
---

# R tidyverse 数据操作：dplyr 管道与数据清洗

tidyverse 是一套风格统一的 R 包集合，其中 `dplyr` 负责数据操作、`tidyr` 负责数据整形。核心思想是**管道串联动词**，让数据处理的每一步都清晰可读。

## 1. 安装与加载

```r
install.packages("tidyverse")   # 一次安装全家桶

library(tidyverse)
# 加载后自动附带: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats, tibble
```

## 2. 管道符 %>%

管道把左边的结果作为右边函数的第一个参数，把"嵌套调用"变成"流水线"：

```r
# 嵌套写法（难读）
round(mean(sqrt(c(1, 2, 3, 4))), 2)

# 管道写法（从左到右读）
c(1, 2, 3, 4) |> sqrt() |> mean() |> round(2)
```

- `%>%`（magrittr）与 `|>`（R 4.1+ 原生）都可用；tidyverse 生态惯用 `%>%`
- 快捷键：RStudio 中 `Ctrl+Shift+M` 插入 `%>%`

## 3. 数据准备

```r
# 创建示例数据：基因表达实验
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

`tibble` 是 tidyverse 的增强数据框：打印更友好、不自动转字符串为因子。

## 4. 核心动词

### 4.1 filter：筛选行

```r
# 表达量大于 20 的样本
expr_data %>%
  filter(expression > 20)

# 多条件：& 且，| 或
expr_data %>%
  filter(condition == "treatment" & expression > 12)

# %in% 匹配集合
expr_data %>%
  filter(gene %in% c("TP53", "MYC"))
```

### 4.2 select：选择列

```r
expr_data %>%
  select(gene, expression)

# 排除列
expr_data %>%
  select(-replicate)

# 辅助选择器
expr_data %>%
  select(starts_with("expr"))      # 列名以 expr 开头
```

### 4.3 mutate：新增/修改列

```r
expr_data %>%
  mutate(
    log2expr = log2(expression),
    group = paste(gene, condition, sep = "_")
  )

# 条件生成列：case_when
expr_data %>%
  mutate(level = case_when(
    expression > 30 ~ "high",
    expression > 10 ~ "medium",
    TRUE            ~ "low"
  ))
```

### 4.4 arrange：排序

```r
expr_data %>%
  arrange(expression)              # 升序

expr_data %>%
  arrange(desc(expression))        # 降序

expr_data %>%
  arrange(gene, desc(expression))  # 多级排序
```

### 4.5 summarise：汇总

```r
expr_data %>%
  summarise(
    mean_expr = mean(expression),
    sd_expr = sd(expression),
    n = n()
  )
```

### 4.6 group_by + summarise：分组汇总（最常用组合）

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

### 4.7 完整管道链

```r
# 每个基因在 treatment 下的平均表达，按降序取前 3
top_genes <- expr_data %>%
  filter(condition == "treatment") %>%
  group_by(gene) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  arrange(desc(mean_expr)) %>%
  slice_max(mean_expr, n = 3)

print(top_genes)
```

## 5. 连接多个表（join）

```r
gene_info <- tibble(
  gene = c("TP53", "BRCA1", "EGFR", "MYC"),
  chromosome = c(17, 17, 7, 8),
  function_desc = c("tumor suppressor", "DNA repair",
                    "receptor kinase", "transcription factor")
)

# 左连接：保留左表所有行
expr_data %>%
  left_join(gene_info, by = "gene")

# 内连接：只保留两边都有的
expr_data %>%
  inner_join(gene_info, by = "gene")

# 右连接 / 全连接
expr_data %>% right_join(gene_info, by = "gene")
expr_data %>% full_join(gene_info, by = "gene")
```

## 6. tidyr：数据整形

### 6.1 pivot_longer：宽表转长表

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

### 6.2 pivot_wider：长表转宽表

```r
long %>%
  pivot_wider(names_from = sample, values_from = expression)
```

### 6.3 缺失值处理

```r
df <- tibble(a = c(1, NA, 3), b = c(NA, 5, 6))

df %>% drop_na()              # 删除含缺失的行
df %>% replace_na(list(a = 0))  # 指定列填充
df %>% fill(a)                # 向下填充
```

## 7. 实战：完整差异分析工作流

```r
# 1. 读取数据
# expr_data <- read_csv("expression.csv")

# 2. 数据清洗
cleaned <- expr_data %>%
  filter(!is.na(expression)) %>%              # 去缺失
  mutate(log2expr = log2(expression)) %>%     # 转换
  left_join(gene_info, by = "gene")           # 注释

# 3. 汇总统计
summary_table <- cleaned %>%
  group_by(gene, condition) %>%
  summarise(
    mean_log2 = mean(log2expr),
    n = n(),
    .groups = "drop"
  )

# 4. 计算处理/对照的差异
diff_table <- summary_table %>%
  pivot_wider(names_from = condition, values_from = mean_log2) %>%
  mutate(log2fc = treatment - control) %>%
  arrange(desc(log2fc))

print(diff_table)
```

## 8. 小结

- **管道**：`%>%` 让数据流从左到右，每一步一个动词
- **dplyr 六大动词**：`filter` / `select` / `mutate` / `arrange` / `summarise` / `group_by`
- **join 家族**：`left_join` / `inner_join` 等按键合并表
- **tidyr**：`pivot_longer` / `pivot_wider` 在宽长表之间转换

下一篇将介绍 ggplot2：用图层语法绘制出版级图表。
