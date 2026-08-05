---
title: "R 语言入门：数据结构、向量化与函数"
date: "2026-08-04"
author: "zorrooz"
tags: ["R语言", "入门", "教程"]
draft: false
description: "R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义"
---

# R 语言入门：数据结构、向量化与函数

R 是统计计算与数据可视化的首选语言。本文讲解 R 的核心：五种数据结构、向量化运算思维、控制流与函数，为后续 tidyverse 与 ggplot2 打下基础。

## 1. 环境准备

### 1.1 安装与 IDE

- 安装 [R](https://cran.r-project.org/)
- 推荐 IDE：[RStudio](https://posit.co/download/rstudio-desktop/)（免费），或 VS Code + R 扩展

### 1.2 基本操作

```r
# 赋值：<- 是 R 的经典赋值符（= 也可用）
x <- 10
x

# 查看帮助
?mean
help("lm")
```

## 2. 基础数据结构

R 的核心数据结构按维度与同质性划分：

| 结构 | 维度 | 元素类型 |
|------|------|----------|
| vector（向量） | 1 | 同质 |
| factor（因子） | 1 | 分类 |
| matrix（矩阵） | 2 | 同质 |
| data.frame（数据框） | 2 | 可异质 |
| list（列表） | N | 可异质 |

### 2.1 向量（vector）

```r
# 创建向量
a <- c(1, 2, 3, 4, 5)            # combine
b <- 1:10                        # 序列
c <- seq(0, 1, by = 0.2)         # 步长序列
d <- rep("A", 3)                 # 重复
e <- sample(1:100, 10)           # 随机抽样

# 类型
typeof(a)      # "double"
is.numeric(a)  # TRUE

# 索引（R 从 1 开始！）
a[1]           # 1
a[c(1, 3)]     # 第 1、3 个元素
a[-2]          # 删除第 2 个
a[a > 3]       # 条件筛选
names(a) <- c("x1", "x2", "x3", "x4", "x5")
a["x1"]        # 按名字取
```

> **R 的索引从 1 开始**，这是与 Python 最大的习惯差异。

### 2.2 因子（factor）

```r
treatment <- factor(c("control", "drug", "control", "drug"),
                    levels = c("control", "drug"))
levels(treatment)          # "control" "drug"
table(treatment)           # 频数统计
```

### 2.3 矩阵（matrix）

```r
m <- matrix(1:12, nrow = 3, ncol = 4)
#      [,1] [,2] [,3] [,4]
# [1,]    1    4    7   10
# [2,]    2    5    8   11
# [3,]    3    6    9   12

m[2, 3]        # 第 2 行第 3 列：8
m[1, ]         # 第 1 行
m[, 2]         # 第 2 列
t(m)           # 转置
m %*% t(m)     # 矩阵乘法
```

### 2.4 数据框（data.frame）

数据框是表格数据的核心结构：

```r
df <- data.frame(
  gene = c("TP53", "BRCA1", "EGFR"),
  expression = c(12.3, 8.9, 45.1),
  chromosome = c(17, 17, 7)
)

# 查看结构
str(df)
summary(df)

# 取列
df$expression          # $ 运算符
df[["expression"]]
df[, "expression"]

# 取行
df[1, ]                # 第一行
df[df$expression > 10, ]

# 增加列
df$log2fc <- log2(df$expression)
```

### 2.5 列表（list）

列表可容纳任意类型：

```r
result <- list(
  name = "分析结果",
  pvalue = 0.003,
  coef = c(1.2, -0.5, 3.1),
  model = lm(expression ~ chromosome, data = df)
)

result$pvalue
result[[3]]            # 数值向量
```

## 3. 向量化运算

**R 的核心思想：对整个向量运算，而非逐元素循环。**

```r
x <- c(1, 2, 3, 4, 5)

x * 2              # 逐元素乘
x + 10
sqrt(x)
log2(x)

# 向量与向量
y <- c(5, 4, 3, 2, 1)
x + y              # 对应位置相加

# 比较运算返回逻辑向量
x > 3              # FALSE FALSE FALSE  TRUE  TRUE

# 统计函数
sum(x); mean(x); median(x); sd(x); range(x)
```

## 4. 控制流

### 4.1 条件

```r
score <- 85

if (score >= 90) {
  grade <- "A"
} else if (score >= 80) {
  grade <- "B"
} else {
  grade <- "C"
}
print(grade)

# 向量化条件：ifelse
values <- c(60, 90, 45, 75)
result <- ifelse(values >= 60, "pass", "fail")
print(result)   # "pass" "pass" "fail" "pass"
```

### 4.2 循环

```r
# for 循环
for (i in 1:5) {
  print(i ^ 2)
}

# 遍历向量元素
genes <- c("TP53", "BRCA1", "EGFR")
for (g in genes) {
  print(paste("分析", g))
}

# while 循环
n <- 0
while (n < 5) {
  n <- n + 1
}
```

### 4.3 尽量避免循环：apply 家族

R 中常用 `apply` 家族替代显式循环：

```r
m <- matrix(1:12, nrow = 3)

apply(m, 1, mean)      # 对每行求均值
apply(m, 2, sum)       # 对每列求和

lapply(list(1:3, 4:6), mean)   # 对列表每个元素操作，返回列表
sapply(list(1:3, 4:6), mean)   # 简化版，返回向量
```

## 5. 函数

### 5.1 定义函数

```r
gc_content <- function(seq) {
  seq <- toupper(seq)
  gc <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
  gc / nchar(seq) * 100
}

gc_content("ATGCCGA")
# [1] 57.14286
```

### 5.2 默认参数与返回值

```r
normalize <- function(x, method = "zscore") {
  if (method == "zscore") {
    (x - mean(x)) / sd(x)
  } else if (method == "minmax") {
    (x - min(x)) / (max(x) - min(x))
  } else {
    stop("未知方法: ", method)
  }
}

normalize(c(1, 2, 3, 4, 5))
normalize(c(1, 2, 3, 4, 5), method = "minmax")
```

### 5.3 匿名函数与函数式编程

```r
# 匿名函数
sapply(1:5, function(i) i ^ 2)

# purrr 风格的 map（需要 tidyverse，下一篇详述）
# purrr::map_dbl(1:5, ~ .x ^ 2)
```

## 6. 实用技巧

```r
# 管道（R 4.1+ 原生管道 |>/ 或 magrittr %>%）
x <- c(1, 2, 3, 4, 5)
x |> mean() |> round(2)

# 缺失值处理
x <- c(1, NA, 3, NA)
is.na(x)
sum(is.na(x))
x[!is.na(x)]
na.omit(x)

# 合并数据框
df1 <- data.frame(id = 1:3, value = c("a", "b", "c"))
df2 <- data.frame(id = 2:4, score = c(90, 85, 88))
merged <- merge(df1, df2, by = "id")   # 类似 SQL join
```

## 7. 小结

- 五种核心结构：向量 / 因子 / 矩阵 / 数据框 / 列表
- **R 索引从 1 开始**；`$` 取列；`[ ]` 取行
- 向量化运算 + `ifelse` + `apply` 家族替代显式循环
- 函数：默认参数、`stop()` 报错、匿名函数

下一篇将介绍 tidyverse：用 `dplyr` 优雅地完成数据操作。
