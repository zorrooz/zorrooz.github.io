---
title: "R语言数据处理示例"
date: "2025-09-19"
author: "zorrooz"
tags: ["R语言", "数据处理", "统计分析", "数据可视化", "生物信息学"]
draft: false
description: "展示 R 语言在生物信息学数据处理中的基本应用"
---

# R语言数据处理示例

## 数据读取与基本操作

```r
# 读取 CSV 文件
data <- read.csv("expression_data.csv")

# 查看数据结构
str(data)
head(data)

# 基本统计
summary(data)
```

## 数据可视化

```r
# 安装并加载 ggplot2
install.packages("ggplot2")
library(ggplot2)

# 创建散点图
ggplot(data, aes(x = Condition, y = Expression)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "基因表达水平比较")
```

## 统计分析

```r
# t检验
t_test_result <- t.test(Expression ~ Condition, data = data)
print(t_test_result)