const n=`---
title: "R语言 ggplot2 数据可视化"
date: "2025-09-19"
author: "zorrooz"
tags: ["R语言", "ggplot2", "数据可视化", "统计图表", "生物信息学"]
draft: false
description: "ggplot2 在生物信息学数据可视化中的应用"
---

# R语言 ggplot2 数据可视化

## 基本语法

\`\`\`r
library(ggplot2)

# 创建散点图
ggplot(data, aes(x = expression, y = condition)) +
  geom_point() +
  theme_minimal() +
  labs(title = "基因表达散点图")
\`\`\`

## 常用图表类型

- 箱线图：\`geom_boxplot()\`
- 柱状图：\`geom_bar()\`
- 线图：\`geom_line()\``;export{n as default};
