const r=`---\r
title: "R语言 ggplot2 数据可视化"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["R语言", "ggplot2", "数据可视化", "统计图表", "生物信息学"]\r
draft: false\r
description: "ggplot2 在生物信息学数据可视化中的应用"\r
---\r
\r
# R语言 ggplot2 数据可视化\r
\r
## 基本语法\r
\r
\`\`\`r\r
library(ggplot2)\r
\r
# 创建散点图\r
ggplot(data, aes(x = expression, y = condition)) +\r
  geom_point() +\r
  theme_minimal() +\r
  labs(title = "基因表达散点图")\r
\`\`\`\r
\r
## 常用图表类型\r
\r
- 箱线图：\`geom_boxplot()\`\r
- 柱状图：\`geom_bar()\`\r
- 线图：\`geom_line()\``;export{r as default};
