const r=`---\r
title: "R语言数据处理示例"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["R语言", "数据处理", "统计分析", "数据可视化", "生物信息学"]\r
draft: false\r
description: "展示 R 语言在生物信息学数据处理中的基本应用"\r
---\r
\r
# R语言数据处理示例\r
\r
## 数据读取与基本操作\r
\r
\`\`\`r\r
# 读取 CSV 文件\r
data <- read.csv("expression_data.csv")\r
\r
# 查看数据结构\r
str(data)\r
head(data)\r
\r
# 基本统计\r
summary(data)\r
\`\`\`\r
\r
## 数据可视化\r
\r
\`\`\`r\r
# 安装并加载 ggplot2\r
install.packages("ggplot2")\r
library(ggplot2)\r
\r
# 创建散点图\r
ggplot(data, aes(x = Condition, y = Expression)) +\r
  geom_boxplot() +\r
  theme_minimal() +\r
  labs(title = "基因表达水平比较")\r
\`\`\`\r
\r
## 统计分析\r
\r
\`\`\`r\r
# t检验\r
t_test_result <- t.test(Expression ~ Condition, data = data)\r
print(t_test_result)`;export{r as default};
