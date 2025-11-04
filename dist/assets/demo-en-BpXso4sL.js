const a=`---\r
title: "R Language Data Processing Example"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["R Language", "Data Processing", "Statistical Analysis", "Data Visualization", "Bioinformatics"]\r
draft: false\r
description: "Demonstrating basic applications of R language in bioinformatics data processing"\r
---\r
\r
# R Language Data Processing Example\r
\r
## Data Reading and Basic Operations\r
\r
\`\`\`r\r
# Read CSV file\r
data <- read.csv("expression_data.csv")\r
\r
# View data structure\r
str(data)\r
head(data)\r
\r
# Basic statistics\r
summary(data)\r
\`\`\`\r
\r
## Data Visualization\r
\r
\`\`\`r\r
# Install and load ggplot2\r
install.packages("ggplot2")\r
library(ggplot2)\r
\r
# Create scatter plot\r
ggplot(data, aes(x = Condition, y = Expression)) +\r
  geom_boxplot() +\r
  theme_minimal() +\r
  labs(title = "Gene Expression Level Comparison")\r
\`\`\`\r
\r
## Statistical Analysis\r
\r
\`\`\`r\r
# t-test\r
t_test_result <- t.test(Expression ~ Condition, data = data)\r
print(t_test_result)`;export{a as default};
