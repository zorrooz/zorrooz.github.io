const a=`---
title: "R Language Data Processing Example"
date: "2025-09-19"
author: "zorrooz"
tags: ["R Language", "Data Processing", "Statistical Analysis", "Data Visualization", "Bioinformatics"]
draft: false
description: "Demonstrating basic applications of R language in bioinformatics data processing"
---

# R Language Data Processing Example

## Data Reading and Basic Operations

\`\`\`r
# Read CSV file
data <- read.csv("expression_data.csv")

# View data structure
str(data)
head(data)

# Basic statistics
summary(data)
\`\`\`

## Data Visualization

\`\`\`r
# Install and load ggplot2
install.packages("ggplot2")
library(ggplot2)

# Create scatter plot
ggplot(data, aes(x = Condition, y = Expression)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Gene Expression Level Comparison")
\`\`\`

## Statistical Analysis

\`\`\`r
# t-test
t_test_result <- t.test(Expression ~ Condition, data = data)
print(t_test_result)
\`\`\``;export{a as default};
