const n=`---\r
title: "R Language ggplot2 Data Visualization"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["R Language", "ggplot2", "Data Visualization", "Statistical Charts", "Bioinformatics"]\r
draft: false\r
description: "Application of ggplot2 in bioinformatics data visualization"\r
---\r
\r
# R Language ggplot2 Data Visualization\r
\r
## Basic Syntax\r
\r
\`\`\`r\r
library(ggplot2)\r
\r
# Create scatter plot\r
ggplot(data, aes(x = expression, y = condition)) +\r
  geom_point() +\r
  theme_minimal() +\r
  labs(title = "Gene Expression Scatter Plot")\r
\`\`\`\r
\r
## Common Chart Types\r
\r
- Box plot: \`geom_boxplot()\`\r
- Bar chart: \`geom_bar()\`\r
- Line chart: \`geom_line()\``;export{n as default};
