const n=`---
title: "R ggplot2 数据可视化：图层语法与常用图表"
date: "2026-08-04"
author: "zorrooz"
tags: ["R语言", "ggplot2", "数据可视化", "教程"]
draft: false
description: "ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表"
---

# R ggplot2 数据可视化：图层语法与常用图表

ggplot2 是 R 最强大的绘图包，基于**图层语法**（Grammar of Graphics）：数据、映射、几何对象、统计变换、坐标、分面、主题，逐层叠加成图。

## 1. 基本原理

\`\`\`r
library(ggplot2)

# 一张图的骨架
ggplot(data = 数据集, aes(x = 变量, y = 变量)) +
  geom_xxx() +          # 几何图层（点/线/柱/箱线...）
  labs(...) +           # 标签
  theme_xxx()           # 主题
\`\`\`

- \`data\`：数据框
- \`aes()\`：美学映射（x、y、color、fill、size、shape）
- \`geom_*\`：几何对象，每层一个

## 2. 内置数据集快速上手

\`\`\`r
# mpg：汽车油耗数据；iris：鸢尾花数据
head(mpg)
head(iris)
\`\`\`

### 2.1 散点图（scatter）

\`\`\`r
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point()

# 按类别着色 + 平滑趋势线
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6)
\`\`\`

### 2.2 箱线图（boxplot）

\`\`\`r
ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5)   # 叠加散点展示分布
\`\`\`

### 2.3 直方图（histogram）

\`\`\`r
ggplot(iris, aes(x = Petal.Length)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white")
\`\`\`

### 2.4 柱状图（bar）

\`\`\`r
# 计数柱状图
ggplot(mpg, aes(x = class)) +
  geom_bar()

# 数值柱状图（需要先聚合）
library(dplyr)
avg_hwy <- mpg %>%
  group_by(class) %>%
  summarise(mean_hwy = mean(hwy), .groups = "drop")

ggplot(avg_hwy, aes(x = reorder(class, mean_hwy), y = mean_hwy)) +
  geom_col(fill = "#3b82f6") +
  coord_flip()          # 横置，长类别名更易读
\`\`\`

### 2.5 折线图（line）

\`\`\`r
# 时间序列示例
df <- data.frame(
  time = seq(0, 20, by = 2),
  value = sin(seq(0, 20, by = 2)) + rnorm(11, 0, 0.1)
)

ggplot(df, aes(x = time, y = value)) +
  geom_line(linewidth = 1, color = "#3b82f6") +
  geom_point(size = 2)
\`\`\`

## 3. 分面（facet）

按变量拆分成多个子图：

\`\`\`r
# 按 drv 分面：一排
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_wrap(~drv)

# 双变量分面：网格
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_grid(cyl ~ drv)
\`\`\`

## 4. 颜色主题

### 4.1 标度控制（scale）

\`\`\`r
# 连续色阶
ggplot(mpg, aes(x = displ, y = hwy, color = year)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red")

# 离散色板
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Set1")     # 内置 RColorBrewer 色板

# 手动指定
scale_color_manual(values = c("#3b82f6", "#10b981", "#f59e0b"))
\`\`\`

### 4.2 主题美化

\`\`\`r
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 2.5)

# 经典白底主题
p + theme_bw()

# 极简主题
p + theme_minimal()

# 自定义细节
p +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",                 # 图例位置
    panel.grid.minor = element_blank(),      # 去掉次要网格线
    plot.title = element_text(face = "bold", size = 16)
  )
\`\`\`

## 5. 标签与导出

\`\`\`r
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 2.5) +
  labs(
    title = "鸢尾花测量数据",
    subtitle = "Sepal 尺寸分布",
    x = "花萼长度 (cm)",
    y = "花萼宽度 (cm)",
    color = "物种"
  ) +
  theme_minimal()

# 保存：ggsave 自动识别扩展名
ggsave("iris_scatter.png", p, width = 8, height = 6, dpi = 300)
ggsave("iris_scatter.pdf", p, width = 8, height = 6)
\`\`\`

> 论文投稿建议导出 **PDF 矢量图**（300 dpi 的 PNG 用于网页）。

## 6. 生物信息学实战：火山图

差异表达分析的标准可视化：

\`\`\`r
# 模拟差异表达结果
set.seed(42)
de <- data.frame(
  gene = paste0("GENE", 1:500),
  log2fc = rnorm(500, 0, 1.8),
  pvalue = runif(500, 0, 1)
)

# 标记显著基因
de <- de %>%
  mutate(significance = case_when(
    abs(log2fc) > 1 & pvalue < 0.05 ~ "up/down",
    pvalue < 0.05 ~ "significant",
    TRUE ~ "ns"
  ))

ggplot(de, aes(x = log2fc, y = -log10(pvalue), color = significance)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("up/down" = "#ef4444", "significant" = "#3b82f6", "ns" = "#9ca3af")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  labs(title = "差异表达火山图", x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()
\`\`\`

## 7. 小结

- 图层语法：\`ggplot(data, aes()) + geom_xxx() + labs() + theme_xxx()\`
- 常用几何：\`geom_point\` / \`geom_boxplot\` / \`geom_histogram\` / \`geom_col\` / \`geom_line\`
- 分组用 \`aes(color/fill)\`，拆图用 \`facet_wrap\` / \`facet_grid\`
- 主题 \`theme_minimal\` + \`labs\` + \`ggsave(dpi=300)\` 即出版级输出

至此 R 教程三部曲完成：入门 → tidyverse → ggplot2。
`;export{n as default};
