---
title: "Markdown 与 LaTeX 公式渲染示例"
date: "2026-07-10"
author: "zorrooz"
tags: ["Markdown", "LaTeX", "KaTeX", "数学公式", "生物信息学"]
draft: false
description: "演示本站 Markdown 渲染能力与 KaTeX 数学公式支持（demo 文章）"
---

# Markdown 与 LaTeX 公式渲染示例

> 本文用于测试博客的 Markdown 渲染能力与 LaTeX 数学公式支持。
> 内容均为演示用途。

## 行内公式

生物信息学中常见的公式可以直接写在正文里，例如样本的覆盖度

$\text{Coverage} = \frac{N_{\text{reads}} \times L_{\text{read}}}{G_{\text{genome}}}$

其中 $N_{\text{reads}}$ 为 reads 总数，$L_{\text{read}}$ 为平均读长，$G_{\text{genome}}$ 为基因组大小。
另一个例子：碱基质量值与错误概率的关系 $Q = -10 \log_{10} P_{\text{error}}$。

## 块级公式

双美元符包裹的公式会单独成行显示：

$$
\text{Phred quality: } Q = -10 \log_{10} P
$$

高斯分布（正态分布）的概率密度函数：

$$
f(x \mid \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)
$$

## 多行公式

使用 `aligned` 环境排版方程组：

$$
\begin{aligned}
\hat{\mu} &= \frac{1}{n} \sum_{i=1}^{n} x_i \\
\hat{\sigma}^2 &= \frac{1}{n-1} \sum_{i=1}^{n} (x_i - \hat{\mu})^2
\end{aligned}
$$

## 矩阵与特殊符号

$$
\mathbf{S} =
\begin{pmatrix}
s_{11} & s_{12} & \cdots & s_{1k} \\
s_{21} & s_{22} & \cdots & s_{2k} \\
\vdots & \vdots & \ddots & \vdots \\
s_{m1} & s_{m2} & \cdots & s_{mk}
\end{pmatrix},
\qquad
\alpha \cdot \beta \cdot \gamma \cdot \delta
$$

## 表格

| 元素 | 符号 | 说明 |
|------|------|------|
| 下标 | `x_i` | $x_i$ |
| 上标 | `x^2` | $x^2$ |
| 分式 | `\frac{a}{b}` | $\frac{a}{b}$ |
| 根式 | `\sqrt{x}` | $\sqrt{x}$ |
| 求和 | `\sum_{i=1}^n` | $\sum_{i=1}^{n} i = \frac{n(n+1)}{2}$ |

## 引用与列表

> KaTeX 支持的公式由 `remark-math` + `rehype-katex` 在构建期渲染，
> 页面加载无需额外 JS 计算。

- [x] 行内公式 `$...$`
- [x] 块级公式 `$$...$$`
- [x] 表格、引用、任务列表
- [ ] 待办事项示例（未完成）

## 代码块

```python
import math

def phred_to_probability(q: int) -> float:
    """将 Phred 质量值转换为错误概率"""
    return 10 ** (-q / 10)

print(phred_to_probability(30))  # 0.001
```
