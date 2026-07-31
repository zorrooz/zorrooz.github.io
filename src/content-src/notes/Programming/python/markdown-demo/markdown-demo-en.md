---
title: "Markdown & LaTeX Formula Rendering Demo"
date: "2026-07-10"
author: "zorrooz"
tags: ["Markdown", "LaTeX", "KaTeX", "Math", "Bioinformatics"]
draft: false
description: "Demo of Markdown rendering and KaTeX math formula support"
---

# Markdown & LaTeX Formula Rendering Demo

> This article tests the blog's Markdown rendering capabilities and LaTeX math support.
> All content is for demonstration purposes.

## Inline Formulas

Common formulas in bioinformatics can be written inline, e.g. sequencing coverage

$\text{Coverage} = \frac{N_{\text{reads}} \times L_{\text{read}}}{G_{\text{genome}}}$

where $N_{\text{reads}}$ is the total number of reads, $L_{\text{read}}$ is the average read length, and $G_{\text{genome}}$ is the genome size.
Another example: the relationship between base quality and error probability $Q = -10 \log_{10} P_{\text{error}}$.

## Display Formulas

Formulas wrapped in double dollar signs are rendered on their own line:

$$
\text{Phred quality: } Q = -10 \log_{10} P
$$

The probability density function of the Gaussian (normal) distribution:

$$
f(x \mid \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)
$$

## Multi-line Formulas

Use the `aligned` environment for equation systems:

$$
\begin{aligned}
\hat{\mu} &= \frac{1}{n} \sum_{i=1}^{n} x_i \\
\hat{\sigma}^2 &= \frac{1}{n-1} \sum_{i=1}^{n} (x_i - \hat{\mu})^2
\end{aligned}
$$

## Matrices & Special Symbols

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

## Tables

| Element | Syntax | Rendered |
|---------|--------|----------|
| Subscript | `x_i` | $x_i$ |
| Superscript | `x^2` | $x^2$ |
| Fraction | `\frac{a}{b}` | $\frac{a}{b}$ |
| Radical | `\sqrt{x}` | $\sqrt{x}$ |
| Summation | `\sum_{i=1}^n` | $\sum_{i=1}^{n} i = \frac{n(n+1)}{2}$ |

## Blockquotes & Lists

> Formulas are rendered at build time by `remark-math` + `rehype-katex`,
> so no client-side JavaScript is needed on page load.

- [x] Inline formula `$...$`
- [x] Display formula `$$...$$`
- [x] Tables, blockquotes, task lists
- [ ] A pending todo item

## Code Block

```python
import math

def phred_to_probability(q: int) -> float:
    """Convert Phred quality score to error probability"""
    return 10 ** (-q / 10)

print(phred_to_probability(30))  # 0.001
```
