---
title: "MaxQuant 分析 DIA 数据参数详解"
date: "2025-09-14"
author: "zorrooz"
tags: ["MaxQuant", "蛋白质组", "DIA", "质谱", "参数配置", "定量分析"]
draft: false
description: "详解 MaxQuant 处理 DIA 蛋白质组数据时的关键参数设置与常见坑点。"
---

# MaxQuant 分析 DIA 数据参数详解

关键参数：

- `Match between runs`: ✅ 启用提升定量覆盖
- `LFQ`: ✅ 标记用于无标定量
- `Enzyme`: Trypsin/P

建议使用默认 FDR 0.01，避免过度放宽。
