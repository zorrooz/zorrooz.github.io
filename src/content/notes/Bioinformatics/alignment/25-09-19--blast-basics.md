---
title: "BLAST 序列比对基础"
date: "2025-09-19"
tags: ["BLAST", "序列比对", "数据库搜索", "同源性", "生物信息学"]
draft: false
description: "BLAST 工具在序列比对和同源性分析中的基本应用"
---

# BLAST 序列比对基础

## BLAST 类型

- **blastn**: 核酸序列比对核酸数据库
- **blastp**: 蛋白序列比对蛋白数据库  
- **blastx**: 核酸序列翻译后比对蛋白数据库
- **tblastn**: 蛋白序列比对翻译的核酸数据库

## 基本命令

```bash
# 本地BLAST搜索
blastp -query protein.fa -db nr -out results.txt