---
title: BioCrawler 项目说明
date: 2025-09-20
tags: [生物信息学, 爬虫, 数据获取]
category: [项目, BioCrawler]
---

# BioCrawler 项目说明

BioCrawler 是一个用于生物数据网站的通用抓取与清洗工具，支持任务编排、断点续抓与结构化导出。

## 功能特点
- 灵活的任务编排与多站点适配
- 丰富的抓取策略（分页、并发、重试、代理池）
- 数据清洗与校验流水线
- 结构化导出（CSV/JSON/SQLite）

## 快速开始
```bash
# 安装依赖与启动
pip install -r requirements.txt
python -m biocrawler run --job ncbi-pubmed
```

## 目录结构
- jobs/ 任务配置
- spiders/ 站点适配
- pipelines/ 数据清洗
- exporters/ 导出器