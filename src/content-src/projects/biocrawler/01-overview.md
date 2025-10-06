---
title: 设计概览与架构
date: 2025-09-21
tags: [架构, 设计, 可扩展性]
category: [项目, BioCrawler]
---

# 设计概览与架构

本文介绍 BioCrawler 的核心模块、依赖关系与可扩展点。

## 核心模块
1. Scheduler：任务调度与重试
2. Fetcher：请求发送与代理池
3. Parser：结构化解析
4. Pipeline：数据清洗
5. Exporter：数据落地

## 可扩展点
- 自定义解析器
- 自定义清洗规则
- 多种导出格式