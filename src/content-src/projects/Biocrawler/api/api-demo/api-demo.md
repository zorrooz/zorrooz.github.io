---
title: "BioCrawler API 接口文档"
date: "2025-09-01"
author: "zorrooz"
tags: ["API", "Python", "爬虫"]
draft: false
description: "一款生物信息学爬虫工具"
---

# BioCrawler API 接口文档

**以下内容皆为虚构！**

## 接口概述

提供 RESTful API 接口，支持生物数据的查询和下载。

## 核心接口

### 1. 序列查询接口

```
GET /api/v1/sequence/{accession_id}
```

**参数说明：**
- `accession_id`: 数据库登录号（如：NM_001301717）

**返回示例：**
```json
{
  "accession": "NM_001301717",
  "sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
  "length": 38,
  "organism": "Homo sapiens"
}
```

### 2. 批量下载接口

```
POST /api/v1/batch
```

**请求体：**
```json
{
  "accessions": ["NM_001301717", "NM_001301718"],
  "format": "fasta"
}
```

## 认证方式

```bash
# 使用 API Key
curl -H "X-API-Key: your_key_here" \
     https://api.biocrawler.example.com/api/v1/sequence/NM_001301717
```

## 错误处理

| 状态码 | 说明 |
|--------|------|
| 200 | 成功 |
| 404 | 数据不存在 |
| 429 | 请求频率限制 |
| 500 | 服务器错误 |

## 使用限制

- 免费版：每小时 100 次请求
- 专业版：每小时 1000 次请求