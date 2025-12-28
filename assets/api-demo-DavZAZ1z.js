const r=`---\r
title: "BioCrawler API 接口文档"\r
date: "2025-09-01"\r
author: "zorrooz"\r
tags: ["API", "Python", "爬虫"]\r
draft: false\r
description: "一款生物信息学爬虫工具"\r
---\r
\r
# BioCrawler API 接口文档\r
\r
**以下内容皆为虚构！**\r
\r
## 接口概述\r
\r
提供 RESTful API 接口，支持生物数据的查询和下载。\r
\r
## 核心接口\r
\r
### 1. 序列查询接口\r
\r
\`\`\`\r
GET /api/v1/sequence/{accession_id}\r
\`\`\`\r
\r
**参数说明：**\r
- \`accession_id\`: 数据库登录号（如：NM_001301717）\r
\r
**返回示例：**\r
\`\`\`json\r
{\r
  "accession": "NM_001301717",\r
  "sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",\r
  "length": 38,\r
  "organism": "Homo sapiens"\r
}\r
\`\`\`\r
\r
### 2. 批量下载接口\r
\r
\`\`\`\r
POST /api/v1/batch\r
\`\`\`\r
\r
**请求体：**\r
\`\`\`json\r
{\r
  "accessions": ["NM_001301717", "NM_001301718"],\r
  "format": "fasta"\r
}\r
\`\`\`\r
\r
## 认证方式\r
\r
\`\`\`bash\r
# 使用 API Key\r
curl -H "X-API-Key: your_key_here" \\\r
     https://api.biocrawler.example.com/api/v1/sequence/NM_001301717\r
\`\`\`\r
\r
## 错误处理\r
\r
| 状态码 | 说明 |\r
|--------|------|\r
| 200 | 成功 |\r
| 404 | 数据不存在 |\r
| 429 | 请求频率限制 |\r
| 500 | 服务器错误 |\r
\r
## 使用限制\r
\r
- 免费版：每小时 100 次请求\r
- 专业版：每小时 1000 次请求`;export{r as default};
