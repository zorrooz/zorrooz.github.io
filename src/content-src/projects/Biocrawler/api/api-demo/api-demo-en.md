---
title: "BioCrawler API Documentation"
date: "2025-09-01"
author: "zorrooz"
tags: ["API", "Python", "Crawler"]
draft: false
description: "A bioinformatics crawler tool"
---

# BioCrawler API Documentation

**The following content is entirely fictional!**

## API Overview

Provides RESTful API interfaces for querying and downloading biological data.

## Core Interfaces

### 1. Sequence Query Interface

```
GET /api/v1/sequence/{accession_id}
```

**Parameter Description:**
- `accession_id`: Database accession number (e.g., NM_001301717)

**Response Example:**
```json
{
  "accession": "NM_001301717",
  "sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
  "length": 38,
  "organism": "Homo sapiens"
}
```

### 2. Batch Download Interface

```
POST /api/v1/batch
```

**Request Body:**
```json
{
  "accessions": ["NM_001301717", "NM_001301718"],
  "format": "fasta"
}
```

## Authentication Method

```bash
# Using API Key
curl -H "X-API-Key: your_key_here" \
     https://api.biocrawler.example.com/api/v1/sequence/NM_001301717
```

## Error Handling

| Status Code | Description |
|--------|------|
| 200 | Success |
| 404 | Data Not Found |
| 429 | Rate Limit Exceeded |
| 500 | Server Error |

## Usage Limits

- Free Tier: 100 requests per hour
- Pro Tier: 1000 requests per hour