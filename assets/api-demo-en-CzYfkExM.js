const r=`---\r
title: "BioCrawler API Documentation"\r
date: "2025-09-01"\r
author: "zorrooz"\r
tags: ["API", "Python", "Crawler"]\r
draft: false\r
description: "A bioinformatics crawler tool"\r
---\r
\r
# BioCrawler API Documentation\r
\r
**The following content is entirely fictional!**\r
\r
## API Overview\r
\r
Provides RESTful API interfaces for querying and downloading biological data.\r
\r
## Core Interfaces\r
\r
### 1. Sequence Query Interface\r
\r
\`\`\`\r
GET /api/v1/sequence/{accession_id}\r
\`\`\`\r
\r
**Parameter Description:**\r
- \`accession_id\`: Database accession number (e.g., NM_001301717)\r
\r
**Response Example:**\r
\`\`\`json\r
{\r
  "accession": "NM_001301717",\r
  "sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",\r
  "length": 38,\r
  "organism": "Homo sapiens"\r
}\r
\`\`\`\r
\r
### 2. Batch Download Interface\r
\r
\`\`\`\r
POST /api/v1/batch\r
\`\`\`\r
\r
**Request Body:**\r
\`\`\`json\r
{\r
  "accessions": ["NM_001301717", "NM_001301718"],\r
  "format": "fasta"\r
}\r
\`\`\`\r
\r
## Authentication Method\r
\r
\`\`\`bash\r
# Using API Key\r
curl -H "X-API-Key: your_key_here" \\\r
     https://api.biocrawler.example.com/api/v1/sequence/NM_001301717\r
\`\`\`\r
\r
## Error Handling\r
\r
| Status Code | Description |\r
|--------|------|\r
| 200 | Success |\r
| 404 | Data Not Found |\r
| 429 | Rate Limit Exceeded |\r
| 500 | Server Error |\r
\r
## Usage Limits\r
\r
- Free Tier: 100 requests per hour\r
- Pro Tier: 1000 requests per hour`;export{r as default};
