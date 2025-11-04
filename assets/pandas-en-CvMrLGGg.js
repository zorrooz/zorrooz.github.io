const n=`---\r
title: "Python Pandas Data Processing Practice"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Python", "Pandas", "Data Processing", "Data Analysis", "Data Cleaning"]\r
draft: false\r
description: "Common operations and techniques of Pandas library in bioinformatics data processing"\r
---\r
\r
# Python Pandas Data Processing Practice\r
\r
## Data Reading\r
\r
\`\`\`python\r
import pandas as pd\r
\r
# Read CSV file\r
df = pd.read_csv('data.csv')\r
\r
# Read Excel file\r
df = pd.read_excel('data.xlsx')\r
\`\`\`\r
\r
## Data Cleaning\r
\r
\`\`\`python\r
# Handle missing values\r
df.fillna(0, inplace=True)\r
\r
# Remove duplicate rows\r
df.drop_duplicates(inplace=True)`;export{n as default};
