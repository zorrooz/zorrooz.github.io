const t=`---\r
title: "Processing FASTQ Quality Statistics with Pandas"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Python", "Pandas", "FASTQ", "Data Cleaning", "Bioinformatics", "Beginner"]\r
draft: false\r
description: "Demonstrating how to quickly analyze sequencing data quality reports using Python Pandas"\r
---\r
\r
# Processing FASTQ Quality Statistics with Pandas\r
\r
Brief example: Reading \`fastqc_data.txt\` and extracting basic statistics:\r
\r
\`\`\`python\r
import pandas as pd\r
df = pd.read_csv("fastqc_data.txt", sep="\\t", comment="#")\r
print(df.head())`;export{t as default};
