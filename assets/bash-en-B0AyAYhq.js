const n=`---
title: "Bash Shell Scripting"
date: "2025-09-19"
author: "zorrooz"
tags: ["Shell", "Bash", "Command Line", "Scripting", "Automation"]
draft: false
description: "Automation Applications of Bash Shell Scripting in Bioinformatics Data Processing"
---

# Bash Shell Scripting

## Basic Syntax

\`\`\`bash
#!/bin/bash

# Variable definition
input_file="data.fastq"
output_dir="results"

# Conditional judgment
if [ -f "$input_file" ]; then
    echo "File exists"
else
    echo "File does not exist"
fi
\`\`\`

## Loop Processing

\`\`\`bash
# Iterate through files
for file in *.fastq; do
    echo "Processing file: $file"
done`;export{n as default};
