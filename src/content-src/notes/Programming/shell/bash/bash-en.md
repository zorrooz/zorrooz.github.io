---
title: "Bash Shell Scripting"
date: "2026-07-06"
author: "zorrooz"
tags: ["Shell", "Bash", "Command Line", "Scripting", "Automation"]
draft: false
description: "Automation Applications of Bash Shell Scripting in Bioinformatics Data Processing"
---

# Bash Shell Scripting

## Basic Syntax

```bash
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
```

## Loop Processing

```bash
# Iterate through files
for file in *.fastq; do
    echo "Processing file: $file"
done

# Read line by line with while
while read -r sample; do
    echo "Processing sample: $sample"
done < sample_list.txt
```

## Functions & Exit Codes

```bash
# Define a wrapper function with error checking
run_alignment() {
    local ref="$1"
    local fq1="$2"
    local fq2="$3"
    bwa mem -t 8 "$ref" "$fq1" "$fq2" || {
        echo "Alignment failed: $fq1 $fq2" >&2
        return 1
    }
}

run_alignment ref.fa r1.fq r2.fq
```

> **Note**: Always use `set -euo pipefail` in scripts to avoid
> undefined variables and silent failures in pipelines.

| Feature | Syntax | Description |
|---------|--------|-------------|
| Variable expansion | `"$var"` | Quote to prevent word splitting |
| Command substitution | `$(cmd)` | Preferred over backticks |
| Arithmetic | `$((a + b))` | Integer arithmetic |
| Parameter expansion | `${var:-default}` | Use default when unset |
