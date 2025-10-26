const n=`---
title: "Bash Shell 脚本编程"
date: "2025-09-19"
author: "zorrooz"
tags: ["Shell", "Bash", "命令行", "脚本编程", "自动化"]
draft: false
description: "Bash Shell 脚本在生物信息学数据处理中的自动化应用"
---

# Bash Shell 脚本编程

## 基本语法

\`\`\`bash
#!/bin/bash

# 变量定义
input_file="data.fastq"
output_dir="results"

# 条件判断
if [ -f "$input_file" ]; then
    echo "文件存在"
else
    echo "文件不存在"
fi
\`\`\`

## 循环处理

\`\`\`bash
# 遍历文件
for file in *.fastq; do
    echo "处理文件: $file"
done`;export{n as default};
