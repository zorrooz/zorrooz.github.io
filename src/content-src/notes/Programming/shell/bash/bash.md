---
title: "Bash Shell 脚本编程"
date: "2026-07-06"
author: "zorrooz"
tags: ["Shell", "命令行", "脚本编程", "自动化"]
draft: false
description: "Bash Shell 脚本在生物信息学数据处理中的自动化应用"
---

# Bash Shell 脚本编程

## 基本语法

```bash
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
```

## 循环处理

```bash
# 遍历文件
for file in *.fastq; do
    echo "处理文件: $file"
done

# while 逐行读取
while read -r sample; do
    echo "开始处理样本: $sample"
done < sample_list.txt
```

## 函数与退出码

```bash
# 定义函数：带检查的比对包装
run_alignment() {
    local ref="$1"
    local fq1="$2"
    local fq2="$3"
    bwa mem -t 8 "$ref" "$fq1" "$fq2" || {
        echo "比对失败: $fq1 $fq2" >&2
        return 1
    }
}

run_alignment ref.fa r1.fq r2.fq
```

> **提示**：脚本中建议始终使用 `set -euo pipefail`，
> 避免未定义变量和管道中静默失败带来的隐患。

| 特性 | 语法 | 说明 |
|------|------|------|
| 变量引用 | `"$var"` | 加引号防止分词 |
| 命令替换 | `$(cmd)` | 推荐代替反引号 |
| 算术运算 | `$((a + b))` | 整数运算 |
| 参数展开 | `${var:-default}` | 未设置时取默认值 |
