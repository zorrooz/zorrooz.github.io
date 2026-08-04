const n=`---
title: "Bash 编程：变量、条件、循环与实用脚本"
date: "2026-08-04"
author: "zorrooz"
tags: ["Bash", "Shell", "脚本编程", "教程"]
draft: false
description: "Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本"
---

# Bash 编程：变量、条件、循环与实用脚本

Bash 脚本把命令行操作组织成可复用的程序，是自动化批处理的核心工具。本文从脚本基础到实战，覆盖变量、条件、循环、函数与常见任务。

## 1. 脚本基础

### 1.1 脚本结构

\`\`\`bash
#!/bin/bash
# 第一行 shebang：指定解释器

# 给脚本执行权限
chmod +x myscript.sh
# 运行
./myscript.sh
# 或用 bash 直接运行（无需执行权限）
bash myscript.sh
\`\`\`

### 1.2 调试

\`\`\`bash
bash -x script.sh    # 跟踪执行，显示每条命令
# 脚本内开启：set -x（跟踪）/ set -e（出错即退出）/ set -u（未定义变量报错）
\`\`\`

生产脚本推荐在开头加：

\`\`\`bash
#!/bin/bash
set -euo pipefail   # 出错即停、未定义变量报错、管道失败检测
\`\`\`

## 2. 变量

### 2.1 定义与引用

\`\`\`bash
name="zorrooz"        # 注意：= 两边不能有空格
echo "Hello, $name"   # 双引号内变量展开
echo 'Hello, $name'   # 单引号内原样输出
echo "长度: \${#name}" # 字符串长度
\`\`\`

### 2.2 命令替换与算术

\`\`\`bash
# 命令替换：$(...)
today=$(date +%F)
echo "今天是 $today"
files=$(ls | wc -l)
echo "共 $files 个文件"

# 算术：$((...))
a=10
b=3
echo $((a + b))    # 13
echo $((a % b))    # 1
\`\`\`

### 2.3 位置参数

\`\`\`bash
#!/bin/bash
echo "脚本名: $0"
echo "第一个参数: $1"
echo "第二个参数: $2"
echo "所有参数: $@"
echo "参数个数: $#"

# 运行：./script.sh file1 file2
\`\`\`

### 2.4 特殊变量

\`\`\`bash
$?    # 上一条命令的退出码（0 成功，非 0 失败）
$$    # 当前进程 PID
$!    # 后台进程 PID
\`\`\`

## 3. 条件判断

### 3.1 if / elif / else

\`\`\`bash
#!/bin/bash
score=$1

if (( score >= 90 )); then
    echo "优秀"
elif (( score >= 60 )); then
    echo "及格"
else
    echo "不及格"
fi
\`\`\`

### 3.2 test 表达式

\`\`\`bash
# 文件判断
[ -f file.txt ] && echo "是普通文件"
[ -d dir ] && echo "是目录"
[ -e path ] && echo "存在"
[ -s file ] && echo "非空文件"

# 字符串比较
[ "$name" = "zorrooz" ] && echo "匹配"
[ -z "$var" ] && echo "变量为空"     # -n 非空

# 数值比较（[ ] 内）
[ "$a" -gt "$b" ] && echo "a > b"    # -eq -ne -gt -ge -lt -le

# 逻辑组合
[ -f "$file" ] && [ -s "$file" ] && echo "非空文件"
[ -f "$file" ] || echo "文件不存在"
\`\`\`

### 3.3 case 分支

\`\`\`bash
#!/bin/bash
case "$1" in
  start)
    echo "启动服务"
    ;;
  stop|kill)
    echo "停止服务"
    ;;
  *)
    echo "用法: $0 {start|stop}"
    exit 1
    ;;
esac
\`\`\`

## 4. 循环

### 4.1 for 循环

\`\`\`bash
# 遍历列表
for fruit in apple banana cherry; do
    echo "处理 $fruit"
done

# 遍历范围
for i in {1..5}; do
    echo "第 $i 次"
done

# 遍历文件（最常用）
for file in *.fasta; do
    echo "处理 $file"
done

# C 风格
for ((i = 0; i < 10; i++)); do
    echo $i
done
\`\`\`

### 4.2 while 循环

\`\`\`bash
# 逐行读取文件
while IFS= read -r line; do
    echo "读取: $line"
done < input.txt

# 计数循环
count=0
while [ $count -lt 5 ]; do
    echo $count
    count=$((count + 1))
done
\`\`\`

### 4.3 break / continue

\`\`\`bash
for i in {1..10}; do
    [ $i -eq 3 ] && continue   # 跳过 3
    [ $i -eq 7 ] && break      # 到 7 终止
    echo $i
done
\`\`\`

## 5. 函数

\`\`\`bash
#!/bin/bash

# 定义函数
greet() {
    local name=$1           # local：局部变量
    echo "Hello, $name!"
    return 0                # 返回值（退出码）
}

# 调用
greet "zorrooz"

# 函数返回数值（注意：不是 return，是 echo 捕获）
add() {
    echo $(( $1 + $2 ))
}
result=$(add 3 4)
echo "3 + 4 = $result"
\`\`\`

## 6. 数组与字符串

\`\`\`bash
# 数组
genes=("TP53" "BRCA1" "EGFR")
echo "\${genes[0]}"          # TP53
echo "\${genes[@]}"          # 所有元素
echo "\${#genes[@]}"         # 长度
genes+=("MYC")              # 追加

# 字符串处理
seq="ATGCCGTAA"
echo "\${seq:0:3}"           # 截取前 3 个字符：ATG
echo "\${seq//T/U}"          # 全部替换 T→U
echo "\${seq/AT/XX}"         # 替换第一个 AT
echo "\${seq^^}"             # 转大写（bash 4+）
\`\`\`

## 7. 实战脚本

### 7.1 批量重命名

\`\`\`bash
#!/bin/bash
# 把所有 .txt 后缀改为 .log
set -euo pipefail

for file in *.txt; do
    [ -e "$file" ] || continue          # 无匹配时跳过
    mv "$file" "\${file%.txt}.log"
    echo "重命名: $file -> \${file%.txt}.log"
done
\`\`\`

### 7.2 生信批处理：批量跑 FASTQ 质控

\`\`\`bash
#!/bin/bash
# 用法: ./qc_pipeline.sh data_dir output_dir
set -euo pipefail

DATA_DIR=\${1:?"请提供数据目录"}
OUT_DIR=\${2:?"请提供输出目录"}
mkdir -p "$OUT_DIR"

for fq in "$DATA_DIR"/*_R1.fastq.gz; do
    base=$(basename "$fq" _R1.fastq.gz)
    echo "=== 处理样本: $base ==="

    fastqc -o "$OUT_DIR" "$fq"
    fastqc -o "$OUT_DIR" "\${DATA_DIR}/\${base}_R2.fastq.gz"
done

echo "全部完成，结果位于 $OUT_DIR"
\`\`\`

### 7.3 日志监控与告警

\`\`\`bash
#!/bin/bash
# 监控日志中 ERROR 数量，超过阈值打印告警
LOG_FILE="app.log"
THRESHOLD=\${1:-10}

count=$(grep -c "ERROR" "$LOG_FILE" || true)
echo "当前 ERROR 数量: $count"

if [ "$count" -gt "$THRESHOLD" ]; then
    echo "⚠️ 错误数量超过阈值 $THRESHOLD！"
    exit 1
fi
\`\`\`

## 8. 常见坑位提醒

\`\`\`bash
# 1. 变量赋值等号两边不能有空格
name = "x"    # 错误！会被当成命令
name="x"      # 正确

# 2. 引号：路径含空格必须加引号
file="my data.txt"
cat "$file"   # 正确
cat $file     # 错误：拆成两个参数

# 3. [ ] 内部两端必须有空格
[ "$a" = "$b" ]   # 正确
["$a"="$b"]       # 错误

# 4. 管道中 set -e 的坑：grep 无匹配返回非 0 会中断
grep "x" file || true    # 显式容忍
\`\`\`

## 9. 小结

- 脚本头：\`#!/bin/bash\` + \`set -euo pipefail\`
- 变量：\`$name\`、\`\${#name}\`、\`$(cmd)\`、\`$((...))\`
- 判断：\`[ 条件 ]\`、\`(( 算术 ))\`、\`case\`
- 循环：\`for file in *.txt\`、\`while read line\`
- 函数：\`local\` 变量、echo 捕获返回值

至此编程方向教程完成：Python ×3、R ×3、Linux ×1、Bash ×1，从零到实战。
`;export{n as default};
