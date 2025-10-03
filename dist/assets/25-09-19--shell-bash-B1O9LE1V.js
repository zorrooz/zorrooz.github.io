const r=`---\r
title: "Bash Shell 脚本编程"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Shell", "Bash", "命令行", "脚本编程", "自动化"]\r
draft: false\r
description: "Bash Shell 脚本在生物信息学数据处理中的自动化应用"\r
---\r
\r
# Bash Shell 脚本编程\r
\r
## 基本语法\r
\r
\`\`\`bash\r
#!/bin/bash\r
\r
# 变量定义\r
input_file="data.fastq"\r
output_dir="results"\r
\r
# 条件判断\r
if [ -f "$input_file" ]; then\r
    echo "文件存在"\r
else\r
    echo "文件不存在"\r
fi\r
\`\`\`\r
\r
## 循环处理\r
\r
\`\`\`bash\r
# 遍历文件\r
for file in *.fastq; do\r
    echo "处理文件: $file"\r
done`;export{r as default};
