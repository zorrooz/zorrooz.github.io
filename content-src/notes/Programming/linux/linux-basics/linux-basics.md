---
title: "Linux 命令行基础：文件系统、权限与文本处理"
date: "2026-08-04"
author: "zorrooz"
tags: ["Linux", "命令行", "教程"]
draft: false
description: "Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理"
---

# Linux 命令行基础：文件系统、权限与文本处理

Linux 命令行（Shell）是生物信息学与科学计算的基本功：服务器操作、流程搭建、批量处理都离不开它。本文覆盖最常用的命令，建议配合 WSL 或远程服务器边学边练。

## 1. 终端与 Shell

```bash
echo $SHELL        # 查看当前 shell，通常为 /bin/bash
whoami             # 当前用户
pwd                # 当前目录（print working directory）
```

## 2. 文件系统导航

```bash
ls                  # 列出文件
ls -l               # 详细信息（权限、大小、时间）
ls -a               # 包含隐藏文件（以 . 开头）
ls -lh              # 人类可读大小

cd /home/user       # 进入目录
cd ..               # 上级目录
cd ~                # 用户主目录
cd -                # 上一个目录

pwd                 # 显示当前位置
```

**路径规则**：`/` 是根目录；`.` 当前目录；`..` 上级目录；`~` 主目录。

## 3. 文件与目录操作

```bash
# 创建
touch file.txt      # 创建空文件
mkdir data          # 创建目录
mkdir -p a/b/c      # 递归创建多级目录

# 复制 / 移动 / 删除
cp file.txt copy.txt
cp -r data/ data_backup/    # 复制目录
mv file.txt newname.txt     # 重命名/移动
rm file.txt                 # 删除文件
rm -r data/                 # 删除目录（危险！）
rm -rf data/                # 强制递归删除（慎用）

# 查看内容
cat file.txt        # 输出全部
less file.txt       # 分页浏览（q 退出，/搜索）
head -n 5 file.txt  # 前 5 行
tail -n 5 file.txt  # 后 5 行
tail -f log.txt     # 实时跟踪（查看日志常用）
```

> **rm -rf 是危险命令**：执行前确认路径，或先 `ls` 验证。

## 4. 通配符与重定向

```bash
# 通配符
ls *.fasta          # 所有 .fasta 文件
ls data_??.txt      # ? 匹配单个字符
ls [abc]*           # 以 a/b/c 开头

# 重定向
ls > list.txt       # 覆盖写入
ls >> list.txt      # 追加写入
ls 2> error.log     # 错误输出重定向
ls 2>/dev/null      # 丢弃错误输出

# 管道：前一个命令的输出作为后一个的输入
ls -l | wc -l                   # 统计文件数
history | grep python           # 在历史命令中搜索
```

## 5. 权限管理

```bash
ls -l
# -rw-r--r-- 1 user group 1234 Aug 4 10:00 file.txt
# ^权限       ^属主 ^属组

chmod 755 script.sh    # rwxr-xr-x：属主读写执行，其他只读执行
chmod +x script.sh     # 添加执行权限（运行脚本必需）
chmod -w file.txt      # 去除写权限

chown user:group file  # 修改属主属组（需要 root）
```

权限数字：`r=4`、`w=2`、`x=1`，三位分别表示属主/属组/其他。

## 6. 文本处理三剑客

### 6.1 grep：搜索

```bash
grep "TP53" genes.txt          # 查找包含 TP53 的行
grep -i "tp53" genes.txt       # 忽略大小写
grep -v "comment" file.txt     # 反向匹配（排除）
grep -c "gene" file.txt        # 计数
grep -n "pattern" file.txt     # 显示行号
grep -r "TODO" src/            # 递归搜索目录

# 管道组合
ps aux | grep python           # 找 python 进程
cat reads.fastq | grep -c "^@" # 统计 FASTQ 中序列条数
```

### 6.2 sed：流编辑

```bash
sed -n '5,10p' file.txt        # 打印第 5-10 行
sed 's/old/new/' file.txt      # 替换每行第一个匹配
sed 's/old/new/g' file.txt     # 全局替换
sed -i 's/old/new/g' file.txt  # 直接修改文件（-i 原地）
sed '/^#/d' config.conf        # 删除注释行
```

### 6.3 awk：按列处理

```bash
awk '{print $1}' file.txt      # 打印第一列
awk -F',' '{print $1, $3}' data.csv   # 按逗号分隔
awk '$3 > 50 {print $1}' data.txt     # 条件过滤
awk '{sum += $2} END {print sum}' data.txt   # 第二列求和
awk 'NR > 1 {print}' file.txt  # 跳过表头
```

> awk 字段分隔符默认空白；`-F','` 指定逗号。

## 7. 进程管理

```bash
ps aux                # 查看所有进程
ps aux | grep python  # 查找特定进程
top                   # 实时监控（q 退出）
htop                  # 增强版监控（更直观）

kill PID              # 终止进程
kill -9 PID           # 强制终止
pkill python          # 按名称终止

# 后台运行
python script.py &    # 后台运行
nohup python script.py > log.txt 2>&1 &   # 脱离终端运行（服务器常用）
jobs                  # 查看后台任务
fg                    # 调回前台
```

> 远程服务器跑长任务用 `nohup ... &` 或 `tmux` / `screen`，避免终端断开导致任务中断。

## 8. 实用组合

```bash
# 统计文件数量与大小
ls -lh | awk '{print $5, $9}'
du -sh data/          # 目录总大小
df -h                 # 磁盘空间

# 压缩与解压
tar -czf archive.tar.gz dir/    # 打包压缩
tar -xzf archive.tar.gz         # 解压
zip -r archive.zip dir/
unzip archive.zip

# 查找文件
find . -name "*.log"            # 按名查找
find . -size +100M              # 按大小查找
which python                    # 命令所在路径
```

## 9. WSL 与开发环境

Windows 下推荐 WSL（Windows Subsystem for Linux）：

```bash
# 安装（管理员 PowerShell）
wsl --install

# 进入 Ubuntu 子系统
wsl

# 在 WSL 中配置开发环境
sudo apt update
sudo apt install build-essential git python3 python3-pip
```

VS Code 连接 WSL：安装 Remote - WSL 扩展后，`code .` 直接打开。

## 10. 小结

- 导航：`cd` / `ls` / `pwd`；操作：`cp` / `mv` / `rm` / `mkdir`
- 三剑客：`grep` 搜索、`sed` 替换、`awk` 按列处理
- 管道 `|` 与重定向 `>` 是组合命令的核心
- 权限 `chmod 755` / `chmod +x`；进程 `ps` / `kill` / `nohup`

下一篇将介绍 Bash 编程：把命令组织成可复用的脚本。
