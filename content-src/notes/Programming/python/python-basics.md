---
title: "Python 编程入门：环境、语法与数据类型"
date: "2026-08-04"
author: "zorrooz"
tags: ["Python", "入门", "教程"]
draft: false
description: "从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码"
---

# Python 编程入门：环境、语法与数据类型

Python 是一门语法简洁、生态丰富的解释型语言，广泛应用于科学计算、数据处理与脚本自动化。本文从零开始，覆盖环境搭建、基础语法、内置数据类型与流程控制，所有示例均可直接运行。

## 1. 环境搭建

### 1.1 安装 Python

访问 [python.org](https://www.python.org/) 下载对应平台的安装包。安装时务必勾选 **Add Python to PATH**。

验证安装：

```bash
python --version
# Python 3.12.x
```

> 在 Linux / WSL 下通常使用 `python3` 命令，Windows 下为 `python`。

### 1.2 虚拟环境（强烈推荐）

虚拟环境为每个项目隔离依赖，避免污染全局环境：

```bash
# 创建虚拟环境
python -m venv .venv

# 激活（Windows）
.venv\Scripts\activate
# 激活（Linux / macOS / WSL）
source .venv/bin/activate
```

激活后提示符前会出现 `(.venv)`。用 [`pip install`](https://pip.pypa.io/en/stable/) 安装依赖：

```bash
pip install numpy pandas
```

将依赖清单导出到文件，方便复现：

```bash
pip freeze > requirements.txt
```

### 1.3 交互式解释器与脚本

```bash
python            # 进入 REPL 交互模式
python script.py  # 运行脚本文件
```

在 REPL 中可直接输入表达式并立即得到结果，适合快速验证想法。

## 2. 第一个程序

```python
print("Hello, Python!")
```

保存为 `hello.py` 后运行：

```bash
python hello.py
# Hello, Python!
```

## 3. 基础语法

### 3.1 注释

```python
# 单行注释

"""
多行注释（实际上是多行字符串，
但常被用作文档说明）
"""
```

### 3.2 变量与赋值

Python 是动态类型语言，变量无需声明类型：

```python
name = "zorrooz"     # 字符串
age = 25             # 整数
height = 1.78        # 浮点数
is_student = True    # 布尔值

# 同时赋值多个变量
a, b, c = 1, 2, 3

# 交换两个变量的值
a, b = b, a
```

变量命名遵循 PEP 8：小写字母 + 下划线（`snake_case`）。

### 3.3 输入输出

```python
name = input("请输入你的名字：")
print("你好，", name)

# f-string 格式化（Python 3.6+，最推荐）
print(f"你好，{name}，今年 {age} 岁")

# format 方法
print("你好，{}，今年 {} 岁".format(name, age))
```

## 4. 内置数据类型

### 4.1 数字

```python
x = 10          # int
y = 3.14        # float
z = 2 + 3j      # complex

# 常用运算
print(7 // 2)   # 整除，得 3
print(7 % 2)    # 取余，得 1
print(2 ** 10)  # 幂运算，得 1024
print(round(3.14159, 2))  # 四舍五入保留 2 位
```

### 4.2 字符串

```python
s = "hello"
t = 'world'

# 拼接与重复
print(s + " " + t)   # hello world
print(s * 3)         # hellohellohello

# 索引与切片（左闭右开）
print(s[0])          # h
print(s[-1])         # o
print(s[1:3])        # el

# 常用方法
print(s.upper())          # HELLO
print(s.replace("l", "L"))  # heLLo
print(len(s))             # 5
print(", ".join(["a", "b", "c"]))  # a, b, c
```

### 4.3 列表（list）

列表是有序、可变的容器：

```python
fruits = ["apple", "banana", "cherry"]

fruits.append("orange")      # 末尾追加
fruits.insert(0, "grape")    # 指定位置插入
fruits.remove("apple")       # 按值删除
popped = fruits.pop()        # 弹出末尾元素

print(fruits[0])             # 索引
print(fruits[-1])            # 最后一个
print(fruits[1:3])           # 切片

# 列表推导式（非常常用）
squares = [x ** 2 for x in range(10)]
print(squares)  # [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
```

### 4.4 元组（tuple）

元组是不可变列表，常用于返回多个值：

```python
point = (3, 4)
x, y = point          # 解包
print(x, y)           # 3 4

def minmax(nums):
    return min(nums), max(nums)

lo, hi = minmax([3, 1, 4, 1, 5])
print(lo, hi)         # 1 5
```

### 4.5 字典（dict）

字典是键值对容器，查询效率高：

```python
person = {"name": "zorrooz", "age": 25}

print(person["name"])            # 取值
print(person.get("city", "未知"))  # 安全的取值，带默认值
person["city"] = "Lanzhou"       # 新增/修改

for key, value in person.items():
    print(f"{key}: {value}")

# 字典推导式
squares = {x: x ** 2 for x in range(5)}
print(squares)  # {0: 0, 1: 1, 2: 4, 3: 9, 4: 16}
```

### 4.6 集合（set）

集合是无序、去重的容器：

```python
a = {1, 2, 3, 3, 3}
print(a)              # {1, 2, 3}，自动去重

b = {2, 3, 4}
print(a & b)          # 交集 {2, 3}
print(a | b)          # 并集 {1, 2, 3, 4}
print(a - b)          # 差集 {1}
```

## 5. 流程控制

### 5.1 条件判断

```python
score = 85

if score >= 90:
    grade = "A"
elif score >= 80:
    grade = "B"
elif score >= 70:
    grade = "C"
else:
    grade = "D"

print(f"等级：{grade}")

# 三元表达式
status = "pass" if score >= 60 else "fail"
```

### 5.2 for 循环

```python
# 遍历可迭代对象
for fruit in ["apple", "banana"]:
    print(fruit)

# range 生成数列
for i in range(5):        # 0 到 4
    print(i)

for i in range(1, 10, 2):  # 1, 3, 5, 7, 9
    print(i)

# enumerate 同时取索引和值
for idx, fruit in enumerate(["a", "b"]):
    print(idx, fruit)
```

### 5.3 while 循环

```python
n = 0
while n < 5:
    print(n)
    n += 1          # 注意：忘记递增会死循环
```

### 5.4 break / continue / else

```python
# break：提前终止循环
for i in range(10):
    if i == 3:
        break
    print(i)        # 输出 0 1 2

# continue：跳过本次迭代
for i in range(5):
    if i == 2:
        continue
    print(i)        # 输出 0 1 3 4

# for-else：循环正常结束（未被 break）时执行
for i in range(3):
    if i == 99:
        break
else:
    print("循环正常结束")
```

## 6. 综合练习

统计一段文本中每个单词的出现次数：

```python
text = "python is fun and python is powerful"

words = text.split()
counter = {}

for word in words:
    counter[word] = counter.get(word, 0) + 1

for word, count in sorted(counter.items(), key=lambda x: -x[1]):
    print(f"{word}: {count}")
```

输出：

```
python: 2
is: 2
fun: 1
and: 1
powerful: 1
```

## 7. 小结

- 用 `venv` 隔离项目依赖，用 `pip freeze` 记录依赖
- 内置类型：`int` / `float` / `str` / `list` / `tuple` / `dict` / `set`
- 列表与字典推导式、f-string 是日常高频语法
- 流程控制：`if/elif/else`、`for`、`while`、`break` / `continue`

下一篇将介绍函数、类与模块，进入真正的工程化编程。
