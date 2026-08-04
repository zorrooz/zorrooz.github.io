const n=`---
title: "Python 进阶：函数、类与模块"
date: "2026-08-04"
author: "zorrooz"
tags: ["Python", "进阶", "教程"]
draft: false
description: "深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包"
---

# Python 进阶：函数、类与模块

掌握基础语法后，本文深入 Python 的函数式与面向对象编程，帮助你写出结构清晰、可复用的代码。

## 1. 函数

### 1.1 定义与调用

\`\`\`python
def greet(name, greeting="Hello"):
    """返回问候语。greeting 是带默认值的参数。"""
    return f"{greeting}, {name}!"

print(greet("zorrooz"))            # Hello, zorrooz!
print(greet("zorrooz", "Hi"))      # Hi, zorrooz!
\`\`\`

### 1.2 参数传递

\`\`\`python
def show(a, b, *args, **kwargs):
    print(f"a={a}, b={b}")
    print(f"位置参数: {args}")      # 元组
    print(f"关键字参数: {kwargs}")  # 字典

show(1, 2, 3, 4, x=5, y=6)
# a=1, b=2
# 位置参数: (3, 4)
# 关键字参数: {'x': 5, 'y': 6}
\`\`\`

- \`*args\`：收集多余的位置参数为元组
- \`**kwargs\`：收集多余的关键字参数为字典

### 1.3 解包传参

\`\`\`python
def add(a, b, c):
    return a + b + c

nums = [1, 2, 3]
print(add(*nums))            # 6，列表解包

data = {"a": 1, "b": 2, "c": 3}
print(add(**data))           # 6，字典解包
\`\`\`

### 1.4 lambda 匿名函数

\`\`\`python
square = lambda x: x ** 2
print(square(5))             # 25

# 与 sorted / map / filter 配合
words = ["banana", "apple", "cherry"]
print(sorted(words, key=lambda w: len(w)))
# ['apple', 'banana', 'cherry']

nums = [1, 2, 3, 4, 5]
print(list(map(lambda x: x * 2, nums)))       # [2, 4, 6, 8, 10]
print(list(filter(lambda x: x % 2 == 0, nums)))  # [2, 4]
\`\`\`

### 1.5 闭包

函数内部定义函数并引用外部变量：

\`\`\`python
def make_counter():
    count = 0
    def counter():
        nonlocal count
        count += 1
        return count
    return counter

c = make_counter()
print(c())   # 1
print(c())   # 2
\`\`\`

\`nonlocal\` 用于修改外层函数中的变量。

## 2. 装饰器

装饰器是"接收函数、返回函数"的高阶函数，用于在不修改原函数的情况下增强其行为：

\`\`\`python
import time

def timer(func):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        print(f"{func.__name__} 耗时 {elapsed:.4f}s")
        return result
    return wrapper

@timer
def slow_sum(n):
    return sum(range(n))

print(slow_sum(10_000_000))
# slow_sum 耗时 0.35xx s
\`\`\`

带参数的装饰器：

\`\`\`python
def repeat(times):
    def decorator(func):
        def wrapper(*args, **kwargs):
            for _ in range(times):
                func(*args, **kwargs)
        return wrapper
    return decorator

@repeat(3)
def hello():
    print("Hi!")

hello()
# Hi! Hi! Hi!
\`\`\`

## 3. 类与面向对象

### 3.1 基本定义

\`\`\`python
class Sequence:
    """表示一条生物序列。"""

    # 类属性：所有实例共享
    alphabet = "ACGT"

    def __init__(self, seq_id, seq):
        """构造方法：初始化实例属性。"""
        self.seq_id = seq_id
        self.seq = seq.upper()

    def length(self):
        return len(self.seq)

    def gc_content(self):
        gc = self.seq.count("G") + self.seq.count("C")
        return gc / len(self.seq) * 100

    def __repr__(self):
        return f"Sequence({self.seq_id}, {self.length()}bp)"

s = Sequence("seq1", "atgcgta")
print(s.length())        # 7
print(f"GC: {s.gc_content():.1f}%")   # GC: 57.1%
print(s)                 # Sequence(seq1, 7bp)
\`\`\`

### 3.2 继承

\`\`\`python
class Protein(Sequence):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(self, seq_id, seq):
        super().__init__(seq_id, seq)   # 调用父类构造方法

    def molecular_weight(self):
        # 简化计算：每个氨基酸约 110 Da
        return self.length() * 110

p = Protein("prot1", "MKWVTFISLL")
print(p.molecular_weight())   # 1210
\`\`\`

### 3.3 魔术方法

常用魔术方法让对象支持内置操作：

\`\`\`python
class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __add__(self, other):      # +
        return Vector(self.x + other.x, self.y + other.y)

    def __eq__(self, other):       # ==
        return self.x == other.x and self.y == other.y

    def __repr__(self):            # print / repr
        return f"Vector({self.x}, {self.y})"

v1 = Vector(1, 2)
v2 = Vector(3, 4)
print(v1 + v2)          # Vector(4, 6)
print(v1 == Vector(1, 2))  # True
\`\`\`

### 3.4 属性（property）

用 \`@property\` 把方法变成属性访问，可加入校验：

\`\`\`python
class Person:
    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not value.strip():
            raise ValueError("名字不能为空")
        self._name = value

p = Person("zorrooz")
p.name = "bio"
print(p.name)
\`\`\`

## 4. 异常处理

\`\`\`python
try:
    num = int(input("输入一个整数："))
    result = 100 / num
except ValueError:
    print("输入的不是整数")
except ZeroDivisionError:
    print("不能除以零")
else:
    print(f"结果: {result}")
finally:
    print("无论是否出错都会执行")
\`\`\`

自定义异常：

\`\`\`python
class InvalidSequenceError(Exception):
    pass

def validate(seq):
    if not set(seq) <= set("ACGT"):
        raise InvalidSequenceError(f"包含非法字符: {seq}")

try:
    validate("ATGXYZ")
except InvalidSequenceError as e:
    print(f"校验失败: {e}")
\`\`\`

## 5. 模块与包

### 5.1 模块导入

\`\`\`python
import math                       # 整个模块
from math import sqrt, pi         # 导入指定名字
import numpy as np                # 别名
from collections import Counter   # 常用：计数

# Counter 示例
from collections import Counter
cnt = Counter("ATGCCGA")
print(cnt)          # Counter({'C': 2, 'G': 2, 'A': 2, 'T': 1})
print(cnt.most_common(2))
\`\`\`

### 5.2 包结构

\`\`\`
myproject/
├── __init__.py        # 标记为包（Python 3.3+ 可省略）
├── utils/
│   ├── __init__.py
│   └── fasta.py       # 定义 read_fasta()
└── main.py
\`\`\`

\`\`\`python
# main.py
from utils.fasta import read_fasta   # 相对包路径导入
\`\`\`

### 5.3 \`if __name__ == "__main__"\`

让脚本既可被导入也可直接运行：

\`\`\`python
def main():
    print("运行主逻辑")

if __name__ == "__main__":
    main()
\`\`\`

## 6. 类型提示（typing）

类型提示提高可读性，配合 IDE 静态检查：

\`\`\`python
from typing import List, Dict, Optional

def count_bases(seq: str) -> Dict[str, int]:
    """返回每种碱基的出现次数。"""
    return {b: seq.count(b) for b in "ACGT"}

def find_motif(seq: str, motif: str) -> Optional[int]:
    idx = seq.find(motif)
    return idx if idx >= 0 else None

print(count_bases("ATGCCGA"))
\`\`\`

## 7. 小结

- \`*args\` / \`**kwargs\`、lambda、闭包、装饰器是函数式编程核心
- 类：\`__init__\`、继承、魔术方法、\`@property\`
- 异常：\`try/except/else/finally\`，自定义异常继承 \`Exception\`
- 模块化：包目录 + \`if __name__ == "__main__"\` 守卫
- 类型提示让代码更可维护

下一篇将介绍 Python 数据处理实战：文件 IO、正则表达式与 NumPy/Pandas。
`;export{n as default};
