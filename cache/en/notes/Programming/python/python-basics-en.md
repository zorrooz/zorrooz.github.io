---
title: "Python Programming Introduction: Environment, Syntax, and Data Types"
date: "2026-08-04"
author: "zorrooz"
tags: ["Python", "Getting Started", "Tutorial"]
draft: false
description: "Learn Python from scratch: environment setup, basic syntax, built-in data types, and flow control, with runnable example code"
---

# Python Programming Introduction: Environment, Syntax, and Data Types

Python is an interpreted language with concise syntax and a rich ecosystem, widely used in scientific computing, data processing, and script automation. This article starts from scratch, covering environment setup, basic syntax, built-in data types, and flow control, with all examples ready to run.

## 1. Environment Setup

### 1.1 Installing Python

Visit [python.org](https://www.python.org/) to download the installer for your platform. Be sure to check **Add Python to PATH** during installation.

Verify the installation:

```bash
python --version
# Python 3.12.x
```

> On Linux / WSL, the command is typically `python3`, while on Windows it is `python`.

### 1.2 Virtual Environments (Highly Recommended)

Virtual environments isolate dependencies for each project, avoiding pollution of the global environment:

```bash
# Create a virtual environment
python -m venv .venv

# Activate (Windows)
.venv\Scripts\activate
# Activate (Linux / macOS / WSL)
source .venv/bin/activate
```

After activation, a `(.venv)` prefix appears in the prompt. Use [`pip install`](https://pip.pypa.io/en/stable/) to install dependencies:

```bash
pip install numpy pandas
```

Export the dependency list to a file for easy reproduction:

```bash
pip freeze > requirements.txt
```

### 1.3 Interactive Interpreter and Scripts

```bash
python            # Enter REPL interactive mode
python script.py  # Run a script file
```

In the REPL, you can directly enter expressions and get immediate results, ideal for quickly testing ideas.

## 2. First Program

```python
print("Hello, Python!")
```

Save it as `hello.py` and run:

```bash
python hello.py
# Hello, Python!
```

## 3. Basic Syntax

### 3.1 Comments

```python
# Single-line comment

"""
Multi-line comment (actually a multi-line string,
but often used for documentation)
"""
```

### 3.2 Variables and Assignment

Python is a dynamically typed language; variables do not require type declarations:

```python
name = "zorrooz"     # String
age = 25             # Integer
height = 1.78        # Float
is_student = True    # Boolean

# Assign multiple variables at once
a, b, c = 1, 2, 3

# Swap two variables
a, b = b, a
```

Variable naming follows PEP 8: lowercase letters plus underscores (`snake_case`).

### 3.3 Input and Output

```python
name = input("Please enter your name: ")
print("Hello,", name)

# f-string formatting (Python 3.6+, most recommended)
print(f"Hello, {name}, you are {age} years old")

# format method
print("Hello, {}, you are {} years old".format(name, age))
```

## 4. Built-in Data Types

### 4.1 Numbers

```python
x = 10          # int
y = 3.14        # float
z = 2 + 3j      # complex

# Common operations
print(7 // 2)   # Integer division, gives 3
print(7 % 2)    # Modulo, gives 1
print(2 ** 10)  # Exponentiation, gives 1024
print(round(3.14159, 2))  # Round to 2 decimal places
```

### 4.2 Strings

```python
s = "hello"
t = 'world'

# Concatenation and repetition
print(s + " " + t)   # hello world
print(s * 3)         # hellohellohello

# Indexing and slicing (left-closed, right-open)
print(s[0])          # h
print(s[-1])         # o
print(s[1:3])        # el

# Common methods
print(s.upper())          # HELLO
print(s.replace("l", "L"))  # heLLo
print(len(s))             # 5
print(", ".join(["a", "b", "c"]))  # a, b, c
```

### 4.3 Lists

Lists are ordered, mutable containers:

```python
fruits = ["apple", "banana", "cherry"]

fruits.append("orange")      # Append at the end
fruits.insert(0, "grape")    # Insert at a specific position
fruits.remove("apple")       # Remove by value
popped = fruits.pop()        # Pop the last element

print(fruits[0])             # Index
print(fruits[-1])            # Last element
print(fruits[1:3])           # Slicing

# List comprehension (very common)
squares = [x ** 2 for x in range(10)]
print(squares)  # [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
```

### 4.4 Tuples

Tuples are immutable lists, often used for returning multiple values:

```python
point = (3, 4)
x, y = point          # Unpacking
print(x, y)           # 3 4

def minmax(nums):
    return min(nums), max(nums)

lo, hi = minmax([3, 1, 4, 1, 5])
print(lo, hi)         # 1 5
```

### 4.5 Dictionaries

Dictionaries are key-value containers with high lookup efficiency:

```python
person = {"name": "zorrooz", "age": 25}

print(person["name"])            # Get value
print(person.get("city", "unknown"))  # Safe get with default
person["city"] = "Lanzhou"       # Add/Modify

for key, value in person.items():
    print(f"{key}: {value}")

# Dictionary comprehension
squares = {x: x ** 2 for x in range(5)}
print(squares)  # {0: 0, 1: 1, 2: 4, 3: 9, 4: 16}
```

### 4.6 Sets

Sets are unordered, deduplicated containers:

```python
a = {1, 2, 3, 3, 3}
print(a)              # {1, 2, 3}, automatically deduplicated

b = {2, 3, 4}
print(a & b)          # Intersection {2, 3}
print(a | b)          # Union {1, 2, 3, 4}
print(a - b)          # Difference {1}
```

## 5. Flow Control

### 5.1 Conditional Statements

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

print(f"Grade: {grade}")

# Ternary expression
status = "pass" if score >= 60 else "fail"
```

### 5.2 for Loops

```python
# Iterate over an iterable
for fruit in ["apple", "banana"]:
    print(fruit)

# range generates a sequence
for i in range(5):        # 0 to 4
    print(i)

for i in range(1, 10, 2):  # 1, 3, 5, 7, 9
    print(i)

# enumerate gets both index and value
for idx, fruit in enumerate(["a", "b"]):
    print(idx, fruit)
```

### 5.3 while Loops

```python
n = 0
while n < 5:
    print(n)
    n += 1          # Note: forgetting to increment causes an infinite loop
```

### 5.4 break / continue / else

```python
# break: terminate the loop early
for i in range(10):
    if i == 3:
        break
    print(i)        # Outputs 0 1 2

# continue: skip the current iteration
for i in range(5):
    if i == 2:
        continue
    print(i)        # Outputs 0 1 3 4

# for-else: executes when the loop completes normally (not broken)
for i in range(3):
    if i == 99:
        break
else:
    print("Loop completed normally")
```

## 6. Comprehensive Exercise

Count the occurrences of each word in a text:

```python
text = "python is fun and python is powerful"

words = text.split()
counter = {}

for word in words:
    counter[word] = counter.get(word, 0) + 1

for word, count in sorted(counter.items(), key=lambda x: -x[1]):
    print(f"{word}: {count}")
```

Output:

```
python: 2
is: 2
fun: 1
and: 1
powerful: 1
```

## 7. Summary

- Use `venv` to isolate project dependencies and `pip freeze` to record them
- Built-in types: `int` / `float` / `str` / `list` / `tuple` / `dict` / `set`
- List and dictionary comprehensions, and f-strings are high-frequency daily syntax
- Flow control: `if/elif/else`, `for`, `while`, `break` / `continue`

The next article will introduce functions, classes, and modules, moving into true engineering programming.