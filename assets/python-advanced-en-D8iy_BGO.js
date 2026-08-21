const n=`---
title: "Python Advanced: Functions, Classes, and Modules"
date: "2026-08-04"
author: "zorrooz"
tags: ["Python", "Advanced", "Tutorial"]
draft: false
description: "Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages"
---

# Python Advanced: Functions, Classes, and Modules

After mastering the basics, this article dives into Python's functional and object-oriented programming, helping you write clean, reusable code.

## 1. Functions

### 1.1 Definition and Call

\`\`\`python
def greet(name, greeting="Hello"):
    """Return a greeting. greeting is a parameter with a default value."""
    return f"{greeting}, {name}!"

print(greet("zorrooz"))            # Hello, zorrooz!
print(greet("zorrooz", "Hi"))      # Hi, zorrooz!
\`\`\`

### 1.2 Parameter Passing

\`\`\`python
def show(a, b, *args, **kwargs):
    print(f"a={a}, b={b}")
    print(f"Positional args: {args}")      # tuple
    print(f"Keyword args: {kwargs}")  # dict

show(1, 2, 3, 4, x=5, y=6)
# a=1, b=2
# Positional args: (3, 4)
# Keyword args: {'x': 5, 'y': 6}
\`\`\`

- \`*args\`: collects extra positional arguments into a tuple
- \`**kwargs\`: collects extra keyword arguments into a dictionary

### 1.3 Unpacking Arguments

\`\`\`python
def add(a, b, c):
    return a + b + c

nums = [1, 2, 3]
print(add(*nums))            # 6, list unpacking

data = {"a": 1, "b": 2, "c": 3}
print(add(**data))           # 6, dict unpacking
\`\`\`

### 1.4 Lambda Anonymous Functions

\`\`\`python
square = lambda x: x ** 2
print(square(5))             # 25

# Combined with sorted / map / filter
words = ["banana", "apple", "cherry"]
print(sorted(words, key=lambda w: len(w)))
# ['apple', 'banana', 'cherry']

nums = [1, 2, 3, 4, 5]
print(list(map(lambda x: x * 2, nums)))       # [2, 4, 6, 8, 10]
print(list(filter(lambda x: x % 2 == 0, nums)))  # [2, 4]
\`\`\`

### 1.5 Closures

Define a function inside another function and reference outer variables:

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

\`nonlocal\` is used to modify variables in the outer function.

## 2. Decorators

A decorator is a higher-order function that "takes a function and returns a function," used to enhance its behavior without modifying the original function:

\`\`\`python
import time

def timer(func):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        print(f"{func.__name__} took {elapsed:.4f}s")
        return result
    return wrapper

@timer
def slow_sum(n):
    return sum(range(n))

print(slow_sum(10_000_000))
# slow_sum took 0.35xx s
\`\`\`

Decorator with parameters:

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

## 3. Classes and Object-Oriented Programming

### 3.1 Basic Definition

\`\`\`python
class Sequence:
    """Represents a biological sequence."""

    # Class attribute: shared by all instances
    alphabet = "ACGT"

    def __init__(self, seq_id, seq):
        """Constructor: initialize instance attributes."""
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

### 3.2 Inheritance

\`\`\`python
class Protein(Sequence):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(self, seq_id, seq):
        super().__init__(seq_id, seq)   # Call the parent class constructor

    def molecular_weight(self):
        # Simplified calculation: each amino acid is approximately 110 Da
        return self.length() * 110

p = Protein("prot1", "MKWVTFISLL")
print(p.molecular_weight())   # 1210
\`\`\`

### 3.3 Magic Methods

Common magic methods allow objects to support built-in operations:

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

### 3.4 Properties (property)

Use @property to turn methods into attribute access, and add validation:

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
            raise ValueError("Name cannot be empty")
        self._name = value

p = Person("zorrooz")
p.name = "bio"
print(p.name)
\`\`\`

## 4. Exception Handling

\`\`\`python
try:
    num = int(input("Enter an integer: "))
    result = 100 / num
except ValueError:
    print("Not an integer")
except ZeroDivisionError:
    print("Cannot divide by zero")
else:
    print(f"Result: {result}")
finally:
    print("This always executes regardless of errors")
\`\`\`

Custom exceptions:

\`\`\`python
class InvalidSequenceError(Exception):
    pass

def validate(seq):
    if not set(seq) <= set("ACGT"):
        raise InvalidSequenceError(f"Contains invalid characters: {seq}")

try:
    validate("ATGXYZ")
except InvalidSequenceError as e:
    print(f"Validation failed: {e}")
\`\`\`

## 5. Modules and Packages

### 5.1 Module Imports

\`\`\`python
import math                       # Entire module
from math import sqrt, pi         # Import specific names
import numpy as np                # Alias
from collections import Counter   # Common: counting

# Counter example
from collections import Counter
cnt = Counter("ATGCCGA")
print(cnt)          # Counter({'C': 2, 'G': 2, 'A': 2, 'T': 1})
print(cnt.most_common(2))
\`\`\`

### 5.2 Package Structure

\`\`\`
myproject/
├── __init__.py        # Marks it as a package (optional in Python 3.3+)
├── utils/
│   ├── __init__.py
│   └── fasta.py       # Defines read_fasta()
└── main.py
\`\`\`

\`\`\`python
# main.py
from utils.fasta import read_fasta   # Relative package import
\`\`\`

### 5.3 \`if __name__ == "__main__"\`

Allows the script to be both imported and run directly:

\`\`\`python
def main():
    print("Running main logic")

if __name__ == "__main__":
    main()
\`\`\`

## 6. Type Hints (typing)

Type hints improve readability and work with IDE static checks:

\`\`\`python
from typing import List, Dict, Optional

def count_bases(seq: str) -> Dict[str, int]:
    """Return the occurrence count of each base."""
    return {b: seq.count(b) for b in "ACGT"}

def find_motif(seq: str, motif: str) -> Optional[int]:
    idx = seq.find(motif)
    return idx if idx >= 0 else None

print(count_bases("ATGCCGA"))
\`\`\`

## 7. Summary

- \`*args\` / \`**kwargs\`, lambda, closures, and decorators are core to functional programming
- Classes: \`__init__\`, inheritance, magic methods, \`@property\`
- Exceptions: \`try/except/else/finally\`, custom exceptions inherit \`Exception\`
- Modularity: package directories + \`if __name__ == "__main__"\` guard
- Type hints make code more maintainable

The next article will cover Python data processing in practice: file I/O, regular expressions, and NumPy/Pandas.`;export{n as default};
