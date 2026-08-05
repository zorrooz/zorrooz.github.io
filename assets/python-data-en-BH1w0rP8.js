const n=`---
title: "Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas"
date: "2026-08-04"
author: "zorrooz"
tags: ["Python","NumPy","Pandas","Data Processing","Tutorial"]
draft: false
description: "A complete practical tutorial on file I/O, regular expressions, and the scientific computing duo (NumPy/Pandas), covering common bioinformatics data scenarios"
---

# Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas

Data processing is one of the most core application scenarios of Python. This article explains file I/O, regular expressions, and the two pillars of scientific computing: NumPy and Pandas, with examples tailored to common bioinformatics tasks.

## 1. File I/O

### 1.1 Text Files

It is recommended to use the \`with\` statement, which automatically manages file closing:

\`\`\`python
# Writing
with open("output.txt", "w", encoding="utf-8") as f:
    f.write("first line\\n")
    f.write("second line\\n")

# Read all
with open("output.txt", "r", encoding="utf-8") as f:
    content = f.read()
print(content)

# Read line by line (recommended for large files)
with open("output.txt", "r", encoding="utf-8") as f:
    for line in f:
        print(line.strip())   # strip removes the newline
\`\`\`

Mode description: \`r\` read, \`w\` write (overwrite), \`a\` append, \`rb\`/\`wb\` binary.

### 1.2 Parsing FASTA Files

\`\`\`python
def read_fasta(path):
    """Read a FASTA file and return a dictionary {id: sequence}."""
    sequences = {}
    current_id = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:].split()[0]   # take the first field as id
                sequences[current_id] = ""
            elif current_id is not None:
                sequences[current_id] += line
    return sequences

seqs = read_fasta("example.fasta")
for seq_id, seq in seqs.items():
    print(seq_id, len(seq))
\`\`\`

### 1.3 CSV Files

Standard library \`csv\` module:

\`\`\`python
import csv

# Reading
with open("data.csv", "r", newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)          # read as dict based on header
    rows = list(reader)
    print(rows[0]["gene"], rows[0]["expression"])

# Writing
with open("out.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["gene", "value"])
    writer.writerows([["TP53", 12.3], ["BRCA1", 8.9]])
\`\`\`

> For big data scenarios, Pandas' \`read_csv\` is recommended (see Section 4).

## 2. Regular Expressions

Regular expressions are used for text pattern matching. Common functions in the \`re\` module:

\`\`\`python
import re

text = "The protein P12345 has 500 amino acids, code: A1B2C3"

# match: match from the beginning
m = re.match(r"The", text)
print(m.group())          # The

# search: search for the first match
m = re.search(r"\\d+", text)
print(m.group())          # 12345 (first digit string)

# findall: return all matches
print(re.findall(r"[A-Z]\\d", text))   # ['P1', 'A1', 'B2', 'C3']

# sub: replace
print(re.sub(r"\\d+", "#", text))
\`\`\`

### 2.1 Common Metacharacters

| Pattern | Meaning | Example |
|------|------|------|
| \`.\` | Any character (except newline) | \`a.c\` matches \`abc\` |
| \`*\` | Previous character 0 or more times | \`ab*c\` matches \`ac\`, \`abc\` |
| \`+\` | Previous character 1 or more times | \`ab+c\` matches \`abc\` |
| \`?\` | Previous character 0 or 1 time | \`colou?r\` matches color/colour |
| \`^\` / \`$\` | Start of line / End of line | \`^ATG\` begins with ATG |
| \`[abc]\` | Character set | \`[ATGC]\` any base |
| \`\\d\` / \`\\w\` / \`\\s\` | Digit / word character / whitespace | |
| \`{n,m}\` | Repeat n to m times | \`\\d{3}\` three digits |
| \`(…)\` | Group capture | \`(ATG){2}\` |

### 2.2 Practical Information Extraction

Extract coordinates and names from gene annotation text:

\`\`\`python
line = 'gene=TP53;gene_id=ENSG00000141510;chromosome=17;start=7668402;end=7687550'

m = re.search(r"gene=(\\w+);.*?start=(\\d+);end=(\\d+)", line)
if m:
    gene, start, end = m.groups()
    print(f"{gene}: {start}-{end} ({int(end) - int(start)} bp)")
# TP53: 7668402-7687550 (19148 bp)
\`\`\`

### 2.3 Compiling Regular Expressions for Performance

\`\`\`python
motif = re.compile(r"^ATG.*(TAA|TAG|TGA)$")   # rough match for a complete ORF
for seq in ["ATGCCCTAA", "ATGGGG"]:
    print(seq, bool(motif.match(seq)))
\`\`\`

## 3. NumPy: The Foundation of Numerical Computing

### 3.1 Array Creation

\`\`\`python
import numpy as np

a = np.array([1, 2, 3, 4])
b = np.zeros((2, 3))          # 2x3 zero matrix
c = np.ones(5)
d = np.arange(0, 1, 0.2)      # [0.  0.2 0.4 0.6 0.8]
e = np.linspace(0, 1, 5)      # 5 equally spaced points from 0 to 1
f = np.random.rand(3, 3)      # uniformly distributed random numbers
\`\`\`

### 3.2 Vectorized Operations (Core Advantage)

\`\`\`python
x = np.array([1.0, 2.0, 3.0])
y = np.array([4.0, 5.0, 6.0])

print(x + y)          # [5. 7. 9.], no loop needed
print(x * y)          # [4. 10. 18.]
print(np.sqrt(x))     # element-wise square root
print(x.mean(), x.std())   # 2.0 0.816...

# Broadcasting: shapes automatically expand
matrix = np.array([[1, 2], [3, 4]])
print(matrix * 10)            # multiply each element by 10
print(matrix + np.array([100, 200]))  # row broadcast
\`\`\`

### 3.3 Indexing and Slicing

\`\`\`python
arr = np.arange(10)
print(arr[2:8:2])        # [2 4 6], start:end:step
print(arr[arr > 5])      # boolean mask [6 7 8 9]

m = np.random.rand(4, 4)
print(m[0, :])           # first row
print(m[:, 1])           # second column
print(m[1:3, 1:3])       # submatrix
\`\`\`

### 3.4 Bioinformatics Example: Normalizing an Expression Matrix

\`\`\`python
# Simulate a gene expression matrix: 4 genes × 5 samples
expr = np.array([
    [12.3, 15.1, 9.8, 20.4, 11.2],
    [3.2, 4.1, 2.9, 5.0, 3.6],
    [80.5, 92.3, 75.1, 100.2, 88.7],
    [45.6, 41.2, 48.9, 39.8, 50.1],
])

# Perform row-wise (per gene) Z-score normalization: (x - mean) / std
means = expr.mean(axis=1, keepdims=True)
stds = expr.std(axis=1, keepdims=True)
z = (expr - means) / stds

print(np.round(z, 2))
\`\`\`

## 4. Pandas: Tabular Data Processing

### 4.1 Series and DataFrame

\`\`\`python
import pandas as pd

# Series: one-dimensional labeled array
s = pd.Series([1, 2, 3], index=["a", "b", "c"])

# DataFrame: two-dimensional table
df = pd.DataFrame({
    "gene": ["TP53", "BRCA1", "EGFR", "MYC"],
    "expression": [12.3, 8.9, 45.1, 33.2],
    "chromosome": [17, 17, 7, 8],
})
print(df)
\`\`\`

### 4.2 Reading and Writing

\`\`\`python
df = pd.read_csv("expression.csv")
df = pd.read_csv("data.tsv", sep="\\t")
df = pd.read_excel("data.xlsx", sheet_name="Sheet1")
df = pd.read_json("data.json")

df.to_csv("out.csv", index=False)
\`\`\`

### 4.3 Inspection and Filtering

\`\`\`python
print(df.head())          # first 5 rows
print(df.info())          # column types and missing values
print(df.describe())      # statistical summary of numeric columns

# Column selection
print(df["gene"])
print(df[["gene", "expression"]])

# Row filtering (boolean indexing)
high = df[df["expression"] > 10]
print(high)

# Multiple conditions
sel = df[(df["expression"] > 10) & (df["chromosome"] == 17)]
print(sel)

# isin filtering
sel = df[df["gene"].isin(["TP53", "MYC"])]
\`\`\`

### 4.4 Grouping and Aggregation

\`\`\`python
# Group by chromosome and compute the mean, count, and max of expression
grouped = df.groupby("chromosome")["expression"].agg(["mean", "count", "max"])
print(grouped)

# Group-wise transformation (e.g., standardization within groups)
df["zscore"] = df.groupby("chromosome")["expression"].transform(
    lambda x: (x - x.mean()) / x.std()
)
\`\`\`

### 4.5 Handling Missing Values

\`\`\`python
df = pd.DataFrame({"a": [1, None, 3], "b": [4, 5, None]})

print(df.isna().sum())        # number of missing values per column
print(df.dropna())            # drop rows with missing values
print(df.fillna(0))           # fill with 0
print(df.fillna(df.mean()))   # fill with the mean
\`\`\`

### 4.6 Practical Example: Differential Expression Gene Screening

\`\`\`python
# Simulate: genes + treatment/control expression + p-value
results = pd.DataFrame({
    "gene": [f"GENE{i}" for i in range(100)],
    "log2fc": np.random.normal(0, 1.5, 100),
    "pvalue": np.random.rand(100),
})

# Filter significantly differentially expressed genes: |log2FC| > 1 and p < 0.05
sig = results[(results["log2fc"].abs() > 1) & (results["pvalue"] < 0.05)]

# Sort by log2FC
sig = sig.sort_values("log2fc", ascending=False)
print(f"Number of significant genes: {len(sig)}")
print(sig.head())
\`\`\`

## 5. Summary

- Files: \`with open\` + line-by-line reading for large files; use the \`csv\` module or Pandas for CSV
- Regex: \`re.search\` / \`findall\` / \`sub\`; write small patterns first, then combine
- NumPy: vectorized operations + broadcasting; avoid Python loops
- Pandas: \`read_csv\` → filtering/grouping/missing-value handling → export; the standard workflow for tabular data

With this, the Python tutorial trilogy is complete: Beginner → Intermediate → Data Processing.`;export{n as default};
