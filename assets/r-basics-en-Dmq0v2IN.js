const n=`---
title: "R Language Primer: Data Structures, Vectorization, and Functions"
date: "2026-08-04"
author: "zorrooz"
tags: ["R Language", "Getting Started", "Tutorial"]
draft: false
description: "Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions"
---

# R Language Primer: Data Structures, Vectorization, and Functions

R is the go-to language for statistical computing and data visualization. This article explains the core of R: five data structures, vectorized operations, control flow, and functions, laying the groundwork for subsequent tidyverse and ggplot2.

## 1. Environment Setup

### 1.1 Installation and IDE

- Install [R](https://cran.r-project.org/)
- Recommended IDE: [RStudio](https://posit.co/download/rstudio-desktop/) (free), or VS Code + R extension

### 1.2 Basic Operations

\`\`\`r
# Assignment: <- is the classic R assignment operator (= also works)
x <- 10
x

# View help
?mean
help("lm")
\`\`\`

## 2. Basic Data Structures

R's core data structures can be classified by dimension and homogeneity:

| Structure | Dimension | Element Type |
|-----------|-----------|--------------|
| vector | 1 | Homogeneous |
| factor | 1 | Categorical |
| matrix | 2 | Homogeneous |
| data.frame | 2 | Heterogeneous allowed |
| list | N | Heterogeneous allowed |

### 2.1 Vectors

\`\`\`r
# Create vectors
a <- c(1, 2, 3, 4, 5)            # combine
b <- 1:10                        # sequence
c <- seq(0, 1, by = 0.2)         # step sequence
d <- rep("A", 3)                 # repeat
e <- sample(1:100, 10)           # random sample

# Types
typeof(a)      # "double"
is.numeric(a)  # TRUE

# Indexing (R starts at 1!)
a[1]           # 1
a[c(1, 3)]     # elements 1 and 3
a[-2]          # remove element 2
a[a > 3]       # conditional filtering
names(a) <- c("x1", "x2", "x3", "x4", "x5")
a["x1"]        # access by name
\`\`\`

> **R indexes from 1**, which is the biggest habit difference from Python.

### 2.2 Factors

\`\`\`r
treatment <- factor(c("control", "drug", "control", "drug"),
                    levels = c("control", "drug"))
levels(treatment)          # "control" "drug"
table(treatment)           # frequency counts
\`\`\`

### 2.3 Matrices

\`\`\`r
m <- matrix(1:12, nrow = 3, ncol = 4)
#      [,1] [,2] [,3] [,4]
# [1,]    1    4    7   10
# [2,]    2    5    8   11
# [3,]    3    6    9   12

m[2, 3]        # row 2, column 3: 8
m[1, ]         # row 1
m[, 2]         # column 2
t(m)           # transpose
m %*% t(m)     # matrix multiplication
\`\`\`

### 2.4 Data Frames

Data frames are the core structure for tabular data:

\`\`\`r
df <- data.frame(
  gene = c("TP53", "BRCA1", "EGFR"),
  expression = c(12.3, 8.9, 45.1),
  chromosome = c(17, 17, 7)
)

# View structure
str(df)
summary(df)

# Extract columns
df$expression          # $ operator
df[["expression"]]
df[, "expression"]

# Extract rows
df[1, ]                # first row
df[df$expression > 10, ]

# Add column
df$log2fc <- log2(df$expression)
\`\`\`

### 2.5 Lists

Lists can hold elements of any type:

\`\`\`r
result <- list(
  name = "Analysis Result",
  pvalue = 0.003,
  coef = c(1.2, -0.5, 3.1),
  model = lm(expression ~ chromosome, data = df)
)

result$pvalue
result[[3]]            # numeric vector
\`\`\`

## 3. Vectorized Operations

**R's core idea: operate on entire vectors, not element-by-element loops.**

\`\`\`r
x <- c(1, 2, 3, 4, 5)

x * 2              # element-wise multiplication
x + 10
sqrt(x)
log2(x)

# vector with vector
y <- c(5, 4, 3, 2, 1)
x + y              # add corresponding positions

# comparison operations return logical vectors
x > 3              # FALSE FALSE FALSE  TRUE  TRUE

# statistical functions
sum(x); mean(x); median(x); sd(x); range(x)
\`\`\`

## 4. Control Flow

### 4.1 Conditionals

\`\`\`r
score <- 85

if (score >= 90) {
  grade <- "A"
} else if (score >= 80) {
  grade <- "B"
} else {
  grade <- "C"
}
print(grade)

# vectorized conditional: ifelse
values <- c(60, 90, 45, 75)
result <- ifelse(values >= 60, "pass", "fail")
print(result)   # "pass" "pass" "fail" "pass"
\`\`\`

### 4.2 Loops

\`\`\`r
# for loop
for (i in 1:5) {
  print(i ^ 2)
}

# iterate over vector elements
genes <- c("TP53", "BRCA1", "EGFR")
for (g in genes) {
  print(paste("Analyzing", g))
}

# while loop
n <- 0
while (n < 5) {
  n <- n + 1
}
\`\`\`

### 4.3 Avoid Loops When Possible: The apply Family

R commonly uses the \`apply\` family instead of explicit loops:

\`\`\`r
m <- matrix(1:12, nrow = 3)

apply(m, 1, mean)      # row-wise means
apply(m, 2, sum)       # column-wise sums

lapply(list(1:3, 4:6), mean)   # apply function to each element of a list, returns a list
sapply(list(1:3, 4:6), mean)   # simplified version, returns a vector
\`\`\`

## 5. Functions

### 5.1 Defining Functions

\`\`\`r
gc_content <- function(seq) {
  seq <- toupper(seq)
  gc <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
  gc / nchar(seq) * 100
}

gc_content("ATGCCGA")
# [1] 57.14286
\`\`\`

### 5.2 Default Arguments and Return Values

\`\`\`r
normalize <- function(x, method = "zscore") {
  if (method == "zscore") {
    (x - mean(x)) / sd(x)
  } else if (method == "minmax") {
    (x - min(x)) / (max(x) - min(x))
  } else {
    stop("Unknown method: ", method)
  }
}

normalize(c(1, 2, 3, 4, 5))
normalize(c(1, 2, 3, 4, 5), method = "minmax")
\`\`\`

### 5.3 Anonymous Functions and Functional Programming

\`\`\`r
# anonymous functions
sapply(1:5, function(i) i ^ 2)

# purrr-style map (requires tidyverse, detailed in the next article)
# purrr::map_dbl(1:5, ~ .x ^ 2)
\`\`\`

## 6. Practical Tips

\`\`\`r
# pipe (native pipe |>/ in R 4.1+ or magrittr %>%)
x <- c(1, 2, 3, 4, 5)
x |> mean() |> round(2)

# handle missing values
x <- c(1, NA, 3, NA)
is.na(x)
sum(is.na(x))
x[!is.na(x)]
na.omit(x)

# merging data frames
df1 <- data.frame(id = 1:3, value = c("a", "b", "c"))
df2 <- data.frame(id = 2:4, score = c(90, 85, 88))
merged <- merge(df1, df2, by = "id")   # similar to SQL join
\`\`\`

## 7. Summary

- Five core structures: vectors / factors / matrices / data frames / lists
- **R indexing starts at 1**; \`$\` for columns; \`[ ]\` for rows
- Vectorized operations + \`ifelse\` + the \`apply\` family replace explicit loops
- Functions: default arguments, \`stop()\` for errors, anonymous functions

The next article will introduce tidyverse: using \`dplyr\` for elegant data manipulation.`;export{n as default};
