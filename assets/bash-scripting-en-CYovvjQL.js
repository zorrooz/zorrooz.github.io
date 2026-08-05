const n=`---
title: "Bash Programming: Variables, Conditions, Loops, and Practical Scripts"
date: "2026-08-04"
author: "zorrooz"
tags: ["Bash","Shell","Scripting","Tutorial"]
draft: false
description: "Complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays and string processing, with practical batch processing and bioinformatics scripts"
---

# Bash Programming: Variables, Conditions, Loops, and Practical Scripts

Bash scripts organize command-line operations into reusable programs and are the core tool for automated batch processing. This article covers everything from script basics to practical applications, including variables, conditions, loops, functions, and common tasks.

## 1. Script Basics

### 1.1 Script Structure

\`\`\`bash
#!/bin/bash
# The first line shebang: specifies the interpreter

# Give the script execute permission
chmod +x myscript.sh
# Run
./myscript.sh
# Or run directly with bash (no execute permission needed)
bash myscript.sh
\`\`\`

### 1.2 Debugging

\`\`\`bash
bash -x script.sh    # Trace execution, showing each command
# Enable inside script: set -x (tracing) / set -e (exit on error) / set -u (error on undefined variable)
\`\`\`

It is recommended to add the following at the beginning of production scripts:

\`\`\`bash
#!/bin/bash
set -euo pipefail   # Exit on error, error on undefined variable, detect pipeline failures
\`\`\`

## 2. Variables

### 2.1 Definition and Reference

\`\`\`bash
name="zorrooz"        # Note: no spaces around =
echo "Hello, $name"   # Variable expansion inside double quotes
echo 'Hello, $name'   # Literal output inside single quotes
echo "Length: \${#name}" # String length
\`\`\`

### 2.2 Command Substitution and Arithmetic

\`\`\`bash
# Command substitution: $(...)
today=$(date +%F)
echo "Today is $today"
files=$(ls | wc -l)
echo "Total $files files"

# Arithmetic: $((...))
a=10
b=3
echo $((a + b))    # 13
echo $((a % b))    # 1
\`\`\`

### 2.3 Positional Parameters

\`\`\`bash
#!/bin/bash
echo "Script name: $0"
echo "First argument: $1"
echo "Second argument: $2"
echo "All arguments: $@"
echo "Number of arguments: $#"

# Run: ./script.sh file1 file2
\`\`\`

### 2.4 Special Variables

\`\`\`bash
$?    # Exit code of the last command (0 success, non-zero failure)
$$    # Current process PID
$!    # Background process PID
\`\`\`

## 3. Conditional Statements

### 3.1 if / elif / else

\`\`\`bash
#!/bin/bash
score=$1

if (( score >= 90 )); then
    echo "Excellent"
elif (( score >= 60 )); then
    echo "Pass"
else
    echo "Fail"
fi
\`\`\`

### 3.2 test Expressions

\`\`\`bash
# File checks
[ -f file.txt ] && echo "Is a regular file"
[ -d dir ] && echo "Is a directory"
[ -e path ] && echo "Exists"
[ -s file ] && echo "Non-empty file"

# String comparison
[ "$name" = "zorrooz" ] && echo "Match"
[ -z "$var" ] && echo "Variable is empty"     # -n is non-empty

# Numeric comparison (inside [ ])
[ "$a" -gt "$b" ] && echo "a > b"    # -eq -ne -gt -ge -lt -le

# Logical combinations
[ -f "$file" ] && [ -s "$file" ] && echo "Non-empty file"
[ -f "$file" ] || echo "File does not exist"
\`\`\`

### 3.3 case Branches

\`\`\`bash
#!/bin/bash
case "$1" in
  start)
    echo "Starting service"
    ;;
  stop|kill)
    echo "Stopping service"
    ;;
  *)
    echo "Usage: $0 {start|stop}"
    exit 1
    ;;
esac
\`\`\`

## 4. Loops

### 4.1 for Loops

\`\`\`bash
# Iterate over a list
for fruit in apple banana cherry; do
    echo "Processing $fruit"
done

# Iterate over a range
for i in {1..5}; do
    echo "Iteration $i"
done

# Iterate over files (most common)
for file in *.fasta; do
    echo "Processing $file"
done

# C-style
for ((i = 0; i < 10; i++)); do
    echo $i
done
\`\`\`

### 4.2 while Loops

\`\`\`bash
# Read a file line by line
while IFS= read -r line; do
    echo "Read: $line"
done < input.txt

# Counting loop
count=0
while [ $count -lt 5 ]; do
    echo $count
    count=$((count + 1))
done
\`\`\`

### 4.3 break / continue

\`\`\`bash
for i in {1..10}; do
    [ $i -eq 3 ] && continue   # Skip 3
    [ $i -eq 7 ] && break      # Terminate at 7
    echo $i
done
\`\`\`

## 5. Functions

\`\`\`bash
#!/bin/bash

# Define a function
greet() {
    local name=$1           # local: local variable
    echo "Hello, $name!"
    return 0                # Return value (exit code)
}

# Call
greet "zorrooz"

# Function returning a value (note: not return, but echo capture)
add() {
    echo $(( $1 + $2 ))
}
result=$(add 3 4)
echo "3 + 4 = $result"
\`\`\`

## 6. Arrays and Strings

\`\`\`bash
# Arrays
genes=("TP53" "BRCA1" "EGFR")
echo "\${genes[0]}"          # TP53
echo "\${genes[@]}"          # All elements
echo "\${#genes[@]}"         # Length
genes+=("MYC")              # Append

# String processing
seq="ATGCCGTAA"
echo "\${seq:0:3}"           # Extract first 3 characters: ATG
echo "\${seq//T/U}"          # Replace all T→U
echo "\${seq/AT/XX}"         # Replace first AT
echo "\${seq^^}"             # Convert to uppercase (bash 4+)
\`\`\`

## 7. Practical Scripts

### 7.1 Batch Rename

\`\`\`bash
#!/bin/bash
# Change all .txt extensions to .log
set -euo pipefail

for file in *.txt; do
    [ -e "$file" ] || continue          # Skip if no match
    mv "$file" "\${file%.txt}.log"
    echo "Renamed: $file -> \${file%.txt}.log"
done
\`\`\`

### 7.2 Bioinformatics Batch Processing: FASTQ Quality Control

\`\`\`bash
#!/bin/bash
# Usage: ./qc_pipeline.sh data_dir output_dir
set -euo pipefail

DATA_DIR=\${1:?"Please provide a data directory"}
OUT_DIR=\${2:?"Please provide an output directory"}
mkdir -p "$OUT_DIR"

for fq in "$DATA_DIR"/*_R1.fastq.gz; do
    base=$(basename "$fq" _R1.fastq.gz)
    echo "=== Processing sample: $base ==="

    fastqc -o "$OUT_DIR" "$fq"
    fastqc -o "$OUT_DIR" "\${DATA_DIR}/\${base}_R2.fastq.gz"
done

echo "All done, results are in $OUT_DIR"
\`\`\`

### 7.3 Log Monitoring and Alerting

\`\`\`bash
#!/bin/bash
# Monitor the number of ERRORs in the log, print an alert if it exceeds the threshold
LOG_FILE="app.log"
THRESHOLD=\${1:-10}

count=$(grep -c "ERROR" "$LOG_FILE" || true)
echo "Current ERROR count: $count"

if [ "$count" -gt "$THRESHOLD" ]; then
    echo "⚠️ Error count exceeds threshold $THRESHOLD!"
    exit 1
fi
\`\`\`

## 8. Common Pitfall Reminders

\`\`\`bash
# 1. No spaces around = in variable assignment
name = "x"    # Wrong! Will be treated as a command
name="x"      # Correct

# 2. Quoting: paths containing spaces must be quoted
file="my data.txt"
cat "$file"   # Correct
cat $file     # Wrong: split into two arguments

# 3. Spaces are required inside [ ]
[ "$a" = "$b" ]   # Correct
["$a"="$b"]       # Wrong

# 4. The pitfall of set -e in pipelines: grep with no match returns non-zero and will abort
grep "x" file || true    # Explicitly tolerate
\`\`\`

## 9. Summary

- Script header: \`#!/bin/bash\` + \`set -euo pipefail\`
- Variables: \`$name\`, \`\${#name}\`, \`$(cmd)\`, \`$((...))\`
- Conditionals: \`[ condition ]\`, \`(( arithmetic ))\`, \`case\`
- Loops: \`for file in *.txt\`, \`while read line\`
- Functions: \`local\` variables, echo capture for return values

This completes the programming direction tutorials: Python ×3, R ×3, Linux ×1, Bash ×1, from zero to practical application.`;export{n as default};
