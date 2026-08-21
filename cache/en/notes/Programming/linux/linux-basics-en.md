---
title: "Linux Command Line Basics: File System, Permissions, and Text Processing"
date: "2026-08-04"
author: "zorrooz"
tags: ["Linux", "Command Line", "Tutorial"]
draft: false
description: "Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management"
---

# Linux Command Line Basics: File System, Permissions, and Text Processing

The Linux command line (Shell) is a fundamental skill for bioinformatics and scientific computing: server operations, workflow construction, and batch processing all depend on it. This article covers the most commonly used commands, and it is recommended to practice while learning using WSL or a remote server.

## 1. Terminal and Shell

```bash
echo $SHELL        # View the current shell, usually /bin/bash
whoami             # Current user
pwd                # Current directory (print working directory)
```

## 2. File System Navigation

```bash
ls                  # List files
ls -l               # Detailed information (permissions, size, time)
ls -a               # Include hidden files (starting with .)
ls -lh              # Human-readable sizes

cd /home/user       # Enter directory
cd ..               # Parent directory
cd ~                # User home directory
cd -                # Previous directory

pwd                 # Display current location
```

**Path rules**: `/` is the root directory; `.` is the current directory; `..` is the parent directory; `~` is the home directory.

## 3. File and Directory Operations

```bash
# Create
touch file.txt      # Create an empty file
mkdir data          # Create a directory
mkdir -p a/b/c      # Recursively create multi-level directories

# Copy / Move / Delete
cp file.txt copy.txt
cp -r data/ data_backup/    # Copy directory
mv file.txt newname.txt     # Rename/move
rm file.txt                 # Delete file
rm -r data/                 # Delete directory (dangerous!)
rm -rf data/                # Force recursive delete (use with caution)

# View contents
cat file.txt        # Output all contents
less file.txt       # Page through (q to quit, / to search)
head -n 5 file.txt  # First 5 lines
tail -n 5 file.txt  # Last 5 lines
tail -f log.txt     # Follow in real time (commonly used for viewing logs)
```

> **`rm -rf` is a dangerous command**: confirm the path before executing, or first `ls` to verify.

## 4. Wildcards and Redirection

```bash
# Wildcards
ls *.fasta          # All .fasta files
ls data_??.txt      # ? matches a single character
ls [abc]*           # Starts with a/b/c

# Redirection
ls > list.txt       # Overwrite
ls >> list.txt      # Append
ls 2> error.log     # Redirect error output
ls 2>/dev/null      # Discard error output

# Pipe: output of the previous command becomes input of the next
ls -l | wc -l                   # Count number of files
history | grep python           # Search history commands
```

## 5. Permission Management

```bash
ls -l
# -rw-r--r-- 1 user group 1234 Aug 4 10:00 file.txt
# ^permissions       ^owner ^group

chmod 755 script.sh    # rwxr-xr-x: owner read/write/execute, others read/execute only
chmod +x script.sh     # Add execute permission (required to run scripts)
chmod -w file.txt      # Remove write permission

chown user:group file  # Change owner and group (requires root)
```

Permission numbers: `r=4`, `w=2`, `x=1`, and the three digits represent owner/group/others respectively.

## 6. The Three Musketeers of Text Processing

### 6.1 grep: Search

```bash
grep "TP53" genes.txt          # Find lines containing TP53
grep -i "tp53" genes.txt       # Ignore case
grep -v "comment" file.txt     # Reverse match (exclude)
grep -c "gene" file.txt        # Count
grep -n "pattern" file.txt     # Show line numbers
grep -r "TODO" src/            # Recursively search directory

# Pipe combinations
ps aux | grep python           # Find python processes
cat reads.fastq | grep -c "^@" # Count number of sequences in FASTQ
```

### 6.2 sed: Stream Editing

```bash
sed -n '5,10p' file.txt        # Print lines 5-10
sed 's/old/new/' file.txt      # Replace first match in each line
sed 's/old/new/g' file.txt     # Global replace
sed -i 's/old/new/g' file.txt  # Modify file in place (-i in-place)
sed '/^#/d' config.conf        # Delete comment lines
```

### 6.3 awk: Column-Based Processing

```bash
awk '{print $1}' file.txt      # Print first column
awk -F',' '{print $1, $3}' data.csv   # Separate by comma
awk '$3 > 50 {print $1}' data.txt     # Conditional filter
awk '{sum += $2} END {print sum}' data.txt   # Sum second column
awk 'NR > 1 {print}' file.txt  # Skip header line
```

> awk field separator defaults to whitespace; `-F','` specifies comma.

## 7. Process Management

```bash
ps aux                # View all processes
ps aux | grep python  # Find specific processes
top                   # Real-time monitoring (q to quit)
htop                  # Enhanced monitoring (more intuitive)

kill PID              # Terminate process
kill -9 PID           # Force kill
pkill python          # Terminate by name

# Background execution
python script.py &    # Run in background
nohup python script.py > log.txt 2>&1 &   # Run detached from terminal (common on servers)
jobs                  # View background jobs
fg                    # Bring to foreground
```

> For long-running tasks on remote servers, use `nohup ... &` or `tmux` / `screen` to avoid task interruption due to terminal disconnection.

## 8. Practical Combinations

```bash
# Count files and sizes
ls -lh | awk '{print $5, $9}'
du -sh data/          # Total size of directory
df -h                 # Disk space

# Compression and decompression
tar -czf archive.tar.gz dir/    # Pack and compress
tar -xzf archive.tar.gz         # Decompress
zip -r archive.zip dir/
unzip archive.zip

# Find files
find . -name "*.log"            # Find by name
find . -size +100M              # Find by size
which python                    # Path of the command
```

## 9. WSL and Development Environment

Under Windows, WSL (Windows Subsystem for Linux) is recommended:

```bash
# Install (in administrator PowerShell)
wsl --install

# Enter the Ubuntu subsystem
wsl

# Configure the development environment in WSL
sudo apt update
sudo apt install build-essential git python3 python3-pip
```

VS Code connects to WSL: after installing the Remote - WSL extension, `code .` opens directly.

## 10. Summary

- Navigation: `cd` / `ls` / `pwd`; operations: `cp` / `mv` / `rm` / `mkdir`
- The three musketeers: `grep` for searching, `sed` for replacing, `awk` for column-based processing
- The pipe `|` and redirection `>` are the core of combining commands
- Permissions: `chmod 755` / `chmod +x`; processes: `ps` / `kill` / `nohup`

The next article will introduce Bash programming: organizing commands into reusable scripts.