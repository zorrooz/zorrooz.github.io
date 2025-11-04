const n=`---\r
title: "Bash Shell Script Programming"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Shell", "Bash", "Command Line", "Script Programming", "Automation"]\r
draft: false\r
description: "Automation applications of Bash Shell scripts in bioinformatics data processing"\r
---\r
\r
# Bash Shell Script Programming\r
\r
## Basic Syntax\r
\r
\`\`\`bash\r
#!/bin/bash\r
\r
# Variable definition\r
input_file="data.fastq"\r
output_dir="results"\r
\r
# Conditional judgment\r
if [ -f "$input_file" ]; then\r
    echo "File exists"\r
else\r
    echo "File does not exist"\r
fi\r
\`\`\`\r
\r
## Loop Processing\r
\r
\`\`\`bash\r
# Iterate through files\r
for file in *.fastq; do\r
    echo "Processing file: $file"\r
done`;export{n as default};
