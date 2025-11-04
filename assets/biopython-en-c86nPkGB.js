const n=`---\r
title: "Biopython Bioinformatics Applications"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Python", "Biopython", "Bioinformatics", "Sequence Analysis", "Genomics"]\r
draft: false\r
description: "Application of Biopython library in sequence analysis and genomic data processing"\r
---\r
\r
# Biopython Bioinformatics Applications\r
\r
## Sequence Operations\r
\r
\`\`\`python\r
from Bio import SeqIO\r
from Bio.Seq import Seq\r
\r
# Read FASTA file\r
for record in SeqIO.parse("sequences.fasta", "fasta"):\r
    print(record.id, len(record.seq))\r
\r
# Sequence translation\r
dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")\r
protein_seq = dna_seq.translate()`;export{n as default};
