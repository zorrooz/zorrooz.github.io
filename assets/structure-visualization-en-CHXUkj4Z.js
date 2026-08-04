const n=`---
title: "Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice"
date: "2026-08-04"
author: "zorrooz"
tags: ["Structure Visualization","PyMOL","ChimeraX","Tutorial"]
draft: false
description: "Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display"
---

# Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice

Structure visualization is an everyday tool in structural biology. Using two major software packages, PyMOL and UCSF ChimeraX, this article covers everything from PDB data retrieval to publication-quality rendering, encompassing all the core operations needed in daily research.

## 1. Data Source: The PDB Database

### 1.1 Search and Download

[wwPDB / RCSB PDB](https://www.rcsb.org/) is the global repository for structural data:

\`\`\`bash
# Command-line download (by PDB ID)
wget https://files.rcsb.org/download/1CRN.pdb
wget https://files.rcsb.org/download/7A0A.pdb   # Cryo-EM structure example
\`\`\`

### 1.2 Core Contents of a PDB File

\`\`\`
HEADER    PLANT SEED PROTEIN             08-JAN-81   1CRN
TITLE     CRAMBIN
ATOM      1  N   THR A   1      11.106  16.700  17.083  1.00 20.14           N
ATOM      2  CA  THR A   1      11.579  17.437  17.166  1.00 20.55           C
...
HELIX    1   1 THR A    1  SER A    8  1                                  8
SHEET    1   A 2 ILE A   9  ARG A  12  0
\`\`\`

- \`ATOM\` lines: atomic coordinates (residue, chain, x/y/z, B-factor, element)
- \`HELIX\` / \`SHEET\`: secondary structure annotations
- \`CONECT\`: covalent connections
- The newer \`mmCIF\` format (\`.cif\`) has become mainstream and carries more complete information

## 2. PyMOL Basics

### 2.1 Launch and Load

\`\`\`bash
pymol 1CRN.pdb          # Load directly
# Or after launching
# File > Open
\`\`\`

### 2.2 Basic Display Commands

\`\`\`python
# Load
fetch 1crn               # Fetch from the network
load 1CRN.pdb

# Display modes
show cartoon             # Cartoon (backbone trace)
show sticks              # Sticks (atomic details)
show surface             # Surface
hide lines               # Hide lines

# Coloring
color cyan               # Overall color
spectrum count           # Rainbow colors (by residue number)
color red, resi 1-10     # Color a specific region
\`\`\`

### 2.3 Selections and Operations

\`\`\`python
# Selection syntax
select helix, ss h                   # All helices
select sheet, ss s                   # All beta sheets
select active_site, resi 15+20+45    # Specific residues
select ligand, resn HEM              # Ligand (heme)

# Object operations
show sticks, active_site
color yellow, active_site
zoom active_site
\`\`\`

### 2.4 Measurements

\`\`\`python
# Distance
distance d1, /1CRN//A/15/CA, /1CRN//A/20/CA

# Angle / dihedral
angle a1, /1CRN//A/15/CA, /1CRN//A/16/CA, /1CRN//A/17/CA
dihedral dh1, /1CRN//A/15/CA, /1CRN//A/16/CA, /1CRN//A/17/CA, /1CRN//A/18/CA
\`\`\`

### 2.5 Interaction Analysis

\`\`\`python
# Hydrogen bonds
distance hbond, (resn HEM), (resn HIS)

# Contact interface
select interface, chain A within 4.5 of chain B
\`\`\`

## 3. PyMOL Script Batch Processing

Write commonly used operations as scripts for reproducible execution:

\`\`\`python
# render_view.py
fetch 1crn
hide everything
show cartoon
color spectrum, ss
set cartoon_transparency, 0.2

# Zoom into the binding site and render
zoom resi 1-20
set ray_shadows, 1
set ray_opaque_background, 0
ray 1200, 900
png 1crn_view.png, dpi=300
\`\`\`

\`\`\`bash
pymol -cq render_view.py    # -c command-line mode, -q quiet
\`\`\`

## 4. UCSF ChimeraX Basics

ChimeraX has a friendly, modern interface (Qt GUI) and supports batch operations via the \`open\` command:

### 4.1 Opening Structures

\`\`\`bash
chimerax 1CRN.pdb
# Or via command
open 1crn
open 7a0a
\`\`\`

### 4.2 Common Commands

\`\`\`python
# Display modes
cartoon
stick
surface

# Coloring
color bychain          # Color by chain
color byhetero         # Distinguish ligands/ions
color byattribute bfactor   # Color by B-factor

# Selections
select :15-20          # Residues 15-20
select /A               # Chain A
select ligand
\`\`\`

### 4.3 Cryo-EM Density Map Display (ChimeraX's Strength)

\`\`\`python
# Open a density map (.mrc)
open map.mrc

# Adjust isosurface level
volume level 0.01
volume #2 level 0.05

# Superpose density map and model
open model.pdb
open map.mrc
volume #2 level 0.02
color #1 cornflowerblue
transparency #2 30

# Local density inspection
volume zone #2 near :45-60 radius 5
\`\`\`

**ChimeraX's \`volume zone\` is the core tool for checking local density quality**: display density around residues to assess whether side chains are discernible.

## 5. Key Points for High-Quality Rendering

### 5.1 Color Principles

- Cartoon coloring: by chain (\`bychain\`) or by domain
- Key residues: use a single highlight color (yellow/red/blue), avoid overusing rainbows
- Ligands: \`byhetero\` automatically assigns element-based colors (C gray, N blue, O red, S yellow)

### 5.2 Lighting and Materials

\`\`\`python
# PyMOL
set ray_shadows, 1
set specular, 0.5
set ambient, 0.3

# ChimeraX
set lightMode full
graphics silhouettes true
\`\`\`

### 5.3 Resolution Output

\`\`\`bash
# PyMOL publication-quality rendering
ray 2400, 1800
png figure.png, dpi=600

# ChimeraX
save figure.png width 2400 height 1800 supersample 3
\`\`\`

## 6. Structure Comparison and Superposition

\`\`\`python
# PyMOL: align two structures
align model2, model1

# ChimeraX: matchmaker
matchmaker #2 to #1
\`\`\`

**RMSD** is the standard metric for measuring structural similarity (typically reported as Cα or all-atom RMSD).

## 7. Interactive Viewing and Sharing

- **Mol* / NGL**: web-based 3D viewing (built into PDBe, RCSB)
- **PyMOL Web**: export sessions as interactive HTML pages
- Commonly used for scientific communication: PyMOL sessions (\`.pse\`), ChimeraX sessions (\`.cxs\`)

## 8. Summary

- PDB is the data source; mmCIF is the new standard
- PyMOL: powerful commands, scriptable batch processing (\`pymol -cq script.py\`)
- ChimeraX: modern GUI + density map handling (\`volume\`, \`volume zone\`) is irreplaceable
- Publication-quality rendering: hide the clutter → highlight the key → high-quality lighting → high-resolution export

The next article will cover atomic modeling: how to build and refine atomic models from density maps.`;export{n as default};
