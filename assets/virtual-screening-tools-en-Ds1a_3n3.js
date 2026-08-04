const n=`---
title: "Mainstream Tools for Virtual Screening"
date: "2026-08-04"
author: "zorrooz"
tags: ["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"]
draft: false
description: "A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies"
---

# Mainstream Tools for Virtual Screening

Virtual screening (VS) uses computation to search for potential active molecules in million-scale compound libraries and is a core approach in the early stage of drug discovery. This article introduces mainstream docking software, the complete workflow, and important considerations.

## 1. Basic Logic of Virtual Screening

\`\`\`
Receptor structure (protein/nucleic acid)
   ↓
① Receptor preparation (add hydrogens, protonation, grid)
   ↓
② Compound library (ZINC / Enamine / self-built library)
   ↓
③ Ligand preparation (3D conformation, protonation, force field parameters)
   ↓
④ Molecular docking (conformational sampling + scoring)
   ↓
⑤ Ranking and selection (Top 100-1000)
   ↓
⑥ Experimental validation (IC50 / activity assay)
\`\`\`

## 2. Mainstream Docking Software

### 2.1 AutoDock Vina (Top Open-Source Choice)

The most popular open-source docking software, offering a good balance between speed and accuracy:

\`\`\`bash
# 1. Receptor and ligand preparation (using MGLTools or scripts)
#    Receptor: pdb → pdbqt (add hydrogens, merge nonpolar hydrogens)
#    Ligand: sdf/mol2 → pdbqt

# 2. Configuration file (conf.txt)
receptor = receptor.pdbqt
ligand = ligand.pdbqt
center_x = 10.5
center_y = 20.3
center_z = -5.8
size_x = 22
size_y = 22
size_z = 22
exhaustiveness = 16
num_modes = 9

# 3. Run
vina --config conf.txt --out results.pdbqt --log log.txt

# 4. Interpretation: the lower the affinity (kcal/mol), the better
#  mode |   affinity | dist from best mode
# -----+------------+----------------------
#    1 |   -9.2     | 0.000
#    2 |   -8.7     | 2.315
\`\`\`

### 2.2 Glide (Schrödinger, Commercial)

The standard in the pharmaceutical industry:

\`\`\`bash
# Prepare workflow (Maestro GUI or command line)
glide -overwrite -adjust \\
  -receptor receptor.maegz \\
  -ligand ligands.maegz \\
  -p glide_docking.inp \\
  -NJOBS 16

# Three precision modes
# HTVS (high-throughput, fast) → SP (standard precision) → XP (high precision, slow)
\`\`\`

**Advantages**: Accurate flexible ligand handling, well-calibrated scoring, and support from the large pharma ecosystem.

### 2.3 DOCK 6 (Open Source)

Classic geometric matching algorithm:

\`\`\`bash
dock6 -i dock.in -o dock.out
# Sphere generation → orientation search → scoring
# Strengths: binding pocket shape matching, initial screening of large libraries
\`\`\`

### 2.4 Other Commonly Used Software

| Software | Features |
|------|------|
| **AutoDock4** | Classic genetic algorithm, high accuracy but slow |
| **rDock** | Open source, suitable for large libraries, easy to script |
| **GOLD** | Genetic algorithm, developed by CCDC |
| **LeDock** | Fast and easy to use (free for academia) |
| **Plants** | Based on Ant Colony Optimization |
| **smina** | Enhanced version of Vina (supports custom scoring) |

## 3. Receptor Preparation (Critical First Step)

### 3.1 Structure Sources

- Crystal/cryo-EM structures (PDB)
- AlphaFold-predicted structures (when no experimental structure is available)
- Note: handle flexibility (side chains) near the binding pocket

### 3.2 Key Preparation Points

\`\`\`bash
# Common tools
# Schrödinger Protein Prep Wizard (commercial)
# ADFRsuite / MGLTools prepare_receptor4.py (free)

python prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt

# Key steps
# 1. Remove water molecules (except conserved waters)
# 2. Add hydrogens, determine protonation states (around pH 7.4)
# 3. Assign bond orders and charges (Amber/CHARMM force fields)
# 4. Define the binding pocket (known active site or cavity detection: fpocket / DoGSiteScorer)
\`\`\`

## 4. Ligand Libraries and Preparation

### 4.1 Compound Library Sources

| Library | Scale | Features |
|----|------|------|
| ZINC20 | Billions | Free, downloadable, with 3D conformations |
| Enamine REAL | > 40 billion | Highly synthesizable, commercially available |
| ChEMBL | Bioactivity data | Known active compounds |
| Self-built library | Custom | Derivative design |

### 4.2 Ligand Preparation

\`\`\`bash
# 3D conformation generation
obabel ligand.smi -O ligand.sdf --gen3d -p 7.4

# Multiple conformations (flexibility is important, especially for ring systems)
obabel ligand.sdf -O conf.sdf --conformer --nconf 50

# Protonation and charges
# Schrödinger LigPrep / Open Babel / RDKit

# RDKit preparation (Python)
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
\`\`\`

## 5. Scoring Functions

### 5.1 Three Types of Scoring Functions

| Type | Representatives | Principle |
|------|------|------|
| Force-field-based | DOCK, AutoDock4 | van der Waals + electrostatics + internal energy |
| Empirical | Glide, Vina | Statistical fitting to experimental data |
| Knowledge-based | GOLD/ASP | Statistics of atom-pair distance distributions |

### 5.2 Limitations of Scoring Functions (Must-Know)

- **Not precise**: correlation coefficients are typically 0.3–0.6; they can only rank, not provide absolute affinities
- Entropy effects, solvent effects, and induced fit are difficult to describe accurately
- **Recommendation**: docking scores + cross-validation with multiple tools + refined calculations using molecular dynamics (MM/PBSA, FEP)

## 6. Complete High-Throughput Screening Workflow (Hands-On HTS)

\`\`\`bash
# 1. Library preparation: split SDF + convert to pdbqt
mkdir -p ligands_pdbqt
obabel library.sdf -O ligands_pdbqt/lig_.pdbqt -m -p 7.4

# 2. Batch docking (smina/Vina multi-file mode)
smina -r receptor.pdbqt -l library.sdf \\
  --center_x 10.5 --center_y 20.3 --center_z -5.8 \\
  --size_x 22 --size_y 22 --size_z 22 \\
  --num_modes 1 --cpu 32 \\
  -o docked.sdf --log scores.txt

# 3. Ranking (by affinity)
sort -k2 -n scores.txt | head -100
\`\`\`

### 6.1 Multi-Stage Funnel Strategy

\`\`\`
Coarse screening (HTVS / Vina, million-scale) → top 5-10%
  ↓
Standard precision (SP / multiple conformations) → top 10-20%
  ↓
High precision (XP / induced fit) → Top 100-1000
  ↓
Consensus scoring (cross-validation with multiple software tools) → Top 50-200
  ↓
Visual inspection (reasonableness of binding modes) → Top 10-50
  ↓
Experimental testing
\`\`\`

### 6.2 Consensus Scoring

Scoring functions from different software tools complement each other; taking intersections/ranking averages can significantly improve enrichment rates:

\`\`\`bash
# Example: take the intersection of Top 100 from Vina + Glide + LeDock
# Or average after normalizing ranks
\`\`\`

## 7. Validation and Advanced Methods

### 7.1 Benchmarking

- **DUD-E**: a standard dataset (containing actives and decoys) used to evaluate a software's ability to discriminate
- Metrics: AUC, enrichment factors (EF1%, EF5%)

### 7.2 Molecular Dynamics Post-Processing

\`\`\`bash
# Docking results → MD simulation (GROMACS/AMBER) → MM/PBSA binding free energy
gmx mdrun -deffnm complex
gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ct traj.xtc ...
\`\`\`

### 7.3 Combining Structure-Based and Ligand-Based Approaches

- **Structure-based** (SBDD): docking + MD + FEP
- **Ligand-based** (LBDD): pharmacophore, QSAR, similarity search
- Combination strategies are the mainstream of modern virtual screening

## 8. Common Pitfalls

| Pitfall | Consequence | Avoidance |
|------|------|------|
| Incorrect pocket definition | Entire library misaligned | Use a known active site or experimental information |
| Rigid receptor | Misses induced fit | Flexible side chains / ensemble docking |
| Incorrect protonation states | Distorted electrostatic scoring | pH-dependent preparation (pH 7.4) |
| Overconfidence in scoring functions | High false-positive rate | Consensus scoring + visual inspection |
| Ignoring synthesizability | Hits cannot be synthesized | Filter PAINS and drug-likeness rules (Lipinski/Veber) |

## 9. Summary

- Three tiers: **Vina** (open source, fast) → **Glide** (commercial standard) → **DOCK/rDock** (large libraries)
- Workflow: receptor preparation → ligand library → docking → multi-stage funnel → experimental validation
- Scoring functions can only rank; consensus scoring + MD refinement improve reliability
- The end point of virtual screening is always experimentation

This completes the computational biology direction: protein design tools + virtual screening tools, with two main threads forming a closed loop of "design-screen" computational drug discovery.`;export{n as default};
