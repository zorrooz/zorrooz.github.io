const n=`---
title: "Mainstream Virtual Screening Tools"
date: "2026-08-04"
author: "zorrooz"
tags: ["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"]
draft: false
description: "A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies"
---

# Mainstream Virtual Screening Tools

Virtual Screening (VS) uses computational methods to search for potential active molecules within million-scale compound libraries, serving as a core approach in the early stages of drug discovery. This article introduces mainstream docking software, complete workflows, and key considerations.

## 1. Basic Logic of Virtual Screening

\`\`\`
Receptor structure (protein/nucleic acid)
   ↓
① Receptor preparation (hydrogen addition, protonation, grid)
   ↓
② Compound library (ZINC / Enamine / custom library)
   ↓
③ Ligand preparation (3D conformation, protonation, force field parameters)
   ↓
④ Molecular docking (conformational sampling + scoring)
   ↓
⑤ Ranking and filtering (Top 100-1000)
   ↓
⑥ Experimental validation (IC50 / activity assay)
\`\`\`

## 2. Mainstream Docking Software

### 2.1 AutoDock Vina (Top Open-Source Choice)

The most popular open-source docking software, offering a good balance between speed and accuracy:

\`\`\`bash
# 1. Receptor and ligand preparation (using MGLTools or scripts)
#    Receptor: pdb → pdbqt (add hydrogens, merge non-polar hydrogens)
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

# 4. Interpretation: lower affinity (kcal/mol) is better
#  mode |   affinity | dist from best mode
# -----+------------+----------------------
#    1 |   -9.2     | 0.000
#    2 |   -8.7     | 2.315
\`\`\`

### 2.2 Glide (Schrödinger, Commercial)

The pharmaceutical industry standard:

\`\`\`bash
# Preparation workflow (Maestro GUI or command line)
glide -overwrite -adjust \\
  -receptor receptor.maegz \\
  -ligand ligands.maegz \\
  -p glide_docking.inp \\
  -NJOBS 16

# Three precision modes
# HTVS (high-throughput, fast) → SP (standard precision) → XP (extra precision, slow)
\`\`\`

**Advantages**: Accurate ligand flexibility handling, well-calibrated scoring functions, and support from a large pharmaceutical ecosystem.

### 2.3 DOCK 6 (Open Source)

Classic geometric matching algorithm (shape matching):

\`\`\`bash
dock6 -i dock.in -o dock.out
# Sphere generation → orientation search → scoring
# Strengths: binding pocket shape matching, primary screening of large libraries
\`\`\`

### 2.4 Other Commonly Used Software

| Software | Features |
|----------|----------|
| **AutoDock4** | Classic genetic algorithm, high accuracy but slow |
| **rDock** | Open source, suitable for large libraries, easy to script |
| **GOLD** | Genetic algorithm, developed by CCDC |
| **LeDock** | Fast and easy to use (free for academia) |
| **Plants** | Based on Ant Colony Optimization |
| **smina** | Enhanced version of Vina (supports custom scoring) |

## 3. Receptor Preparation (Critical First Step)

### 3.1 Structure Sources

- Crystal/Cryo-EM structures (PDB)
- AlphaFold predicted structures (when no experimental structure is available)
- Note: Handle flexibility (side chains) near the binding pocket

### 3.2 Preparation Key Points

\`\`\`bash
# Common tools
# Schrödinger Protein Prep Wizard (commercial)
# ADFRsuite / MGLTools prepare_receptor4.py (free)

python prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt

# Key steps
# 1. Remove water molecules (except conserved waters)
# 2. Add hydrogens, determine protonation states (near pH 7.4)
# 3. Assign bond orders and charges (Amber/CHARMM force field)
# 4. Define binding pocket (known active site or cavity detection: fpocket / DoGSiteScorer)
\`\`\`

## 4. Ligand Libraries and Preparation

### 4.1 Compound Library Sources

| Library | Size | Features |
|---------|------|----------|
| ZINC20 | Billions | Free, downloadable, includes 3D conformations |
| Enamine REAL | > 40 billion | Highly synthesizable, commercially available |
| ChEMBL | Activity data | Known active compounds |
| Custom library | Tailored | Derivative design |

### 4.2 Ligand Preparation

\`\`\`bash
# 3D conformation generation
obabel ligand.smi -O ligand.sdf --gen3d -p 7.4

# Multiple conformations (important for flexibility, especially ring systems)
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

| Type | Representative | Principle |
|------|----------------|-----------|
| Force field-based | DOCK, AutoDock4 | van der Waals + electrostatics + internal energy |
| Empirical | Glide, Vina | Statistical fitting to experimental data |
| Knowledge-based | GOLD/ASP | Statistical analysis of atom pair distance distributions |

### 5.2 Limitations of Scoring Functions (Must Know)

- **Imprecise**: Correlation coefficients are typically 0.3–0.6; they can only rank, not provide absolute affinities
- Entropic effects, solvation effects, and induced fit are difficult to describe accurately
- **Recommendation**: Docking scores + cross-validation with multiple tools + molecular dynamics refinement (MM/PBSA, FEP)

## 6. Complete High-Throughput Screening Workflow (HTS in Practice)

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
Coarse screening (HTVS / Vina, million-scale) → Top 5-10%
  ↓
Standard precision (SP / multiple conformations) → Top 10-20%
  ↓
High precision (XP / induced fit) → Top 100-1000
  ↓
Consensus scoring (cross-validation with multiple software) → Top 50-200
  ↓
Visual inspection (binding mode plausibility) → Top 10-50
  ↓
Experimental testing
\`\`\`

### 6.2 Consensus Scoring

Scoring functions from different software are complementary; taking the intersection or averaging ranks can significantly improve enrichment rates:

\`\`\`bash
# Example: intersection of Top 100 from Vina + Glide + LeDock
# Or average after rank normalization
\`\`\`

## 7. Validation and Advanced Methods

### 7.1 Benchmark Testing

- **DUD-E**: Standard dataset (containing actives and decoys) to evaluate a software's discriminating ability
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
- Combined strategies represent the mainstream of modern virtual screening

## 8. Common Pitfalls

| Pitfall | Consequence | Avoidance |
|---------|-------------|-----------|
| Incorrect pocket definition | Entire library mispositioned | Use known active sites or experimental information |
| Rigid receptor | Missed induced fit | Flexible side chains / ensemble docking |
| Incorrect protonation states | Distorted electrostatic scoring | pH-dependent preparation (pH 7.4) |
| Overconfidence in scoring functions | High false positive rate | Consensus scoring + visual inspection |
| Ignoring synthesizability | Screened compounds cannot be synthesized | Filter PAINS, drug-likeness rules (Lipinski/Veber) |

## 9. Summary

- Three tiers: **Vina** (open-source, fast) → **Glide** (commercial standard) → **DOCK/rDock** (large libraries)
- Workflow: receptor preparation → ligand library → docking → multi-stage funnel → experimental validation
- Scoring functions can only rank; consensus scoring + MD refinement improves reliability
- The endpoint of virtual screening is always experimentation

This completes the computational biology track: protein design tools + virtual screening tools, with the two main threads forming a closed-loop "design-screen" computational drug discovery pipeline.`;export{n as default};
