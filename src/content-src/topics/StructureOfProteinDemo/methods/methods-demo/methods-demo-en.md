---
title: "Experimental Methods - Protein Structure Determination"
date: "2025-07-01"
author: "zorrooz"
tags: ["Methods", "Cryo-EM", "Sample Preparation", "Data Collection"]
draft: false
description: "Demonstration of experimental methods for protein structure determination"
---

# Experimental Methods

**The following content is entirely fictional!**

## 1. Sample Preparation

### Protein Expression and Purification
```bash
# Expression system
Host: E. coli BL21(DE3)
Induction conditions: 0.5 mM IPTG, 18°C, 16 hours

# Purification steps
1. Ni-NTA affinity chromatography
2. Gel filtration chromatography (Superdex 200)
3. Concentrate to 3 mg/mL
```

### Negative Stain EM Screening
- Stain: 2% Uranyl acetate
- Data acquisition: 200 kV, 50,000x magnification
- Screening criteria: Monodisperse, intact particles

## 2. Cryo-EM Data Collection

### Sample Preparation
- Apply 3 μL sample onto Quantifoil R1.2/1.3 grid
- Blot with filter paper, rapid freezing in liquid nitrogen
- Working temperature: -180°C

### Acquisition Parameters
| Parameter | Value |
|------|-----|
| Acceleration voltage | 300 kV |
| Pixel size | 1.5 Å |
| Electron dose | 50 e⁻/Å² |
| Number of frames | 50 frames |

## 3. Data Processing

### Motion Correction
- Software: MotionCor2
- Corrects inter-frame motion and drift

### Particle Picking
- Method: Template matching
- Software: RELION 3.1
- Number of particles: ~50,000

### 2D Classification
- Number of classes: 20
- Selection criteria: Clear features

## 4. 3D Reconstruction

### Initial Model
- Method: Random conical tilt reconstruction
- Resolution: ~8 Å

### Refined Reconstruction
- Iterations: 20 rounds
- Resolution: 3.0 Å
- Resolution criterion: FSC 0.143

## 5. Model Building and Refinement

### Initial Model Building
- Software: Coot
- Method: Automated modeling + manual correction

### Structure Refinement
- Software: Phenix.refine
- Refinement items:
  - Atomic coordinates
  - B-factors
  - Occupancies

### Quality Validation
- Ramachandran plot: 96.5% favored regions
- Geometry deviations: Bond length 0.01 Å, bond angle 1.2°
- MolProbity score: 1.8