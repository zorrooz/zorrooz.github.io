const n=`---
title: "Cryo-EM Single Particle Analysis: Full Data Processing Workflow"
date: "2026-08-04"
author: "zorrooz"
tags: ["Cryo-Electron Microscopy", "cryo-EM", "Data Processing", "RELION"]
draft: false
description: "From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment"
---

# Cryo-EM Single Particle Analysis: Full Data Processing Workflow

Cryo-EM single particle analysis (SPA) is a core technique for determining near-atomic resolution structures of biological macromolecules. This article outlines the complete pipeline from micrographs to atomic models, introducing the purpose, common software, and key parameters for each step.

## 1. Workflow Overview

\`\`\`
Sample preparation (cryo-grid preparation)
   ↓
Data acquisition (Titan Krios / Glacios, etc.)
   ↓
① Preprocessing: motion correction + CTF estimation
   ↓
② Particle picking
   ↓
③ 2D classification (removing bad particles)
   ↓
④ Initial model (Ab initio)
   ↓
⑤ 3D classification (heterogeneity analysis)
   ↓
⑥ 3D refinement
   ↓
⑦ Resolution assessment (FSC curve)
   ↓
⑧ Atomic model building and refinement
\`\`\`

## 2. Sample Preparation and Data Acquisition

### 2.1 Cryo-grid Preparation

- Apply the protein sample (typically 0.5–5 mg/mL) to a holey carbon grid
- Use a Vitrobot / Leica GP or similar device to rapidly plunge into liquid ethane, forming **vitreous ice**
- Ideal ice thickness: slightly thinner than the particle diameter; ice that is too thick reduces signal-to-noise, while ice that is too thin may deform particles

### 2.2 Data Acquisition

- 300 kV electron microscope (Titan Krios) + direct electron detector (Gatan K3 / Falcon 4)
- **Super-resolution mode** + dose fractionation (total dose about 40–60 e⁻/Å²)
- Each movie typically contains 30–50 frames for motion correction
- Target: 10,000–20,000 micrographs (high-resolution structures usually require millions of particles)

## 3. Preprocessing

### 3.1 Motion Correction

Sample drift caused by the electron beam blurs images and must be corrected frame by frame:

\`\`\`bash
# Motion correction in RELION 4
relion_run_motioncorr \\
  --i Micrographs/ \\
  --o MotionCorr/ \\
  --j 8 \\
  --pipeline_control MotionCorr/job001/
\`\`\`

Common software: **MotionCor2** (classic), **built-in RELION**, **cryoSPARC Patch Motion Correction**.

### 3.2 CTF Estimation

The CTF (Contrast Transfer Function) describes the phase flipping and defocus effects of electron microscope imaging:

\`\`\`bash
# CTFFIND4 (classic command line)
ctffind --input_micrograph MotionCorr/mic1.mrc \\
        --output_diag mic1_diag.mrc \\
        --defU 20000 --defV 20000

# Gctf (faster GPU version)
Gctf --input 'MotionCorr/*.mrc' --apix 1.0
\`\`\`

Common software: **CTFFIND4**, **Gctf**, **cryoSPARC Patch CTF**.

**Key criteria**: The resolution cutoff is usually taken as 3–4 Å; discard micrographs with excessive astigmatism or severe drift.

## 4. Particle Picking

### 4.1 Traditional Methods

- **Template-based**: use 2D class averages as templates and perform correlation searches
- Software: auto-picking in RELION, Template Picker in cryoSPARC

### 4.2 Deep Learning Approaches (currently mainstream)

\`\`\`bash
# Topaz: CNN-based particle picking
topaz train --model-path model.pkl --epochs 20 \\
  --num-workers 4 --train-data particles.toml

topaz extract --model-path model.pkl \\
  --micrographs micrographs.toml \\
  --output-prefix picked --radius 100
\`\`\`

Mainstream tools: **Topaz**, **crYOLO** (very fast, supports real-time picking), **cryoSPARC Blob Picker**.

> In modern workflows, a blob/template picker is often used first to quickly select a subset for 2D classification, and clean class averages are then used to train Topaz/crYOLO models, significantly improving picking quality for complex samples (membrane proteins, small proteins).

## 5. 2D Classification: Particle Screening

2D classification clusters particles by projection orientation while **removing ice contamination, mis-picked particles, and denatured particles**:

\`\`\`bash
# RELION 2D classification
relion_refine --o Class2D/ \\
  --particle_images particles.star \\
  --ctf --K 100 --iter 25 --tau2_fudge 4 \\
  --flatten_solvent --zero_mask --pool 30 \\
  --per_particle_ctf --j 8
\`\`\`

**Quality criteria**:
- Class averages should show clear secondary structure (α-helices, β-sheets)
- Retain 50–80% of the total particles as "good"
- High-quality class averages are a prerequisite for subsequent 3D reconstruction

## 6. 3D Reconstruction: Initial Model and Classification

### 6.1 Initial Model (Ab initio)

No reference model is required; reconstruction starts from random initial structures:

\`\`\`bash
# cryoSPARC
cryosparc_abinitio ... --num_classes 3
\`\`\`

- Commonly used: **cryoSPARC Ab-initio** (fast and robust), RELION 3D initial model
- Generate 2–4 classes and select the class with the clearest structural features as a reference

### 6.2 3D Classification: Handling Conformational Heterogeneity

Biological samples often exist in multiple conformations (open/closed, ligand-bound/unbound):

\`\`\`bash
# RELION 3D classification
relion_refine --o Class3D/ \\
  --ref initial.mrc --K 4 --iter 40 \\
  --tau2_fudge 8 --particle_diameter 200
\`\`\`

**Strategy**: First use coarse classification (K=4–6) to inspect the conformational distribution, then sub-classify the conformation of interest (K=2–3). Ultimately, each 3D class corresponds to one conformation.

## 7. 3D Refinement

### 7.1 Standard Refinement

\`\`\`bash
relion_refine --o Refine3D/ \\
  --ref class3d.mrc --particle_diameter 200 \\
  --auto_mask --solvent_mask \\
  --ctf --per_particle_ctf --j 8
\`\`\`

Key parameters:
- Start with **C1 symmetry**; for symmetric molecules, specify C/D/O symmetry (e.g., C3, D7)
- A **solvent mask** must be used; otherwise the FSC is artificially high
- Use local refinement to handle flexible regions

### 7.2 Bayesian Polishing

Performs trajectory fitting and weighting at the movie-frame level for particles, improving high-frequency information:

\`\`\`bash
relion_polish --i Refine3D/run_data.star \\
  --o Polish/ \\
  --model_type 2 --nr_iter 5
\`\`\`

### 7.3 CTF Refinement

Refines per-particle defocus, astigmatism, and beam tilt:

\`\`\`bash
relion_ctf_refine --i Polish/particles.star \\
  --o CtfRefine/ \\
  --fit_defocus --fit_magnification --fit_tilt
\`\`\`

## 8. Resolution Assessment: FSC Curve

**FSC (Fourier Shell Correlation)**: Particles are randomly split into two half-sets (half-maps), reconstructed separately, and the correlation coefficient between the two maps is computed:

- **FSC = 0.143 criterion**: gold-standard resolution determination (after half-map correction)
- **FSC = 0.5**: conservative estimate (uncorrected)

\`\`\`bash
relion_postprocess --i Refine3D/run_half1_class001.mrc \\
  --mask mask.mrc --autob_masking \\
  --angpix 1.0
\`\`\`

**Quality checklist**:
| Metric | Good standard |
|------|----------|
| Global resolution (FSC 0.143) | ≤ 3.5 Å (high resolution) |
| Number of particles | ≥ 100,000 (hundreds of thousands typically needed for 3 Å) |
| Orientation distribution | Uniform (no preferred orientation) |
| Density continuity | Side-chain density clearly visible (< 3 Å) |

## 9. Comparison of Common Software Stacks

| Step | RELION (command line) | cryoSPARC (GUI) |
|------|----------------------|------------------|
| Motion correction | Built-in / MotionCor2 | Patch Motion |
| CTF | CTFFIND4 / Gctf | Patch CTF |
| Picking | Template / Topaz integration | Blob / Template / Topaz |
| 2D/3D | Stable and reliable | Fast, interactive |
| Non-uniform refinement | — | **NU-refinement** (powerful for flexible macromolecules) |

> **Workflow recommendation**: Use cryoSPARC for preprocessing and initial reconstruction, and RELION for high-precision refinement and polishing; alternatively, use cryoSPARC for the entire workflow (more user-friendly for beginners).

## 10. Summary

- Preprocessing (motion correction + CTF) determines the upper limit of data quality
- Deep learning-based particle picking (Topaz/crYOLO) significantly improves performance on complex samples
- 2D classification screens particles; 3D classification handles conformational heterogeneity
- Three key refinement steps: solvent mask + Bayesian polishing + CTF refinement
- FSC = 0.143 is the gold standard for resolution

Future articles will cover structure visualization (PyMOL/ChimeraX) and atomic modeling (Coot/phenix), turning density maps into atomic models.`;export{n as default};
