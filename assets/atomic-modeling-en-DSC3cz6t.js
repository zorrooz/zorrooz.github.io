const n=`---
title: "Atomic Modeling and Refinement"
date: "2026-08-04"
author: "zorrooz"
tags: ["Atomic Modeling", "Coot", "phenix", "Structure Refinement", "Tutorial"]
draft: false
description: "Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics"
---

# Atomic Modeling and Refinement

After obtaining a high-quality density map, the next step is to build and refine an atomic model, ultimately yielding publishable, analyzable coordinate files. This article covers the complete workflow from automatic modeling to quality assessment.

## 1. Workflow Overview

\`\`\`
Density map (.mrc)
   ↓
① Initial model (automatic modeling / homologous model / AF prediction)
   ↓
② Model placement into density (docking / rigid body)
   ↓
③ Iterative refinement: Coot manual correction ↔ phenix automatic refinement
   ↓
④ Quality assessment (MolProbity / EMRinger / Q-score)
   ↓
⑤ Validation and publication (PDB deposition)
\`\`\`

## 2. Initial Model Sources

### 2.1 Automatic Modeling (Mainstream in the Deep Learning Era)

- **ModelAngelo**: Automatic modeling based on the AlphaFold architecture, building models directly from density maps
- **DeepTracer**: Fast fully automatic modeling (GPU)
- **phenix.map_to_model**: Classic automatic backbone/side-chain building

\`\`\`bash
# ModelAngelo usage example
model_angelo --halfmaps half1.mrc half2.mrc \\
  --output-dir model_angelo_out \\
  --num-workers 8

# phenix.map_to_model
phenix.map_to_model map.mrc resolution=3.0 \\
  output_model=auto_model.pdb
\`\`\`

### 2.2 Homology Modeling and Molecular Replacement

- If a homologous structure exists: use **phenix.dock_in_map** to place it into the density
- AlphaFold-predicted models: directly used as initial models (commonly used in cryo-EM)
- At low resolution (> 5 Å), place rigid bodies first, then gradually refine

## 3. Coot: Manual Modeling and Correction

Coot is the industry standard for manual model correction (GUI tool).

### 3.1 Common Operations

\`\`\`
File > Open Coordinates...   # Load model
File > Open Map...           # Load density map (FSC correlation calculated automatically)

# Basic corrections
1. Backbone adjustment: drag Cα (Real Space Refine Zone)
2. Side-chain orientation: Rotamers
3. Missing residues: Add Terminal Residue / Build
4. Delete extra: Delete Residue Range
5. Mutate residues: Mutate & Autofit
6. Ligand building: Ligand Builder / restraints
\`\`\`

### 3.2 Key Concept: Real Space Refine

\`\`\`
Calculate > Real Space Refine Zone
\`\`\`

- Select a density region, and Coot optimizes the geometry and density fit of that region
- Repeat: fix → refine → check → fix again
- After each refinement round, return to Coot to check problematic regions (Ramachandran outliers, clashes)

### 3.3 Ligand Modeling

\`\`\`bash
# Generate ligand restraints
phenix.elbow ligand.smiles    # Generate chemical restraints from SMILES
# Or
phenix.ligand_idealization ligand.cif
\`\`\`

Use \`Ligand > Ligand Builder\` in Coot for manual building, then Real Space Refine.

## 4. phenix Automatic Refinement

### 4.1 Basic Refinement

\`\`\`bash
phenix.real_space_refine model.pdb \\
  map.mrc \\
  resolution=3.0 \\
  output.prefix=refined \\
  nproc=8
\`\`\`

### 4.2 Refinement Strategies

\`\`\`bash
# Multiple rounds of iteration: global first, then local
phenix.real_space_refine model.pdb map.mrc \\
  resolution=3.0 \\
  strategy=individual_sites+individual_adp \\
  restraints=ligand.cif \\
  nproc=8
\`\`\`

Commonly used strategy combinations:
- \`rigid_body\`: rigid-body refinement first at low resolution
- \`individual_sites\`: per-atom coordinate refinement (high resolution)
- \`individual_adp\`: B-factor refinement
- \`torsion\`: torsion-angle refinement (commonly used for proteins)
- \`local_grid\`: local grid search (flexible regions)

### 4.3 Automatic Correction Loop (NCS and Water Molecules)

\`\`\`bash
# Find water molecules (at high resolution)
phenix.find_peaks_holes map.mrc map.mrc resolution=3.0

# NCS restraints (significantly improves quality when symmetric subunits are present)
phenix.real_space_refine ... strategy=individual_sites+individual_adp ncs_search.enabled=True
\`\`\`

## 5. Model Quality Assessment

### 5.1 Geometric Quality: MolProbity

\`\`\`bash
phenix.molprobity refined.pdb
\`\`\`

Core metrics:

| Metric | Good standard (3 Å) | Meaning |
|------|----------------|------|
| Ramachandran favored | > 95% | Backbone dihedral angle validity |
| Ramachandran outliers | < 0.5% | Outlier backbone |
| Rotamer outliers | < 1% | Abnormal side-chain conformations |
| Clashscore | < 5 | Atomic collisions |
| MolProbity score | < 2 | Composite geometry score |

### 5.2 Density Fit Quality

- **EMRinger**: assesses side-chain-to-density fit (> 2 indicates good)
- **Q-score**: assesses local density quality (1.0 = perfect, 0.5 = acceptable)
- **CC_mask**: model-density correlation coefficient (> 0.8 good)
- **map-model FSC**: FSC comparison between the model and the two half-maps

\`\`\`bash
phenix.map_model_cc refined.pdb map.mrc
phenix.em_ringer refined.pdb map.mrc
\`\`\`

### 5.3 Overfitting Detection

**Free component (CC-free)**: correlation computed independently from the half-map used in refinement. Must be monitored during refinement:

\`\`\`
Before refinement: CC_work ≈ CC_free
After refinement: CC_work significantly > CC_free → sign of overfitting
\`\`\`

**Standard practice**: use only one half-map for refinement (or full map + independent half-map validation), and use the other half-map to compute CC_free.

## 6. Complete Refinement Workflow Example

\`\`\`bash
# 1. Automatic modeling
model_angelo --halfmaps half1.mrc half2.mrc --output-dir auto

# 2. Initial refinement
phenix.real_space_refine auto/model.pdb map.mrc \\
  resolution=3.1 strategy=rigid_body+individual_sites \\
  output.prefix=round1

# 3. Coot manual correction (interactive)
coot round1_refined.pdb map.mrc

# 4. Second refinement round (add B-factors)
phenix.real_space_refine coot_out.pdb map.mrc \\
  resolution=3.1 strategy=individual_sites+individual_adp \\
  output.prefix=round2

# 5. Quality assessment
phenix.molprobity round2_refined.pdb
phenix.map_model_cc round2_refined.pdb map.mrc
phenix.em_ringer round2_refined.pdb map.mrc

# 6. Iterate until all metrics pass
\`\`\`

## 7. PDB Deposition

Before publication, prepare:

1. Coordinate file (.pdb / .cif)
2. Density maps (half-maps + full map, .mrc)
3. Validation reports (MolProbity, EMRinger, etc.)
4. FSC curves and mask information

Submit via [wwPDB deposition](https://deposit.wwpdb.org/) to obtain PDB ID and EMDB ID.

## 8. Common Problems

| Problem | Cause | Solution |
|------|------|------|
| Blurred side-chain density | Insufficient resolution / incorrect B-factors | Increase resolution, check orientation distribution |
| Many Ramachandran outliers | Incorrect backbone building | Rebuild the region in Coot based on density |
| Overfitting | Excessive refinement | Monitor with CC_free, reduce iterations |
| Poor ligand fit | Incomplete restraints | Regenerate with phenix.elbow |
| Poor local region quality | Flexibility/conformational heterogeneity | Separate conformations with 3D classification, local refinement |

## 9. Summary

- Three sources for model building: automatic modeling (ModelAngelo/DeepTracer), homology/AF prediction, de novo building
- Alternating iterations of Coot manual correction and phenix automatic refinement is the standard workflow
- Triple quality check: geometry (MolProbity) + fit (EMRinger/Q-score) + overfitting (CC_free)
- Sub-3 Å maps should aim for side-chain-level confidence; ligands and modifications require separate validation

This completes the structural biology pipeline: data processing → technical review → visualization → atomic modeling, forming a complete loop from experiment to model.`;export{n as default};
