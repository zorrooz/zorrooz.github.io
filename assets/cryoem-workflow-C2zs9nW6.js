const n=`---
title: "冷冻电镜单颗粒分析：数据处理全流程"
date: "2026-08-04"
author: "zorrooz"
tags: ["冷冻电镜", "cryo-EM", "数据处理", "RELION"]
draft: false
description: "从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估"
---

# 冷冻电镜单颗粒分析：数据处理全流程

冷冻电镜单颗粒分析（Single Particle Analysis, SPA）是解析生物大分子近原子分辨率结构的核心技术。本文梳理从显微照片到原子模型的完整流程，介绍各步骤的目的、常用软件与关键参数。

## 1. 流程总览

\`\`\`
样品制备（冷冻制样）
   ↓
数据采集（Titan Krios / Glacios 等）
   ↓
① 预处理：运动校正 + CTF 估计
   ↓
② 颗粒挑选（Particle Picking）
   ↓
③ 2D 分类（筛除坏颗粒）
   ↓
④ 初始模型（Ab-initio）
   ↓
⑤ 3D 分类（异质性分析）
   ↓
⑥ 3D 精修（Refinement）
   ↓
⑦ 分辨率评估（FSC 曲线）
   ↓
⑧ 原子模型搭建与精修
\`\`\`

## 2. 样品制备与数据采集

### 2.1 冷冻制样

- 将蛋白样品（浓度通常 0.5–5 mg/mL）施加到载网（holey carbon grid）
- 用 Vitrobot / Leica GP 等设备快速浸入液乙烷，形成**玻璃态冰**
- 理想冰厚：略薄于颗粒直径；冰太厚信噪比差，太薄颗粒易变形

### 2.2 数据采集

- 300 kV 电镜（Titan Krios）+ 直接电子探测器（Gatan K3 / Falcon 4）
- **超分辨率模式** + 剂量分配（total dose 约 40–60 e⁻/Å²）
- 每张照片（movie）通常 30–50 帧，用于运动校正
- 目标：10,000–20,000 张照片（高分辨率结构通常需要百万级颗粒）

## 3. 预处理（Preprocessing）

### 3.1 运动校正（Motion Correction）

电子束引起的样品漂移会模糊图像，需要按帧校正：

\`\`\`bash
# RELION 4 中的 motion correction
relion_run_motioncorr \\
  --i Micrographs/ \\
  --o MotionCorr/ \\
  --j 8 \\
  --pipeline_control MotionCorr/job001/
\`\`\`

常用软件：**MotionCor2**（经典）、**RELION 内置**、**cryoSPARC Patch Motion Correction**。

### 3.2 CTF 估计

CTF（Contrast Transfer Function，衬度传递函数）描述了电子显微镜成像的相位翻转与欠焦效应：

\`\`\`bash
# CTFFIND4（经典命令行）
ctffind --input_micrograph MotionCorr/mic1.mrc \\
        --output_diag mic1_diag.mrc \\
        --defU 20000 --defV 20000

# Gctf（更快的 GPU 版本）
Gctf --input 'MotionCorr/*.mrc' --apix 1.0
\`\`\`

常用软件：**CTFFIND4**、**Gctf**、**cryoSPARC Patch CTF**。

**关键判据**：分辨率截断（Resolution limit）通常取 3–4 Å；剔除 astigmatism 过大或漂移严重的照片。

## 4. 颗粒挑选（Particle Picking）

### 4.1 传统方法

- **模板匹配**（Template-based）：用 2D 类平均做模板，相关性搜索
- 软件：RELION 的 auto-picking、cryoSPARC 的 Template Picker

### 4.2 深度学习方法（当前主流）

\`\`\`bash
# Topaz：基于 CNN 的颗粒挑选
topaz train --model-path model.pkl --epochs 20 \\
  --num-workers 4 --train-data particles.toml

topaz extract --model-path model.pkl \\
  --micrographs micrographs.toml \\
  --output-prefix picked --radius 100
\`\`\`

主流工具：**Topaz**、**crYOLO**（速度极快，支持实时挑选）、**cryoSPARC Blob Picker**。

> 现代流程常先用 Blob/模板快速挑一批做 2D 分类，再用干净类平均训练 Topaz/crYOLO 模型，显著提升复杂样品（膜蛋白、小蛋白）的挑选质量。

## 5. 2D 分类：筛选颗粒

2D 分类将颗粒按投影方向聚类，同时**剔除冰污染、错误挑选、变性颗粒**：

\`\`\`bash
# RELION 2D classification
relion_refine --o Class2D/ \\
  --particle_images particles.star \\
  --ctf --K 100 --iter 25 --tau2_fudge 4 \\
  --flatten_solvent --zero_mask --pool 30 \\
  --per_particle_ctf --j 8
\`\`\`

**质量判据**：
- 类平均应显示清晰的二级结构（α 螺旋、β 片层）
- 保留占总量 50–80% 的"好"颗粒
- 好的类平均是后续 3D 重构质量的前提

## 6. 3D 重构：初始模型与分类

### 6.1 初始模型（Ab-initio）

无需参考模型，从随机起始重构：

\`\`\`bash
# cryoSPARC
cryosparc_abinitio ... --num_classes 3
\`\`\`

- 常用：**cryoSPARC Ab-initio**（速度快、稳健）、RELION 3D initial model
- 生成 2–4 个类，选择结构特征最清晰的类作为参考

### 6.2 3D 分类：处理构象异质性

生物样品常存在多个构象（开放/闭合、配体结合/未结合）：

\`\`\`bash
# RELION 3D classification
relion_refine --o Class3D/ \\
  --ref initial.mrc --K 4 --iter 40 \\
  --tau2_fudge 8 --particle_diameter 200
\`\`\`

**策略**：先粗分（K=4–6）看构象分布，再对感兴趣构象细分（K=2–3），最终每个 3D 类对应一个构象。

## 7. 3D 精修（Refinement）

### 7.1 常规精修

\`\`\`bash
relion_refine --o Refine3D/ \\
  --ref class3d.mrc --particle_diameter 200 \\
  --auto_mask --solvent_mask \\
  --ctf --per_particle_ctf --j 8
\`\`\`

关键参数：
- **C1 对称性**起步，对称分子可指定 C/D/O 对称（如 C3、D7）
- **溶剂掩膜**（solvent mask）必须使用，否则 FSC 虚高
- 局部精修（local refinement）处理柔性区域

### 7.2 贝叶斯抛光（Bayesian Polishing）

对 movie 帧级颗粒做轨迹拟合与加权，提升高频信息：

\`\`\`bash
relion_polish --i Refine3D/run_data.star \\
  --o Polish/ \\
  --model_type 2 --nr_iter 5
\`\`\`

### 7.3 各向异性校正（CTF refinement）

精修每个颗粒的欠焦、像散、束倾斜：

\`\`\`bash
relion_ctf_refine --i Polish/particles.star \\
  --o CtfRefine/ \\
  --fit_defocus --fit_magnification --fit_tilt
\`\`\`

## 8. 分辨率评估：FSC 曲线

**FSC（Fourier Shell Correlation）**：把颗粒随机分成两个半集（half-maps），分别重构后比较两图的相关系数：

- **FSC = 0.143 准则**：金标准分辨率判定（半图校正后）
- **FSC = 0.5**：保守估计（未校正）

\`\`\`bash
relion_postprocess --i Refine3D/run_half1_class001.mrc \\
  --mask mask.mrc --autob_masking \\
  --angpix 1.0
\`\`\`

**质量检查表**：
| 指标 | 良好标准 |
|------|----------|
| 全局分辨率（FSC 0.143） | ≤ 3.5 Å（高分辨） |
| 颗粒数 | ≥ 10 万（3 Å 级通常需数十万） |
| 方向分布 | 均匀（无 Preferred orientation） |
| 密度连续性 | 侧链密度清晰可见（< 3 Å） |

## 9. 常用软件栈对比

| 环节 | RELION（命令行） | cryoSPARC（GUI） |
|------|------------------|------------------|
| 运动校正 | 内置 / MotionCor2 | Patch Motion |
| CTF | CTFFIND4 / Gctf | Patch CTF |
| 挑选 | 模板 / Topaz 集成 | Blob / Template / Topaz |
| 2D/3D | 稳定可靠 | 快、交互好 |
| 非均匀精修 | — | **NU-refinement**（柔性大分子利器） |

> **工作流建议**：cryoSPARC 做预处理与初始重构，RELION 做高精度精修与抛光；或全流程 cryoSPARC（对新手更友好）。

## 10. 小结

- 预处理（运动校正 + CTF）决定数据质量上限
- 深度学习颗粒挑选（Topaz/crYOLO）显著提升复杂样品表现
- 2D 分类筛颗粒，3D 分类处理构象异质性
- 精修三件套：溶剂掩膜 + 贝叶斯抛光 + CTF refinement
- FSC = 0.143 是分辨率金标准

后续文章将介绍结构可视化（PyMOL/ChimeraX）与原子建模（Coot/phenix），把密度图变成原子模型。
`;export{n as default};
