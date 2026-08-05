---
title: "虚拟筛选主流工具"
date: "2026-08-04"
author: "zorrooz"
tags: ["虚拟筛选", "分子对接", "AutoDock Vina", "计算生物学"]
draft: false
description: "虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略"
---

# 虚拟筛选主流工具

虚拟筛选（Virtual Screening, VS）通过计算在百万级化合物库中寻找潜在活性分子，是药物发现早期阶段的核心手段。本文介绍主流对接软件、完整流程与注意事项。

## 1. 虚拟筛选的基本逻辑

```
受体结构（蛋白/核酸）
   ↓
① 受体准备（加氢、质子化、网格）
   ↓
② 化合物库（ZINC / Enamine / 自建库）
   ↓
③ 配体准备（3D 构象、质子化、力场参数）
   ↓
④ 分子对接（构象采样 + 打分）
   ↓
⑤ 排序筛选（Top 100-1000）
   ↓
⑥ 实验验证（IC50 / 活性测定）
```

## 2. 主流对接软件

### 2.1 AutoDock Vina（开源首选）

最流行的开源对接软件，速度与精度平衡良好：

```bash
# 1. 受体与配体准备（用 MGLTools 或脚本）
#    受体：pdb → pdbqt（加氢、合并非极性氢）
#    配体：sdf/mol2 → pdbqt

# 2. 配置文件（conf.txt）
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

# 3. 运行
vina --config conf.txt --out results.pdbqt --log log.txt

# 4. 解读：affinity（kcal/mol）越低越好
#  mode |   affinity | dist from best mode
# -----+------------+----------------------
#    1 |   -9.2     | 0.000
#    2 |   -8.7     | 2.315
```

### 2.2 Glide（Schrödinger，商业）

制药行业标准：

```bash
# 准备工作流（Maestro GUI 或命令行）
glide -overwrite -adjust \
  -receptor receptor.maegz \
  -ligand ligands.maegz \
  -p glide_docking.inp \
  -NJOBS 16

# 三种精度模式
# HTVS（高通量，快）→ SP（标准精度）→ XP（高精度，慢）
```

**优势**：精确的配体柔性处理、完善的打分校准、大型药企生态支持。

### 2.3 DOCK 6（开源）

经典几何匹配算法（shape matching）：

```bash
dock6 -i dock.in -o dock.out
# 球集（sphere）生成 → 取向搜索 → 打分
# 擅长：结合口袋形状匹配、大库初筛
```

### 2.4 其他常用软件

| 软件 | 特点 |
|------|------|
| **AutoDock4** | 经典遗传算法，精度高速度慢 |
| **rDock** | 开源、适合大库、易脚本化 |
| **GOLD** | 遗传算法，CCDC 出品 |
| **LeDock** | 快速、易用（学术免费） |
| **Plants** | 基于 Ant Colony 优化 |
| **smina** | Vina 的增强版（支持自定义打分） |

## 3. 受体准备（关键第一步）

### 3.1 结构来源

- 晶体/冷冻电镜结构（PDB）
- AlphaFold 预测结构（无实验结构时）
- 注意：结合口袋附近的柔性（侧链）处理

### 3.2 准备要点

```bash
# 常用工具
# Schrödinger Protein Prep Wizard（商业）
# ADFRsuite / MGLTools prepare_receptor4.py（免费）

python prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt

# 关键步骤
# 1. 去除水分子（保守水除外）
# 2. 加氢、确定质子化状态（pH 7.4 附近）
# 3. 分配键序与电荷（Amber/CHARMM 力场）
# 4. 定义结合口袋（已知活性位点 or 空腔检测：fpocket / DoGSiteScorer）
```

## 4. 配体库与准备

### 4.1 化合物库来源

| 库 | 规模 | 特点 |
|----|------|------|
| ZINC20 | 数十亿 | 免费、可下载、带 3D 构象 |
| Enamine REAL | > 400 亿 | 可合成性高、商购 |
| ChEMBL | 活性数据 | 已知活性化合物 |
| 自建库 | 定制 | 衍生物设计 |

### 4.2 配体准备

```bash
# 3D 构象生成
obabel ligand.smi -O ligand.sdf --gen3d -p 7.4

# 多构象（柔性重要，尤其环状体系）
obabel ligand.sdf -O conf.sdf --conformer --nconf 50

# 质子化与电荷
# Schrödinger LigPrep / Open Babel / RDKit

# RDKit 准备（Python）
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
```

## 5. 打分函数（Scoring Functions）

### 5.1 三类打分函数

| 类型 | 代表 | 原理 |
|------|------|------|
| 力场型 | DOCK、AutoDock4 | 范德华+静电+内能 |
| 经验型 | Glide、Vina | 统计拟合实验数据 |
| 知识型 | GOLD/ASP | 原子对距离分布统计 |

### 5.2 打分函数的局限（必须知道）

- **不精确**：相关系数通常 0.3–0.6，只能排序不能给绝对亲和力
- 熵效应、溶剂效应、诱导拟合难以准确描述
- **建议**：对接分数 + 多工具交叉验证 + 分子动力学（MM/PBSA、FEP）精算

## 6. 完整高通量筛选流程（HTS 实操）

```bash
# 1. 库准备：SDF 拆分 + 转 pdbqt
mkdir -p ligands_pdbqt
obabel library.sdf -O ligands_pdbqt/lig_.pdbqt -m -p 7.4

# 2. 批量对接（smina/Vina 多文件模式）
smina -r receptor.pdbqt -l library.sdf \
  --center_x 10.5 --center_y 20.3 --center_z -5.8 \
  --size_x 22 --size_y 22 --size_z 22 \
  --num_modes 1 --cpu 32 \
  -o docked.sdf --log scores.txt

# 3. 排序（按 affinity）
sort -k2 -n scores.txt | head -100
```

### 6.1 多阶段漏斗策略

```
粗筛（HTVS / Vina，百万级）→ 前 5-10%
  ↓
标准精度（SP / 多构象）→ 前 10-20%
  ↓
高精度（XP / 诱导拟合）→ Top 100-1000
  ↓
共识打分（多个软件交叉）→ Top 50-200
  ↓
视觉检查（结合模式合理性）→ Top 10-50
  ↓
实验测试
```

### 6.2 共识打分（Consensus Scoring）

不同软件打分函数互补，取交集/排名平均可显著提升富集率：

```bash
# 示例：Vina + Glide + LeDock 三软件 Top 100 取交集
# 或按排名归一化后取平均
```

## 7. 验证与高级方法

### 7.1 基准测试（Benchmark）

- **DUD-E**：标准数据集（含活性物与诱饵 decoys），评估软件区分能力
- 指标：AUC、富集因子（EF1%、EF5%）

### 7.2 分子动力学后处理

```bash
# 对接结果 → MD 模拟（GROMACS/AMBER）→ MM/PBSA 结合自由能
gmx mdrun -deffnm complex
gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ct traj.xtc ...
```

### 7.3 结构基与配体基结合

- **基于结构**（SBDD）：对接 + MD + FEP
- **基于配体**（LBDD）：药效团、QSAR、相似性搜索
- 组合策略是现代虚拟筛选主流

## 8. 常见陷阱

| 陷阱 | 后果 | 规避 |
|------|------|------|
| 口袋定义错误 | 全库错位 | 用已知活性位点或实验信息 |
| 受体刚性 | 错过诱导拟合 | 柔性侧链 / 集合对接（ensemble docking） |
| 质子化状态错误 | 静电打分失真 | pH 相关准备（pH 7.4） |
| 打分函数过信 | 假阳性率高 | 共识打分 + 视觉检查 |
| 忽略可合成性 | 筛出无法合成 | 过滤 PAINS、类药性规则（Lipinski/Veber） |

## 9. 小结

- 三梯队：**Vina**（开源快速）→ **Glide**（商业标准）→ **DOCK/rDock**（大库）
- 流程：受体准备 → 配体库 → 对接 → 多阶段漏斗 → 实验验证
- 打分函数只能排序；共识打分 + MD 精算提升可信度
- 虚拟筛选的终点永远是实验

至此计算生物学方向完成：蛋白质设计工具 + 虚拟筛选工具，两条主线构成"设计-筛选"的计算药物发现闭环。
