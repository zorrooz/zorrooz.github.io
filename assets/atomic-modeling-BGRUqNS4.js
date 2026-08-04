const n=`---
title: "原子建模与精修"
date: "2026-08-04"
author: "zorrooz"
tags: ["原子建模", "Coot", "phenix", "结构精修", "教程"]
draft: false
description: "冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标"
---

# 原子建模与精修

拿到高质量密度图后，下一步是搭建原子模型并精修，最终获得可发布、可分析的坐标文件。本文覆盖从自动建模到质量评估的完整工作流。

## 1. 工作流总览

\`\`\`
密度图（.mrc）
   ↓
① 初始模型（自动建模 / 同源模型 / AF 预测）
   ↓
② 模型放入密度（docking / rigid body）
   ↓
③ 迭代精修：Coot 手动修正 ↔ phenix 自动精修
   ↓
④ 质量评估（MolProbity / EMRinger / Q-score）
   ↓
⑤ 验证与发布（PDB deposition）
\`\`\`

## 2. 初始模型来源

### 2.1 自动建模（深度学习时代的主流）

- **ModelAngelo**：基于 AlphaFold 架构的自动建模，直接从密度图搭模型
- **DeepTracer**：快速全自动建模（GPU）
- **phenix.map_to_model**：经典自动建主链/侧链

\`\`\`bash
# ModelAngelo 使用示例
model_angelo --halfmaps half1.mrc half2.mrc \\
  --output-dir model_angelo_out \\
  --num-workers 8

# phenix.map_to_model
phenix.map_to_model map.mrc resolution=3.0 \\
  output_model=auto_model.pdb
\`\`\`

### 2.2 同源建模与分子替换

- 有同源结构：用 **phenix.dock_in_map** 放入密度
- AlphaFold 预测模型：直接作为初始模型（cryo-EM 中常用）
- 低分辨率（> 5 Å）时先放刚性体，再逐步细化

## 3. Coot：手动建模与修正

Coot 是手动模型修正的行业标准（GUI 工具）。

### 3.1 常用操作

\`\`\`
File > Open Coordinates...   # 载入模型
File > Open Map...           # 载入密度图（自动算 FSC 相关）

# 基本修正
1. 主链调整：拖拽 Cα（Real Space Refine 区域）
2. 侧链方向：Rotamers（旋转异构体）
3. 缺失残基：Add Terminal Residue / Build
4. 删除多余：Delete Residue Range
5. 突变残基：Mutate & Autofit
6. 配体搭建：Ligand Builder / restraints
\`\`\`

### 3.2 关键理念：Real Space Refine

\`\`\`
Calculate > Real Space Refine Zone
\`\`\`

- 选中密度区域，Coot 优化该区域几何与密度拟合
- 重复：修正 → 精修 → 检查 → 再修正
- 每轮精修后回到 Coot 检查不良区域（Ramachandran 离群点、clash）

### 3.3 配体建模

\`\`\`bash
# 生成配体 restraints
phenix.elbow ligand.smiles    # 从 SMILES 生成化学约束
# 或
phenix.ligand_idealization ligand.cif
\`\`\`

Coot 中 \`Ligand > Ligand Builder\` 手动搭建，再 Real Space Refine。

## 4. phenix 自动精修

### 4.1 基本精修

\`\`\`bash
phenix.real_space_refine model.pdb \\
  map.mrc \\
  resolution=3.0 \\
  output.prefix=refined \\
  nproc=8
\`\`\`

### 4.2 精修策略

\`\`\`bash
# 多轮迭代：先整体后局部
phenix.real_space_refine model.pdb map.mrc \\
  resolution=3.0 \\
  strategy=individual_sites+individual_adp \\
  restraints=ligand.cif \\
  nproc=8
\`\`\`

常用 strategy 组合：
- \`rigid_body\`：低分辨率先做刚性体
- \`individual_sites\`：逐原子坐标精修（高分辨率）
- \`individual_adp\`：B 因子精修
- \`torsion\`：扭转角精修（蛋白质常用）
- \`local_grid\`：局部网格搜索（柔性区域）

### 4.3 自动修正循环（NCS 与水分子）

\`\`\`bash
# 找水分子（高分辨时）
phenix.find_peaks_holes map.mrc map.mrc resolution=3.0

# NCS 约束（有对称亚基时显著提升质量）
phenix.real_space_refine ... strategy=individual_sites+individual_adp ncs_search.enabled=True
\`\`\`

## 5. 模型质量评估

### 5.1 几何质量：MolProbity

\`\`\`bash
phenix.molprobity refined.pdb
\`\`\`

核心指标：

| 指标 | 良好标准（3 Å） | 含义 |
|------|----------------|------|
| Ramachandran  favored | > 95% | 主链二面角合理性 |
| Ramachandran outliers | < 0.5% | 离群主链 |
| Rotamer outliers | < 1% | 侧链构象异常 |
| Clashscore | < 5 | 原子碰撞 |
| MolProbity score | < 2 | 综合几何分数 |

### 5.2 密度拟合质量

- **EMRinger**：评估侧链与密度的拟合（> 2 表示良好）
- **Q-score**：评估局部密度质量（1.0 = 完美，0.5 = 可接受）
- **CC_mask**：模型-密度相关系数（> 0.8 良好）
- **map-model FSC**：模型与两张半图的 FSC 对比

\`\`\`bash
phenix.map_model_cc refined.pdb map.mrc
phenix.em_ringer refined.pdb map.mrc
\`\`\`

### 5.3 过拟合检测

**Free component（CC-free）**：与精修中使用的半图独立计算的相关性。精修时必须监控：

\`\`\`
精修前 CC_work ≈ CC_free
精修后 CC_work 明显 > CC_free → 过拟合信号
\`\`\`

**规范做法**：精修只用一张半图（或全图 + 独立半图验证），用另一张半图计算 CC_free。

## 6. 完整精修流程示例

\`\`\`bash
# 1. 自动建模
model_angelo --halfmaps half1.mrc half2.mrc --output-dir auto

# 2. 初始精修
phenix.real_space_refine auto/model.pdb map.mrc \\
  resolution=3.1 strategy=rigid_body+individual_sites \\
  output.prefix=round1

# 3. Coot 手动修正（交互）
coot round1_refined.pdb map.mrc

# 4. 第二轮精修（加 B 因子）
phenix.real_space_refine coot_out.pdb map.mrc \\
  resolution=3.1 strategy=individual_sites+individual_adp \\
  output.prefix=round2

# 5. 质量评估
phenix.molprobity round2_refined.pdb
phenix.map_model_cc round2_refined.pdb map.mrc
phenix.em_ringer round2_refined.pdb map.mrc

# 6. 迭代直到所有指标达标
\`\`\`

## 7. 提交 PDB

发布前需准备：

1. 坐标文件（.pdb / .cif）
2. 密度图（half-maps + full map，.mrc）
3. 验证报告（MolProbity、EMRinger 等）
4. FSC 曲线与 mask 信息

通过 [wwPDB deposition](https://deposit.wwpdb.org/) 提交，获取 PDB ID 与 EMDB ID。

## 8. 常见问题

| 问题 | 原因 | 解决 |
|------|------|------|
| 侧链密度模糊 | 分辨率不足 / B 因子错误 | 提高分辨率、检查方向分布 |
| Ramachandran 离群多 | 主链建错 | Coot 中按密度重建该区域 |
| 过拟合 | 精修过度 | 用 CC_free 监控、减少迭代 |
| 配体拟合差 | restraints 不完整 | phenix.elbow 重新生成 |
| 局部区域质量差 | 柔性/构象异质性 | 3D 分类分离构象、局部精修 |

## 9. 小结

- 建模三来源：自动建模（ModelAngelo/DeepTracer）、同源/AF 预测、从头搭建
- Coot 手动修正与 phenix 自动精修**交替迭代**是标准节奏
- 质量三重检查：几何（MolProbity）+ 拟合（EMRinger/Q-score）+ 过拟合（CC_free）
- 亚 3 Å 应追求侧链级可信度，配体与修饰需单独验证

至此结构生物学方向完成：数据处理流程 → 技术综述 → 可视化 → 原子建模，构成从实验到模型的完整闭环。
`;export{n as default};
