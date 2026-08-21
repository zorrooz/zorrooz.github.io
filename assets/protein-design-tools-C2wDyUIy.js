const n=`---
title: 蛋白质设计主流工具
date: 2026-08-04
author: zorrooz
tags:
  - 蛋白质设计
  - Rosetta
  - AlphaFold
  - RFdiffusion
  - 计算生物学
draft: false
description: 蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流
---

# 蛋白质设计主流工具

蛋白质设计（Protein Design）是计算生物学的重要方向：给定目标功能，设计出能折叠成特定结构的氨基酸序列。本文梳理从经典能量函数到生成式 AI 的主流工具谱系。

## 1. 设计问题的两种范式

### 1.1 序列设计（Sequence Design）

给定目标结构（骨架），寻找最稳定的氨基酸序列：

\`\`\`
结构（backbone）→ 序列设计 → 氨基酸序列 → 实验验证
\`\`\`

### 1.2 结构设计（De novo Design）

从零设计新结构/新功能：

\`\`\`
功能需求 → 结构生成 → 序列设计 → 实验验证
\`\`\`

## 2. Rosetta：能量函数范式

Rosetta 是结构预测与设计领域的"老牌帝国"，基于物理与知识混合的能量函数。

### 2.1 核心能量函数

- **REF2015**：经典全原子能量函数（范德华、氢键、溶剂化、静电、二面角偏好）
- 打分 = 各项能量加权和，越低越好
- Rosetta 的"设计"本质是**能量最小化**：在序列/结构空间搜索低能态

### 2.2 常用协议

\`\`\`bash
# 序列设计（FixBB：固定骨架设计）
rosetta_scripts.default.linuxgccrelease \\
  -s input.pdb \\
  -parser:protocol fixbb.xml \\
  -out:prefix designed_

# fixbb.xml 核心片段
# <TaskOperations>
#   <RestrictToRepacking name="repack"/>
#   <OperateOnResidueSubset name="no_design_core" selector="..." />
# </TaskOperations>
# <RosettaScripts>
#   <PackRotamersMover name="pack" task_operations="repack"/>
# </RosettaScripts>
\`\`\`

### 2.3 Rosetta 的经典应用

- **酶设计**：催化活性位点支架设计（如 Kemp 消除酶）
- **蛋白-蛋白界面设计**：结合特异性改造
- **binder 设计**：高亲和力结合蛋白
- **热稳定性改造**：计算突变提高 Tm

### 2.4 优点与局限

| 优点 | 局限 |
|------|------|
| 物理可解释、可精细调控 | 计算量大、需要经验调参 |
| 支持任意骨架与修饰 | 对柔性区域预测有限 |
| 生态成熟（文档、教程多） | 学习曲线陡峭 |

## 3. AlphaFold：预测即设计辅助

### 3.1 AlphaFold2（2021）

- 端到端深度学习：序列 → 结构（MSA + 注意力机制）
- 精度：CASP14 原子级精度（GDT > 90）
- 对设计的价值：**快速验证设计序列的可折叠性**

### 3.2 AlphaFold3（2024）

- 统一框架预测蛋白质、核酸、小分子、离子复合物
- 引入扩散模型与配对表示
- 意义：设计靶标（binder 与配体的复合物结构）可直接预测

\`\`\`bash
# ColabFold：AlphaFold2 的轻量实现（GPU 友好）
colabfold_batch input.fasta out_dir \\
  --num-recycle 3 --model-type alphafold2_multimer_v3
\`\`\`

### 3.3 关键用法：验证设计

\`\`\`python
# 设计序列的可信度指标
# pLDDT > 90：高置信
# PAE < 5：结构域间相对位置可信
# 设计验证流程：序列 → AF2 预测 → 与目标结构比对（RMSD）
\`\`\`

## 4. 生成式设计：RFdiffusion 与 ProteinMPNN

### 4.1 RFdiffusion（2023，Baker Lab）

基于**去噪扩散模型**的结构生成器：

- 输入：功能基序（motif）、形状约束、结合靶点
- 输出：全新蛋白质骨架（backbone）
- 应用：binder 设计、对称寡聚体、酶活性位点支架

\`\`\`bash
# RFdiffusion 示例：设计结合蛋白
run_inference.py \\
  scaffoldguided.scaffoldguided=True \\
  scaff_loader.contigs_map.json=contigs.json \\
  inference.output_prefix=outputs/design \\
  inference.num_designs=100 \\
  potentials.guidance_scale=3.0
\`\`\`

### 4.2 ProteinMPNN（2022）

**序列设计神经网络**：给定骨架 → 序列，与 RFdiffusion 形成"生成-编码"闭环：

\`\`\`bash
# ProteinMPNN 序列设计
python protein_mpnn_run.py \\
  --pdbpath scaffold.pdb \\
  --out_folder mpnn_out \\
  --num_seq_per_target 8 \\
  --batch_size 4
\`\`\`

**标准工作流（Baker Lab 范式）**：

\`\`\`
① RFdiffusion：生成骨架
   ↓
② ProteinMPNN：骨架 → 多组候选序列（采样温度 0.1-0.3）
   ↓
③ AF2 反向验证：序列 → 预测结构 → 与骨架 RMSD 比对
   ↓
④ 筛选高 pLDDT / 低 RMSD 候选
   ↓
⑤ 实验表达验证（酵母/噬菌体展示、ELISA）
\`\`\`

### 4.3 其他生成式工具

- **Chroma**：另一个扩散模型（生成 + 序列设计一体）
- **Genie**：更快的蛋白质扩散模型
- **ESMFold**（Meta，2022）：无 MSA 的极速预测（每序列 < 1 秒），适合大规模设计筛选
- **ESM3**：多模态生成（序列+结构+功能）

## 5. binder 设计专题

Binder（结合蛋白）设计是当前最活跃的应用场景：

### 5.1 经典流程（结合现有工具）

\`\`\`bash
# 1. 靶点表面热点分析（Rosetta）
# 2. 热区引导的 binder 骨架生成（RFdiffusion / 或 hotspot 约束）
# 3. ProteinMPNN 序列设计
# 4. AF2 复合物预测验证（binder + target 对接）
# 5. 实验筛选（酵母展示、FACS 分选）
\`\`\`

### 5.2 关键考量

- 结合热点残基（hotspot）的确定：实验突变扫描（alanine scan）或计算（Rosetta ddG）
- 结合亲和力预测：\`FoldX\`、\`Rosetta interface_ddg\`、AF2 PAE
- 多轮"设计-验证"迭代：实验反馈回到计算

## 6. 工具选择速查表

| 任务 | 首选工具 | 备选 |
|------|----------|------|
| 结构预测 | AlphaFold2/3（ColabFold） | ESMFold（速度） |
| 骨架生成 | RFdiffusion | Chroma |
| 序列设计 | ProteinMPNN | Rosetta FixBB |
| 热稳定性设计 | Rosetta（ddg） | FoldX |
| 结合亲和力评估 | Rosetta interface_ddg | FoldX / AF2 PAE |
| 界面设计 | Rosetta 协议 | RFdiffusion + MPNN |
| 大规模筛选 | ESMFold + ProteinMPNN | — |

## 7. 实验验证闭环

计算设计必须回到实验：

\`\`\`
表达（E. coli / 酵母）→ 纯化 → 表征
  → SEC / DSC（稳定性）
  → 表面等离子共振 SPR / BLI（亲和力）
  → 结构验证（晶体 / 冷冻电镜 / AlphaFold 对照）
\`\`\`

**设计成功率参考**：文献报道 AF2 引导的 RFdiffusion binder 设计，实验阳性率可达 10–20%（远超传统方法），但每个成功案例背后都有多轮迭代。

## 8. 小结

- **经典范式**：Rosetta 能量函数（可解释、可调控）
- **预测范式**：AlphaFold 系列（验证、复合物预测）
- **生成范式**：RFdiffusion + ProteinMPNN + AF2 验证闭环（当前主流）
- 设计成功的衡量标准始终是**实验验证**

下一篇将介绍虚拟筛选的主流工具：把小分子与靶点对接的完整工具链。
`;export{n as default};
