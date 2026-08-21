---
title: "生物大分子结构可视化：PyMOL 与 ChimeraX 实战"
date: "2026-08-04"
author: "zorrooz"
tags: ["结构可视化", "PyMOL", "ChimeraX", "教程"]
draft: false
description: "用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示"
---

# 生物大分子结构可视化：PyMOL 与 ChimeraX 实战

结构可视化是结构生物学的日常工具。本文以 PyMOL 与 UCSF ChimeraX 两大主流软件为例，从 PDB 数据获取到出版级渲染，覆盖日常科研所需的全部核心操作。

## 1. 数据来源：PDB 数据库

### 1.1 检索与下载

[wwPDB / RCSB PDB](https://www.rcsb.org/) 是全球结构数据仓库：

```bash
# 命令行下载（按 PDB ID）
wget https://files.rcsb.org/download/1CRN.pdb
wget https://files.rcsb.org/download/7A0A.pdb   # 冷冻电镜结构示例
```

### 1.2 PDB 文件的核心内容

```
HEADER    PLANT SEED PROTEIN             08-JAN-81   1CRN
TITLE     CRAMBIN
ATOM      1  N   THR A   1      11.106  16.700  17.083  1.00 20.14           N
ATOM      2  CA  THR A   1      11.579  17.437  17.166  1.00 20.55           C
...
HELIX    1   1 THR A    1  SER A    8  1                                  8
SHEET    1   A 2 ILE A   9  ARG A  12  0
```

- `ATOM` 行：原子坐标（残基、链、x/y/z、B 因子、元素）
- `HELIX` / `SHEET`：二级结构注释
- `CONECT`：共价连接
- 新格式 `mmCIF`（.cif）已成为主流，信息更完整

## 2. PyMOL 基础

### 2.1 启动与加载

```bash
pymol 1CRN.pdb          # 直接加载
# 或启动后
# File > Open
```

### 2.2 基础显示命令

```python
# 载入
fetch 1crn               # 从网络获取
load 1CRN.pdb

# 显示模式
show cartoon             # 卡通（主链走向）
show sticks              # 棍状（原子细节）
show surface             # 表面
hide lines               # 隐藏线条

# 着色
color cyan               # 整体着色
spectrum count           # 彩虹色（按残基序号）
color red, resi 1-10     # 指定区域着色
```

### 2.3 选择与操作

```python
# 选择语法
select helix, ss h                   # 所有螺旋
select sheet, ss s                   # 所有片层
select active_site, resi 15+20+45    # 指定残基
select ligand, resn HEM              # 配体（血红素）

# 对象操作
show sticks, active_site
color yellow, active_site
zoom active_site
```

### 2.4 测量

```python
# 距离
distance d1, /1CRN//A/15/CA, /1CRN//A/20/CA

# 角度 / 二面角
angle a1, /1CRN//A/15/CA, /1CRN//A/16/CA, /1CRN//A/17/CA
dihedral dh1, /1CRN//A/15/CA, /1CRN//A/16/CA, /1CRN//A/17/CA, /1CRN//A/18/CA
```

### 2.5 相互作用分析

```python
# 氢键
distance hbond, (resn HEM), (resn HIS)

# 接触面
select interface, chain A within 4.5 of chain B
```

## 3. PyMOL 脚本批处理

把常用操作写成脚本，可重复执行：

```python
# render_view.py
fetch 1crn
hide everything
show cartoon
color spectrum, ss
set cartoon_transparency, 0.2

# 放大结合位点并渲染
zoom resi 1-20
set ray_shadows, 1
set ray_opaque_background, 0
ray 1200, 900
png 1crn_view.png, dpi=300
```

```bash
pymol -cq render_view.py    # -c 命令行模式，-q 静默
```

## 4. UCSF ChimeraX 基础

ChimeraX 界面友好、现代（Qt GUI），且支持 `open` 命令批量操作：

### 4.1 打开结构

```bash
chimerax 1CRN.pdb
# 或命令
open 1crn
open 7a0a
```

### 4.2 常用命令

```python
# 显示模式
cartoon
stick
surface

# 配色
color bychain          # 按链着色
color byhetero         # 配体/离子区分
color byattribute bfactor   # 按 B 因子

# 选择
select :15-20          # 残基 15-20
select /A               # A 链
select ligand
```

### 4.3 冷冻电镜密度图展示（ChimeraX 强项）

```python
# 打开密度图（.mrc）
open map.mrc

# 调整等值面水平
volume level 0.01
volume #2 level 0.05

# 密度图与模型叠加
open model.pdb
open map.mrc
volume #2 level 0.02
color #1 cornflowerblue
transparency #2 30

# 区域密度查看（局部）
volume zone #2 near :45-60 radius 5
```

**ChimeraX 的 `volume zone` 是检查局部密度质量的核心工具**：在残基周围显示密度，判断侧链是否可辨。

## 5. 高质量渲染要点

### 5.1 色彩原则

- 卡通着色：按链（bychain）或按结构域
- 关键残基：统一高亮色（黄/红/蓝），避免彩虹滥用
- 配体：`byhetero` 自动区分元素色（C 灰、N 蓝、O 红、S 黄）

### 5.2 光照与材质

```python
# PyMOL
set ray_shadows, 1
set specular, 0.5
set ambient, 0.3

# ChimeraX
set lightMode full
graphics silhouettes true
```

### 5.3 分辨率输出

```bash
# PyMOL 出版级渲染
ray 2400, 1800
png figure.png, dpi=600

# ChimeraX
save figure.png width 2400 height 1800 supersample 3
```

## 6. 结构比对与叠合

```python
# PyMOL：对齐两个结构
align model2, model1

# ChimeraX：matchmaker
matchmaker #2 to #1
```

**RMSD** 是衡量结构相似度的标准指标（通常报告 Cα 或全原子 RMSD）。

## 7. 交互式查看与分享

- **Mol* / NGL**：网页端 3D 查看（PDBe、RCSB 内置）
- **PyMOL Web**：把会话导出为 HTML 交互页面
- 科研交流常用：PyMOL session（.pse）、ChimeraX session（.cxs）

## 8. 小结

- PDB 是数据源头，mmCIF 是新标准
- PyMOL：命令强大、脚本化批处理（`pymol -cq script.py`）
- ChimeraX：现代 GUI + 密度图处理（`volume`、`volume zone`）无可替代
- 出版级渲染：隐藏杂项 → 突出关键 → 高质量光线 → 高分辨率导出

下一篇将介绍原子建模：如何从密度图搭建并精修原子模型。
