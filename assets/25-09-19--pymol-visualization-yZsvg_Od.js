const r=`---\r
title: "PyMOL 分子结构可视化"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["PyMOL", "分子可视化", "蛋白质结构", "3D结构", "结构生物学"]\r
draft: false\r
description: "PyMOL 在蛋白质三维结构可视化和分析中的应用"\r
---\r
\r
# PyMOL 分子结构可视化\r
\r
## 基本操作\r
\r
\`\`\`python\r
# 加载蛋白质结构\r
load protein.pdb\r
\r
# 显示卡通图\r
show cartoon\r
\r
# 着色\r
color red, chain A\r
color blue, chain B\r
\`\`\`\r
\r
## 常用功能\r
\r
- 结构比对\r
- 表面电荷显示\r
- 活性位点标记`;export{r as default};
