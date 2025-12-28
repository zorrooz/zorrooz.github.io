const r=`---\r
title: "实验方法 - 蛋白质结构解析"\r
date: "2025-07-01"\r
author: "zorrooz"\r
tags: ["Methods", "Cryo-EM", "Sample Preparation", "Data Collection"]\r
draft: false\r
description: "蛋白质结构解析实验方法演示"\r
---\r
\r
# 实验方法\r
\r
**以下内容皆为虚构！**\r
\r
## 1. 样品制备\r
\r
### 蛋白质表达与纯化\r
\`\`\`bash\r
# 表达系统\r
宿主：E. coli BL21(DE3)\r
诱导条件：0.5 mM IPTG, 18°C, 16小时\r
\r
# 纯化步骤\r
1. Ni-NTA 亲和层析\r
2. 凝胶过滤层析 (Superdex 200)\r
3. 浓缩至 3 mg/mL\r
\`\`\`\r
\r
### 负染电镜筛选\r
- 染色剂：2% 醋酸铀酰\r
- 数据采集：200 kV, 50,000x 放大倍数\r
- 筛选标准：单分散、完整颗粒\r
\r
## 2. 冷冻电镜数据采集\r
\r
### 样品制备\r
- 3 μL 样品滴加在 Quantifoil R1.2/1.3 网格\r
- 滤纸吸干，液氮快速冷冻\r
- 工作温度：-180°C\r
\r
### 采集参数\r
| 参数 | 值 |\r
|------|-----|\r
| 加速电压 | 300 kV |\r
| 像素尺寸 | 1.5 Å |\r
| 电子剂量 | 50 e⁻/Å² |\r
| 帧数 | 50 frames |\r
\r
## 3. 数据处理\r
\r
### 运动矫正\r
- 软件：MotionCor2\r
- 修正帧间运动和漂移\r
\r
### 颗粒挑选\r
- 方法：模板匹配\r
- 软件：RELION 3.1\r
- 颗粒数量：~50,000\r
\r
### 二维分类\r
- 类别数：20\r
- 选择标准：清晰特征\r
\r
## 4. 三维重建\r
\r
### 初始模型\r
- 方法：随机圆锥重构\r
- 分辨率：~8 Å\r
\r
### 精细重建\r
- 迭代：20 轮\r
- 分辨率：3.0 Å\r
- 分辨率标准：FSC 0.143\r
\r
## 5. 模型构建与优化\r
\r
### 初始模型搭建\r
- 软件：Coot\r
- 方法：自动建模 + 手动修正\r
\r
### 结构优化\r
- 软件：Phenix.refine\r
- 优化内容：\r
  - 原子坐标\r
  - 温度因子\r
  - 占有率\r
\r
### 质量验证\r
- Ramachandran plot：96.5% 优势区域\r
- 几何偏差：键长 0.01 Å，键角 1.2°\r
- MolProbity score：1.8`;export{r as default};
