const r=`---\r
title: "Python 机器学习入门"\r
date: "2025-09-19"\r
author: "zorrooz"\r
tags: ["Python", "机器学习", "scikit-learn", "分类", "回归"]\r
draft: false\r
description: "使用 scikit-learn 进行生物数据机器学习分析"\r
---\r
\r
# Python 机器学习入门\r
\r
## 数据预处理\r
\r
\`\`\`python\r
from sklearn.preprocessing import StandardScaler\r
from sklearn.model_selection import train_test_split\r
\r
# 数据标准化\r
scaler = StandardScaler()\r
X_scaled = scaler.fit_transform(X)\r
\r
# 数据集划分\r
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)\r
\`\`\`\r
\r
## 常用算法\r
\r
- 逻辑回归：\`LogisticRegression\`\r
- 随机森林：\`RandomForestClassifier\`\r
- 支持向量机：\`SVC\``;export{r as default};
