---
title: "Python 机器学习入门"
date: "2025-09-19"
tags: ["Python", "机器学习", "scikit-learn", "分类", "回归"]
draft: false
description: "使用 scikit-learn 进行生物数据机器学习分析"
---

# Python 机器学习入门

## 数据预处理

```python
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# 数据标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 数据集划分
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
```

## 常用算法

- 逻辑回归：`LogisticRegression`
- 随机森林：`RandomForestClassifier`
- 支持向量机：`SVC`