const r=`---\r
title: "BioCrawler 安装指南"\r
date: "2025-08-01"\r
author: "zorrooz"\r
tags: ["Installation", "Python", "Setup"]\r
draft: false\r
description: "BioCrawler 安装和配置演示"\r
---\r
\r
# BioCrawler 安装指南\r
\r
**以下内容皆为虚构！**\r
\r
## 环境要求\r
\r
- Python 3.8+\r
- pip 20.0+\r
- 操作系统：Windows 10+, macOS 10.15+, Linux\r
\r
\r
## 快速安装\r
\r
### 1. 使用 pip 安装\r
\r
\`\`\`bash\r
pip install biocrawler\r
\`\`\`\r
\r
### 2. 从源码安装\r
\r
\`\`\`bash\r
git clone https://github.com/example/biocrawler.git\r
cd biocrawler\r
pip install -e .\r
\`\`\`\r
\r
## 配置设置\r
\r
### 创建配置文件\r
\r
\`\`\`bash\r
# 初始化配置\r
biocrawler init\r
\r
# 编辑配置文件\r
nano ~/.biocrawler/config.yaml\r
\`\`\`\r
\r
### 配置示例\r
\r
\`\`\`yaml\r
api:\r
  key: "your_api_key_here"\r
  timeout: 30\r
\r
database:\r
  cache_dir: "~/.biocrawler/cache"\r
  max_cache_size: 10GB\r
\r
download:\r
  threads: 4\r
  retries: 3\r
\`\`\`\r
\r
## 验证安装\r
\r
\`\`\`bash\r
# 测试连接\r
biocrawler test\r
\r
# 查看版本\r
biocrawler --version\r
\r
# 显示帮助\r
biocrawler --help\r
\`\`\`\r
\r
## 常见问题\r
\r
**Q: 安装失败怎么办？**\r
A: 确保 Python 版本 >= 3.8，并尝试更新 pip\r
\r
**Q: 如何更新？**\r
A: \`pip install --upgrade biocrawler\`\r
\r
## Docker 安装\r
\r
\`\`\`bash\r
docker pull biocrawler/biocrawler:latest\r
docker run -it biocrawler/biocrawler:latest`;export{r as default};
