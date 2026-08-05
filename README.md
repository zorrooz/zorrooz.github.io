# gblog

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/zorrooz/zorrooz.github.io)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)
[![GitHub stars](https://img.shields.io/github/stars/zorrooz/gblog?style=social)](https://github.com/zorrooz/zorrooz.github.io)
[![English Version](https://img.shields.io/badge/English-Version-blue)](README_EN.md)


[English](README_EN.md) | [中文](README.md)

gblog是基于 Vue 3 + Vite 构建的现代化个人博客系统，支持中英文双语、主题切换、Markdown 文章管理等功能。

## ✨ 项目介绍

gblog 是一个功能完善的静态博客系统，专为个人知识管理和内容展示设计。项目采用 Vue 3 组合式 API 开发，使用 Vite 作为构建工具，支持响应式设计和现代化的用户体验。

### 🎯 核心特性

- 双语支持：内置中文和英文双语切换功能
- 主题切换：支持亮色/暗色/跟随系统三种模式
- Markdown 渲染：完整的 Markdown 支持，包含代码高亮和数学公式
- 内容复制：文章与表格一键复制（表格为 TSV 格式，可直接粘贴到 Excel）
- 内容组织：支持分类、资源、关于页等多维度内容管理
- 自动翻译：集成 AI 翻译工具，支持内容自动中英翻译
- 响应式设计：适配各种屏幕尺寸
- 快速加载：基于 Vite 的快速构建和热更新

### 🛠️ 技术栈

- 前端框架：Vue 3 (Composition API)
- 构建工具：Vite 7
- 状态管理：Pinia
- 路由管理：Vue Router 4
- UI 框架：Bootstrap 5
- 主题实现：Sass
- Markdown 处理：unified生态（remark和rehype）
- 国际化：vue-i18n
- 部署：GitHub Pages

### 📁 项目结构

内容数据独立在 `data` 分支（本地 worktree 挂到 `../blog-data`），`main` 分支只含代码：

```
src/  
├── views/               # 页面组件
├── components/          # 可复用组件
├── stores/              # Pinia 状态管理
├── utils/               # 工具函数和生成器
└── router/              # 路由配置

# 数据分支（../blog-data，data 分支）
content-src/             # 第一层 src：纯手写中文源（仅这里编辑，严格无英文）
│   ├── notes/           # 笔记文章（唯一含 markdown 的目录）  
│   ├── projects/        # 项目元数据（纯 yaml，无文章）  
│   ├── topics/          # 课题元数据（纯 yaml，无文章）  
│   ├── categories.yaml  # 分类定义  
│   ├── about.yaml       # 关于页面  
│   └── resources.yaml   # 资源页面  
cache/                   # 第二层 cache：机器维护的持久态（入库）
│   ├── en/              #   英文翻译层（镜像 content-src，-en 后缀身份）
│   ├── tag-mapping.json #   zh→en 标签映射
│   └── .translate-state.json # 翻译增量状态
content/                 # 第三层 final：生成 JSON + html（可再生，不入库）
```

### 🚀 快速开始

### 🌍 环境

- **Node.js**: `>=23.6.0`（Node 原生 type stripping 直跑 `.ts` 生成器）

### 💻 安装

1. 克隆仓库并挂数据 worktree（与仓库同级目录 `../blog-data`）

   ```
   git clone https://github.com/zorrooz/zorrooz.github.io.git  
   cd zorrooz.github.io
   git worktree add ../blog-data data
   ```

2. 安装依赖

   ```
   npm install
   ```

3. 启动开发服务器
   
   ```
   npm run dev
   ```

   开发服务器将默认运行在 http://localhost:5173 。
4. 构建生产版本

   ```
   npm run build
   ```

5. 预览生产构建

   ```
   npm run preview
   ```

6. 发布（编辑 `../blog-data/content-src/**` 后）

   提交并推送数据分支，GitHub Actions 自动重新构建并部署到 GitHub Pages：

   ```
   npm run deploy
   ```

### 📖 详细用法

gblog 采用纯 Markdown + YAML 元数据的内容管理方式。系统通过构建时生成器将源文件转换为优化的 JSON 文件，实现快速的内容加载和导航。

**工作流程：**

1. 编写阶段：在数据分支的 `content-src/` 目录（`../blog-data/content-src/`）创建 Markdown 文件和 YAML 配置
2. 构建阶段：`prebuild`脚本自动运行生成器，将内容转换为 JSON（生成到 `../blog-data/content/`）
3. 运行阶段：Vue 组件加载 JSON 文件，渲染页面内容
4. 发布阶段：`npm run deploy` 提交并推送数据分支，CI 自动重新构建部署

#### 分类内容管理

项目支持三种内容类型，在`content-src/categories.yaml`中定义：

1. 笔记 (notes)

   用于学习记录和技术文档：

   ```
   - name: "分类标识符"  
     title: "显示标题"  
     desc: "分类描述说明"  
     date: "创建日期"  
     categories:  
       子分类键: "子分类显示名称"  
       另一个子分类键: "子分类显示名称"
       ...
   ```

2. 项目 (projects)

   用于代码项目和实战案例：

   ```
   - name: "项目标识符"  
     title: "项目显示名称"  
     desc: "项目功能描述"  
     github: "GitHub仓库地址"  
     date: "项目创建日期"  
     categories:  
       子分类键: "子分类显示名称"  
       另一个子分类键: "子分类显示名称"
       ...
   ```

3. 课题 (topics)

   用于研究课题和学术内容：
   
   ```
   - name: "课题标识符"  
     title: "课题显示名称"  
     desc: "研究内容描述"  
     doi: "学术标识符"  
     date: "课题结束日期"  
     categories:  
       子分类键: "子分类显示名称"  
       另一个子分类键: "子分类显示名称"
       ...
   ```

#### Markdown 文件编写

本项目采用纯 Markdown，支持标准的 Markdown 语法，请在`content-src/categories.yaml`中编辑元信息。

文件组织结构：

```
content-src/  
├── categories.yaml          # 分类定义文件  
├── notes/                   # 笔记目录（唯一含 markdown 的目录）  
│   ├── 分类标识符/  
│   │   ├── 子分类/  
│   │   │   └── 文章名/
│   │   │       └── 文章名.md  
│   │   └── 其他子分类/  
│   │       └── 文章名/
│   │           └── 文章名.md  
│   └── 其他分类/  
```

> 注：**项目（projects）与课题（topics）为纯 yaml 元数据**，在 `categories.yaml` 中定义（github/doi/url/language/license 等字段），不再有 markdown 文章；分类页卡片直接外链到 GitHub / DOI。

工作流程：

1. 在数据分支的 `content-src`（`../blog-data/content-src/`）目录下创建或编辑 Markdown 文件
2. 编写文章内容
3. 运行`npm run dev`查看效果（dev 服务器自动重新生成内容并热刷新）
4. 发布：`npm run deploy`

#### 资源页面配置

资源页面采用三级层次结构：分类 → 子分类 → 资源项，在`content-src/resources.yaml`中定义：

```
# 顶级分类  
- title: "分类名称"  
  children:  
    # 子分类  
    - title: "子分类名称"  
      items:  
        # 具体资源  
        - name: "资源名称"  
          url: "https://example.com"  
          desc: "资源描述"
```

#### 关于页面配置

关于页面采用三段式结构：个人介绍 + 内容区块 + 联系方式，在`content-src/about.yaml`中定义：

```
# 个人简介（必填）  
introduction: "这里是你的自我介绍，支持多行文本"  
  
# 内容区块（可选，可多个）  
section:  
  - title: "区块标题"  
    items:  
      - item: "项目名称或职位"  
        desc: "详细描述或说明"  
      - item: "另一个项目"  
        desc: "补充信息"  
  
# 联系方式（可选，可多个）  
contacts:  
  - label: "联系类型"  
    value: "显示文本"  
    link: "链接地址"  
    icon: "Font Awesome 图标类"
```

#### 使用翻译工具

在开始翻译前，请先完成API Key配置：

请在`scripts`目录下新建`llmConfig.ts`文件（该文件被 gitignore，不会入库），并填入以下内容：

```
export default {
  url: 'your_api_endpoint_url',  // 替换为实际的 API 端点 URL
  apikey: 'your_api_key_here',   // 替换为实际的 API 密钥
  model: 'your_model_name',      // 替换为使用的模型名称
};
```

配置完成后，可使用以下命令进行翻译操作：

- 增量翻译（推荐）

  ```
  npm run translate
  ```

- 仅翻译新文件
  
  ```
  npm run translate -- translate --new
  ```

- 强制重新翻译所有文件

  ```
  npm run translate -- translate --force
  ```

- 检查翻译状态

  ```
  npm run translate -- status
  ```

- 更多用法
  
  使用`npm run translate -- help`查看详细用法。

翻译完成后，请手动构建完成内容更新：

```
npm run prebuild
```

#### 内容生成器

一般情况下，内容生成器会在构建时自动运行，无需手动干预`runAllGenerators`。
如需手动更新所有内容，可运行：

```
npm run prebuild
```

生成器会依次处理：

1. 生成笔记、项目、课题索引
2. 分类结构生成（依赖步骤1）
3. 文章列表生成（依赖笔记信息）
4. 标签云生成（依赖文章信息）

## ⚙️ 致谢

本项目由 Zorrooz 个人开发，采用了多个 LLM 及 AI IDE 辅助编码，感谢以下产品提供的服务：

- [deepseek](https://chat.deepseek.com/)
- [千问](https://www.qianwen.com/)
- [豆包](https://www.doubao.com/chat/)
- [Trae](https://www.trae.cn/)
- [CodeBuddy](https://www.codebuddy.ai/)
- [Visual Studio Code (GitHub Copilot)](https://code.visualstudio.com/)


## 📄 许可证

版权所有 (c) 2025-2026 Zorrooz。保留所有权利。

本项目采用 MIT 许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

---
**如果本项目对你有帮助，请考虑 Star ⭐ 本项目，感谢！**
