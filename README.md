# gblog

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)
[![GitHub stars](https://img.shields.io/github/stars/zorrooz/zorrooz.github.io?style=social)](https://github.com/zorrooz/zorrooz.github.io)
[![English Version](https://img.shields.io/badge/English-Version-blue)](README_EN.md)

[English](README_EN.md) | [中文](README.md)

gblog 是基于 Vue 3 + Vite 7 构建的现代化双语静态博客系统，支持中英文切换、亮暗主题、Markdown 文章管理、全文检索与 PWA 离线访问。

## ✨ 项目介绍

gblog 采用 **代码/数据双分支** 架构：`main` 分支只含代码，所有内容数据独立存放在 `data` 分支（本地经 git worktree 挂到 `../blog-data`）。内容以 Markdown + YAML 源文件撰写，构建期由 Node 生成器转换为 JSON/HTML 产物，前端以 SSG 预渲染 + 运行时加载的方式呈现，兼顾 SEO 与首屏速度。

### 🎯 核心特性

- **双语内容**：`/zh`、`/en` 语言前缀路由，内容按「中文手写源 + AI 翻译镜像」双层管理
- **主题切换**：亮色 / 暗色 / 跟随系统三种模式（Bootstrap 5 主题变量 + Sass）
- **SSR / SSG 预渲染**：vite-ssg 构建期真实渲染每个页面，利于 SEO
- **全文检索**：MiniSearch 客户端搜索，索引由构建期预生成（`search-index.json`）
- **Markdown 渲染**：GFM + 代码高亮（highlight.js）+ 数学公式（KaTeX）
- **内容复制**：文章一键复制纯文本；表格复制为 TSV，可直接粘贴到 Excel
- **完整阅读体验**：三栏文章布局（目录树 + 正文 + OnThisPage）、TOC 抽屉、阅读时长、标签筛选、回到顶部
- **PWA 离线**：vite-plugin-pwa 自动生成 Service Worker，内容产物预缓存
- **自动翻译**：DeepSeek API 增量翻译（内容哈希判定，零成本跳过未变化文件）
- **响应式设计**：桌面/平板/移动端统一边距与断点

### 🛠️ 技术栈

| 领域 | 技术 |
|------|------|
| 前端框架 | Vue 3（Composition API + `<script>` 混合） |
| 构建工具 | Vite 7（ESM） |
| SSR / SSG | vite-ssg（构建期逐页预渲染） |
| 状态管理 | Pinia |
| 路由 | Vue Router 4（history 模式，语言前缀 + 旧 URL 重定向） |
| UI 框架 | Bootstrap 5 + Sass |
| Markdown | unified（remark-parse → remark-gfm → remark-breaks → remark-math → remark-rehype → rehype-highlight → rehype-katex → rehype-stringify） |
| 国际化 | vue-i18n（App 层） + 内容层双语镜像 |
| 搜索 | MiniSearch（索引构建期预生成） |
| PWA | vite-plugin-pwa（autoUpdate + generateSW） |
| AI 翻译 | OpenAI SDK + DeepSeek（增量翻译器 + 标签映射器） |
| 部署 | GitHub Actions + GitHub Pages |
| Node | `>=23.6.0`（原生 type stripping 直跑 `.ts` 脚本） |

### 📁 项目结构

```
# 代码分支（main）
src/                          # 浏览器应用（含 SSR；不含任何 Node 构建脚本）
├── main.ts                   # Entry: ViteSSG → Pinia → Router → i18n（内联 v-reveal 指令）
├── App.vue                   # AppHeader + router-view + AppFooter
├── router/index.ts           # /{zh,en} × 5 路由 + 旧 URL 重定向
├── i18n/                     # vue-i18n 实例 + locales/{zh-CN,en-US}.ts
├── stores/                   # theme.ts / locale.ts（Pinia）
├── config.ts                 # SITE 常量、locale 映射、主题模式
├── types.ts                  # 与 content/*.json 产物一一对应的领域类型
├── views/                    # Home / Category / Resource / About / Article
├── components/
│   ├── layout/               # AppHeader, NavActions, AppFooter, PostList,
│   │                         # RenderMarkdown, NavigationTree, OnThisPage
│   └── widgets/              # BackToTop, SearchModal, TocDrawer
├── composables/              # useFloatingButton / useLocalizedContent
├── utils/                    # 浏览器运行时工具（不 import scripts/）
│   ├── contentLoader.ts      # import.meta.glob 强类型加载（@data 别名）
│   ├── navigation.ts         # 语言前缀路径、语言切换、tag 跳转、文章路径转换
│   ├── scroll.ts             # 回到顶部 + 弹层滚动锁定
│   ├── clipboard.ts          # 复制（Clipboard API + 回退）
│   └── readingTime.ts        # 阅读时长估算
└── assets/                   # 字体 OTF / styles / avatar

scripts/                      # Node 内容工具链（构建期运行，不进浏览器）
├── dataConfig.ts             # 数据目录唯一接入点（支持 GBLOG_DATA_DIR）
├── markdownProcessor.ts      # unified pipeline: md → HTML
├── runAllGenerators.ts       # 生成器编排（顺序见下）
├── llmConfig.ts              # DeepSeek API 配置（gitignore，勿提交）
├── generators/               # core/（共享 IO）+ 11 个生成器
├── translator/               # AI 增量翻译（CLI：translate/status/help）
└── tagMerger/                # zh→en 标签映射补齐 + 一致性修复（CLI）

vite/contentDevPlugin.ts      # dev 插件：监听数据目录 → 重跑生成器 → full-reload

public/                       # favicon.png / apple-touch-icon.png / icon-512.png（PWA 图标）

# 数据分支（data，本地挂 ../blog-data）
content-src/                  # 第一层 src：纯手写中文源（仅这里编辑）
│   ├── categories.yaml       # 分类 + notes/projects/topics 全部元数据
│   ├── about.yaml            # 关于页
│   ├── resources.yaml        # 资源页
│   └── notes/<分类>/<子分类>/<文章>/<文章>.md
cache/                        # 第二层 cache：机器维护的持久态（入库）
│   ├── en/                   #   英文翻译层（镜像 content-src，-en 身份后缀）
│   ├── tag-mapping.json      #   zh→en 标签映射
│   └── .translate-state.json #   翻译增量状态（源相对路径 → 内容 hash）
content/                      # 第三层 final：生成 JSON + html（可再生，不入库）
```

## 🚀 快速开始

### 🌍 环境要求

- **Node.js** `>=23.6.0`（Node 原生 type stripping 直跑 `.ts` 生成器）

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

### 常用命令

| 命令 | 说明 |
|------|------|
| `npm run dev` | 开发服务器 http://localhost:5173（监听数据目录，自动重生成内容并热刷新） |
| `npm run prebuild` | 仅运行内容生成器 |
| `npm run build` | 生成器 + vite-ssg 构建（SSG 预渲染 + PWA） |
| `npm run preview` | 预览生产构建 |
| `npm run lint` | ESLint（自动修复） |
| `npm run typecheck` | vue-tsc + tsc 双项目类型检查 |
| `npm run format` | Prettier（src/ scripts/ vite/） |
| `npm run data:translate` | AI 增量翻译（写 cache/en） |
| `npm run data:tag-merge` | zh→en 标签映射增量补齐 |
| `npm run data:pack` | 打包文章（assets/ 引用图 → 文章同目录；无参数 = 递归全部，幂等） |
| `npm run data:deploy` | 提交并推送数据分支变更，触发 CI 重新构建部署 |
| `npm run data:publish` | 仅推送数据分支，不提交本地变更 |

> 命令规范：数据分支操作统一 `data:` 前缀；`translate`/`tagmerge`/`pack`/`deploy` 为旧名兼容别名。

### 发布流程

1. 编辑 `../blog-data/content-src/**`（只写中文源）
2. 本地验证：`npm run prebuild`（必要时先 `npm run data:translate` 生成英文）
3. 发布：`npm run data:deploy`（自动 `add -A → commit → push origin data`）
4. GitHub Actions 自动重新构建并部署到 GitHub Pages

> 也可用 `npm run data:publish` 只推送数据分支而不提交本地变更。

## 📖 详细用法

### 内容三层模型

- **第一层 `content-src/`**：纯手写中文源，**唯一人工编辑区**，严格不含机器产物
- **第二层 `cache/`**：机器维护的持久态（英文翻译镜像、标签映射、翻译状态），**入库留证**
- **第三层 `content/`**：可再生的最终产物（JSON + HTML），gitignore 不入库

CI 构建时设置 `GBLOG_NO_TRANSLATE=1` 跳过翻译与标签补齐，直接使用已入库的 `cache/` 产物。

### 分类内容管理

三种内容类型均在 `content-src/categories.yaml` 中定义：

1. **笔记（notes）**——学习记录与技术文档，有 Markdown 文章

   ```yaml
   - name: "分类标识符"
     title: "显示标题"
     desc: "分类描述说明"
     date: "创建日期"
     categories:
       子分类键: "子分类显示名称"
       ...
   ```

2. **项目（projects）**——纯 yaml 元数据，无文章页，卡片外链 GitHub

   ```yaml
   - name: "项目标识符"
     title: "项目显示名称"
     desc: "项目功能描述"
     github: "GitHub 仓库地址"
     date: "创建日期"
     categories:
       子分类键: "子分类显示名称"
       ...
   ```

3. **课题（topics）**——纯 yaml 元数据，无文章页，卡片外链 DOI

   ```yaml
   - name: "课题标识符"
     title: "课题显示名称"
     desc: "研究内容描述"
     doi: "学术标识符"
     date: "结束日期"
     categories:
       子分类键: "子分类显示名称"
       ...
   ```

> 项目/课题可携带 `github` / `doi` / `url` / `status` / `language` / `license` / `journal` / `year` / `authors` 等完整元数据字段。

### Markdown 文章编写

文章目录结构（目录名必须匹配 `categories.yaml` 中 notes 段的 `name` 字段）：

```
content-src/notes/
├── <分类标识符>/
│   ├── <子分类>/
│   │   └── <文章名>/
│   │       └── <文章名>.md
│   └── <其他子分类>/
│       └── <文章名>/
│           └── <文章名>.md
└── <其他分类>/
```

每篇文章自带 YAML frontmatter：

```yaml
---
title: 文章标题
date: 2025-01-01
author: zorrooz
tags: [标签1, 标签2]
draft: false
description: 摘要（用于列表卡片与搜索）
---
正文 Markdown 内容……
```

> **中英 frontmatter `tags` 必须一一对应**（同一篇文章 zh/en 标签数量与语义对称），否则首页标签云双语标签数不一致。

### 双语机制

- 英文**绝不手写**：`npm run data:translate` 将中文源增量翻译到 `cache/en/` 镜像（`-en` 后缀是内容身份，是 URL 的一部分）
- 翻译器按**内容 hash** 判定增量：源文件未变化即跳过，命中零成本
- 标签经 `cache/tag-mapping.json` 映射表查表翻译，缺失映射由 `npm run data:tag-merge` 增量补齐
- 生成器按 locale 取源：`srcDirFor(locale)`（zh→content-src，en→cache/en），输出恒为 `content/*` 与 `content/*-en*`

#### 翻译工具

翻译前需在 `scripts/llmConfig.ts` 配置 API（该文件被 gitignore，不会入库）：

```ts
export default {
  url: 'https://api.deepseek.com',
  apikey: 'your_api_key_here',
  model: 'your_model_name',
  // thinking: true, // 可选：显式开启 DeepSeek 思考模式（默认关闭，省 token 且 temperature 生效）
}
```

常用命令（`--` 之后的参数传给 CLI）：

```bash
npm run data:translate                                # 增量翻译（默认 content-src，推荐）
npm run data:translate -- translate --new             # 仅翻译新文件
npm run data:translate -- translate --force           # 强制重新翻译所有文件
npm run data:translate -- status <目录>               # 检查翻译状态
npm run data:translate -- help                        # 查看完整用法
npm run data:tag-merge                                # 补齐 zh→en 标签映射
npm run data:pack -- notes/Programming/bash/bash-scripting   # 打包单篇文章
npm run data:pack                                     # 递归打包全部文章
```

### 资源页面配置

`content-src/resources.yaml` 采用三级层次：分类 → 子分类 → 资源项：

```yaml
- title: "分类名称"
  children:
    - title: "子分类名称"
      items:
        - name: "资源名称"
          url: "https://example.com"
          desc: "资源描述"
```

### 关于页面配置

`content-src/about.yaml` 三段式结构：个人介绍 + 内容区块 + 联系方式：

```yaml
introduction: "自我介绍，支持多行文本"
section:
  - title: "区块标题"
    items:
      - item: "项目名称或职位"
        desc: "详细描述或说明"
contacts:
  - label: "联系类型"
    value: "显示文本"
    link: "链接地址"
    icon: "Font Awesome 图标类"
```

### 内容生成器

`npm run prebuild`（或 `npm run build` 时自动）按严格顺序运行 11 个生成器（`scripts/generators/`，编排见 `scripts/runAllGenerators.ts`）：

```
1.  generateNotes      笔记索引（notes/**/*.md → notes.json）
2.  generateProjects   项目元数据（categories.yaml → projects.json）
3.  generateTopics     课题元数据（categories.yaml → topics.json）
4.  generateCategories 分类结构（YAML + 前三步产物 → categories.json）
5.  generatePosts      文章列表（notes + categories → posts.json）
6.  generateTags       标签云（posts → tags.json）
    （1-6 先 zh-CN，翻译/标签补齐后重复 en-US）
7.  增量翻译           内容 hash 判定，写 cache/en（GBLOG_NO_TRANSLATE=1 关闭）
8.  tagMerger          标签映射补齐 + 以 zh 为基准修复 -en 文件 tags
9.  generateHtml       文章 md → content/html/**（预渲染 HTML）
10. generateSitemap    public/sitemap.xml + robots.txt
11. generateSearchIndex categories + html → search-index{,-en}.json
```

产物约定：除 html/sitemap 外均为 `content/{名称}.json`，英文追加 `-en` 后缀；`generateAbout` / `generateResources` 无依赖、独立运行。

### 搜索

构建期由 `generateSearchIndex` 将文章元数据与预渲染 HTML 清洗为纯文本索引（`content/search-index{,-en}.json`），运行时由 SearchModal 用 MiniSearch 全文检索，支持中英双语索引切换。

### 路由与旧链接

- 所有页面带语言前缀：`/zh/`、`/en/`（文章为 `/zh/article/notes/...`）
- 旧的无前缀 URL（`/article/...`、`/category` 等）自动 302 重定向到首选语言前缀
- 部署到 GitHub Pages 时构建 `404.html` 兜底，PWA 离线可访问

## ⚙️ CI / 部署

`.github/workflows/build.yml`：`main` 与 `data` 分支双检出 → Node 24 安装依赖 → 设置 `GBLOG_DATA_DIR` 与 `GBLOG_NO_TRANSLATE=1` → `prebuild + build` → peaceiris 部署到 `gh-pages` 分支。

## 🤝 致谢

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
