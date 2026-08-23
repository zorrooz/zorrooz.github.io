# gblog

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)

[English](README_EN.md) | [中文](README.md)

基于 Vue 3 + Vite 7 的双语静态博客系统：Markdown + YAML 写作，构建期生成 JSON/HTML，SSG 预渲染 + PWA 离线。代码与内容分离双分支，Obsidian 即创作前端。

## ✨ 特性

- **双语**：`/zh` `/en` 语言前缀路由；中文手写源 + AI 翻译镜像（`cache/en`），增量翻译按内容 hash 跳过未变更
- **写作即 Obsidian**：vault 指向 `content-src`，模板 / 集中配图 / 边写边预览；`data:pack` 一键导出自包含文章包
- **阅读体验**：三栏文章布局（目录树 + 正文 + OnThisPage）、移动端 TOC 抽屉、阅读进度与时长、标签筛选、文章/代码块/表格复制
- **检索与渲染**：MiniSearch 全文检索（索引构建期预生成）、GFM + 代码高亮 + KaTeX 公式
- **主题**：亮 / 暗 / 跟随系统，reduce-motion 感知
- **PWA**：Service Worker 预缓存内容产物，离线可访问

## 🧰 技术栈

| 层 | 选型 |
|---|---|
| 框架 | Vue 3（`<script setup>`）+ Pinia + Vue Router |
| 构建 | Vite 7 + vite-ssg（SSG 预渲染）+ vite-plugin-pwa |
| 样式 | Bootstrap 5（栅格/工具类）+ 全局设计令牌 SCSS |
| 内容 | remark/rehype 管线（GFM、高亮、KaTeX）、js-yaml |
| 检索 | MiniSearch（CJK 二元分词） |
| 工具链 | Node 原生 type stripping 直跑 `.ts`（≥23.6），无构建编译层 |

## 🚀 快速开始

### 环境

- Node.js `>=23.6.0`

### 安装

```bash
git clone https://github.com/zorrooz/zorrooz.github.io.git
cd zorrooz.github.io
git worktree add ../blog-data data   # 数据分支挂到同级目录（代码/内容分离）
npm install
npm run dev                          # http://localhost:5173
```

### 常用命令

| 命令 | 说明 |
|------|------|
| `npm run dev` | 开发服务器 :5173（监听数据目录，自动重生成内容 + 热刷新） |
| `npm run build` | 生成器 + vite-ssg 构建（SSG 预渲染 + PWA） |
| `npm run prebuild` | 仅运行内容生成器（本地验证内容） |
| `npm run preview` | 预览生产构建 |
| `npm run lint` / `npm run typecheck` / `npm run format` | ESLint / 双项目类型检查 / Prettier |
| `npm run data:translate` | AI 增量翻译（写 `cache/en`） |
| `npm run data:tag-merge` | zh→en 标签映射增量补齐 |
| `npm run data:pack` | 导出文章为自包含包（md + 引用图 → `exports/`；`--out`/`--zip`） |
| `npm run data:deploy` | 提交并推送数据分支变更，触发 CI 重新构建部署 |
| `npm run data:publish` | 仅推送数据分支（不提交本地变更） |

> 数据分支操作统一 `data:` 前缀；`translate`/`tagmerge`/`pack`/`deploy` 为旧名兼容别名。

## 📝 内容写作

数据在 `../blog-data`，三层模型：

```
content-src/   第一层 · 纯手写中文源（唯一人工编辑区）
cache/         第二层 · 机器维护（英文镜像 / 标签映射 / 翻译状态，入库）
content/       第三层 · 可再生的最终产物（JSON + HTML，gitignore）
```

- **文章**：`content-src/notes/<分类>/<子分类>/<slug>.md`（扁平，一篇一 md），目录名匹配 `categories.yaml`
- **frontmatter**：`title` / `date` / `author` / `tags` / `draft` / `description`
- **配图**：统一放 `content-src/assets/`，引用 `assets/xx.png`（站点直接解析，无需打包）；分享/迁移时 `data:pack` 导出
- **英文**：绝不手写，`npm run data:translate` 生成；中英 frontmatter `tags` 必须一一对应

### 发布

```bash
npm run prebuild        # 本地验证
npm run data:deploy     # 提交推送 data 分支 → GitHub Actions 自动构建部署
```

## 🏗️ 项目结构

```
src/                        浏览器应用（Vue 3 + Pinia + vue-i18n；含 SSR）
├── main.ts                 应用入口：ViteSSG → Pinia → Router → i18n → v-reveal 注册
├── App.vue                 布局骨架 + SEO head + 全局快捷键（/ 打开搜索）
├── directives/reveal.ts    v-reveal 滚动入场指令
├── router/                 /{zh,en} × 5 路由 + 旧 URL 重定向 + scrollBehavior
├── config.ts               站点配置单一来源（locale 映射 / 主题模式 / 路由前缀等）
├── i18n/                   vue-i18n 实例 + locale 单一来源 + 类型 schema 校验
├── stores/                 Pinia：theme（主题）/ locale（语言）
├── views/                  页面：Home / Category / Resource / About / Article / NotFound
├── components/
│   ├── layout/             布局骨架：AppHeader / AppFooter / NavActions
│   ├── docs/               文档内容域：RenderMarkdown / OnThisPage / NavigationTree / TreeNode
│   ├── widgets/            浮动交互件：BackToTop / SearchModal / TocDrawer / FloatingButton
│   └── common/             通用件：PostList / PostCard / Pagination / IconButton / ModalOverlay / EmptyState
├── composables/            useLocalizedContent / usePageMeta / useSearch / useTagNavigation /
│                           useCopyFeedback / useFloatingButton / useReadingProgress /
│                           useStickySidebars / useRoutePagination
├── utils/                  纯工具：contentLoader / navigation / scroll / markdownImages /
│                           markdownDom / toc / articles / pagination / clipboard 等
└── assets/                 字体 OTF / styles / avatar.*

scripts/                    Node 内容工具链（原生 .ts 直跑）
├── dataConfig.ts           数据目录唯一接入点（GBLOG_DATA_DIR 可覆盖）
├── runAllGenerators.ts     生成器编排（locale 步骤表 + 依赖顺序）
├── lib/                    共享层：fs / text / cli / frontmatter / llm / tags / metadata…
├── generators/             11 个生成器 + markdownProcessor（md→HTML 渲染管线）
└── tools/                  独立 CLI 工具（每工具一目录：lib + *Cli.ts）
    ├── packArticle/        文章自包含包导出
    ├── translator/         AI 增量翻译
    └── tagMerger/          zh→en 标签映射

vite/contentDevPlugin.ts    dev 插件：监听数据目录 → 重跑生成器 → full-reload
.github/workflows/          CI：main + data 双检出 → prebuild + build → 部署 gh-pages
```

## 🔄 工作原理

```
写作（Obsidian → content-src/*.md,*.yaml）
        │ dev: contentDevPlugin 监听变更      CI: GBLOG_NO_TRANSLATE=1 跳过翻译
        ▼
生成器 scripts/generators/*
        ├─ 增量翻译（DeepSeek）→ cache/en 镜像
        ├─ notes/posts/categories/tags JSON
        ├─ markdownProcessor 预渲染 HTML
        └─ sitemap / robots.txt / search-index
        ▼
浏览器 src/
        ├─ SSG：vite-ssg 按 /zh、/en 逐页预渲染 + PWA 预缓存
        └─ 运行时：路由切换仅换数据，语言切换保持阅读位置
```

## 📄 许可证

MIT © Zorrooz — 见 [LICENSE](LICENSE)。
