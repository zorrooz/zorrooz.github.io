# AGENTS.md — gblog

A static personal blog system built with Vue 3 + Vite 7 + Bootstrap 5. Markdown + YAML source files are processed at build-time into JSON metadata, loaded by a Vue SPA at runtime. Bilingual (zh-CN / en-US), theme switching, GitHub Pages deployment.

## Quick Start

| Command | Description |
|---------|-------------|
| `npm run dev` | Dev server at http://localhost:5173 |
| `npm run build` | Generators → `vite build` |
| `npm run prebuild` | Content generators only |
| `npm run preview` | Preview production build |
| `npm run lint` | ESLint with auto-fix |
| `npm run format` | Prettier on `src/` |
| `npm run translate` | AI translation CLI (incremental) |
| `npm run deploy` | Build + gh-pages deploy |

**Always lint and build after changes.**

## Documents (约定文档索引)

| 文档 | 说明 |
|------|------|
| `docs/modernization-plan.md` | **架构现代化执行计划**（Phase 0-6 分阶段执行；每阶段完成需更新其状态表） |
| `.qoder/repowiki/` | **代码索引（Qoder 生成）**：项目 Wiki 文档树，按需加载上下文 |

**注意**：架构改造期间（Phase 0-6 进行中），下文 Architecture 描述的是现状基线；涉及改造的文件与命令以 `docs/modernization-plan.md` 为准，每阶段完成后同步更新本文件。

### Repowiki 代码索引（按需加载）

`.qoder/repowiki/` 是 Qoder 自动生成的代码索引，含模块级知识卡与完整文档树。**需要了解某模块/子系统时，先按此加载上下文，再动手改代码**：

```
.qoder/repowiki/
├── knowledge/zh/_index.yaml        # 模块索引（modules: name → scope 文件清单/子模块）
├── knowledge/zh/<模块>/_module.yaml# 模块元数据（scope、依赖、关联）
├── knowledge/zh/<模块>/*.md        # 模块知识卡：概览 / 架构设计 / 技术栈 / 应用结构 / 开发规范
├── zh/content/                     # 完整文档树：架构、内容系统、样式、指南、API 参考等
└── zh/meta/repowiki-metadata.json  # 文档间关系（parent-child 等）
```

使用规则：
1. 定位模块：先读 `knowledge/zh/_index.yaml`，按 `dir_name`/`scope` 找到目标模块
2. 加载上下文：读该模块知识卡与 `zh/content/` 对应章节，再读涉及的源码
3. 修改代码前，按需加载受影响模块的文档，避免凭猜测改动
4. 注意：该索引为**生成产物，可能滞后于代码**；内容与源码冲突时以源码为准，并在完成后可触发重新生成

## Architecture

```
src/
├── main.js                   # Entry: createApp → Pinia → Router → i18n
├── App.vue                   # <AppHeader> + <router-view> + <AppFooter>
├── router/index.js           # 5 hash-mode routes
├── stores/
│   ├── app.js                # Pinia: theme + locale state
│   ├── i18n.js               # vue-i18n (zh-CN / en-US)
│   ├── locales/{zh-CN,en-US}.js
│   └── styles/global.scss    # Bootstrap SCSS + CSS custom properties
├── views/
│   ├── Home.vue              # ProfileCard + TagCloud + PostList
│   ├── Category.vue          # Notes/Projects/Topics cards with stats
│   ├── Resource.vue          # Hierarchical resource catalog
│   ├── About.vue             # Timeline + sections + contacts
│   └── Article.vue           # NavigationTree + RenderMarkdown + OnThisPage
├── components/
│   ├── layout/               # AppHeader, AppFooter, PostList, ProfileCard, etc.
│   └── widgets/              # BackToTop, TocDrawer
├── content/                  # GENERATED JSON (do not edit)
├── content-src/              # SOURCE: YAML + Markdown (edit here)
│   ├── {categories,about,resources}.yaml
│   ├── notes/                # articles under category/subcategory/
│   ├── projects/
│   └── topics/
├── utils/
│   ├── contentLoader.js      # Runtime: import.meta.glob for JSON & MD
│   ├── markdownProcessor.js  # unified pipeline: remark → rehype
│   ├── runAllGenerators.js   # Build orchestration
│   ├── generators/           # 8 build-time Node.js scripts
│   └── translator/           # AI translation via DeepSeek API
└── assets/
    ├── fonts/                # Agave-Regular, SourceHanSansSC
    └── icons/                # gblog.svg, theme/lang toggle PNGs
```

## Framework

| Layer | Tech |
|-------|------|
| UI | Vue 3 Composition API (`<script>`, Options API mix) |
| Build | Vite 7 (ESM, `"type": "module"`) |
| State | Pinia |
| Routing | Vue Router 4, hash mode (`createWebHashHistory`) |
| i18n | vue-i18n (locale in localStorage) |
| Styling | Bootstrap 5 SCSS + CSS custom properties |
| Markdown | unified → remark-parse → remark-gfm → remark-breaks → remark-math → remark-rehype → rehype-highlight → rehype-katex → rehype-stringify |
| Highlight | highlight.js (CDN, theme per `data-bs-theme`) |
| Icons | Font Awesome 6.4 (CDN in index.html) |
| Deployment | gh-pages |
| Node | `^20.19.0 \|\| >=22.12.0` |

## Content Pipeline (Critical Order)

```
1. generateNotes     — reads content-src/notes/**/*.md → notes.json
2. generateProjects  — reads content-src/projects/**/*.md → projects.json
3. generateTopics    — reads content-src/topics/**/*.md → topics.json
4. generateCategories— YAML + notes/projects/topics → categories.json
5. generatePosts     — notes.json + categories.json → posts.json
6. generateTags      — posts.json → tags.json
(1-6 for zh-CN, then repeat for en-US)

generateAbout and generateResources run independently at import time in runAllGenerators.js.
```

## Content Source Organization

```
content-src/{notes|projects|topics}/<category>/<subcategory>/<article-dir>/<article-dir>.md
```

- Directory name under `notes/`/`projects/`/`topics/` must match a `name` field in `categories.yaml`.
- No YAML frontmatter in .md — metadata is in `categories.yaml`.

### URL Mapping

File: `notes/Omics/genomics/bwa/bwa.md`
URL: `/#/article/notes/Omics/genomics/bwa/bwa`
Route: `{ name: 'Article', params: { path: ['notes', 'Omics', 'genomics', 'bwa', 'bwa'] } }`

## Internationalization (Dual Layer)

**Layer 1 — UI:** `src/stores/locales/{zh-CN,en-US}.js` via vue-i18n. Controls nav, buttons, labels.

**Layer 2 — Content:** `-en` suffix pattern:
- Chinese: `article.md`, `categories.json`
- English: `article-en.md`, `categories-en.json`
- Auto-fallback in `contentLoader.js`: if English file missing, loads Chinese.

## Theme System

Three modes cycle: `auto` → `light` → `dark`. Persisted in `localStorage.theme`.
Mechanism: `document.documentElement.setAttribute('data-bs-theme', mode)`
CSS custom properties in `global.scss` for both themes. Highlight.js theme swapped via MutationObserver.

## Key Conventions

### Code Style
- No semicolons, single quotes, 2-space indent, LF, 100-char width
- Vue: `<script>` (Options API with `setup()` return)
- No unnecessary comments
- `@/` alias for `src/`

### Imports
- ES modules only (no `require`)
- Named exports preferred over default

### Content Types
- **Notes** (学习): technical learning articles, organized by topic
- **Projects** (项目): software documentation, GitHub-linked
- **Topics** (课题): research case studies, DOI-linked

### Adding Content
1. Create directory under `content-src/{notes|projects|topics}/`
2. Create `.md` file
3. Add category in `categories.yaml`
4. Run `npm run prebuild`
5. English: `article-en.md` or `npm run translate`

### Utilities
- `contentLoader.js` provides: `loadPosts()`, `loadCategories()`, `loadNotes()`, `loadTags()`, `loadAbout()`, `loadResources()`, `loadMarkdownContent()`
- `markdownProcessor.js` export: `renderMarkdown(markdown)` → HTML string
- `appStore`: `toggleTheme()`, `toggleLanguage()`, `initTheme()`, `initLocale()`

## ESLint / Prettier
- Flat config (`eslint.config.js`) with Vue + Prettier
- `browser` globals for `*.{js,vue}`, `node` globals for `src/utils/**/*.js`
- Ignored: `dist/`, `dist-ssr/`, `coverage/`
