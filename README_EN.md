# gblog

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)

[English](README_EN.md) | [中文](README.md)

A bilingual static blog built with Vue 3 + Vite 7: write in Markdown + YAML, generate JSON/HTML at build time, SSG pre-rendered with PWA offline support. Code and content live in separate branches; Obsidian is the authoring frontend.

## ✨ Features

- **Bilingual**: `/zh` `/en` locale-prefixed routing; hand-written Chinese sources + AI-translated mirror (`cache/en`) with hash-based incremental translation
- **Obsidian authoring**: point a vault at `content-src` — templates, central `assets/`, live preview; `data:pack` exports self-contained article bundles in one command
- **Reading experience**: three-column article layout (tree nav + content + on-this-page), mobile TOC drawer, reading progress & time, tag filtering, copy article / code block / table
- **Search & rendering**: MiniSearch full-text search (index pre-built at build time), GFM + syntax highlighting + KaTeX math
- **Theming**: light / dark / system, with reduce-motion support
- **PWA**: service worker pre-caches content artifacts for offline access

## 🧰 Tech Stack

| Layer | Choice |
|---|---|
| Framework | Vue 3 (`<script setup>`) + Pinia + Vue Router |
| Build | Vite 7 + vite-ssg (SSG) + vite-plugin-pwa |
| Styling | Bootstrap 5 (grid/utilities) + global design-token SCSS |
| Content | remark/rehype pipeline (GFM, highlighting, KaTeX), js-yaml |
| Search | MiniSearch (CJK bigram tokenization) |
| Toolchain | Node native type stripping runs `.ts` directly (≥23.6), no compile layer |

## 🚀 Getting Started

### Prerequisites

- Node.js `>=23.6.0`

### Install

```bash
git clone https://github.com/zorrooz/zorrooz.github.io.git
cd zorrooz.github.io
git worktree add ../blog-data data   # mount the data branch as a sibling worktree (code/content split)
npm install
npm run dev                          # http://localhost:5173
```

### Commands

| Command | Description |
|------|------|
| `npm run dev` | Dev server :5173 (watches the data dir, regenerates content + hot reload) |
| `npm run build` | Generators + vite-ssg build (SSG + PWA) |
| `npm run prebuild` | Run content generators only (local content check) |
| `npm run preview` | Preview the production build |
| `npm run lint` / `npm run typecheck` / `npm run format` | ESLint / dual-project type check / Prettier |
| `npm run data:translate` | AI incremental translation (writes `cache/en`) |
| `npm run data:tag-merge` | Incremental zh→en tag mapping sync |
| `npm run data:pack` | Export articles as self-contained bundles (md + images → `exports/`; `--out`/`--zip`) |
| `npm run data:deploy` | Commit & push the data branch, triggering CI build & deploy |
| `npm run data:publish` | Push the data branch only (no local commit) |

> Data-branch commands use the `data:` prefix; `translate`/`tagmerge`/`pack`/`deploy` are legacy aliases.

## 📝 Writing Content

Data lives in `../blog-data`, organized in three layers:

```
content-src/   Layer 1 · hand-written Chinese sources (the only human-edited area)
cache/         Layer 2 · machine-maintained state (English mirror / tag mapping, committed)
content/       Layer 3 · regenerable final artifacts (JSON + HTML, gitignored)
```

- **Articles**: `content-src/notes/<category>/<subcategory>/<slug>.md` (flat, one md per post); directory names must match `categories.yaml`
- **Frontmatter**: `title` / `date` / `author` / `tags` / `draft` / `description`
- **Images**: keep them in `content-src/assets/` and reference as `assets/xx.png` (resolved by the site, no bundling needed); use `data:pack` to export for sharing
- **English**: never hand-written — generated via `npm run data:translate`; frontmatter `tags` must map 1:1 between locales

### Publishing

```bash
npm run prebuild        # verify locally
npm run data:deploy     # commit & push data branch → GitHub Actions builds and deploys
```

## 🏗️ Project Structure

```
src/                        Browser app (Vue 3 + Pinia + vue-i18n; includes SSR)
├── main.ts                 Entry: ViteSSG → Pinia → Router → i18n → v-reveal registration
├── App.vue                 Layout shell + SEO head + global shortcut (/ opens search)
├── directives/reveal.ts    v-reveal scroll-in directive
├── router/                 /{zh,en} × 5 routes + legacy URL redirects + scrollBehavior
├── config.ts               Single source of site config (locale maps / theme modes / route prefixes…)
├── i18n/                   vue-i18n instance + locale single source + typed schema validation
├── stores/                 Pinia stores: theme / locale
├── views/                  Pages: Home / Category / Resource / About / Article / NotFound
├── components/
│   ├── layout/             Layout shell: AppHeader / AppFooter / NavActions
│   ├── docs/               Document domain: RenderMarkdown / OnThisPage / NavigationTree / TreeNode
│   ├── widgets/            Floating widgets: BackToTop / SearchModal / TocDrawer / FloatingButton
│   └── common/             Shared UI: PostList / PostCard / Pagination / IconButton / ModalOverlay / EmptyState
├── composables/            useLocalizedContent / usePageMeta / useSearch / useTagNavigation /
│                           useCopyFeedback / useFloatingButton / useReadingProgress /
│                           useStickySidebars / useRoutePagination
├── utils/                  Pure utilities: contentLoader / navigation / scroll / markdownImages /
│                           markdownDom / toc / articles / pagination / clipboard etc.
└── assets/                 Font OTFs / styles / avatar.*

scripts/                    Node content toolchain (native .ts execution)
├── dataConfig.ts           Single entry point for data dirs (overridable via GBLOG_DATA_DIR)
├── runAllGenerators.ts     Generator orchestration (locale step table + dependency order)
├── lib/                    Shared layer: fs / text / cli / frontmatter / llm / tags / metadata…
├── generators/             11 generators + markdownProcessor (md→HTML pipeline)
└── tools/                  Standalone CLI tools (one directory per tool: lib + *Cli.ts)
    ├── packArticle/        Self-contained article bundle exporter
    ├── translator/         AI incremental translation
    └── tagMerger/          zh→en tag mapping

vite/contentDevPlugin.ts    Dev plugin: watches the data dir → reruns generators → full-reload
.github/workflows/          CI: dual checkout of main + data → prebuild + build → deploy gh-pages
```

## 🔄 How It Works

```
Authoring (Obsidian → content-src/*.md,*.yaml)
        │ dev: contentDevPlugin watches changes       CI: GBLOG_NO_TRANSLATE=1 skips translation
        ▼
Generators scripts/generators/*
        ├─ Incremental translation (DeepSeek) → cache/en mirror
        ├─ notes/posts/categories/tags JSON
        ├─ markdownProcessor pre-rendered HTML
        └─ sitemap / robots.txt / search-index
        ▼
Browser src/
        ├─ SSG: vite-ssg pre-renders per /zh, /en page + PWA precaching
        └─ Runtime: route changes only swap data; language switch keeps reading position
```

## 📄 License

MIT © Zorrooz — see [LICENSE](LICENSE).
