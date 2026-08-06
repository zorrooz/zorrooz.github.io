# gblog

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)
[![GitHub stars](https://img.shields.io/github/stars/zorrooz/zorrooz.github.io?style=social)](https://github.com/zorrooz/zorrooz.github.io)
[![中文版](https://img.shields.io/badge/中文-Version-blue)](README.md)

[English](README_EN.md) | [中文](README.md)

gblog is a modern bilingual static blog system built with Vue 3 + Vite 7, featuring English/Chinese switching, light/dark themes, Markdown article management, full-text search, and PWA offline support.

## ✨ Project Introduction

gblog uses a **code/data dual-branch** architecture: the `main` branch contains code only, while all content lives on the `data` branch (attached locally via a git worktree at `../blog-data`). Content is written as Markdown + YAML sources, converted to JSON/HTML artifacts by Node generators at build time, and presented through SSG pre-rendering with runtime loading — balancing SEO and first-paint speed.

### 🎯 Core Features

- **Bilingual Content**: `/zh` and `/en` locale-prefixed routing; content managed as "hand-written Chinese source + AI-translated mirror"
- **Theme Switching**: light / dark / follows-system modes (Bootstrap 5 theme variables + Sass)
- **SSR / SSG Pre-rendering**: vite-ssg renders every page at build time for SEO
- **Full-text Search**: MiniSearch client-side search, index pre-generated at build time (`search-index.json`)
- **Markdown Rendering**: GFM + syntax highlighting (highlight.js) + math formulas (KaTeX)
- **Content Copying**: one-click article text copy; tables copied as TSV, paste-ready for Excel
- **Complete Reading Experience**: three-column article layout (tree nav + content + on-this-page), TOC drawer, reading time, tag filtering, back-to-top
- **PWA Offline**: vite-plugin-pwa auto-generates a Service Worker and precaches content artifacts
- **Auto Translation**: DeepSeek API incremental translation keyed on content hashes (untouched files cost nothing to skip)
- **Responsive Design**: unified margins and breakpoints across desktop / tablet / mobile

### 🛠️ Tech Stack

| Area | Technology |
|------|------------|
| Frontend | Vue 3 (Composition API + `<script>` mix) |
| Build tool | Vite 7 (ESM) |
| SSR / SSG | vite-ssg (per-page pre-rendering) |
| State | Pinia |
| Routing | Vue Router 4 (history mode, locale prefix + legacy redirects) |
| UI | Bootstrap 5 + Sass |
| Markdown | unified (remark-parse → remark-gfm → remark-math → rehype-highlight → rehype-katex) |
| i18n | vue-i18n (app layer) + mirrored content layer |
| Search | MiniSearch (index pre-generated at build time) |
| PWA | vite-plugin-pwa (autoUpdate + generateSW) |
| AI Translation | OpenAI SDK + DeepSeek (incremental translator + tag mapper) |
| Deployment | GitHub Actions + GitHub Pages |
| Node | `>=23.6.0` (native type stripping runs `.ts` scripts directly) |

### 📁 Project Structure

```
# Code branch (main)
src/                          # Browser app (incl. SSR; no Node build scripts)
├── main.ts                   # Entry: ViteSSG → Pinia → Router → i18n
├── App.vue                   # AppHeader + router-view + AppFooter
├── router/index.ts           # /{zh,en} × 5 routes + legacy URL redirects
├── i18n/                     # vue-i18n instance + locales/{zh-CN,en-US}.ts
├── stores/                   # theme.ts / locale.ts (Pinia)
├── config/site.ts            # SITE constants, locale maps, theme modes
├── views/                    # Home / Category / Resource / About / Article
├── components/
│   ├── layout/               # AppHeader, NavActions, AppFooter, PostList,
│   │                         # RenderMarkdown, NavigationTree, OnThisPage
│   └── widgets/              # BackToTop, SearchModal, TocDrawer
├── composables/              # useFloatingButton / useLocalizedContent
├── utils/                    # Browser runtime helpers (do not import scripts/)
│   ├── contentLoader.ts      # Typed import.meta.glob loading (@data alias)
│   ├── navigation.ts         # Locale paths, language switch, tag nav, article paths
│   ├── scroll.ts             # Back-to-top + overlay scroll locking
│   ├── clipboard.ts          # Copy (Clipboard API + fallback)
│   ├── readingTime.ts        # Reading time estimate
│   └── reveal.ts             # v-reveal scroll-in directive
├── types/content.ts          # Domain types matching content/*.json artifacts
└── assets/                   # Fonts / styles / avatar

scripts/                      # Node content toolchain (build-time, never bundled)
  ├── dataConfig.ts           # Single entry for data dirs (supports GBLOG_DATA_DIR)
  ├── markdownProcessor.ts    # unified pipeline: md → HTML
  ├── runAllGenerators.ts     # Generator orchestration
  ├── llmConfig.ts            # DeepSeek API config (gitignored, never committed)
  ├── generators/             # core/ (shared IO) + 11 generators
  ├── translator/             # AI incremental translation (CLI: translate/status/help)
  └── tagMerger/              # zh→en tag mapping + consistency fix (CLI)

vite/contentDevPlugin.ts      # Dev plugin: watches data dir → rerun generators → full-reload

# Data branch (data, worktree at ../blog-data)
content-src/                  # Layer 1: hand-written Chinese sources (edit only here)
  ├── categories.yaml         # categories + all notes/projects/topics metadata
  ├── about.yaml              # about page
  ├── resources.yaml          # resources page
  └── notes/<cat>/<sub>/<article>/<article>.md
cache/                        # Layer 2: machine-maintained persistent state (committed)
  ├── en/                     #   English translation layer (mirrors content-src, -en identity)
  ├── tag-mapping.json        #   zh→en tag mapping
  └── .translate-state.json   #   incremental translation state (source → content hash)
content/                      # Layer 3: generated JSON + html (regenerable, not committed)
```

## 🚀 Getting Started

### 🌍 Environment

- **Node.js** `>=23.6.0` (native type stripping runs `.ts` generators directly)

### 💻 Installation

1. Clone the repository and attach the data worktree (sibling dir `../blog-data`)

   ```
   git clone https://github.com/zorrooz/zorrooz.github.io.git
   cd zorrooz.github.io
   git worktree add ../blog-data data
   ```

2. Install dependencies

   ```
   npm install
   ```

### Common Commands

| Command | Description |
|---------|-------------|
| `npm run dev` | Dev server at http://localhost:5173 (watches data dir; regenerates content + HMR) |
| `npm run prebuild` | Content generators only |
| `npm run build` | Generators + vite-ssg build (SSG pre-render + PWA) |
| `npm run preview` | Preview the production build |
| `npm run lint` | ESLint (auto-fix) |
| `npm run typecheck` | vue-tsc + tsc dual-project type check |
| `npm run format` | Prettier (`src/ scripts/ vite/`) |
| `npm run translate` | AI incremental translation (writes `cache/en`) |
| `npm run tagmerge` | Incrementally complete zh→en tag mappings |
| `npm run deploy` | Commit and push the data branch; CI rebuilds and deploys |

### Publish Flow

1. Edit `../blog-data/content-src/**` (Chinese sources only)
2. Verify locally: `npm run prebuild` (run `npm run translate` first if new content needs English)
3. Publish: `npm run deploy` (auto `add -A → commit → push origin data`)
4. GitHub Actions rebuilds and deploys to GitHub Pages

> `npm run data:publish` pushes the data branch only, without committing local changes.

## 📖 Detailed Usage

### Three-Layer Content Model

- **Layer 1 `content-src/`**: hand-written Chinese sources, the only human-edited area, strictly machine-free
- **Layer 2 `cache/`**: machine-maintained persistent state (English mirror, tag mapping, translation state), committed as evidence
- **Layer 3 `content/`**: regenerable final artifacts (JSON + HTML), gitignored

CI builds with `GBLOG_NO_TRANSLATE=1` to skip translation/tag-sync and reuse the committed `cache/` artifacts.

### Category Content Management

Three content types are defined in `content-src/categories.yaml`:

1. **Notes** — learning records & technical articles, with Markdown files

   ```yaml
   - name: "category_identifier"
     title: "Display Title"
     desc: "Category description"
     date: "creation_date"
     categories:
       subcategory_key: "Subcategory Display Name"
       ...
   ```

2. **Projects** — YAML metadata only, no article pages; cards link to GitHub

   ```yaml
   - name: "project_identifier"
     title: "Project Display Name"
     desc: "Project description"
     github: "GitHub repo URL"
     date: "creation_date"
     categories:
       subcategory_key: "Subcategory Display Name"
       ...
   ```

3. **Topics** — YAML metadata only, no article pages; cards link to DOI

   ```yaml
   - name: "topic_identifier"
     title: "Topic Display Name"
     desc: "Research description"
     doi: "academic_identifier"
     date: "end_date"
     categories:
       subcategory_key: "Subcategory Display Name"
       ...
   ```

> Projects/topics may carry full metadata fields: `github` / `doi` / `url` / `status` / `language` / `license` / `journal` / `year` / `authors`, etc.

### Writing Markdown Articles

Article directory structure (directory names must match the `name` fields in the `notes` section of `categories.yaml`):

```
content-src/notes/
├── <category_identifier>/
│   ├── <subcategory>/
│   │   └── <article_name>/
│   │       └── <article_name>.md
│   └── <other_subcategory>/
│       └── <article_name>/
│           └── <article_name>.md
└── <other_category>/
```

Each article carries YAML frontmatter:

```yaml
---
title: Article title
date: 2025-01-01
author: zorrooz
tags: [tag1, tag2]
draft: false
description: Summary (used in list cards and search)
---
Markdown body here…
```

> **zh/en frontmatter `tags` must correspond one-to-one** (equal count & semantics for the same article), otherwise the homepage tag cloud will show mismatched counts between languages.

### Bilingual Mechanism

- English is **never hand-written**: `npm run translate` incrementally translates Chinese sources into the `cache/en/` mirror (`-en` suffix is the content identity and part of the URL)
- Incremental detection is keyed on **content hash**: untouched sources are skipped, costing zero tokens
- Tags are looked up in `cache/tag-mapping.json`; missing mappings are filled incrementally by `npm run tagmerge`
- Generators pick sources per locale: `srcDirFor(locale)` (zh→content-src, en→cache/en); outputs are always `content/*` and `content/*-en*`

#### Translation Tool

Before translating, configure the API in `scripts/llmConfig.ts` (gitignored, never committed):

```ts
export default {
  url: 'https://api.deepseek.com',
  apikey: 'your_api_key_here',
  model: 'your_model_name',
}
```

Common commands (args after `--` are passed to the CLI):

```bash
npm run translate                                   # incremental (default content-src, recommended)
npm run translate -- translate --new                # translate new files only
npm run translate -- translate --force              # force re-translate everything
npm run translate -- status <directory>             # check translation status
npm run translate -- help                           # full usage
npm run tagmerge                                    # fill in zh→en tag mappings
```

### Resources Page Configuration

`content-src/resources.yaml` uses a three-level hierarchy: Category → Subcategory → Resource Item:

```yaml
- title: "Category Name"
  children:
    - title: "Subcategory Name"
      items:
        - name: "Resource Name"
          url: "https://example.com"
          desc: "Resource description"
```

### About Page Configuration

`content-src/about.yaml` uses a three-section structure: introduction + content blocks + contact info:

```yaml
introduction: "Self-introduction, multi-line supported"
section:
  - title: "Block Title"
    items:
      - item: "Project name or role"
        desc: "Detailed description"
contacts:
  - label: "Contact Type"
    value: "Display text"
    link: "Link address"
    icon: "Font Awesome icon class"
```

### Content Generators

`npm run prebuild` (or automatically during `npm run build`) runs 11 generators in a strict order (`scripts/generators/`; orchestrated in `scripts/runAllGenerators.ts`):

```
1.  generateNotes      post index (notes/**/*.md → notes.json)
2.  generateProjects   metadata (categories.yaml → projects.json)
3.  generateTopics     metadata (categories.yaml → topics.json)
4.  generateCategories category tree (YAML + outputs 1-3 → categories.json)
5.  generatePosts      post list (notes + categories → posts.json)
6.  generateTags       tag cloud (posts → tags.json)
    (1-6 run for zh-CN first, then repeated for en-US after translate tag sync)
7.  Incremental translation   content-hash keyed, writes cache/en (GBLOG_NO_TRANSLATE=1 disables)
8.  tagMerger           fill tag mappings + fix -en tag consistency
9.  generateHtml        article md → content/html/** (pre-rendered HTML)
10. generateSitemap     public/sitemap.xml + robots.txt
11. generateSearchIndex categories + html → search-index{,-en}.json
```

Artifact convention: everything but html/sitemap lands in `content/{name}.json`; English gets a `-en` suffix. `generateAbout` / `generateResources` run independently (no dependencies).

### Search

At build time `generateSearchIndex` strips article metadata + pre-rendered HTML to plain text (`content/search-index{,-en}.json`); at runtime the SearchModal performs full-text search via MiniSearch, honoring the current locale's index.

### Routing & Legacy Links

- Every page is locale-prefixed: `/zh/`, `/en/` (articles: `/zh/article/notes/...`)
- Legacy unprefixed URLs (`/article/...`, `/category`, etc.) are 302-redirected to the preferred locale prefix
- A `404.html` fallback is built for GitHub Pages; the PWA is accessible offline

## ⚙️ CI / Deployment

`.github/workflows/build.yml`: checks out `main` + `data` branches → Node 24 install → sets `GBLOG_DATA_DIR` + `GBLOG_NO_TRANSLATE=1` → `prebuild + build` → peaceiris deploys to the `gh-pages` branch.

## 🤝 Acknowledgments

This project was developed personally by Zorrooz, with assistance from multiple LLMs and AI IDEs. Thanks to the following services:

- [deepseek](https://chat.deepseek.com/)
- [千问](https://www.qianwen.com/)
- [豆包](https://www.doubao.com/chat/)
- [Trae](https://www.trae.cn/)
- [CodeBuddy](https://www.codebuddy.ai/)
- [Visual Studio Code (GitHub Copilot)](https://code.visualstudio.com/)

## 📄 License

Copyright (c) 2025-2026 Zorrooz. All rights reserved.

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
**If this project is helpful to you, please consider giving it a star ⭐. Thank you!**