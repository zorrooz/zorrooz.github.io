# gblog

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)

[English](README_EN.md) | [中文](README.md)

A bilingual static blog built with Vue 3 + Vite 7: write in Markdown + YAML, generate JSON/HTML at build time, SSG pre-rendered with PWA offline support. Code and content live in separate branches; Obsidian is the authoring frontend.

## ✨ Features

- **Bilingual**: `/zh` `/en` locale-prefixed routing; hand-written Chinese sources + AI-translated mirror (`cache/en`) with incremental translation that skips unchanged files at zero cost
- **Obsidian authoring**: point a vault at `content-src` — templates, central `assets/`, live preview; `data:pack` exports self-contained article bundles in one command
- **Reading experience**: three-column article layout (tree nav + content + on-this-page), TOC drawer, reading time, tag filtering, copy article/table as TSV
- **Search & rendering**: MiniSearch full-text search (index pre-generated at build time), GFM + syntax highlighting + KaTeX
- **Theme & polish**: light/dark/system, ambient header light + micro-grain + blur entrance, reduce-motion aware
- **PWA**: Service Worker precaches content artifacts for offline access

## 🚀 Quick Start

### Requirements

- Node.js `>=23.6.0` (native type stripping runs `.ts` scripts directly)

### Install

```bash
git clone https://github.com/zorrooz/zorrooz.github.io.git
cd zorrooz.github.io
git worktree add ../blog-data data   # attach data branch as a sibling (code/content separation)
npm install
```

### Commands

| Command | Description |
|---------|-------------|
| `npm run dev` | Dev server :5173 (watches data dir, regenerates content + HMR) |
| `npm run build` | Generators + vite-ssg build (SSG pre-rendering + PWA) |
| `npm run prebuild` | Content generators only (local content validation) |
| `npm run preview` | Preview production build |
| `npm run lint` / `npm run typecheck` / `npm run format` | ESLint / dual-project typecheck / Prettier |
| `npm run data:translate` | AI incremental translation (writes `cache/en`) |
| `npm run data:tag-merge` | Incremental zh→en tag mapping |
| `npm run data:pack` | Export articles as self-contained bundles (md + images → `exports/`; `--out`/`--zip`) |
| `npm run data:deploy` | Commit & push data branch changes; CI rebuilds and deploys |
| `npm run data:publish` | Push data branch only (no local commit) |

> Data-branch commands use the `data:` prefix; `translate`/`tagmerge`/`pack`/`deploy` are legacy aliases.

## 📝 Writing Content

Data lives in `../blog-data` with a three-layer model:

```
content-src/   layer 1 · hand-written Chinese sources (only human-edited area)
cache/         layer 2 · machine-maintained (English mirror / tag mapping / translate state, committed)
content/       layer 3 · regenerable final artifacts (JSON + HTML, gitignored)
```

- **Articles**: `content-src/notes/<category>/<subcategory>/<slug>.md` (flat, one md per article); directory names must match `categories.yaml`
- **Frontmatter**: `title` / `date` / `author` / `tags` / `draft` / `description`
- **Images**: put them in `content-src/assets/`, reference as `assets/xx.png` (the site resolves them directly, no bundling); use `data:pack` to export for sharing/migration
- **English**: never hand-written — `npm run data:translate` generates it; zh/en frontmatter `tags` must correspond one-to-one

### Publishing

```bash
npm run prebuild        # validate locally
npm run data:deploy     # commit & push data branch → GitHub Actions builds and deploys
```

## 🏗️ Project Structure

```
src/            browser app (Vue 3 + Pinia + vue-i18n; incl. SSR, no Node build scripts)
scripts/        Node content toolchain: lib/ (shared), generators/ (11 generators),
                translator/ (AI translation), tagMerger/ (tag mapping), packArticle.ts (export)
vite/           contentDevPlugin (dev watcher → regenerate → full-reload)
.github/        CI: checkout main + data → prebuild + build → deploy gh-pages
```

## 📄 License

MIT © Zorrooz — see [LICENSE](LICENSE).
