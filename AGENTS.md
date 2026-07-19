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
