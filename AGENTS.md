# AGENTS.md — gblog

A static personal blog system built with Vue 3 + Vite 7 + Bootstrap 5. Markdown + YAML source files are processed at build-time into JSON metadata, loaded by a Vue SPA at runtime. Bilingual (zh-CN / en-US), GitHub Pages deployment.

**仓库分成三分支**：
- `main` — 纯代码（不含任何 content）
- `data` — 纯数据（content-src + cache），本地经 git worktree 挂到同级目录 `../blog-data`
- `gh-pages` — 构建产物（GitHub Actions 自动生成，勿手动提交）

## Quick Start

> **前置**：首次克隆后需挂数据 worktree —— `git worktree add ../blog-data data`（与仓库同级目录）。

| Command | Description |
|---------|-------------|
| `npm run dev` | Dev server at http://localhost:5173（数据目录由 contentDevPlugin 监听） |
| `npm run build` | Generators → vite-ssg build（SSG 预渲染 + PWA） |
| `npm run prebuild` | Content generators only |
| `npm run preview` | Preview production build |
| `npm run lint` | ESLint with auto-fix |
| `npm run typecheck` | vue-tsc (app) + tsc (node) 双项目类型检查 |
| `npm run format` | Prettier on `src/ scripts/ vite/` |
| `npm run data:translate` | AI translation CLI（增量，写 cache/en） |
| `npm run data:tag-merge` | zh→en 标签映射增量补齐 CLI |
| `npm run data:pack` | 打包文章（assets/ 引用图 → 文章同目录并改写引用；无参数 = 递归全部，幂等） |
| `npm run data:deploy` | 提交数据分支变更并推送（add -A → commit → push origin data，触发 CI） |
| `npm run data:publish` | 仅推送数据分支（`git -C ../blog-data push origin data`，触发 CI） |

> **命令规范**：数据分支操作统一 `data:` 前缀（kebab-case）；`translate`/`tagmerge`/`pack`/`deploy`
> 为旧名兼容别名（内部转调 `data:*`）。构建/质量命令用生态标准名（dev/build/prebuild/lint/typecheck/format/preview）。

**Always lint and build after changes.**
**发布流程**：编辑 `../blog-data/content-src/**` → `npm run prebuild` 本地验证 → `npm run data:deploy`（自动提交 data 分支）→ GitHub Actions 重新构建并部署到 gh-pages。

## Documents (约定文档索引)

| 文档 | 说明 |
|------|------|
| `.qoder/repowiki/` | **代码索引（Qoder 生成）**：项目 Wiki 文档树，按需加载上下文 |

### Repowiki 代码索引（按需加载）

`.qoder/repowiki/` 是 Qoder 自动生成的代码索引，含模块级知识卡与完整文档树。**需要了解某模块/子系统时，先按此加载上下文，再动手改代码**：

```
.qoder/repowiki/
├── knowledge/zh/_index.yaml        # 模块索引（modules: name → scope 文件清单/子模块）
├── knowledge/zh/<模块>/_module.yaml# 模块元数据（scope、依赖、关联）
├── knowledge/zh/<模块>/*.md        # 模块知识卡：概览 / 架构设计 / 技术栈 / 应用结构 / 开发规范
├── zh/content/                     # 完整文档树：架构、内容系统、指南、API 参考等
└── zh/meta/repowiki-metadata.json  # 文档间关系（parent-child 等）
```

使用规则：
1. 定位模块：先读 `knowledge/zh/_index.yaml`，按 `dir_name`/`scope` 找到目标模块
2. 加载上下文：读该模块知识卡与 `zh/content/` 对应章节，再读涉及的源码
3. 修改代码前，按需加载受影响模块的文档，避免凭猜测改动
4. 注意：该索引为**生成产物，可能滞后于代码**；内容与源码冲突时以源码为准，并在完成后可触发重新生成

## Architecture

```
src/                        # 浏览器应用（含 SSR；不包含任何 Node 构建脚本）
├── main.ts                    # Entry: ViteSSG → Pinia → Router → i18n；内联 v-reveal 指令定义
├── App.vue                   # <AppHeader> + <router-view> + <AppFooter>
├── router/index.ts           # 语言前缀路由（/{zh,en} × 5 路由）+ 旧 URL 重定向
├── i18n/
│   ├── index.ts              # vue-i18n 实例（createI18n，zh-CN/en-US）
│   └── locales/{zh-CN,en-US}.ts   # 也被 scripts/generators/generateCategories.ts 消费（跨层引用）
├── stores/                   # Pinia only
│   ├── theme.ts              # useThemeStore（theme 状态 + toggleTheme/initTheme）
│   └── locale.ts             # useLocaleStore（locale 状态）
├── config.ts                 # SITE、locale 映射（LOCALE_MAP/SEGMENT_OF/HTML_LANG）、THEME_MODES
├── views/
│   ├── Home.vue              # hero（greeting + name + stats）+ 标签云（模板内）+ PostList
│   ├── Category.vue          # Notes/Projects/Topics cards with stats
│   ├── Resource.vue          # Hierarchical resource catalog（侧栏 + 卡片）
│   ├── About.vue             # 头部（左介绍 + 右名片） + stats + 时间线 + sections
│   └── Article.vue           # NavigationTree + RenderMarkdown + OnThisPage
├── components/
│   ├── layout/               # AppHeader, NavActions, AppFooter, PostList, RenderMarkdown, NavigationTree, TreeNode, OnThisPage
│   ├── widgets/              # BackToTop, SearchModal, TocDrawer, FloatingButton（BackToTop/TocDrawer 共享浮动按钮）
│   └── common/               # IconButton（统一 stroke SVG 图标）、ModalOverlay（弹层基座）、PostCard、Pagination、EmptyState
├── composables/
│   ├── useLocalizedContent.ts# 内容页统一加载模式（首次加载 + locale 重载 + 异常回退）
│   ├── usePageMeta.ts        # 统一页面 SEO title（i18n 键 + BLOG_TITLE）
│   ├── useSearch.ts          # MiniSearch 全文检索（CJK 分词 + locale 切换重建索引）
│   ├── useTagNavigation.ts   # 标签筛选跳转（Home/PostList/Article 共用）
│   ├── useCopyFeedback.ts    # 复制成功/失败反馈（1.2s 复原）
│   └── useFloatingButton.ts  # 浮动按钮拖拽 + 底沿广播协议（FLOATING_BASE_EVENT 常量）
├── utils/                    # 浏览器运行时纯工具（不 import 任何 scripts/ 代码）
│   ├── contentLoader.ts      # Runtime: import.meta.glob for JSON & MD（@data 别名），强类型 load*
│   ├── navigation.ts         # 站内导航域：toLocalePath/switchLocale/goToTag/toArticle/normalizeArticleKey
│   ├── articles.ts           # flattenCategoryArticles（分类 → 文章列表展开）
│   ├── icons.ts              # 操作图标 path 数据单一来源（IconButton/iconSvg 共用）
│   ├── pagination.ts         # getVisiblePages 窗口化分页纯函数
│   ├── format.ts / tags.ts / url.ts  # 数字缩写 / 标签云 / URL·DOI 规范化纯函数
│   ├── scroll.ts             # scrollToTop/scrollToHeading（reduce-motion 感知）+ 弹层滚动锁定
│   ├── clipboard.ts          # copyText（Promise<boolean>，Clipboard API + execCommand 回退）
│   └── readingTime.ts        # 阅读时长估算
├── types.ts                  # 领域类型：与 content/*.json 产物一一对应（Post/Note/Tag/Category*/Resource/About）
└── assets/
    ├── fonts/                # 完整字体 OTF（思源宋体×2、思源黑体、Agave）
    ├── styles/               # global.scss + highlight/
    └── avatar.*              # 关于页头像（任意格式，存在即自动使用）
```

```
scripts/                    # Node 内容工具链（不被浏览器打包；Node 原生 type stripping 直跑 .ts）
├── lib/                     # 共享工具层：llm（配置/客户端/调用骨架）、fs（walk/walkAsync/JSON·YAML IO）、text、cli、frontmatter、tags、tagMapping、yamlEntries
├── dataConfig.ts           # 数据目录统一配置（唯一接入点，支持 GBLOG_DATA_DIR；EN_SUFFIX 常量）
├── markdownProcessor.ts    # unified pipeline: remark → rehype
├── runAllGenerators.ts     # Build orchestration（locale 步骤表 + 依赖顺序）
├── packArticle.ts          # 打包文章 CLI（assets 图 → 文章同目录）
├── llmConfig.ts            # DeepSeek API 配置；被 gitignore 忽略（API key，勿提交）
├── generators/             # core/（兼容 barrel，实现已迁移 scripts/lib/）+ 11 个生成器
├── translator/             # AI translation via DeepSeek API（CLI）
└── tagMerger/              # zh→en 标签映射增量工具（CLI）

vite/
└── contentDevPlugin.ts     # Vite dev 插件：监听数据目录变更→runAllGenerators→full-reload
```

### 数据分支（data）目录结构

代码 `src/` 不包含任何 content。内容数据全部位于数据分支（`../blog-data`），按**三层模型**组织：

```
content-src/                 # 第一层 src — 纯手写中文源（严格无机器产物，仅这里编辑）
│   ├── {categories,about,resources}.yaml
│   ├── assets/              #   配图集中目录（Obsidian 附件默认位置；发布前 npm run data:pack 打包进文章）
│   └── notes/<cat>/<sub>/<slug>.md     # 扁平：一篇文章一个 md（slug = 文件名，无同名目录）
cache/                       # 第二层 cache — 机器维护的持久态（入库），全部由工具生成
│   ├── en/                  #   英文翻译层：镜像 content-src 结构，文件名沿用 -en 身份后缀
│   │   │                       （categories-en.yaml、notes/.../<slug>-en.md）
│   │   └── …                #   -en 后缀是「内容身份」（URL 的一部分），cache/en 目录表示机器层
│   ├── tag-mapping.json     #   zh→en 标签映射 / 标签合并映射
│   └── .translate-state.json#   翻译增量状态（源相对路径 → 内容 hash）
content/                      # 第三层 final — 生成 JSON + html（可再生，gitignore，不入库）
```

- **中文源**：只写 `content-src/`。**英文**：绝不手写，`npm run data:translate` 生成到 `cache/en/`（LLM 非确定 + 有成本，故入库存证，CI 靠 `GBLOG_NO_TRANSLATE=1` 跳过生成）。
- 生成器按 locale 取源：`srcDirFor(locale)`（zh→content-src，en→cache/en），输出恒为 `content/*` 与 `content/*-en*`。
- 中英实体按镜像相对路径配对：`cache/en/notes/.../<slug>-en.md` ↔ `content-src/notes/.../<slug>.md`。

## Framework

| Layer | Tech |
|-------|------|
| Vue | Vue 3 Composition API (`<script>`, Options API mix) |
| Build | Vite 7 (ESM, `"type": "module"`) |
| SSR / SSG | vite-ssg（构建期预渲染每页；`concurrency: 1` 串行保 locale 不被并发覆盖） |
| State | Pinia |
| Routing | Vue Router 4（history 模式，由 vite-ssg 创建；`/{zh,en}` 语言前缀 + 旧 URL 重定向） |
| PWA | vite-plugin-pwa（`autoUpdate` + `generateSW`，precache 内容产物） |
| i18n | vue-i18n (locale in localStorage) |
| Search | MiniSearch（客户端全文检索，索引由 generateSearchIndex 预生成） |
| Markdown | unified → remark-parse → remark-gfm → remark-breaks → remark-math → remark-rehype → rehype-highlight → rehype-katex → rehype-stringify |
| Deployment | gh-pages |
| Node | `>=23.6.0`（Node 原生 type stripping 直跑 `.ts` 脚本，最低 22.6 需 `--experimental-strip-types`） |

## Content Pipeline (Critical Order)

```
1. generateNotes     — reads srcDirFor(locale)/notes/**/*.md → content/notes.json
2. generateProjects  — reads `projects` section of categories.yaml → content/projects.json
3. generateTopics    — reads `topics` section of categories.yaml → content/topics.json
4. generateCategories— YAML + notes/projects/topics → content/categories.json
5. generatePosts     — notes.json + categories.json → content/posts.json
6. generateTags      — posts.json → content/tags.json
(1-6 for zh-CN, then repeat for en-US)

4.5 incremental translate — 内容 hash 增量翻译（仅缺失/变化的 zh 源 → cache/en 镜像 -en 文件；tag 经 zh→en 映射查表翻译；GBLOG_NO_TRANSLATE=1 关闭）
4.6 tagMerger mapping     — ensureTagTranslation 增量补齐 zh→en 标签映射（缓存于 cache/tag-mapping.json，命中 0 token）
4.7 tag consistency fix   — 以 zh 文件为基准自动重写 cache/en 的 -en 文件 tags（解决中英标签数量/名称不一致，失败仅告警）
8.5 tags-consistency check— 构建产物级中英标签一致性校验（输出 [OK]/[Warn]）
9.  generateHtml     — md → `content/html/**`（markdownProcessor 预渲染文章 HTML）
10. generateSitemap  — `public/sitemap.xml` + `public/robots.txt`（gitignore，CI 自动生成）
11. generateSearchIndex — categories.json + content/html → `content/search-index{,-en}.json`

generateAbout and generateResources run independently at import time in runAllGenerators.ts.
```

> 所有路径基于数据分支（`../blog-data`，CI 中为 `$GBLOG_DATA_DIR`）；第一层 content-src 纯手写，第二层 cache（含 en/ 翻译层）入库存证，第三层 content/ 可再生产物不入库。

## Content Source Organization

**只有 notes 有 markdown 文章**；projects/topics 为纯 yaml 元数据（无 md 文件、无文章页，卡片外链到 GitHub/DOI）。

```
content-src/notes/<category>/<subcategory>/<slug>.md        # 扁平：一篇文章一个 md
cache/en/notes/<category>/<subcategory>/<slug>-en.md       # 英文镜像
```

- notes 目录名必须匹配 `categories.yaml` 中 `notes` 段的 `name` 字段；子分类目录名匹配其 `categories` 映射 key。
- **图片**：写作期统一放 `content-src/assets/`（引用 `assets/xx.png`，Obsidian 最短路径格式，站点已支持解析）；发布前 `npm run data:pack -- notes/<cat>/<sub>/<slug>`（或 `npm run data:pack` 递归全部）将引用图片复制到文章同目录并改写引用（同步改写 cache/en 镜像，幂等）。
- 每篇文章 md 自带 **YAML frontmatter**：`title` / `date` / `author` / `tags` / `draft` / `description`（zh/en 各自文件分别定义）
- `categories.yaml` 定义分类层级（`name`/`title`/`desc`/`categories` 子分类映射）；projects/topics 的全部元数据（`github`/`doi`/`url`/`status`/`language`/`license`/`journal`/`year`/`authors` 等）也在此定义
- **中英 frontmatter `tags` 必须一一对应**（同一篇文章 zh/en 标签数量与语义对称），否则首页标签云双语标签数不一致（如 zh `Shell` / en 误写 `Shell`+`Bash`）

### URL Mapping

路由为 history 模式（`/zh`/`/en` 语言前缀，vite-ssg 构建期真实预渲染每页）：

```
File:  notes/Omics/genomics/bwa.md
URL:   /zh/article/notes/Omics/genomics/bwa        (en → /en/article/notes/Omics/genomics/bwa-en)
Route: { path: '/:locale/article/:path*', name: 'zh-Article'|'en-Article', props: true }
       (zh-*/en-* 前缀命名，其余为 zh-Home/en-Home 等)
旧 URL /article/...（无前缀）→ 302 重定向到 preferredLocaleSegment() 前缀
旧「每文一目录」URL（.../bwa/bwa 或 .../bwa/bwa-en）→ Article 路由 beforeEnter 自动重定向到新 URL
```

## Internationalization (Dual Layer)

**Layer 1 — App:** `src/i18n/locales/{zh-CN,en-US}.ts` via vue-i18n. Controls nav, buttons, labels.
语言包结构以 zh-CN 为基准（`src/i18n/schema.ts` 的 `AppMessages`），en-US 用 `satisfies` 做编译期键校验；
SEO title 统一走 `usePageMeta`（i18n meta.* 键 + `BLOG_TITLE`）。

**Layer 2 — Content:** locale 源目录分层（`srcDirFor(locale)`）+ `-en` 后缀身份：
- Chinese: `content-src/article.md` → 输出 `content/categories.json`
- English: `cache/en/article-en.md` → 输出 `content/categories-en.json`
- Auto-fallback in `contentLoader.ts`: if English file missing, loads Chinese.

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
1. 创建 `../blog-data/content-src/notes/<cat>/<sub>/<slug>.md`（slug 英文小写连字符，含 YAML frontmatter；Obsidian 中可用 `templates/new-post.md` 模板）
2. 配图先放 `content-src/assets/`，发布前 `npm run data:pack -- notes/<cat>/<sub>/<slug>` 打包
3. 在 `categories.yaml` 添加分类
4. 运行 `npm run prebuild` 本地验证
5. 英文：`npm run data:translate`（生成 cache/en，不手写）
6. 发布：`npm run data:deploy`（提交并推送 data 分支 → CI 自动构建部署）

### Utilities
- `contentLoader.ts` provides: `loadPosts()`, `loadCategories()`, `loadNotes()`, `loadTags()`, `loadAbout()`, `loadResources()`, `loadMarkdownContent()`；locale 判定统一走 `src/locale.ts` 的 `currentLocale()`
- `navigation.ts`：`toLocalePath()`、`switchLocale()`、`goToTag()`、`toArticle()`、`normalizeArticleKey()`（去 .md/-en 后比较）、`articlePathFromUrl()`、`toArticleRoutePath()`、`joinRoutePathParam()`
- `articles.ts`：`flattenCategoryArticles()`（分类 → 文章列表展开）；`icons.ts`：操作图标 path 单一来源（`ICON_MARKS`/`iconSvg`）
- `pagination.ts`：`getVisiblePages()` 窗口化分页纯函数；`format.ts`/`tags.ts`/`url.ts`：数字缩写/标签云/URL·DOI 规范化
- `scroll.ts`：`scrollToTop()`、`scrollToHeading(el, offset)`（reduce-motion 感知 + 锁延迟）、`isScrollLocked()`、`lock/unlockScrollOverflow`、`lock/unlockScrollPosition`；`clipboard.ts`：`copyText(): Promise<boolean>`
- composables：`useLocalizedContent`（内容页加载模式）、`usePageMeta`（SEO title）、`useSearch`（MiniSearch + locale 重建）、`useTagNavigation`、`useCopyFeedback`、`useFloatingButton`（`FLOATING_BASE_EVENT` 常量）
- `markdownProcessor.ts` export: `renderMarkdown(markdown)` → HTML string（构建期专用，客户端勿 import）
- `useThemeStore`: `theme: ThemeMode`, `toggleTheme()`, `initTheme()`（**必须以 localStorage 为准覆盖**：SSR 预渲染内嵌的 pinia state 恒为 `auto`，hydrate 会覆盖客户端初始值，`initTheme` 里重新读 `initialTheme()` 才能持久化）；`useLocaleStore`: `locale`, `setLocale()`, `initLocale()`

## Layout & Shape Conventions（布局与形状规则）

### 内容宽度

- 内容页（首页/分类/资源/关于）：`.page-section` max-width **1120px**，左右 padding 48px（`<992px` 32px，`<768px` 20px）
- 文章页：Bootstrap `container` max-width **1280px**，正文 `.article-content` max-width 820px 居中，左右 sidebar 各约 220px
- **移动端水平边距全局统一 20px**（header / footer / 内容页 / 文章页 gutter），不得各页面自定
- **显示缩放断点**：`<1200px`（含 175% 缩放等效视口）header `.container` 放宽为 `max-width: 100%` 贴边；footer 在大缩放下内容居中排列
- 原因：文章页需容纳三栏（导航/正文/TOC）故放宽；其余页面单栏 1120px 保证行宽

### 圆角分层（controls 6 / buttons 8 / cards 12 / panels 16 / pills full）

| 层级 | 值 | 用途 |
|------|-----|------|
| `--radius-sm` 6px | 树形导航项、输入框、kbd、search-result | 小控件 |
| `--radius-btn` 8px | 文本/数字按钮（page-btn、icon-btn、chip-close） | 按钮 |
| `--radius` 12px | 卡片（cat-card、res-card、offcanvas-card、post 卡片） | 容器 |
| `--radius-lg` 16px | 大面板（search-panel）、头像块 | 强调容器 |
| `--radius-pill` 9999px | 标签 pill、圆形图标按钮 | 全圆角 |

**规则**：圆形（`border-radius: 50%`）仅用于图标操作按钮（footer 外链、BackToTop、TocDrawer、cat-card/res-card 外链、page-btn--nav）；方形按钮一律 8px。同一容器内不可混用两种按钮形状。

### 页面 header 模式

- **列表页标题**（首页 posts-header、分类页 section header）：26px serif `section-title` + hairline（首页 posts-header 无 indicator；分类页 section header 用 tint 图标块）+ hairline 下边框
- **页面级 H1**（分类/资源/关于页头、文章标题）：`article-title` 风格（clamp(32-46px) serif，`text-wrap: balance`）
- 规则：一页只出现一种 header 组合；列表分组标题用「小号大写标签 + hairline」或「带 indicator 的 section-title」，不重复堆叠 H1 样式

### 侧栏导航与 TOC（VitePress/D3 风格）

- **分组标题**（NavigationTree/资源页侧栏）：13px 加粗正文（非大写小标签），组间用 hairline 分隔
- **条目**：13px 纯文字链接（无圆角底块），hover 变主色；激活项 = 主色 + 左侧 2px 主色指示条（**不加粗**）
- **可折叠目录**：`tree-item--folder`（button）+ 右侧 `fa-chevron-right` caret（展开时 rotate 90°）；**默认全部展开**，仅记录用户手动折叠（`collapsedIds`），路由切换不重置
- **子列表导轨**：`.tree-sublist`/`.res-group__items` 用 `border-left: 1px solid var(--line)` + `padding-left: 16px`；激活指示条 `left: -17px` 压在导轨线上
- **TOC（OnThisPage）**：标题 en `On this page` / zh `本页目录`；`.otp-content` 左侧发丝导轨 + `.otp-marker`（2px 主色条，`top` 平滑过渡跟随激活项，公式 `link.offsetTop + 36`）；链接 13px **换行显示**（不截断省略），hover/激活变 `--fg`，不做圆角底块
- 分类页计数（`.cat-section__count`）为纯文本 meta（等宽数字），不用 pill

### 卡片与 hover（首页 editorial 语言）

- 卡片（cat-card / res-card / about-cell / article-nav-item / Bootstrap `.card`）：hover = 仅边框变主色微光，**无背景变化、无阴影**；浮动层（抽屉/弹窗/浮动按钮）才用阴影
- 首页 post-item：无边框卡片，行间 hairline 分隔 + 序号（等宽灰）+ hover 标题变主色、箭头 `translateX(3px)` 微位移（基准风格，其他页面向其看齐）

### 动效与无障碍

- **v-reveal 入场**（`main.ts` 指令 + `global.scss` keyframes）：`opacity 0→1 + translateY(14px)→0 + blur(5px)→0`，`--reveal-delay` 做 stagger；`prefers-reduced-motion` 用户由指令直接跳过（不挂类）。
- **skip link**：`App.vue` 的 `.skip-link`（Tab 聚焦时显示，`main#main-content` 为目标）；新增全局可访问入口时保留该约定。
- 所有交互元素：hover/active/focus-visible 三态齐全；过渡 140-220ms（`--dur-fast`/`--dur-base`）；滚动平滑走 `scrollToHeading`（reduce-motion 感知）。

### 404

- 语言前缀下的未知路径 → `NotFound.vue`（品牌 404：巨型 serif 数字 + 说明 + 返回首页 pill）；无前缀未知路径先补 locale 前缀再落入 404；`dist/404.html` 由构建 onFinished 复制 index.html 生成
- **返回按钮 `flex: none` 不可移除**：App.vue 的 `.main-content > * > * { flex: 1 }` 链式布局会拉伸 404 短页面的最后一个子元素（曾出现 200px 高的大椭圆 pill）

### 字体

思源宋体（标题）、思源黑体（时间线中文/代码块中文注释）与 Agave（等宽）均以**完整字体文件**直接引用（`src/assets/fonts/*.otf`，`@font-face` + `font-display: swap`），**无子集化、无任何构建脚本**。新增内容不影响字体渲染。

### 图标

- **操作图标**（header search/theme/language/menu、复制按钮）：内联 stroke SVG——`viewBox 24`、`stroke-width 1.75`、`fill none`、`stroke currentColor`、渲染 18px，颜色随 `--fg-3` → hover `--primary` 自动变化。**禁止 PNG 位图图标**。
- **内容/装饰图标**（footer 外链、时间线、section 卡片）：FontAwesome（`fa-solid` / `fa-brands`），类名由数据驱动（`about.yaml` 的 `icon` 字段、页面模板硬编码类名），不依赖任何本地图标文件。

### 品牌与头像

- 品牌为**文字 Logo**：`zorrooz's blog`（`.app-nav__brand`，`--font-serif` 宋体），无图形 logo。
- favicon 与 PWA 图标：`public/{favicon,apple-touch-icon,icon-512}.png`，由 `index.html`（favicon/apple-touch + og:image）与 `vite.config.ts` 的 PWA manifest（64/180/512 + maskable）引用。
- 关于页头像：`src/assets/avatar.*`（任意图片格式），`About.vue` 通过 `import.meta.glob('../../assets/avatar.*')` 自动接入（有图用图、无图回退字母 monogram）；头像框为**黑色线框**（`1.5px solid var(--ink)`，无底色）。

### 复制功能

- **复制文章**：Article.vue 头部按钮，克隆 `.markdown-body` 清理辅助元素后取 `innerText`（标题 + 正文纯文本）。
- **复制表格**：RenderMarkdown 对每个 `<table>` 包 `.table-block-wrapper` + header + 复制按钮，复制为 **TSV**（`\t` 连接单元格，可直接粘贴 Excel）。

## Color System — 三色映射规范

gblog 仅用 **3 色族**：蓝色（主色）、灰色（中性色）、黑/白（基底）。所有语义 token 必须从这三族派生，禁止硬编码 rgba。

### 色板

| Token | 亮色值 | 暗色值 | 用途 |
|-------|--------|--------|------|
| `blue-600` | `#3b82f6` | `#6cb2ff` | **主色**：链接、active、按钮、焦点环、进度条 |
| `blue-500` | `#5aa0f8` | `#4d9fff` | 辅助蓝：选中背景、轻量强调（非交互主色） |
| `blue-200` | `#e6f1ff` | `#1e3352` | 强 tint：hover 背景加深 |
| `blue-100` | `#f2f8ff` | `#16263d` | 弱 tint：标签/卡片 hover 底、selection 底 |
| `gray-50` | `#f5f6f8` | `#1a202b` | surface-2：次级表面、灰底 |
| `gray-100` | `#e4e7ec` | `#232936` | 分割线 |
| `gray-200` | `#d0d5dd` | `#323a4a` | 强分割线、序号 |
| `gray-400` | `#6f7684` | `#8e96a6` | 弱文字（meta、caption） |
| `gray-500` | `#4d5563` | `#a2a9b8` | 次文字 |
| `ink` | `#16191f` | `#e7e9ee` | 正文/标题（近黑/近白） |

### 语义映射规则

```
--primary        = blue-600   (#3b82f6)   主交互色（链接、active、按钮填充）
--primary-hover  = blue-700   (#2563eb)   hover 加深（比 primary 暗一档）
--primary-active = blue-800   (#1d4ed8)   active/click 最暗
--primary-light  = blue-500   (#5aa0f8)   轻量强调（选中背景、非核心高亮）
--tint           = blue-100   浅底（hover 背景、tag 底、selection 底）
--tint-strong    = blue-200   中浅底（hover 背景加深）
```

### 映射禁忌

1. **主色反转**：`--primary` 必须是最深的交互蓝（blue-600），不能是浅蓝（blue-500）。hover 在 primary 基础上再暗一档。
2. **硬编码 rgba**：禁止 `rgba(90, 160, 248, x)` 或 `rgba(59, 130, 246, x)`，必须用 `var(--primary)` + opacity 或 `color-mix()`。
3. **灰色跳跃**：灰色 scale 必须等距递进（50→100→200→400→500），不能跳级。
4. **暗色 primary 亮度反转**：暗色模式下 primary 仍是最亮蓝（blue-600 = `#6cb2ff`），hover 更亮（`#93c5fd`），保持「hover = 更强」语义。

### 组件颜色分配

| 组件区域 | 颜色 Token | 说明 |
|----------|-----------|------|
| 正文/标题 | `--fg`（ink） | 最高对比 |
| 次要文字 | `--fg-2`（gray-500） | 描述、摘要 |
| 弱文字/meta | `--fg-3`（gray-400） | 日期、分类、caption |
| 链接默认 | `--primary`（blue-600） | 正文内链 |
| 链接 hover | `--primary-hover` | 加深一档 |
| Active 导航 | `--primary` + `font-weight: 600` | 当前页高亮 |
| 按钮填充 | `--primary` | CTA、分页 active |
| 按钮 hover | `--primary-hover` | 深一档 |
| 标签默认 | `--surface-2` + `--fg-2` | 灰底灰字 |
| 标签 hover | `--tint` + `--primary` | 蓝底蓝字 |
| 卡片边框 | `--line`（gray-100） | 默认分割 |
| 卡片 hover 边框 | `rgba(var(--primary), 0.3)` | 蓝色微光 |
| 分割线 | `--line` | hairline |
| 表面背景 | `--surface`（白/`#141821`） | 卡片、面板底 |
| 次级表面 | `--surface-2`（gray-50） | 灰底区域 |
| 代码行内 | `--primary-light`（blue-500） | `color` 属性 |
| 代码块背景 | `--app-markdown-code-bg` | 灰底 |
| 引用块边框 | `--primary`（3px left border） | 标识性 |
| 焦点环 | `--primary`（2px outline） | 无障碍 |
| 选中背景 | `--tint`（blue-100 + 16% opacity） | 文字选中 |

### 添加新组件时

1. 文字色：优先用 `--fg` / `--fg-2` / `--fg-3`，不自造色值
2. 交互色：链接/按钮/active 一律用 `--primary`
3. 背景层次：`--surface` → `--surface-2` → `--tint`，三级递进
4. 边框：`--line`（默认） → `--line-strong`（强调）
5. hover 效果：边框用 `rgba(color, 0.3)`，背景用 `--tint`

## ESLint / Prettier
- Flat config (`eslint.config.js`) with Vue + Prettier
- TS parsing: `@typescript-eslint/parser` for `**/*.ts` + vue files with `lang="ts"`
- `browser` globals for `*.{js,vue}`, `node` globals for `scripts/**/*.{js,ts}` + `vite/**/*.ts`
- Ignored: `dist/`, `dist-ssr/`, `coverage/`
