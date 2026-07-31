# gblog 架构现代化执行计划

> 状态：**执行中（Phase 0-1 已完成）**
> 制定日期：2026-07-31
> 路由来源：AGENTS.md → 本文件（约定文档）
> 执行方式：由用户切换到 build 模式，按 Phase 顺序逐个执行；每个 Phase 完成后更新本文件顶部的状态表

## 阶段状态表

| Phase | 内容 | 状态 |
|-------|------|------|
| Phase 0 | 基建：dev 热更新 + SSR 安全 + 高亮本地化 | ✅ 已完成 |
| Phase 1 | Markdown 渲染移入构建期（bundle 瘦身） | ✅ 已完成 |
| Phase 2 | SSG 预渲染 + history 路由（核心升级） | ⬜ 待执行 |
| Phase 3 | 内容层整理（生成器共享 core + sitemap） | ⬜ 待执行 |
| Phase 4 | 组件现代化（`<script setup>` 迁移） | ⬜ 待执行 |
| Phase 5 | TypeScript 渐进迁移 | ⬜ 待执行 |
| Phase 6 | 体验增强（PWA / 搜索 / 资产本地化，可选） | ⬜ 待执行 |

---

## 0. 目标与约束（红线，任何 Phase 不得违反）

1. **现有功能与样式不能被破坏**（允许小幅度提升）：
   - 全部页面路由、中英切换、主题三态切换（auto/light/dark）、文章 prev/next、TOC（OnThisPage + TocDrawer）、代码复制按钮、heading 锚点、sticky 双栏、移动端断点布局，行为必须与改造前一致。
   - 不修改 `content-src/` 任何源文件；不修改 `global.scss` 的 CSS 变量体系。
2. **文件组织与代码可读性不能下降**：
   - 生成器一律放 `src/utils/generators/`；生成物放 `src/content/`；加载逻辑收口在 `src/utils/contentLoader.js`；构建工具插件放 `src/utils/` 下。
   - 沿用既有代码风格：无分号、单引号、2 空格缩进、ESM、命名导出优先、无多余注释。
3. **不引入额外后端**：一切增强都必须落在构建期或纯静态运行时。
4. **每阶段验收流程（强制）**：
   ```
   npm run lint
   npm run build
   npm run preview
   ```
   然后浏览器手动回归清单：
   - [ ] 首页文章列表 / 分类页 / 资源页 / 关于页正常
   - [ ] 主题 auto→light→dark 循环正常，代码高亮随之切换
   - [ ] 中英切换正常（UI + 文章内容 + 回退逻辑）
   - [ ] 文章页：正文渲染、标题锚点、代码复制、TOC、prev/next、sticky 侧栏
   - [ ] 窄屏（≤991px）侧栏降级为流式布局
   - [ ] 直接刷新深层链接不 404（Phase 2 起）
5. 每个 Phase 完成后：更新本文件「阶段状态表」+ 在 git commit message 中标注 Phase 编号。

---

## 1. 现状事实（证据，2026-07-31 核查）

| 事实 | 证据位置 |
|------|----------|
| 纯 SPA + hash 路由，无预渲染 | `src/router/index.js:9`（`createWebHashHistory`） |
| Markdown 在浏览器端用 unified 全链路渲染 | `src/utils/markdownProcessor.js:51-64` |
| 高亮主题走 CDN + MutationObserver 动态换样式 | `src/utils/markdownProcessor.js:15-50` |
| 内容加载：JSON 为 eager glob、markdown 为 lazy raw glob | `src/utils/contentLoader.js:17-22` |
| 生成器为 8 个手写 Node 脚本，dev 模式不监听 content-src | `src/utils/runAllGenerators.js`，`npm run prebuild` 手动触发 |
| 组件为 Options API + setup() 混搭 | `src/views/Article.vue:70-377` |
| SSR 不安全：顶层访问 localStorage/document | `src/stores/i18n.js:5-6`、`src/stores/theme.js`、`src/stores/app.js:8,11` |
| `src/content/` 生成 JSON 已入库（16 个文件） | `git ls-files src/content` |
| 内容规模：zh 13 篇 + en 13 篇，content-src 无图片 | 2026-07-31 统计 |
| Font Awesome 6.4 走 CDN | `index.html:8-11` |
| 部署：gh-pages 推送 dist；user site → base `/` | `package.json` `deploy` 脚本 |

---

## 2. 总体路线图

```
Phase 0（基建）→ Phase 1（构建期渲染）→ Phase 2（SSG/路由）→ Phase 3（内容层）→ Phase 4（组件）→ Phase 5（TS）→ Phase 6（可选）
```

依赖关系：Phase 1 是 Phase 2 的前提（预渲染需同步拿到文章 HTML，无需 Suspense/asyncData）；Phase 0b（SSR 安全）是 Phase 2 的前提。

---

## Phase 0 — 基建：dev 热更新 + SSR 安全 + 高亮本地化

**目标**：开发体验提升；为 Phase 2 的 SSR 预渲染扫清 DOM/localStorage 依赖；消除运行时 CDN 高亮主题。

**执行记录（2026-07-31）**：
- 0a 偏差：chokidar v4 已移除 glob 支持，`contentDev()` 改用 `watch('src/content-src')` 目录递归监听（新增/删除文件仍走 `server.restart()`）
- 0b 附项：`runAllGenerators.js` 重构为导出 `runAllGenerators()` 供插件复用（`isDirectRun` 守卫保留 CLI 直跑）；generateAbout/generateResources 加同款守卫
- 附项：修复基线 lint 22 处错误（未用变量/保留键 `_otp*`/`_lockedScrollY` 改名、可选 catch、稀疏数组、vite.config.js 补 node globals、`Article.vue:80` 未用 glob 删除）——验收流程强制 lint 通过所需，均为机械改动

### 0a. Dev 内容热更新插件
- 新建 `src/utils/contentDevPlugin.js`：导出 Vite 插件 `contentDev()`
  - `configureServer(server)`：启动时跑一次 `runAllGenerators()`；用 `chokidar`（vite 自带依赖，无需新增）监听 `src/content-src/**`，防抖 300ms 后重跑生成器并 `server.ws.send('full-reload')`
- `vite.config.js`：
  - `plugins` 增加 `contentDev()`（dev 专属，内部用 `apply: 'serve'` 或 `configureServer` 天然隔离）
  - `server: { watch: { ignored: ['**/src/content/**'] } }`（防止生成器写文件触发 reload 死循环）
- 注意：生成器可能产出新 HTML 文件（Phase 1 后），lazy glob 在 full-reload 后会重新扫描；若新文件未被拾取，回退方案为 `server.restart()`

### 0b. SSR 安全改造（行为不变，纯守卫）
- `src/stores/i18n.js:5-6`：`localStorage` / `document.documentElement.lang` 改为 `typeof window` 守卫的惰性初始化
- `src/stores/theme.js`：`setTheme` 内访问 `document/window/localStorage` 加守卫；SSR 下直接返回
- `src/stores/app.js:8,11`：theme/locale 初始值改为惰性读取（守卫后取值）
- 核查其余视图/组件顶层是否存在 `window/document/localStorage`（已知 `Article.vue:100` 已有守卫）

### 0c. highlight.js 主题本地化
- 从 `node_modules/highlight.js/styles/` 复制：
  - `github.css` → `src/assets/styles/highlight/github.css`（原样，作为默认 light）
  - `github-dark-dimmed.css` → `src/assets/styles/highlight/github-dark-dimmed.css`（**所有选择器加 `[data-bs-theme='dark']` 前缀**）
- `src/main.js` 在 `global.scss` 之后引入两文件
- 删除 `src/utils/markdownProcessor.js:15-50` 的 `loadHighlightStyle` + MutationObserver 块；文件顶部注释改为「构建期专用（Node 环境），勿在客户端 import」
- 保留 `RenderMarkdown.vue` 中 `.markdown-body .hljs { background: var(--app-markdown-code-bg) }` 覆盖规则

**验证**：build + preview；黑暗/亮色主题下代码块配色与改造前一致；禁用网络后主题切换仍生效（无 CDN 请求）。

---

## Phase 1 — Markdown 渲染移入构建期

**目标**：unified/remark/rehype/lowlight 从客户端 bundle 中移除（大幅瘦身）；文章 HTML 成为静态生成物，为 Phase 2 预渲染铺路。

### 1a. 新生成器 `src/utils/generators/generateHtml.js`
- 复用 `src/utils/markdownProcessor.js` 的 `renderMarkdown()`（Node 端运行 unified 全链路：GFM、数学、KaTeX、代码高亮）
- 文件枚举规则与 `generateNotes.js:145-151` 一致：zh 取全部 `.md`（排除 `-en.md`）、en 仅取 `-en.md`
- 输出：`src/content/html/<相对路径>.html`（如 `html/notes/Omics/genomics/bwa/bwa.html`、`html/notes/.../bwa-en.html`）
- `.gitignore` 追加 `src/content/html/`（不进库；JSON 入库策略不变）
- `src/utils/runAllGenerators.js` 在 tags 之后追加 zh/en 两轮 `generateHtml`

### 1b. 加载层切换 `src/utils/contentLoader.js`
- `markdownModules` glob 改为 `../content/html/**/*.html`（`query: '?raw'`、lazy）
- 新增 `loadHtmlContent(filePath)`：`-en` 优先 + 中文回退（回退语义与现 `loadMarkdownContent` 完全一致）
- 删除旧 `loadMarkdownContent()`（含三路 possiblePaths 逻辑，净简化）

### 1c. 渲染组件 `src/components/layout/RenderMarkdown.vue`
- 移除 `import { renderMarkdown } from '@/utils/markdownProcessor'`
- props 语义改为接收**已渲染 HTML**（prop 名可保留 `rawMarkdown` 以减少改动面，实现内部处理 HTML 字符串）
- `rewriteImageLinks`：删除 `![..](..)` markdown 级分支，保留 `<img src>` 分支（对 HTML 字符串同样机制；`@` 解析与 passthrough 规则原样保留）
- `enhanceHeadings` / `enhanceCodeBlocks` / 复制按钮逻辑不动（纯 DOM 增强，与是否预渲染无关）

### 1d. `src/views/Article.vue` 适配
- `loadArticleContent()` 改为同步 map 查询（HTML 全部进 bundle，无需 `await`；签名保留 async 以兼容调用方）
- 阅读时长：`readingMinutes` 改为对 HTML 去标签文本 `text.length / 800`（与现指标同源，允许的微小差异）
- 删除 `Article.vue:80` 重复声明的 md glob

**验证**：build 产物对比 bundle 体积显著下降；全部 13 篇中文文章正文渲染逐篇与改造前对照（含代码块、表格、数学公式、标题层级）。

---

## Phase 1 — Markdown 渲染移入构建期

**目标**：unified/remark/rehype/lowlight 从客户端 bundle 中移除（大幅瘦身）；文章 HTML 成为静态生成物，为 Phase 2 预渲染铺路。

**执行记录（2026-07-31）**：
- 1a-1d 按计划完成；主 bundle 从 1,047 kB → 302 kB（gzip 323 → 98 kB），KaTeX 字体文件不再打包
- 源文章无数学公式（`$` 均为 shell 变量），katex 差异未触发；`katex.min.css` 保留（防未来公式内容）
- `loadHtmlContent()` 回退语义与旧 `loadMarkdownContent` 对齐（en 优先 `-en.html` 再中文，zh 按路径原样）

---

## Phase 2 — SSG 预渲染 + history 路由（核心升级）

**目标**：真实 URL（SEO/分享友好）、每篇文章独立静态 HTML、首屏加速。

### 2a. 依赖与 spike 验证（先验证再动代码）
- 安装 `vite-ssg` + `@unhead/vue`
- 先做最小 spike：临时最小 `main.js` 跑通 `vite-ssg build`，确认与 Vite 7 / Vue 3.5 兼容（peer deps + 实际构建）
- **回退预案（若 vite-ssg 不兼容）**：自写 `src/entry-server.js`（约 80 行），用 Vite `ssrLoadModule` + `@vue/server-renderer` 逐路由 `renderToString`；方案 B 仍为纯静态，不引入后端

### 2b. 路由切换 `src/router/index.js`
- `createWebHashHistory` → `createWebHistory(import.meta.env.BASE_URL)`
- 路由表结构不变（`/article/:path*` 通配保留）

### 2c. 入口改造 `src/main.js`
- 改为 `export const createApp = ViteSSG(App, { routes }, ({ app, router, initialState, isClient }) => {...})` 模式
- Pinia：按 vite-ssg 官方文档 initialState 模式（SSR 时 `initialState.pinia = pinia.state.value`，客户端恢复）
- i18n：SSR 固定 `zh-CN`（`isClient` 时再按 localStorage 恢复）——预渲染中文版，英文仍运行时切换（已确认取舍，双路径 i18n 留待 Phase 6）
- `package.json`：`build` 改为 `vite-ssg build`（`prebuild` 生命周期钩子不变，生成器照跑）

### 2d. 动态文章路由枚举
- 在 `src/main.js` 导出 `includedRoutes(paths)`（或 vite.config `ssgOptions.includedRoutes`）：读取生成期 `src/content/categories.json`（zh），展开全部文章路径 → `/article/<type>/<group>/...`
- 文章 HTML 已在 Phase 1 生成 → 预渲染时 Article.vue 同步取得内容，无需 Suspense/asyncData

### 2e. SEO head（顺带收益）
- `@unhead/vue` `useHead`：首页 / 分类页 / 文章页设置 `title` + `description`（来自 JSON 元数据）；`index.html` 补基础 meta（description、og: 基础项）

### 2f. GH Pages 兼容
- 构建后复制 `index.html` → `dist/404.html`（history 模式 SPA fallback；user site base `/` 无需改）
- 实现方式：`ssgOptions.onFinished` 回调或 postbuild 脚本（优先 ssgOptions 内完成，不新增顶层脚本）

**验证**：
- `npm run build` 后 `dist/article/**/*.html` 存在且正文内容已内联（非空壳）
- `npm run preview` 直接刷新深链接（如 `/article/notes/Omics/genomics/bwa/bwa`）不 404
- 全部导航链路、语言切换、主题切换回归

---

## Phase 3 — 内容层整理（不改行为，只提可读性）

- 新建 `src/utils/generators/core/` 共享模块，从 `generateNotes.js` 抽取：
  - `walk`（目录遍历）、`normalizeTags`、`parseFrontMatterAndBody`、`markdownToPlain`、`countWordsSmart`
- 核查 `generateProjects.js` / `generateTopics.js` 中的重复实现并收敛引用（先确认行为一致再替换，逐文件对照输出 JSON）
- 新增 `src/utils/generators/generateSitemap.js`：读 `posts.json`（zh）→ 写 `public/sitemap.xml` + `public/robots.txt`（纯静态资源，vite 自动拷贝）
- `runAllGenerators.js` 追加 sitemap 步骤

**验证**：生成 JSON 与改造前 diff 完全一致（除时间戳类字段）；build + preview 正常。

---

## Phase 4 — 组件现代化（机械迁移，逐个验证）

- 全部视图/组件迁移 `<script setup>`；`Article.vue` 优先（700 行最大，先啃硬骨头）
- 迁移语义对照：props → `defineProps`、emits → `defineEmits`、computed/watch/lifecycle 一一对应；`setup()` 返回 `t/locale` 的混搭写法移除（改 `useI18n` 直接用）
- Vue 3.5 特性落地：`useTemplateRef`（替换手动 `this.$refs`）、`onWatcherCleanup`（清 resize/scroll 监听）
- **不做** unplugin-vue-router（手动路由表已足够清晰，避免过度工具化）

**验证**：每迁移一个组件即 `npm run lint` + 浏览器回归该组件涉及的页面；全部完成后完整回归清单。

---

## Phase 5 — TypeScript 渐进迁移

- 新增 `tsconfig.json`（`allowJs: true` + `checkJs: true` 起步）、`vue-tsc` 检查脚本（并入验证流，不阻塞）
- 迁移顺序（每步可独立合入、随时暂停）：
  1. `src/utils/generators/**`（纯 Node、无 Vue 依赖，风险最低、收益最大）
  2. `src/utils/contentLoader.js` / `markdownProcessor.js`
  3. `src/stores/**`、`src/router/**`
  4. 视图/组件逐文件 `lang="ts"`
- JS/TS 混存不阻塞；每步保持 lint + build 通过

**验证**：`vue-tsc --noEmit` 无新增错误；功能回归。

---

## Phase 6 — 体验增强（可选，视前几阶段节奏决定是否执行）

1. **PWA 离线**：`vite-plugin-pwa`，预缓存内容 HTML + JSON，离线可读全部文章
2. **站内全文搜索**：构建期生成 minisearch 索引（跟随 generatePosts），运行时零后端模糊搜索
3. **资产本地化**：Font Awesome 由 CDN 换本地子集（先统计实际用到的 icon 再裁剪）
4. **双语言路径化** `/zh/` `/en/`（SEO 完整化；改动面最大，最后做）

---

## 附：风险与回退预案汇总

| 风险 | 位置 | 预案 |
|------|------|------|
| vite-ssg 与 Vite 7 不兼容 | Phase 2a | 先 spike 验证；失败则自写 `entry-server.js`（方案 B） |
| 预渲染后文章内容为空壳 | Phase 2 | 检查 Article.vue 内容加载是否同步化（Phase 1d 前提）；`includedRoutes` 路径枚举与 categories.json 对照 |
| dev 下新 HTML 文件不被 glob 拾取 | Phase 0a | full-reload 不生效时改 `server.restart()` |
| 生成器重构引入输出差异 | Phase 3 | 重构前后 JSON 逐字段 diff |
| 组件迁移引入行为回归 | Phase 4 | 单组件迁移 + 即时回归，发现即回滚该文件 |
| `.gitignore` 变更影响仓库状态 | Phase 1a | 仅追加 `src/content/html/`，不动既有 JSON 追踪 |
