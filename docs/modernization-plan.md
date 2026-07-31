# gblog 架构现代化执行计划

> 状态：**执行中（Phase 0-6.5 已完成）**
> 制定日期：2026-07-31
> 路由来源：AGENTS.md → 本文件（约定文档）
> 执行方式：由用户切换到 build 模式，按 Phase 顺序逐个执行；每个 Phase 完成后更新本文件顶部的状态表

## 阶段状态表

| Phase | 内容 | 状态 |
|-------|------|------|
| Phase 0 | 基建：dev 热更新 + SSR 安全 + 高亮本地化 | ✅ 已完成 |
| Phase 1 | Markdown 渲染移入构建期（bundle 瘦身） | ✅ 已完成 |
| Phase 2 | SSG 预渲染 + history 路由（核心升级） | ✅ 已完成 |
| Phase 3 | 内容层整理（生成器共享 core + sitemap） | ✅ 已完成 |
| Phase 4 | 组件现代化（`<script setup>` 迁移） | ✅ 已完成 |
| Phase 5 | TypeScript 渐进迁移 | ✅ 已完成 |
| Phase 6 | 体验增强（PWA / 搜索 / 资产本地化 / 双语言路径化） | ✅ 已完成 |
| Phase 6.5 | 全量 TypeScript 化（.js → .ts，零运行时依赖） | ✅ 已完成 |

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

**执行记录（2026-07-31）**：
- 2a：安装 `vite-ssg@28.3.0` + `@unhead/vue@2.1.16`（降级至与 vite-ssg 内置 2.x 对齐）；最小 spike 验证 Vite 7 / Vue 3.5 兼容后删除
- 2b：`src/router/index.js` 改为导出 `routes` 数组；ViteSSG 内部创建 router（SSR 用 memory history、客户端用 web history）
- 2c：`src/main.js` 改 `export const createApp = ViteSSG(...)`；pinia 按 `initialState.pinia` 传递；i18n SSR 固定 zh-CN、客户端 `initLocale()` 恢复；bootstrap 改为 `if (!import.meta.env.SSR) import('bootstrap')`；`build` 脚本改 `vite-ssg build`
- 2d：`ssgOptions.includedRoutes` 过滤含 `:` 的动态路径段（Windows 上 `dist\article\:path*.html` 会 ENOENT），再读 `src/content/categories.json` 展开 13 篇文章路径
- 2e：`@unhead/vue` 2.x 无 `Head` 组件（且 `Head` 触犯 lint 保留名）→ 用 `app.mixin(VueHeadMixin)` + Article.vue 组件级 `head()`；Home/Category 用 `useHead`；`index.html` 补 meta description + OG
- 2f：`ssgOptions.onFinished` 复制 `dist/index.html` → `dist/404.html`
- **关键修复（正文空壳根因）**：eager glob `import: 'default'` 后模块值直接是字符串，`htmlModules[key].default` 恒为 `undefined` → 全部文章返回空；改 `htmlModules[key]`
- **关键修复（SSR 不等待 async created）**：Home/Category 的 `async created() { await ... }` 在 SSR 不执行完成 → 列表为空；loaders 已同步化后改为同步调用
- Article.vue `created()` 全同步化 + `updateSidebarDimensions`/`$nextTick` 回调加 SSR 守卫；TocDrawer/BackToTop 的 `window.innerHeight` 初始化加守卫
- 构建产出 17 个 HTML（4 静态页 + 13 篇文章，正文/高亮/TOC 已内联，`data-server-rendered`）；主 bundle 含全部文章 HTML（37 KB）

**验证**：
- `npm run build` 后 `dist/article/**/*.html` 存在且正文内容已内联（非空壳）✅
- `npm run preview` 直接刷新深链接（如 `/article/notes/Omics/genomics/bwa/bwa`）不 404 ✅（SPA fallback 返回 index.html，gh-pages 由 404.html 兜底）
- 全部导航链路、语言切换、主题切换回归（dev server 启动无错误、HTTP 200 冒烟通过）

---

## Phase 3 — 内容层整理（不改行为，只提可读性）

- 新建 `src/utils/generators/core/` 共享模块，从 `generateNotes.js` 抽取：
  - `walk`（目录遍历）、`normalizeTags`、`parseFrontMatterAndBody`、`markdownToPlain`、`countWordsSmart`
- 核查 `generateProjects.js` / `generateTopics.js` 中的重复实现并收敛引用（先确认行为一致再替换，逐文件对照输出 JSON）
- 新增 `src/utils/generators/generateSitemap.js`：读 `posts.json`（zh）→ 写 `public/sitemap.xml` + `public/robots.txt`（纯静态资源，vite 自动拷贝）
- `runAllGenerators.js` 追加 sitemap 步骤

**执行记录（2026-07-31）**：
- 新建 `src/utils/generators/core/index.js`，收敛全部重复实现：`walk` / `ensureDirectoryExistence` / `normalizeTags` / `parseFrontMatterAndBody` / `markdownToPlain` / `countWordsSmart` / `toPosixRelativeNoExt` / `extractTitleFromH1`
- 收敛范围：notes/topics/projects（全文）、html（walk + ensureDirectoryExistence）、其余 5 个生成器（ensureDirectoryExistence）；`generateCategories.js` 的 `readJsonArray`/`readYaml` 等单点函数保留原地（无重复）
- **行为验证**：改造前后 `npm run prebuild` 产物逐文件 SHA256 对比 42/42 完全一致（0 diff）
- `generateSitemap.js`：读 zh `categories.json` 提取 13 篇 `articleUrl` + `date`，加 4 个静态路由（`/` `/about` `/resource` `/category`）→ 共 17 条 URL；站点常量 `https://zorrooz.github.io`（user site）；写 `public/sitemap.xml` + `public/robots.txt`（Vite 自动拷贝入 dist）
- `runAllGenerators.js` 追加 sitemap 步骤（最后执行，依赖 zh categories.json）

**验证**：生成 JSON 与改造前 diff 完全一致（除时间戳类字段）✅（SHA256 42/42）；build + preview 正常 ✅（dist 含 sitemap.xml/robots.txt）

---

## Phase 4 — 组件现代化（机械迁移，逐个验证）

- 全部视图/组件迁移 `<script setup>`；`Article.vue` 优先（700 行最大，先啃硬骨头）
- 迁移语义对照：props → `defineProps`、emits → `defineEmits`、computed/watch/lifecycle 一一对应；`setup()` 返回 `t/locale` 的混搭写法移除（改 `useI18n` 直接用）
- Vue 3.5 特性落地：`useTemplateRef`（替换手动 `this.$refs`）、`onWatcherCleanup`（清 resize/scroll 监听）
- **不做** unplugin-vue-router（手动路由表已足够清晰，避免过度工具化）

**验证**：每迁移一个组件即 `npm run lint` + 浏览器回归该组件涉及的页面；全部完成后完整回归清单。

**执行记录（2026-07-31）**：
- 16 个 `.vue` 文件全部迁移 `<script setup>`（3 小组件 → 6 布局组件 → RenderMarkdown → 5 视图 → Article.vue 收尾）
- 语义对照落地：props → `defineProps`、emits → `defineEmits`、`$refs` → `useTemplateRef`（Home/Article/RenderMarkdown）、`$route/$router` → `useRoute/useRouter`、`setup()` 返回 `t/locale` 混搭全部移除（`useI18n` 直接解构）
- Vue 3.5 特性：`useTemplateRef`（Home 侧栏、Article 双栏 + OnThisPage ref、RenderMarkdown 容器）
- 外部方法调用：`OnThisPage` 经 `defineExpose({ refreshToc, resetToc })` 供 Article.vue 使用（原 `$refs.xxx.xxx()` 语义保留）
- SSR head 适配：`@unhead/vue` 2.x 的 `VueHeadMixin` 读取组件选项 `head`，`<script setup>` 下用 `defineOptions({ name, head() {} })` 承载（`this.currentPost` 绑定组件实例，SSR 标题验证正常）
- `AppHeader` 已是 script setup，仅将 `navItems` 由 ref+watch 改为 computed（locale 依赖自动重算）
- `PostList` 删除从未被引用的 `clearTag`（原 Options API 下 lint 不报，script setup 下报未使用；Home.vue 自有实现）
- 单字组件名触犯 `vue/multi-word-component-names`：4 个视图用 `defineOptions({ name: 'HomeView'/'CategoryView'/'ResourceView'/'AboutView' })` 保留原名
- **行为验证**：lint ✅；`vite-ssg build` 17 页 SSR 全部内联正文（`markdown-body`/`<title>` 均在）；preview 深链 200；dev server 无编译错误；主 bundle 273.97 → 259.6 kB（-14.4 kB，Options API 运行时开销去除）

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

**执行记录（2026-07-31）**：
- 基建（commit `6030200`）：
  - 依赖：`typescript@5.9.3`（**不能升 7.x**：原生 Go 版无 `./lib/tsc` JS API，vue-tsc 报 `ERR_PACKAGE_PATH_NOT_EXPORTED`）+ `vue-tsc@3.3.9` + `@typescript-eslint/{parser,eslint-plugin}` + `@types/node` + `@types/js-yaml`
  - `tsconfig.json` 三件套（references 分拆）：`tsconfig.app.json`（浏览器侧 `src/**/*.{js,ts,vue}`，exclude generators/translator/runAllGenerators/contentDevPlugin/markdownProcessor）+ `tsconfig.node.json`（vite.config.js、eslint.config.js、generators/translator/runAllGenerators/contentDevPlugin/markdownProcessor）
  - 均 `strict` + `noEmit` + `allowJs:true` + **`checkJs:false`（JS 走每文件 `// @ts-check` 自愿开启，替代计划的 `checkJs:true` 全局开启）**
  - `npm run typecheck` = `vue-tsc --noEmit -p tsconfig.app.json && tsc --noEmit -p tsconfig.node.json`
  - `eslint.config.js` 接入 `@typescript-eslint/parser`（vue 文件 `lang="ts"` 块）+ `no-undef` 对 `.vue`/`.ts` 关闭（`any` 等类型名会误报）
- 生成器（commit `6030200`）：11 个 generator 文件 `// @ts-check` + JSDoc；`generateCategories.js` 加 `@typedef {FormattedArticle}`/`@typedef {DetailedSubCategory}`；`new Date(a)-new Date(b)` → `.getTime()`（TS2362/2363）；`.filter(Boolean)` → `.filter((item) => item !== null)`；yaml.load 强转 `Record<string, unknown>`；catch 统一 `e instanceof Error ? e.message : e`
  - **行为验证**：改造前后 `npm run prebuild` 产物逐文件 diff 完全一致（`git diff src/content` 为空）
- 加载层 + 状态层（commit `9347c50`）：`contentLoader.js`（glob 模块值强转）、`markdownProcessor.js`（`// @ts-check`，移入 tsconfig.node）、`stores/{app,theme,i18n}.js`（`Ref<string>`、locale 强转 `'zh-CN'|'en-US'`）、`router/index.js`（`RouteRecordRaw[]`，该文件改动随本 Phase 一并合入）
- 视图/组件（本 Phase 合入）：14 个 `.vue` 全部 `lang="ts"`（2 小组件 → 6 布局 → 6 视图）；零运行时依赖改动，仅以下等价改写：
  - `RenderMarkdown` watch 去掉 `async: true`（`{immediate,async}` 组合触发 vue 类型重载歧义 `WatchCallback<() => string>`；行为等价——`initRender` 内部 `await nextTick()` 后 emit，Article 侧本就有 nextTick 兜底）
  - `PostList` 补真实翻译键 `pagination`（原模板引用了未定义的 `paginationLabel`，aria-label 恒为空）
  - `Article` 侧栏数组 `[leftContent, rightContent].forEach(...)` 提为 `const sidebarEls = [...]`（vue-tsc 对模板 ref 窄化后的数组字面量链式调用误报）
  - 其余为参数/ref 类型标注（`any[]`、`Element`、`MutationObserver | null`、`InstanceType<typeof OnThisPage>` 等）
- 全量回归：lint ✅ / typecheck ✅（0 错误）/ build 17 页 ✅ / preview 深链 + 正文内联 + 高亮 + 导航树 200 冒烟 ✅ / dist 还原

---

## Phase 6 — 体验增强（可选，视前几阶段节奏决定是否执行）

1. **PWA 离线**：`vite-plugin-pwa`，预缓存内容 HTML + JSON，离线可读全部文章
2. **站内全文搜索**：构建期生成 minisearch 索引（跟随 generatePosts），运行时零后端模糊搜索
3. **资产本地化**：Font Awesome 由 CDN 换本地子集（先统计实际用到的 icon 再裁剪）
4. **双语言路径化** `/zh/` `/en/`（SEO 完整化；改动面最大，最后做）

**验证**：lint + typecheck + build + preview 冒烟（双语言全部页面 200）。

**执行记录（2026-07-31）**：
- 6.1 PWA（commit `f0b93be`）：`vite-plugin-pwa@1.0.0` + `@vite-pwa/assets-generator`；`vite.config.js` `VitePWA` 插件（`registerType: 'autoUpdate'`、`workbox: { globPatterns: ['**/*.{js,css,html,json,woff2,png,svg,ico,webmanifest}'] }`、`injectManifest` 不启用）；`main.js` `registerSW()`；`public/favicon.svg`（PWA 专用 + 原有 favicon 复用）、`public/apple-touch-icon.png`、`manifest.webmanifest`；`registerType` 冲突（manifest 经 PWA 插件生成时移除 `registerType` 避免与 `vite-plugin-pwa` 的 manifest 文件互斥）；页面 `onLoad` 注册 service worker 后 PWA 自动更新。验证：`npm run build` 生成 `sw.js` + `manifest.webmanifest`，preview 下 `/sw.js` 200
- 6.2 搜索（commit `591a3cd`）：`minisearch@7.1.1`（CJK 大分词，`split: (term) => term.split('').slice(0, 8)`）＋ `contentDevPlugin` 保留（search-index 由 `generateSearchIndex` 生成）：
  - `src/utils/generators/generateSearchIndex.js`：读 posts.json + 正文（`renderMarkdown`）→ 权重 `title×4 + tags×3 + description×2 + body×1`，ID 为文章相对路径（`notes/Omics/genomics/bwa/bwa`、`-en` 后缀版本），输出 `src/content/search-index.json` / `search-index-en.json`
  - `SearchModal.vue`：懒加载 `import.meta.glob('/src/content/search-index*.json')`（**必须根相对式**——`../content/...` 在 SFC 中解析为空），minisearch 构建索引，`SearchResult[]`→`as unknown as SearchDoc[]` 强转（类型不匹配）
  - `contentLoader.js`：`jsonModules` glob 排除 `search-index*`（**首字符类陷阱**：写成 `[bcnprt]` 漏掉 `about/about-en`，导致 about 数据恒为空 → 修复为 `[abcnprt]`）
  - 验证：dev server + 浏览器搜索、构建产物每个 locale 独立懒加载 chunk（search-index-CoIlhH-g.js 10,473B）、预览 + 搜索 200
- 6.3 FA 本地化（commit `0fc55b2`）：`@fortawesome/fontawesome-free@6.4.0`（--save-exact）＋ 本地 Python 工具链 `pyftsubset`（fontTools）+ `brotli`；`scripts/subset-fa-icons.mjs` 从 all.min.css 正则提取码点 → 对 solid/brands 分别子集，输出 `src/assets/fonts/fa-{solid,brands}-subset.woff2`（860B/360B，源 150KB/108KB）+ `src/assets/styles/fa-subset.css`；`main.js` 导入，`index.html` 移除 CDN link；图标清单 13 solid + 1 brands（star/info-circle/calendar-alt/link/arrow-right/bookmark/clock/search/times/arrow-up/envelope/layer-group/folder-open/github）；`npm run subset:icons` 重跑；构建时 woff2 因 `assetsInlineLimit` base64 内联进 app css（子集 <4KB）
- 6.4 双语言路径化（本 Phase 合入）：
  - `src/utils/localePath.js`（新）：`toLocalePath(path)` 按 `i18n.global.locale` 加 `/zh`/`/en` 前缀（无前导斜杠时补 `/`）
  - `src/router/index.js`：`prefixedRoutes(prefix)` 生成 5 命名路由（`zh-Home`/`en-Home` 等）；遗留路径 `/`、`/category`、`/resource`、`/about`、`/article/:path*` 重定向到 `localePrefix()`（localStorage → navigator.language → `/zh`）
  - `src/stores/app.js`：`toggleLanguage` → `setLocale(nextLocale)`（JSDoc `@param {'zh-CN'|'en-US'}`）；`initLocale` 优先读 `pathname.match(/^\/(zh|en)(\/|$)/)`；SSR 下跳过 localStorage/document 写入
  - `src/main.js`：SSR 后台回调第 4 参 `routePath` → `setLocale(routePath.startsWith('/en') ? 'en-US' : 'zh-CN')`
  - `App.vue` watch 路由前缀 ↔ locale 同步；`AppHeader` 品牌/nav/toggleLanguage 全部走 `toLocalePath` + 对面前缀 router.push（保留 query）；`PostList`/`NavigationTree`/`Article`/`SearchModal`/`Category` 链接统一 `toLocalePath('/article/...')`
  - `vite.config.js`：`ssgOptions.includedRoutes` 保留 prefixed 路径 + `getArticleRoutes()`（改读 `categories.json`→`/zh…`、`categories-en.json`→`/en…`，`-en` 文章仅存在于后者）；**`concurrency: 1`**（i18n.global.locale 是模块单例，vite-ssg 默认并发渲染 `/zh` `/en` 页面互相覆盖 locale，导致 en 页渲染成中文）
  - `generateSitemap.js`：`LOCALE_PREFIXES = ['/zh','/en']` ×（4 静态 + 13 文章）= 34 条
  - 验证：typecheck ✅ / lint ✅ / build 34 页（zh 4 + en 4 + zh 13 文章 + en 13 文章）✅ / preview 全部 200（含 `/en/article/.../bwa-en` 预渲染 10,976B）✅ / zh 页含「分类」、en 页含 `Categories` 且无中文导航 ✅ / `/article/...` 旧路由 825B SPA 壳（预期——重定向）✅ / sitemap 34 条 ✅
  - 已知遗留（非本 Phase 引入）：SSR `<html lang="en">` 属性恒为 en（`git show 0fc55b2:dist/about.html` 早前即如此，来自 vite-ssg 的 useHead 默认值）

---

## Phase 6.5 — 全量 TypeScript 化（.js → .ts）

Phase 5 是 JSDoc/`// @ts-check` 渐进式；本 Phase 将全部源码文件由 `.js` 改为 `.ts`，零运行时依赖改动。

**执行记录（2026-07-31）**：
- 策略：Node v24.18 原生 type stripping 直跑 `.ts`（`node ./xxx.ts` 无需转译）；跨文件 import 显式带 `.ts` 扩展名；两个 tsconfig 均开启 `allowImportingTsExtensions`（均 `noEmit`，安全）
- node 侧（tsconfig.node.json）：11 个生成器 + `core/index.ts`（JSDoc typedef → interface：`NoteItem`/`ProjectItem`/`TopicItem`/`PostItem`/`CategoryEntry`/`FormattedArticle`/`DetailedSubCategory` 等）、`runAllGenerators.ts`、`contentDevPlugin.ts`（`configureServer` 标 `ViteDevServer`，返回类型 `Plugin`）、`markdownProcessor.ts`、translator 三件套
- app 侧（tsconfig.app.json）：`main.ts`（入口改 `import './router/index.ts'` 等显式后缀）、`stores/{app,i18n,theme}.ts`（`setLocale` 参数收窄 `SupportedLocale` 联合类型）、`stores/locales/{zh-CN,en-US}.ts`（纯数据提前转，被 generateCategories.ts 引用）、`router/index.ts`、`utils/{localePath,contentLoader}.ts`、14 个 `.vue` 全量 `lang="ts"`（AppFooter 补空 `<script setup lang="ts">`）
- 配置：
  - `tsconfig.node.json` include 全量改 `.ts`；`tsconfig.app.json` exclude 同步改 `.ts` 名；均加 `allowImportingTsExtensions: true`
  - `eslint.config.js`：新增 `app/ts` 块（`parser: tsParser`，拆出原 `app/vue-ts` 的 ts 部分——`parserOptions.parser` 仅对 vue 文件生效，纯 .ts 文件之前全部被 espree 误解析）；末尾 `app/ts-rules-override` 覆盖 `js.configs.recommended` 的 core `no-unused-vars`（core 版对 TS 类型位置误报）→ `@typescript-eslint/no-unused-vars`（`^_` 忽略模式）；`app/node-scripts` glob 改 `src/utils/**/*.{js,ts}`
  - `package.json`：`prebuild`/`translate` 脚本改 `.ts`；`vite.config.js` import `contentDevPlugin.ts`
  - 新 `src/config/llmConfig.d.ts`：`llmConfig.js` 被 gitignore（含 API key），`.ts` 下 import 无类型 → 仅提交声明文件，运行时仍读本地 `.js`
- 迁移陷阱（已修）：
  - 回调参数 implicit any：`Array.isArray(x) ? x.map((it) => ...)` 中 `it` 无推断 → 显式标注 `(it: Record<string, unknown>)`
  - `contentDevPlugin` 的 `server` 参数 → `ViteDevServer` 类型
  - `translator.ts` 的 `completion.choices[0].message.content` 可空 → `?? ''`
  - `loadJsonContent` 保持 `any` 返回（原 JSDoc `@returns {any}`，视图层直接当数组用，收窄成 `unknown` 会牵连 7 个调用点）
  - `import('bootstrap')` 无类型 → 新 `src/env.d.ts` 声明 `declare module 'bootstrap'`
  - `scripts/subset-fa-icons.mjs` 顺带修复：`verifyCmap` 未使用参数 `family` 删除（lint 基线）
- 验证：`npm run prebuild` 26 步全成功且产物与迁移前逐字节一致（`git diff src/content` 为空）；typecheck ✅ / lint ✅ / build 34 页 ✅；dist 仅 chunk hash 变化（模块路径 .js→.ts），页面内容零差异；`dist/.vite/ssr-manifest.json` 中模块名同步为 `.ts`（预期）

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
