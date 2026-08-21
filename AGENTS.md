# AGENTS.md — gblog

Vue 3 + Vite 7 双语静态博客（zh-CN / en-US），Markdown + YAML 写作，构建期生成 JSON/HTML，vite-ssg SSG + GitHub Pages。代码与内容分离：

- `main` — 纯代码（不含任何 content）
- `data` — 纯数据（`content-src` + `cache`），本地经 git worktree 挂到同级 `../blog-data`
- `gh-pages` — 构建产物（GitHub Actions 自动生成，勿手动提交）

> **前置**：首次克隆需 `git worktree add ../blog-data data`（与仓库同级）。

## 命令

| 命令 | 说明 |
|------|------|
| `npm run dev` | Dev :5173（contentDevPlugin 监听数据目录，变更→重跑生成器→full-reload） |
| `npm run build` | Generators → vite-ssg build（SSG + PWA） |
| `npm run prebuild` | 仅内容生成器（本地验证内容用） |
| `npm run preview` | 预览生产构建 |
| `npm run lint` / `npm run typecheck` / `npm run format` | ESLint --fix / vue-tsc + tsc / Prettier |
| `npm run data:translate` | AI 增量翻译（写 cache/en） |
| `npm run data:tag-merge` | zh→en 标签映射增量补齐 |
| `npm run data:pack` | 导出文章自包含包（md + 引用图 → `exports/`；`--out`/`--zip`；只读源） |
| `npm run data:deploy` | `git -C ../blog-data add -A && commit && push origin data`（触发 CI） |
| `npm run data:publish` | 仅 `git -C ../blog-data push origin data` |

> 规范：数据分支操作统一 `data:` 前缀（kebab-case）；`translate`/`tagmerge`/`pack`/`deploy` 为旧名别名。
> **Always lint and build after changes.** 发布流程：编辑 content-src → `npm run prebuild` 验证 → `npm run data:deploy` → CI 部署。

## 架构

```
src/                      # 浏览器应用（含 SSR；不含 Node 构建脚本）
├── main.ts               # ViteSSG → Pinia → Router → i18n；内联 v-reveal 指令
├── App.vue               # AppHeader + router-view + AppFooter（含 skip-link → main#main-content）
├── router/index.ts       # /{zh,en} × 5 路由 + 旧 URL 重定向
├── locale.ts             # currentLocale()/persistLocale() 单一来源
├── config.ts             # SITE/BLOG_TITLE、LOCALE_MAP、THEME_MODES、ARTICLE_ROUTE_PREFIX、HEADER_OFFSET
├── i18n/                 # vue-i18n；schema.ts（AppMessages = typeof zhCN），en-US 用 satisfies 编译期校验
├── stores/               # theme.ts（ThemeMode + initTheme）、locale.ts
├── views/                # Home / Category / Resource / About / Article
├── components/           # layout/（AppHeader/NavigationTree/OnThisPage/RenderMarkdown…）
│                         # widgets/（BackToTop/SearchModal/TocDrawer/FloatingButton）
│                         # common/（IconButton/ModalOverlay/PostCard/Pagination/EmptyState）
├── composables/          # useLocalizedContent / usePageMeta / useSearch / useTagNavigation / useCopyFeedback / useFloatingButton
├── utils/                # 浏览器运行时纯工具（不 import 任何 scripts/ 代码）
└── assets/               # 字体 OTF / styles / avatar.*
```

```
scripts/                  # Node 工具链（Node >=23.6 原生 type stripping 直跑 .ts）
├── dataConfig.ts         # 数据目录唯一接入点（GBLOG_DATA_DIR；contentSrcDir/enSrcDir/contentDir/cacheDir）
├── runAllGenerators.ts   # 生成器编排（locale 步骤表 + 依赖顺序；vite dev 插件复用）
├── llmConfig.ts          # DeepSeek API 配置（gitignore，勿提交）
├── lib/                  # 共享工具层：llm / fs / text / cli / frontmatter（含 sanitizeFrontmatter）/ tags
│                         # / tagMapping / yamlEntries / yamlJson（YAML→JSON 骨架）/ metadata（projects/topics）
├── generators/           # 11 个生成器 + markdownProcessor（md→HTML 渲染管线）+ core/（兼容 barrel）
└── tools/                # 独立 CLI 内容工具：packArticle（导出）/ translator（AI 翻译）/ tagMerger（标签映射）

vite/contentDevPlugin.ts  # dev 插件：监听数据目录 → runAllGenerators → full-reload
```

## 数据模型（三层）

```
content-src/    第一层 src   — 纯手写中文源（严格无机器产物，仅这里编辑）
│   ├── {categories,about,resources}.yaml
│   ├── assets/              — 配图集中目录（Obsidian 附件默认位置；站点直接解析，无需打包）
│   └── notes/<cat>/<sub>/<slug>.md   — 扁平：一篇文章一个 md
cache/          第二层 cache — 机器维护的持久态（入库）：en/（英文镜像，-en 身份后缀）、
│                               tag-mapping.json、.translate-state.json（源相对路径 → 内容 hash）
content/        第三层 final — 生成 JSON + html（可再生，gitignore，不入库）
```

- **中文源**：只写 content-src；**英文**：绝不手写，`data:translate` 生成到 cache/en/（LLM 非确定 + 有成本，入库存证，CI 用 `GBLOG_NO_TRANSLATE=1` 跳过）
- 生成器按 locale 取源：`srcDirFor(locale)`（zh→content-src，en→cache/en），输出恒为 `content/*` 与 `content/*-en*`
- 中英实体按镜像相对路径配对：`cache/en/notes/.../<slug>-en.md` ↔ `content-src/notes/.../<slug>.md`

## 内容管线（Critical Order）

```
1.  generateNotes      — srcDirFor(locale)/notes/**/*.md → notes.json
2.  generateProjects   — categories.yaml `projects` 段 → projects.json
3.  generateTopics     — categories.yaml `topics` 段 → topics.json
4.  generateCategories — YAML + notes/projects/topics → categories.json
5.  generatePosts      — notes.json + categories.json → posts.json
6.  generateTags       — posts.json → tags.json
(1-6 先 zh-CN，重复 en-US)
4.5 增量翻译    — 内容 hash 判定（仅缺失/变化的 zh 源 → cache/en；GBLOG_NO_TRANSLATE=1 关闭）
4.6 tagMerger   — ensureTagTranslation 增量补齐 zh→en 映射（cache/tag-mapping.json 命中 0 token）
4.7 一致性修复  — 以 zh 为基准重写 -en 文件 tags（失败仅告警）
8.5 一致性校验  — 构建产物级中英标签一致性校验（[OK]/[Warn]）
9.  generateHtml        — md → content/html/**（markdownProcessor 预渲染）
10. generateSitemap     — public/sitemap.xml + robots.txt（gitignore，CI 生成）
11. generateSearchIndex — categories.json + content/html → content/search-index{,-en}.json

generateAbout / generateResources 无依赖，runAllGenerators 中独立运行。
```

## 内容源约定

- **只有 notes 有 markdown**；projects/topics 为纯 yaml（卡片外链 GitHub/DOI，无文章页）
- notes 目录名必须匹配 `categories.yaml` notes 段的 `name`；子分类目录名匹配 `categories` 映射 key
- **图片**：写作期放 content-src/assets/（引用 `assets/xx.png`，Obsidian 最短路径，站点已支持解析，**无需打包**）；分享/迁移用 `data:pack` 导出
- 每篇 md 带 YAML frontmatter：`title` / `date` / `author` / `tags` / `draft` / `description`
- **中英 frontmatter `tags` 必须一一对应**（数量与语义对称），否则标签云双语数量不一致

### URL 映射

```
File: notes/Omics/genomics/bwa.md
URL:  /zh/article/notes/Omics/genomics/bwa     (en → /en/article/notes/Omics/genomics/bwa-en)
Route: { path: '/:locale/article/:path*', name: 'zh-Article'|'en-Article', props: true }
旧 URL /article/...（无前缀）→ 302 到 preferredLocaleSegment()
旧「每文一目录」URL（.../bwa/bwa 或 .../bwa/bwa-en）→ beforeEnter 自动重定向
```

## 国际化（双层）

1. **App**：vue-i18n（zh-CN/en-US），schema.ts `AppMessages = typeof zhCN` 基准，en-US `satisfies` 编译期校验；SEO title 走 usePageMeta（i18n meta.* 键 + BLOG_TITLE）
2. **Content**：locale 源目录分层 + `-en` 身份后缀；`contentLoader.ts` 英文缺失时自动回退中文

## 代码风格

- No semicolons, single quotes, 2-space indent, LF, 100-char width
- Vue `<script setup>`（Options API mix）；无多余注释；`@/` 别名
- ES modules only；Named exports 优先
- **生成器必须保持与产物逐字节一致**（改结构前先跑 prebuild 对比）

## 布局与形状规则（硬性）

- **内容宽度**：内容页 `.page-section` 1120px（pad 48/32/20）；文章页 container 1280px、正文 `.article-content` **780px** 居中、左右 sidebar 各 **240px**；移动端水平边距全局统一 20px
- **标题层级**：文章标题 `.article-title` clamp(30-40px)；正文 H1 32 / H2 24 / H3 20 / H4 18 / H5 17 / H6 14（移动端 ≤768：24/21/19，≤576：22/19/18；正文 16px）
- **圆角分层**：控件 6 / 按钮 8 / 卡片 12 / 面板 16 / 胶囊 full；圆形（50%）仅用于图标操作按钮，方形按钮一律 8px，同一容器不混用
- **侧栏导航/TOC**：分组标题 13px 加粗正文；条目 13px 纯文字链接（无圆角底块），hover 变主色，激活 = 主色 + 左侧 2px 指示条（不加粗）；可折叠目录 caret 旋转；`.tree-sublist` 导轨 `border-left` + padding-left 16px，激活条 `left: -17px`；TOC 链接 13px 换行显示（不截断）
- **卡片 hover**：仅边框变主色微光（无背景变化、无阴影）；浮动层（抽屉/弹窗/浮动按钮）才用阴影
- **动效**：v-reveal 入场 opacity+translateY(14px)+blur(5px)，`--reveal-delay` stagger，reduce-motion 用户跳过；过渡 140-220ms（`--dur-fast`/`--dur-base`）
- **404**：品牌 404（巨型 serif 数字 + 返回首页 pill）；返回按钮 **`flex: none` 不可移除**（App.vue `.main-content > * > * { flex:1 }` 会拉伸短页面最后子元素）
- **图标**：操作图标 = 内联 stroke SVG（viewBox 24、stroke-width 1.75、fill none、stroke currentColor、18px）；**禁止 PNG 位图图标**；装饰图标用 FontAwesome（数据驱动）
- **字体**：思源宋体（标题）/ 思源黑体 / Agave（等宽）完整文件引用，**无子集化、无构建脚本**
- **主题持久化**：`initTheme()` 必须重读 `initialTheme()` 覆盖 SSR 内嵌 pinia state（恒为 auto，hydrate 会覆盖 localStorage）
- **语言切换不跳位**：zh↔en 切换必须保持当前阅读位置——文章页按 `restoreScrollY` 恢复滚动（双 nextTick 等布局稳定），资源页按索引路径 `[catIdx, subIdx]` 恢复选中子分类（中英 title 不同，不能按名字匹配）
- **复制**：文章复制取 `.markdown-body` innerText（清辅助元素）；表格复制为 TSV（\t 分隔，可贴 Excel）

## 颜色系统（三色族：蓝/灰/黑）

- 色板：`blue-600 #3b82f6`（主色，暗色 `#6cb2ff`）/ `blue-500` / `blue-200` / `blue-100`；`gray-50/100/200/400/500`；`ink #16191f`（暗色 `#e7e9ee`）
- 语义映射：`--primary` = blue-600（**必须是最深交互蓝**）；hover 再暗一档（`--primary-hover`）；active 最暗；`--primary-light` = blue-500；`--tint` = blue-100；`--tint-strong` = blue-200
- **禁忌**：① 禁止主色反转（primary 不能是浅蓝）② 禁止硬编码 rgba（用 `var(--primary)` + opacity 或 color-mix）③ 灰色禁止跳级 ④ 暗色 primary 保持最亮蓝，hover 更亮
- 组件：正文 `--fg` / 次要 `--fg-2` / meta `--fg-3`；链接 `--primary` → hover `--primary-hover`；按钮填充 primary；标签默认 `--surface-2`+`--fg-2`，hover `--tint`+`--primary`；卡片边框 `--line`；分割线 `--line`；选中背景 `--tint`
- 新组件：文字色用 fg 族，交互色用 primary，背景 `--surface`→`--surface-2`→`--tint` 三级，边框 `--line`，hover 边框 `rgba(color,0.3)` / 背景 `--tint`

## ESLint / Prettier

- Flat config；`@typescript-eslint/parser` 处理 `**/*.ts` + vue `lang="ts"`；browser globals 用于 `*.{js,vue}`，node globals 用于 `scripts/**` + `vite/**`；忽略 `dist/`、`dist-ssr/`、`coverage/`
