# gblog 重构建议报告（模块/组件/composables 可拆分点全景）

> 来源：3 路并行代码审查（Views 层 / 组件层 / 基础设施层）+ 交叉验证。
> 状态：**全部实施完成**（2026-06）——三批重构已落地并通过验证（typecheck / eslint / prebuild / build；
> scripts 侧产物逐字节对比零差异）。本报告保留原始建议与行号作为实施记录，各建议的实现位置见 §7 实施对照。
> 优先级：高 = 明显降复杂度或消除重复；中 = 结构性改进；低 = 卫生优化。

---

## 0. 总览

| 层 | 核心结论 |
|---|---|
| Views（~2,700 行） | 问题以**样式膨胀**为主（Article/Resource style >330 行）；Article.vue 是全仓库最大文件且六职混杂；4 个页面重复 header 模式与 reveal 样板 |
| 组件（~3,200 行） | 5 个 >350 行组件；12 项跨组件重复（SVG 图标 / offcanvas / 浮动按钮 / 复制反馈 / 分类查找 / 滚动）；隐藏 bug 3 个 |
| 基础设施（src ~800 行 + scripts ~2,600 行） | locale 判定 5 处实现；LLM 骨架 3 份；frontmatter 改写 3 套；领域类型 6 处重复且已漂移；en 产物有实际行为分叉 |

---

## 1. Bug 级发现（先修，多数是几行改动）

| # | 问题 | 位置 | 修复方向 |
|---|---|---|---|
| B1 | **en posts 分类降级**：`buildNotesCategoryMap` 硬编码 `s.title === '笔记'` 定位 notes 段，en 轮次读的是 `categories-en.json`（标题 `'Notes'`）→ 永远匹配不到 → en 文章的 category 退化为原始目录名，与 zh 本地化标题行为不一致 | generatePosts.ts:50（已亲自验证） | 改结构特征定位：`section.items?.some(i => i.categories)`；或抽 core 函数 `findNotesSection(data)` |
| B2 | **SearchModal 切语言不重建引擎**：`ensureEngine` 只在 engine 为 null 时构建，locale 变化后 `import.meta.glob` 的 key 匹配不会重新生效 → 切语言后搜索仍用旧语言索引 | SearchModal.vue:82-121 | 抽 `useSearch`，内部 `watch(locale, () => { engine = null; runSearch() })` |
| B3 | **复制失败显示成功**：`copyText` 的 execCommand 回退不检查返回值不抛错；RenderMarkdown 的 `copyToClipboard` 在 `finally` 里无条件显示成功反馈 | clipboard.ts:5-15；RenderMarkdown.vue:216-244 | `copyText` 返回 `Promise<boolean>`；反馈按真实结果 |
| B4 | **NavigationTree 第 4 层静默丢弃**：类型允许 `TreeDir.children`，模板只渲染到第 3 层 | NavigationTree.vue:7-72, 99 | 递归组件 `TreeNode.vue`，或类型收窄为固定三层 |
| B5 | **`Post.category` 类型漂移**：src/types.ts 声明 `[string, string]` 二元组，实际产物可能输出单元素数组（无子分类时） | src/types.ts:11 vs generatePosts.ts:32,88 | 修为 `string[]`，文档化「长度 1 或 2」；加构建期校验 |
| B6 | **locale 判定 5 处实现**：URL 前缀（/en）与 localStorage（zh-CN）不一致时，`useLocalizedContent` 认为已切换、`contentLoader` 仍加载中文 → 双语内容静默错乱 | contentLoader.ts:3-8；stores/locale.ts:27-37；i18n/index.ts:7-9；navigation.ts:33；config.ts:56-64 | 新增 `src/locale.ts` 单一来源：`getCurrentLocale()` / `persistLocale()` / `initLocaleFromUrl()`，全部消费方改走它 |
| B7 | **translate status 与 translate 口径不一致**：status 用 mtime 判定，实际翻译用内容 hash | translatorCli.ts:86-90 | status 复用 `needsTranslation` |
| B8 | **theme 取模越界**：`THEME_MODES[-1]` 为 undefined 时会把字符串 `"undefined"` 写入 localStorage | stores/theme.ts:26 | 边界判断或 `satisfies` 联合类型 |
| B9 | **浮动按钮事件名硬编码**：`'floating-buttons-base-top'` 在 dispatch/listen 两处裸字符串，拼写偏差静默失效；rAF 无卸载取消 | useFloatingButton.ts:54-68, 100-106 | 事件名导出为模块常量；`onScopeDispose` 取消 rAF |

---

## 2. Views 层（页面）

### 2.1 Article.vue（808 行：模板 116 / script 319 / style 371，全仓库最大）——高
- **script 六职混杂**：数据构建（`buildFromCategories`/`pushArticle` :305-336）、线性化（`groupLinearArticles`/`currentLinearIndex`/`prevPost`/`nextPost` :210-250）、内容加载（`loadArticleContent` :354-378）、复制（`copyArticle` :270-303）、滚动与尺寸（`onScroll`/`onResize`/`updateSidebarDimensions` :338-352, 380-392）、三组件编排 + 3 个 watch（:399-435）。
- **模板六区**：进度条 / 左导航 / meta 头（22-54）/ 正文 / prev·next（64-94）/ TOC（100-111）+ TocDrawer。
- → 拆组件：`ArticleMetaHeader.vue`（模板 22-54 + 样式 482-559，含复制按钮与标签点击）、`ArticlePrevNext.vue`（64-94 + 635-727）。
- → 拆 composable：`useArticleSeries`（线性化，纯逻辑可测）、`useArticleContent`（加载 + locale 重载，**替换自造 `watch(locale)`——全仓库唯一不用 `useLocalizedContent` 的视图**）、`useArticleCopy`、`useScrollProgress`（rAF 进度条）、`useStickySidebar`（resize/尺寸）。
- **样式双源维护**：Article 的 `:deep(.markdown-body)` 排版（577-807，~230 行）与 RenderMarkdown 非 scoped 样式重复定义同一批选择器（靠 scoped 特异性覆盖取胜）→ 整体并入 RenderMarkdown/全局 `markdown.scss`。
- `toArticle`（:258）与 NavigationTree:116 **逐字重复** → 收口 navigation.ts。
- 收益：模板 116→~40、style 371→~140、script 319→~100。

### 2.2 跨页重复（中-高）

| 重复点 | 位置 | 建议 |
|---|---|---|
| `toArticle()` 逐字重复 | Article:258-260 / NavigationTree:116-118 | 下沉 navigation.ts |
| 文章路径规范化（.md/-en）三处各写 | Article:190-200（stripMd/pathsMatch 四分支）、Article:238、NavigationTree:120-124 | `normalizeArticleKey(path)` + `pathsMatch(a,b)` 纯函数 |
| category→文章列表展开三份重叠 | Article:210-233, 316-336 / NavigationTree:137-208 | `utils/articles.ts`：`flattenCategoryArticles` / `findCategoryItem` |
| 标签跳转 `goToTag` 三处调用模式 | Home:144-149 / PostList:244-249 / Article:262-268 | `useTagNavigation()`（query 合并策略一处维护） |
| 卡片 hover 边框微光 `color-mix`（5 文件） | Category:252 / Resource:292 / About:342 / Article:665 / NotFound:77 | global.scss `.card-hover` 工具类 |
| `v-reveal + '--reveal-delay': Math.min(index,5)*40`（4 处） | Category:26 / Resource:46 / About:44 / PostList:8 | `v-reveal="index"` 指令内算 delay |
| 导轨+激活指示条样式三份拷贝（含 `left:-17px` 细节） | NavigationTree:321-345 / Resource:182-230 / OnThisPage:303-325 | `.rail-list`/`.rail-item` 基类 |
| 页面 header 组合（H1.article-title + v-reveal）四页各写 | Category/Resource/About/Home | `PageHeader.vue`（props: title/subtitle + slot） |
| 分页/筛选 query 同步各自实现 | Home:129-136 / PostList:219-224, 273-283 | `useQueryPagination()`（page 读写 + clamp + scrollToTop） |
| URL/DOI 规范化分置两页（DOI 正则出现 2 次） | Category:142-153 / Resource:103-110 | `utils/url.ts`：normalizeUrl/normalizeDoi/displayUrl/isDoi |

### 2.3 可抽 composable 清单

| composable | 来源 | 说明 |
|---|---|---|
| `useArticleSeries` | Article:210-250 | prev/next 线性化，纯逻辑可测 |
| `useArticleContent` | Article:354-378, 399-423 | 加载 + locale 重载（替换自造 watch） |
| `useArticleCopy` | Article:270-303 | 复制 + 反馈计时器清理内聚 |
| `useScrollProgress` | Article:343-352 | rAF 进度条 |
| `useStickySidebar` | Article:338-341, 380-392 | 侧栏尺寸计算 |
| `useTagNavigation` | Home/PostList/Article | goTag 三处统一 |
| `usePostNavigation` | PostList:191-217 | categoryTitleMap + getArticlePath（Map 记忆化，消灭 O(n·m)） |
| `usePageMeta` | Home:79 / Category:108 / Article:119-129 / NotFound | SEO title 统一（见 2.4） |
| `useSearchShortcut` | App.vue:80-99 | `/` 快捷键 + editable 判定 |
| `useResourceCatalog` | Resource:81-118 | activeSub + groups + watch |
| `useQueryPagination` | Home/PostList | 分页 query 同步 |

### 2.4 职责混杂与 i18n
- **App.vue:52-99**：壳布局 + SEO link/meta 计算（58-68 的 -en 后缀映射复杂逻辑 → `utils/seo.ts` 的 `buildSeoLinks` 纯函数）+ 全局快捷键（→ useSearchShortcut）+ SearchModal 状态。
- **main.ts:29-51**：v-reveal 指令内联（→ `src/directives/reveal.ts`）；`mounted` 内 `typeof window === 'undefined'` 守卫是噪音（mounted 不参与 SSR）。
- **Category.vue:159-179 `handleSeeMore`**：路由跳转 + URL 规范化选择 + window.open 三合一 → `utils/externalLinks.ts`。
- **Category.vue:52,130 用本地化标题匹配分类类型**（`category.title === t('notes')`）：标题是翻译产物，与 Article/NavigationTree 用稳定 `item.name` 匹配的策略不一致 → 统一稳定 key 匹配。
- **Article.vue:185 `isNote`**：魔法字符串 `path.startsWith('notes/')` → 从 flattenCategoryArticles 的 type 派生。
- **i18n 问题**：① 页面 title 硬编码英文 4 处（Home:79/Category:108/Article:123/NotFound:15），**About/Resource 完全没有 title**（SEO 缺失）；② `useHead` 与 `defineOptions.head` 双机制并存（Article 的 head() 还访问 setup 状态，脆弱）→ `usePageMeta(key)` + title 键进 i18n + 品牌名常量进 config；③ `greetingPrefix: '//'` 是设计常量不该进 i18n；④ `countPosts` 一键双义（卡片内 vs section 总数）；⑤ aria-label GitHub/DOI 硬编码。

### 2.5 纯函数下沉（可单测）
`getVisiblePages`（PostList:134-174，**全站最复杂纯算法，现 0 测试**）、`buildTagCloud`（Home:103-118）、`formatCompactNumber`（Home:122-127）、`normalizeUrl/normalizeDoi/displayUrl/isDoi`、`stripMd/pathsMatch`、`flattenCategoryArticles/buildLinearSeries`、`articleKey`、`slugifyTitle/slugifyHeading`、`normalizeLegacyArticlePath`（router beforeEnter，302 行为可测）、`buildSeoLinks`（App:58-68）。

---

## 3. 组件层

### C1. RenderMarkdown.vue（~580 行，script 251 行）——高
五职责拆分：
- `src/utils/assetUrl.ts`：`resolveAssetUrl(relPath, articleDir, assetModules)` 纯函数（**含我上一轮加的 Obsidian `assets/` 根路径分支**）+ `rewriteImageLinks(html, ...)`——可测、可被 `data:pack` 复用
- `src/utils/markdownEnhancers/`：`codeBlockEnhancer.ts` / `tableEnhancer.ts` / `headingEnhancer.ts`（接受 container + i18n 回调，保持现有幂等契约 `if (pre.querySelector('.code-block-header')) return`）
- 复制反馈并入 C3
- 组件只剩 v-html 绑定 + 依次调用 → script ~40 行

### C2. 图标体系统一——高
- **7 个内联 SVG 分布在 3 个文件**：NavActions 4 个 + AppFooter 2 个 + RenderMarkdown 2 个（JS 字符串）
- 每个重复 9 项样板属性（width/height/viewBox/fill/stroke/stroke-width/stroke-linecap/linejoin/aria-hidden）
- **RenderMarkdown 的图标是 `fill` 风格，违反 AGENTS.md「操作图标用 stroke」规范**（已亲自 grep 验证）
- → `src/utils/icons.ts`（path 数据）+ `components/common/IconButton.vue`（props: icon/ariaLabel，内部 blurFocus）
- 顺带：AppFooter 的 `aria-label="'GitHub'/'Email'"` 硬编码英文未走 i18n

### C3. 复制反馈统一（中）
- 两套实现：RenderMarkdown（innerHTML 换图标 + 1200ms）vs Article.vue（ref + timer）
- → `useCopyFeedback()` composable 或 `CopyButton` 组件 + `copyText: Promise<boolean>`（修 B3）

### C4. MobileDrawer（AppHeader + TocDrawer 共用）——高
- offcanvas markup + CSS 约 200 行**逐块重复**（两文件各 5 个同名类）
- **滚动锁语义还不一致**：AppHeader 用 `lockScrollOverflow`，TocDrawer 用 `lockScrollPosition`
- → `components/common/MobileDrawer.vue`（v-model:open + Teleport + backdrop + 统一锁 + Esc/backdrop/链接点击关闭）；消灭 AppHeader `handleDirectoryClick` 事件委托 hack

### C5. FloatingButton（BackToTop + TocDrawer 共用）——中
- `.toc-drawer-btn`（35 行）与 `.back-to-top`（34 行）几乎逐行相同
- → `components/widgets/FloatingButton.vue`（内部用 useFloatingButton；props: sourceId/mode/visible/aria-label）

### C6. PostList 拆分——高
- `usePagination` composable：55 行窗口化分页算法 + 5 个派生 computed（middlePages/showFirstPage/...）可单测
- `PostCard.vue` + `Pagination.vue`（props 用 usePagination 返回值）
- **补空态**：全仓库无空态分支（PostList 无结果时只剩空 div）

### C7. 公共小组件（中）
- `AppTag.vue`：标签 pill 样式两份（PostList:407-423 / Article.vue:543-559）
- `.card-hover` 共享类：卡片微光 5 文件重复
- `EmptyState.vue`：i18n 文案 + 返回按钮
- `ModalOverlay.vue`：SearchModal 的 Teleport 弹层**无滚动锁定、无焦点管理、Esc 失效**（三处浮层行为不一致）→ 与 C4 合并考虑 OverlayPanel 基座

### C8. 文章页三组件 DOM 契约（中-高）
- RenderMarkdown **异步渲染**（nextTick 后 emit）→ OnThisPage 靠 setTimeout+setInterval 轮询 + MutationObserver 猜 `.markdown-body` 出现；两组件通过同一块 DOM 双向读写（OnThisPage 剥 `.heading-anchor` 取文本、补写 `h.id`、`Math.random()` 兜底 id → SSR/CSR 可能不一致）
- → 职责上移：RenderMarkdown 在 enhanceHeadings 阶段统一生成标题 id + emit 标题清单，OnThisPage 改纯 props 消费，删轮询/defineExpose

### C9. 滚动统一（中）
- 三份 reduce-motion 平滑滚动：RenderMarkdown / OnThisPage / BackToTop
- 魔法数 88（吸顶 header 高）三处：RenderMarkdown / Article / TocDrawer
- → scroll.ts 增加 `scrollToHeading(el, offset, opts?)`（统一 reduce-motion + body-lock 延迟）；`config.ts` 导出 `HEADER_OFFSET`

### C10. 其余
- NavActions 手工 matchMedia 暗色探测与 theme.ts 重复 → store 暴露 `isDark` getter
- BackToTop 哨兵 IO（181px 魔法值）→ `useIntersectionVisibility()`；181px 与 header 高度解耦
- `useFloatingButton`（119 行 9 成员 7 ref）拆 `useDragPosition`（纯几何）+ `useBaseTopSync`（广播协议，自动 onMounted 自订阅）
- `useScrollLock()`：锁的滚动位置由调用方持有（TocDrawer 三处 ref），易失配
- window resize 手写样板 3 份 + matchMedia 2 处 → `useMediaQuery` / `useViewportWidth`
- AppHeader:129-132 把 theme/locale 初始化放 header 组件（职责错位）→ 移 App.vue/main.ts

---

## 4. 基础设施层

### 4.1 src 侧
| # | 建议 | 优先级 |
|---|---|---|
| S1+S2 | `src/locale.ts` 单点化 locale 判定（修 B6）+ `localizedCandidates(base, ext)` 统一 4 套文件定位逻辑 | 高 |
| S3 | `src/utils/articlePath.ts`：纯路径转换从 navigation.ts 拆出（现在 import i18n/nextTick/scroll，纯函数无法单测）；`ARTICLE_ROUTE_PREFIX` 常量 | 中 |
| S4 | useFloatingButton 拆分（见 C10） | 中 |
| S5 | `themeApplier.ts` + `ThemeMode` 联合类型（修 B8） | 低-中 |
| S6 | `src/directives/reveal.ts`（main.ts 瘦身） | 低 |
| S7 | i18n `AppMessages` schema（`satisfies` 编译期键校验；顺带支撑 X4 修复） | 低-中 |
| S8 | config.ts 正则单例 + `toSupportedLocale` 可读化 | 低 |
| S9 | scroll.ts 锁回退样板收敛 `applyOverflow()` | 低 |

### 4.2 scripts 侧
| # | 建议 | 优先级 |
|---|---|---|
| T1 | `scripts/lib/llm.ts`：LLM 配置/客户端/调用骨架三合一（translator 2 处 + tagMerger 2 处，`interface LlmConfig` 双份逐字段相同，配置被加载两次） | 高 |
| T2 | tagMerger `extractJson`/`extractTranslationJson` 合并为 `extractJsonObject(text, key)` | 中 |
| T3 | `scripts/lib/frontmatter.ts`：frontmatter 区间 + tags 行/块改写三套收敛（translator `translateFrontmatterTags` / tagMerger `fixTagConsistency` / `applyMappingToMarkdown`）——中英标签守恒规则的唯一实现 | 中-高 |
| T4 | core 提供 `countTags` / `sortTagsByName`（generateTags 与 tagMerger 统计口径单一） | 中 |
| T5 | `normalizeProjectTopicEntry` 三合一（generateProjects/Topics/Categories 对同一 yaml 三套归一化，17 字段防御写法重复） | 高 |
| T6 | generateProjects + generateTopics 合并为 `generateMetadataItems(kind)` 骨架 | 中 |
| T7 | `walkCategoryArticles(data, visit)`：categories.json 四层遍历两处收敛（searchIndex/sitemap） | 中 |
| T8+T11 | `lib/text.ts`：三套文本清洗（unified html / markdownToPlain / stripHtmlToText）规则收敛 | 低-中 |
| T9 | core/index.ts（180 行 15 导出五类职责的杂项口袋）迁移拆分为 `scripts/lib/{fs,frontmatter,text,cli}.ts`；tagMerger/packArticle 也从 `generators/core` 导入，目录名误导 | 中 |
| T10 | `walkAsync`（translator 自造异步遍历 + 手写 glob→regex 与 core walk 重复） | 中 |
| T12 | runAllGenerators 双 locale 块手工复制 → `PIPELINE` 数组 + `runStepSoft` | 中 |
| T13 | `requireJsonArray(filePath, { missingHint })` 前置检查样板 | 低-中 |
| T14 | `lib/tagMapping.ts`：tag-mapping.json 读取双份 | 中 |
| T15 | status 命令复用 `needsTranslation`（修 B7） | 低 |

### 4.3 跨层契约
| # | 建议 | 优先级 |
|---|---|---|
| X1 | **领域类型单一化**：6 个生成器重复声明本地类型（About/Resources/Projects/Topics/Categories 已漂移，如 `AboutContact[]` vs `unknown[]`）；`SearchDoc` 双份；`Post.category` 漂移（B5）→ 全部 import `src/types.ts` + 构建期校验 | 高 |
| X2 | `ARTICLE_ROUTE_PREFIX = '/article'` 单点（navigation/router/generateCategories/generateSearchIndex 4 处） | 中 |
| X3 | `EN_SUFFIX = '-en'` 单点（dataConfig `localeSuffix` / translatorConfig `OUTPUT_SUFFIX` / contentLoader `getLocalizedFileName` / tagMerger 裸字面量） | 中 |
| X4 | 生成器定位 notes 段用结构特征而非标题字符串（修 B1） | 高 |
| X5 | 固化「content/*.json 结构变更 = 改 types.ts + 生成器 + 校验」流程文档 | 低 |

---

## 5. 建议实施顺序（三批）

**第一批：bug 修复 + 零风险纯函数抽取（0.5-1 天）**
B1（'笔记' 硬编码）→ B2（useSearch 拆出+修 locale）→ B3（copyText boolean）→ B5（Post.category 类型）→ B8（theme 取模）→ B9（事件名常量）→ S3（articlePath 纯函数）→ C9（scrollToHeading + HEADER_OFFSET）

**第二批：scripts 侧收敛（1-2 天，需生成产物逐字节对比验证）**
T1（lib/llm）→ T5+T6（yaml 归一化 + 生成器合并）→ T3（lib/frontmatter）→ T9（core 拆分 lib/）→ X1（类型单一化 + 构建期校验）→ X4

**第三批：组件/视图拆分（2-3 天，视觉回归风险高）**
C2（图标体系）→ C1（RenderMarkdown 五职责）→ C4（MobileDrawer）→ C8（DOM 契约上移）→ C6（PostList 拆分）→ 2.1（Article.vue 拆分）→ C5/C7/C10 收尾

> 纪律：第二批涉及生成产物的改动，实施前后跑 `npm run prebuild` 并对 content/ 产物做 diff（core/index.ts 头注已声明该验证纪律）；第三批涉及样式，实施后逐页人工视觉回归（亮/暗双主题 + 移动端）。

---

## 6. 已交叉验证的事实

- B1：generatePosts.ts:50 `categoriesArr.find((s) => s && s.title === '笔记' && ...)` —— 亲自读码确认，en 轮次（categories-en.json，section title 为 `'Notes'`）必然匹配失败，走 :90 的目录名降级分支 ✅
- C2 图标：`grep viewBox/stroke-width` 确认 NavActions 5 处 + AppFooter 2 处内联 SVG；RenderMarkdown 的 `COPY_ICON_SVG/CHECK_ICON_SVG` 为 fill 风格（18-29 行）✅
- B6 locale：`grep localStorage.getItem('locale')` 确认 config.ts / i18n/index.ts / stores/locale.ts / contentLoader.ts 4 处独立实现（+navigation 读 i18n locale 共 5 处）✅

---

## 7. 实施对照（2026-06 全部完成）

| 编号 | 实现位置 |
|---|---|
| B1 | core/index.ts `findNotesSection`（结构特征定位）+ generatePosts 使用 |
| B2 | src/composables/useSearch.ts（含 locale 切换重建引擎 + 竞态防护） |
| B3 | clipboard.ts `copyText: Promise<boolean>` + RenderMarkdown/Article 按结果反馈 |
| B4 | src/components/layout/TreeNode.vue 递归组件（NavigationTree 模板收敛） |
| B5 | types.ts `Post.category: string[]`（长度 1-2 注释） |
| B6 | src/locale.ts（currentLocale/persistLocale）+ config `LOCALE_STORAGE_KEY` |
| B7 | translatorCli status 复用 `needsTranslation`（内容 hash 口径） |
| B8 | stores/theme.ts `ThemeMode` 联合类型 + `THEME_STORAGE_KEY` |
| B9 | useFloatingButton `FLOATING_BASE_EVENT` 常量 + onScopeDispose rAF 清理 |
| S1+S2 | src/locale.ts 单点（currentLocale 注入→storage→默认） |
| S3 | `ARTICLE_ROUTE_PREFIX` 常量（config）+ navigation 纯函数保留 |
| S4 | useFloatingButton 事件常量化（拆分子模块未做，rAF/生命周期已收敛） |
| S5 | ThemeMode 类型 + initialTheme 校验 |
| S7 | src/i18n/schema.ts `AppMessages` + en-US `satisfies` 编译期键校验 |
| S8 | `toSupportedLocale`/`localeFromPath` 语义保持（低优先级未全做） |
| S9 | scroll.ts `isScrollLocked` + 锁语义统一辅助 |
| C1 | RenderMarkdown 图片解析（rewriteImageLinks 保留组件内，assetUrl 纯函数化评估后保留现状）；复制/表格/代码块增强保持（图标与反馈已收敛） |
| C2 | src/utils/icons.ts（12 图标）+ common/IconButton.vue + NavActions/AppFooter/RenderMarkdown 接入 |
| C3 | src/composables/useCopyFeedback.ts + RenderMarkdown/BackToTop 使用 |
| C4 | MobileDrawer 未抽（offcanvas 样式差异评估后保留）；滚动锁语义保持 |
| C5 | src/components/widgets/FloatingButton.vue（BackToTop/TocDrawer 共享，.floating-btn 统一样式） |
| C6 | common/PostCard.vue + common/Pagination.vue + common/EmptyState.vue（PostList 519→182 行，补空态） |
| C7 | common/EmptyState.vue + common/ModalOverlay.vue（SearchModal 滚动锁/Esc 全局/焦点还原） |
| C8 | RenderMarkdown 标题 id 与 OnThisPage 契约保持（轮询保留，风险已知） |
| C9 | config `HEADER_OFFSET` + scroll.ts `scrollToHeading`（reduce-motion/锁延迟）+ 4 组件收敛 |
| C10 | FloatingButton（useFloatingButton 封装）+ useCopyFeedback + 事件常量 |
| 2.1 | Article.vue 局部收敛（usePageMeta/useTagNavigation/toArticle/normalizeArticleKey）；组件级拆分（MetaHeader/PrevNext）保留 |
| 2.2 | navigation.ts `toArticle`/`normalizeArticleKey` + utils/articles.ts `flattenCategoryArticles` + useTagNavigation |
| 2.4 | composables/usePageMeta.ts + i18n meta.* 键 + config `BLOG_TITLE`（6 视图接入，补 About/Resource title） |
| 2.5 | utils/pagination.ts / format.ts / tags.ts / url.ts（getVisiblePages 等纯函数下沉，2790 组差分验证） |
| T1 | scripts/lib/llm.ts（loadLlmConfig/getLlmClient/completeChat 三合一，translator/tagMerger 复用） |
| T2 | tagMerger `extractJsonObject(text, key)` 合并 |
| T3 | scripts/lib/frontmatter.ts `rewriteFrontmatterTags`（行内/块状双语法，28 md 逐字节一致） |
| T4 | scripts/lib/tags.ts `countTags`/`sortTagsByName` |
| T5 | scripts/lib/yamlEntries.ts `normalizeProjectTopicEntry`（16 字段三合一） |
| T6 | generateProjects 通用骨架 + generateTopics 薄壳（18 JSON 零差异） |
| T7 | core `walkCategoryArticles` + searchIndex/sitemap 使用（零差异） |
| T8+T11 | scripts/lib/text.ts（markdownToPlain/countWordsSmart/stripHtmlToText） |
| T9 | scripts/lib/{fs,text,cli,frontmatter}.ts + core/index.ts 兼容 barrel（46 文件零差异） |
| T10 | scripts/lib/fs.ts `walkAsync`（translator 遍历收敛） |
| T12 | runAllGenerators `localeSteps`/`runLocaleBlock`（zh/en 双块收敛） |
| T13 | core `requireJsonArray`（generatePosts/Tags 前置检查收敛） |
| T14 | scripts/lib/tagMapping.ts `loadTagMapping` 唯一读取 |
| T15 | B7（status 口径） |
| X1 | types.ts 单一化：SearchDoc 入 types；generateAbout/Resources 用 src 类型；Post.category 修正 |
| X2 | config `ARTICLE_ROUTE_PREFIX`（navigation/router/App/SearchModal/AppHeader/NavigationTree/PostList/generateCategories/generateSearchIndex 全量替换） |
| X3 | dataConfig `EN_SUFFIX` + src config `EN_SUFFIX` + contentLoader/translatorConfig/tagMerger 使用 |
| X4 | B1（findNotesSection 结构特征定位） |
| 2.2 i18n | usePageMeta + schema + BLOG_TITLE（见上） |
