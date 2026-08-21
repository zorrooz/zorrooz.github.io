# gblog × Obsidian 本地整合适配计划

> 目标：把 Obsidian 变成 gblog 的**本地创作前端**，保留现有构建/发布管线不动，
> 让「在 Obsidian 里写作 → 本地秒级预览 → git 自动同步 → CI 自动部署」成为日常闭环。
>
> 状态：**P0/P1 已实施**（2026-06）——扁平化结构、assets 集中图片、`npm run data:pack` 打包命令、
> 旧 URL 重定向、gitignore、模板均已落地；P2（callout/高亮语法插件）待做。
> 适用范围：`gblog`（代码仓库）+ `blog-data`（数据仓库）

---

## 0. 结论先行（TL;DR）

- **推荐方案**：Obsidian 仅作为**编辑器**，vault 直接指向数据分支的 `content-src/`；
  站点管线（Vue + 生成器 + GitHub Actions）**完全不动**。这是投入产出比最高的路线。
- **零代码改动即可用的部分**（Phase 0）：打开 vault、忽略 `.obsidian/`、配置 Obsidian
  链接格式与附件目录 —— 当天就能用，因为本仓库内容模型与 Obsidian 天然兼容（见 §1）。
- **必须的小代码改动只有 2 处**：`contentDevPlugin` watcher 忽略 `.obsidian/`（防重生成风暴），
  `blog-data/.gitignore` 增加 `.obsidian/`（防 `npm run deploy` 误提交机器配置）。
- **可选增强**（Phase 2）：在 `markdownProcessor.ts` 的 unified 管道里加 remark 插件，
  补上 Obsidian 特有语法（callout / `==高亮==` / `%%注释%%`）的渲染；wikilink 按需二选一。
- **不做**：不引入 Quartz / Obsidian Publish / 自定义 Obsidian 插件重度开发。

---

## 1. 现状盘点：哪些已天然兼容，哪些是差距

### 1.1 已兼容（好消息，几乎不用改）

| 站点现状 | Obsidian 对应 | 兼容性 |
|---|---|---|
| 文章 = `content-src/notes/**/*.md` + YAML frontmatter | Obsidian 笔记 = md + Properties（即 frontmatter） | ✅ 文件模型一致 |
| frontmatter：`title/date/author/tags/draft/description` | Properties 的 text/date/list/checkbox 类型 | ✅ 语义一一对应；Obsidian date 类型输出 `YYYY-MM-DD` 与现状一致 |
| `remark-gfm`：表格/任务列表/脚注 | GFM 全支持 | ✅ |
| `remark-breaks`：单换行 = `<br>` | Obsidian 默认「严格换行关闭」= 单换行渲染 `<br>` | ✅ 行为恰好一致 |
| `remark-math` + `rehype-katex`：`$...$`/`$$...$$` | Obsidian 原生 LaTeX | ✅ |
| 相对路径图片（`![](fig.png)`） | Obsidian 附件（同目录模式） | ✅ `RenderMarkdown.vue` 已用 `import.meta.glob` 收集 `content-src/cache/en` 图片并在开发/构建期解析 |
| `npm run dev` 数据热更新 | Obsidian 保存即写盘 | ✅ `contentDevPlugin` 已实现「变更 → runAllGenerators → full-reload」 |
| 文章正文 = 构建期预渲染 HTML（`generateHtml` → `content/html`） | — | ✅ 语法适配只需改 `markdownProcessor.ts` 一个点，搜索索引（`generateSearchIndex` 消费同一 HTML）自动同步受益 |

### 1.2 差距（需要适配的部分）

| # | 差距 | 影响 | 工作量 |
|---|---|---|---|
| G1 | `[[wikilink]]` / `![[embed]]` 无法解析（remark 无此插件） | 链接以纯文本显示 | 中（或零，见 §5.4） |
| G2 | `> [!note]` callout 渲染为普通引用块（`[!note]` 字面残留） | 样式丢失 | 小 |
| G3 | `==高亮==`、`%%注释%%` 无处理 | 高亮以字面 `==` 显示；注释外泄 | 小 |
| G4 | vault 若指向 `content-src`，`.obsidian/` 变更会触发 dev 重生成风暴 | 开发体验差、无谓 CPU | 极小 |
| G5 | `.obsidian/` 未被 data 分支 gitignore，`npm run deploy` 会误提交机器配置 | 污染「纯手写源」层 | 极小 |
| G6 | 正文 `#tag` 行内标签不会被收集（站点只读 frontmatter tags） | 标签遗漏 | 约定规避 + lint |

---

## 2. 资料调研摘要

调研了「Obsidian vault 作为静态博客内容源」的主流实践，结论归纳：

1. **Obsidian 官方定位**：笔记是可移植的 Markdown 文件，内部链接/嵌入/附件均有标准文件语义
   （[Obsidian Internal links 官方文档](https://obsidian.md/help/links)）。frontmatter（Properties）
   与静态站点生成器约定同源，这是「vault 即内容源」可行性的基础。
2. **vault 直出静态站是成熟路线**：Quartz 是典型代表，直接把 vault 文件夹发布为数字花园
   （[Hivebook 对 Quartz 的架构介绍](https://hivebook.wiki/wiki/quartz-static-site-generator-for-digital-gardens)，
   [dchua 用 Quartz + Fyra 发布 vault 的实战记录](https://dchua.com/posts/2026-06-14-publishing-your-obsidian-vault-with-quartz-and-fyra/)）。
   但这要求**换掉**现有 Vue 站，不在本计划考虑范围——我们只借鉴其「vault 组织 + 语法映射」经验。
3. **git 是 Obsidian 与站点之间的事实同步通道**：obsidian-git 插件支持定时自动 commit/push，
   配置项成熟（[obsidian-git 设置文档（DeepWiki）](https://deepwiki.com/Vinzent03/obsidian-git/6-settings-and-configuration)、
   [obsidian-git 仓库](https://github.com/Vinzent03/obsidian-git)）。中文社区也有现成的
   [Obsidian Git + 博客发布工作流](https://deepfog.top/posts/obsidian-git-workflow/) 可参考。
4. **Obsidian 语法在 remark/unified 生态有现成插件**：
   - wikilink：[`@flowershow/remark-wiki-link`](https://www.npmjs.com/package/@flowershow/remark-wiki-link)
     （维护较活跃的 fork）
   - 综合（callout/高亮/注释/嵌入等）：[`remark-obsidian-md`](https://www.npmjs.com/package/remark-obsidian-md)
     —— 一个包覆盖 G2/G3 的大部分，值得优先评估
5. **vault 作为纯文本 git 目录也利于 AI 代理写作**：如
   [markdown-vault-mcp](https://github.com/mikesusz/markdown-vault-mcp) 一类的 MCP 服务可让
   智能体直接读写 vault（本仓库已重度使用 Qoder/opencode 等代理，可作为远期加分项）。

---

## 3. 目标架构

```
┌──────────────────────────── 本地（作者机） ────────────────────────────┐
│                                                                        │
│  Obsidian（编辑器）                                                     │
│  vault = ../blog-data/content-src        ┌──────────────────────────┐  │
│  ├── notes/<Category>/<sub>/<article>/   │  npm run dev              │  │
│  │     <article>.md  ← 写作            │  contentDevPlugin:        │  │
│  ├── categories.yaml / about.yaml /      │  保存 → runAllGenerators │  │
│  │   resources.yaml                       │  → full-reload 秒级预览  │  │
│  └── .obsidian/（gitignore，不进库）      └──────────────────────────┘  │
│                                                                        │
│  obsidian-git 插件 ──auto commit/push──▶ blog-data（data 分支, git）    │
└────────────────────────────────────────────────────────────────────────┘
                                │ push origin data
                                ▼
┌──────────────────────────── CI（GitHub Actions）──────────────────────┐
│  npm run prebuild（GBLOG_NO_TRANSLATE=1 跳过翻译）→ vite-ssg build      │
│  → 部署 gh-pages  ◀── 现状流程，零改动                                  │
└────────────────────────────────────────────────────────────────────────┘
```

不变式：**三层数据模型不动**（content-src 手写 / cache 机器 / content 产物）；
**英文机器层（cache/en）不进 vault 创作视野**，只由 `npm run translate` 生成。

---

## 4. 方案对比与选择

| 方案 | 说明 | 优点 | 缺点 | 结论 |
|---|---|---|---|---|
| **A. Obsidian 当编辑器（推荐）** | vault = content-src，管线原样 | 零迁移、当天可用、保持站点定制性 | 需少量语法适配（§5-6） | ✅ |
| B. 引入 Quartz 等 Obsidian 原生 SSG | vault 直接发布，弃用 Vue 站 | Obsidian 体验最顺 | 推翻现有站（设计系统/双语/搜索全丢），重构成本高 | ❌ |
| C. 自研 Obsidian 插件深度集成 | 发布命令、状态视图等 | 体验最定制 | 维护一个插件 + 随 Obsidian API 升级，投入产出比低 | ⏸ 远期可选 |

选 A，理由：本仓库的「内容即 md + frontmatter」模型与 Obsidian 几乎零摩擦（§1.1），
差异只在少数 Obsidian 方言语法，用 remark 插件补齐即可，**不换引擎**。

---

## 5. Vault 布局与配置规范

### 5.1 vault 根目录：`content-src`（推荐）

- 所见即「可发布内容」：notes/ 树 + 4 个 yaml，与站点 URL 结构 1:1 对应，Graph view 也干净。
- `cache/`（英文层、tag 映射）与 `content/`（可再生产物）不进入创作视野，符合三层模型语义。
- **备选**：vault 根 = `blog-data` 全库 —— 好处是可顺带浏览/手修 `cache/en` 英文镜像；
  代价是看到大量机器产物（可用 Obsidian「排除文件」功能隐藏 `content/`、`cache/tag-mapping.json` 等）。
  若选此方案，`.obsidian/` 位于 blog-data 根，不在 dev watcher 监听范围内（`contentDevPlugin`
  只监听 `content-src` + `cache/en`），G4 自动消失；但仍需 G5 的 gitignore。

### 5.2 `.obsidian/` 处理（G5，必须）

- data 分支 `.gitignore` 增加：
  ```gitignore
  .obsidian/
  .trash/
  ```
- 不提交 `.obsidian`，但**在本仓库提供一份模板**（如 `vault/` 目录存 app.json /
  community-plugins.json 等最小配置 + 一份 README 说明），新机器一键拷入，避免各机配置漂移。
- 若想共享插件列表，可在模板里固定：Templater、obsidian-git、Linter（可选）。

### 5.3 附件策略（已实施：assets 集中 + 打包命令）

- **写作期**：图片统一放 `content-src/assets/`（vault 根 `assets/` 文件夹）。Obsidian 设置：
  「附件默认存放位置」→ 指定文件夹 → `assets`。
- 引用格式：`![示意图](assets/fig-1.png)`（Obsidian 默认「最短路径」格式即如此；
  `RenderMarkdown.rewriteImageLinks` 已支持 vault 根相对路径解析，见 §8）。
- **发布前**：`npm run data:pack -- notes/<cat>/<sub>/<slug>`（或 `npm run data:pack` 递归全部）把引用图复制到文章同目录并改写引用
  （同步改写 cache/en 镜像引用，不复制文件；幂等；同名冲突自动加 slug 前缀）。
- 命名规范：英文小写连字符（`fig-1.png`、`pipeline.png`），避免中文文件名。

### 5.4 链接格式（G1 的零成本解法）

- Obsidian 设置：**「新链接格式 = 基于最短路径的 Markdown 链接」**（`Files & Links → New link format
  → Shortest path when possible`，并**关闭** `Use [[Wikilinks]]`）。
- 效果：`[[bash-scripting]]` → 写作时插入标准 Markdown 链接，**管线零改动**，站点点击即跳转。
  配合 Obsidian「重命名时自动更新内部链接」，移动/改名文章链接自动维护。
- 代价：丢失 wikilink 的「按文件名解引用」便利（需引用时输入相对路径，Obsidian 有补全提示）。
- **若要保留 wikilink 体验**，见 §6.1 的增强路径（构建期解析）。

---

## 6. 语法兼容层详细设计（Phase 2 核心）

统一在 `scripts/markdownProcessor.ts` 的 unified 管道中处理（注意插件顺序：
remark 系插件必须加在 `remarkRehype` **之前**；rehype 系在它之后）。
`generateHtml` 与 `generateSearchIndex` 共用此 processor，因此**只改一处，全站生效**。

| Obsidian 语法 | 处理策略 | 实现 | 优先级 |
|---|---|---|---|
| `[[slug]]` wikilink | 见 §6.1，二选一 | — | P1 |
| `![[img.png]]` 图片嵌入 | 若采用 Markdown 链接格式，Obsidian 自动写成 `![](...)`，无需处理；若保留 wikilink 格式，构建期把 `![[x]]` 转为 `<img>`（解析同 §6.1） | remark 插件 | P1 |
| `![[other-note]]` 笔记嵌入 | **v1 不支持**：降级为普通链接（转换时提醒）；完整嵌入（构建期 inline 目标 HTML）复杂度高、收益低 | — | P3（不做） |
| `> [!note]/[!tip]/[!warning]…` callout | 渲染为带类型图标与边框的 callout 容器，样式走站点设计系统（`--primary`/`--tint` 等 token，遵循三色规范） | 优先评估 `remark-obsidian-md`；不满意则自写 ~30 行 rehype 插件（`blockquote` 首行匹配 `[!type]` → 改写为 `<div class="callout callout--type">`） | P1 |
| `==text==` 高亮 | 转 `<mark>`，样式 `--tint` 底 + `--ink` 字 | `remark-obsidian-md` 或自写 ~10 行 inline 插件 | P2 |
| `%%注释%%` | 构建期剥离（不进入 HTML/搜索索引/字数统计） | remark 插件（删除文本节点） | P1 |
| 正文 `#tag` 行内标签 | **约定规避**：标签一律写 frontmatter `tags:`；正文行内标签由 lint 报警 | lint 脚本 | P2 |
| 任务列表 `- [ ]`、脚注、表格 | remark-gfm 已支持 | — | ✅ 已有 |
| Mermaid 代码块 | Obsidian 原生渲染，站点无 | v1 降级为代码块展示；远期可引 mermaid.js 客户端渲染 | P3 |

### 6.1 wikilink 二选一决策

**路径 ①（推荐先做）：Obsidian 用 Markdown 链接格式** —— 零代码、零风险，唯一损失是 wikilink
输入体验。适合内容以「文章间互链不多」为特征的博客场景。

**路径 ②（保留 wikilink）：构建期解析** —— 引入 `@flowershow/remark-wiki-link`，配置自定义
`hrefTemplate`：构建时扫描 `notes.json`（或源目录）建立 `basename（slug）→ relativePath` 索引，
`[[slug]]` → `/zh/article/<relativePath>`；**歧义检测**（同 slug 不同分类 → 构建告警，按
categories.yaml 层级消歧或报错）。注意点：
- 解析器同时要处理 `[[slug|别名]]` 的显示文本；
- `remark-wiki-link` 默认把链接渲染为 `<span>`，需用 `hrefTemplate` 或自写 resolver 输出 `<a>`；
- 英文镜像（`cache/en`）解析同一索引，URL 由运行时 locale 前缀决定，构建期只写相对路径。

建议：**先走路径 ① 上线**；若写作中频繁感觉需要 wikilink，再花 1 天切路径 ②（两个方向都预留，
不冲突——路径 ② 的插件对普通 Markdown 链接无副作用）。

---

## 7. 同步与发布工作流（效率核心）

### 7.1 日常闭环（推荐）

```
Obsidian 写作（自动保存）
   │
   ├─▶ [本地] npm run dev：contentDevPlugin 检测变更 → 重生成 → 浏览器热刷新（秒级预览）
   │
   └─▶ [同步] obsidian-git 插件：定时（如每 5 分钟）或手动
            commit + push origin data（相当于自动执行 npm run data:deploy 的提交部分）
                 │
                 ▼
        [CI] GitHub Actions：prebuild → build → 部署 gh-pages（现状，零改动）
```

- 发布节奏策略：**日常小改走 obsidian-git 自动推送**（data 分支每次推送都触发 CI，约 1-2 分钟上线）；
  **批量/敏感改动**（改 categories.yaml 结构、删文章、调 frontmatter）走
  `npm run prebuild` 本地验证 + `npm run data:deploy` 手动推送。
- `npm run data:translate`（英文生成）仍为本地手动步骤，不进 vault 自动流程；需要发布英文时先跑一次。

### 7.2 新文章 SOP（Templater 模板）

- 提供 Templater 模板 `new-note.md`：按「目录即契约」生成
  `notes/<Category>/<subcategory>/<slug>/<slug>.md`，frontmatter 脚手架：
  ```yaml
  ---
  title: "{{title}}"
  date: "{{date:YYYY-MM-DD}}"
  author: "zorrooz"
  tags: []
  draft: true
  description: ""
  ---
  ```
- `draft: true` 默认（站点已支持 draft 过滤）→ 写完改 `false` 即发布，天然防半成品上线。

### 7.3 内容 lint（可选，P2）

新增 `scripts/lintContent.ts`（并入 prebuild 或独立 CLI）：缺失 title/date、`notes/` 目录名不在
`categories.yaml`、draft 未配对（zh draft=true 但 en 已生成）、zh/en tags 不对称、正文行内 `#tag`、
`%%注释%%` 残留等 —— 输出 `[Warn]` 不阻断构建，但让 Obsidian 侧的错误在 dev 终端即时可见。

---

## 8. 具体改动清单

### 代码仓库（gblog）

| 文件 | 改动 | 阶段 |
|---|---|---|
| `vite/contentDevPlugin.ts` | watcher 加 `ignored: [/(^|[/\\])\.(obsidian|trash)([/\\]|$)/]`（G4） | P1 |
| `scripts/markdownProcessor.ts` | 挂 remark 插件：注释剥离、callout、高亮、（可选 wikilink）（G2/G3/G1） | P2 |
| `scripts/lintContent.ts`（新增） | 内容 lint（§7.3） | P2 |
| `vault/`（新增，模板目录） | `.obsidian` 最小配置模板 + 说明（§5.2） | P1 |
| `AGENTS.md` | 增补「Obsidian 工作流」章节（vault 位置、设置项、SOP、lint 说明） | P3 |
| `package.json` | （如走路径②）加 `@flowershow/remark-wiki-link` 或 `remark-obsidian-md` | P2 |

### 数据仓库（blog-data）

| 文件 | 改动 | 阶段 |
|---|---|---|
| `.gitignore` | 增加 `.obsidian/`、`.trash/`（G5，**必须**，否则 deploy 误提交） | P0 |
| `content-src/notes/**/*-en.md`（历史遗留） | 清理或明确忽略（zh 过滤器已排除，仅观感问题） | P3（可选） |

### 前端样式（如做 callout）

- `src/assets/styles/global.scss`：`.callout` 组件样式，严格使用 `--primary/--tint/--line/--fg` 等
  现有 token（遵循三色映射规范），圆角 12px 卡片 + 左侧 3px `--primary` 标识条（与引用块语义一致）。

---

## 9. 分阶段路线图

| Phase | 内容 | 工作量 | 验收标准 |
|---|---|---|---|
| **P0 零代码启用** | vault 打开、data 分支 gitignore、Obsidian 设置（Markdown 链接格式/附件同目录/排除文件）、dev 预览验证 | 0.5 天 | 在 Obsidian 改一篇文章 → 浏览器秒级刷新；`git status` 无 `.obsidian` |
| **P1 体验硬化** | watcher 忽略 `.obsidian`；Templater 新文章模板；vault 模板目录 | 0.5-1 天 | Obsidian 连续操作不再触发重生成风暴；新文章 1 分钟出稿 |
| **P2 语法兼容** | markdownProcessor 插件：注释剥离 → callout → 高亮；评估 remark-obsidian-md；lint 脚本 | 2-3 天 | 旧文章渲染零回归（对比 prebuild 前后 HTML）；callout/高亮样式符合设计系统 |
| **P3 工作流固化** | obsidian-git 配置、AGENTS.md 文档、清理历史 `-en.md` 遗留 | 1 天 | 断网写作 → 联网自动推送 → CI 部署 → 站点更新，全程无手动 git 命令 |
| **P4 远期可选** | wikilink 路径②、Mermaid、自定义 Obsidian 发布插件、MCP 智能体写作 | 按需 | 写作体验最大化（不阻塞主线） |

---

## 10. 风险与「不做」清单

**风险与对策**
1. **英文层被手改**：vault 只含 `content-src`，`cache/en` 不在创作视野；AGENTS.md 明确
   「英文只能站点审阅，源文件只动 zh」；lint 可加「-en 文件不应出现在 content-src」检查。
2. **`npm run data:deploy` 误提交 `.obsidian/`**：P0 必须 gitignore，模板目录提供一致配置。
3. **Obsidian 写文件触发 dev 重生成**：watcher 忽略规则 + 300ms 防抖已存在；Obsidian 自动保存
   频率远低于触发阈值。
4. **新增/删除文件需 dev server 重启**：现有 `contentDevPlugin` 已自动 `server.restart()`，无感。
5. **wikilink 歧义**：若走路径②，构建期检测同 slug 不同路径并告警。

**明确不做**
- 不引入 Quartz / Hugo / Astro 等 Obsidian 原生 SSG 替换现有站（§4 方案 B）。
- 不购买/接入 Obsidian Publish。
- 不把 `cache/en`、`content/` 纳入写作流程。
- 不支持笔记级 `![[embed]]`（P3 预留，不做）。
- 不写自定义 Obsidian 插件（P4 前）。

---

## 11. 参考资料

- [Obsidian 官方：Internal links（内部链接/嵌入/附件语义）](https://obsidian.md/help/links)
- [obsidian-git 插件设置文档（DeepWiki）](https://deepwiki.com/Vinzent03/obsidian-git/6-settings-and-configuration)
- [obsidian-git 仓库（Vinzent03）](https://github.com/Vinzent03/obsidian-git)
- [Obsidian Git + 博客发布工作流（中文实践）](https://deepfog.top/posts/obsidian-git-workflow/)
- [Quartz 静态站生成器架构介绍（Hivebook Wiki）](https://hivebook.wiki/wiki/quartz-static-site-generator-for-digital-gardens)
- [用 Quartz + Fyra 发布 Obsidian Vault 的实战记录](https://dchua.com/posts/2026-06-14-publishing-your-obsidian-vault-with-quartz-and-fyra/)
- [@flowershow/remark-wiki-link（wikilink remark 插件）](https://www.npmjs.com/package/@flowershow/remark-wiki-link)
- [remark-obsidian-md（callout/高亮/注释等综合插件）](https://www.npmjs.com/package/remark-obsidian-md)
- [markdown-vault-mcp（vault 供 AI 代理读写）](https://github.com/mikesusz/markdown-vault-mcp)
- [publishmd（Markdown 发布预处理管线，思路参考）](https://pypi.org/project/publishmd/)
