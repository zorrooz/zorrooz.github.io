# src/utils/ 工具文档

这个目录包含了项目中使用的各种工具函数和脚本，用于内容加载、Markdown 处理、数据生成和翻译等功能。

## 工具概述

### 1. contentLoader.js
统一的内容加载工具函数，根据当前语言动态加载对应的 JSON 或 MD 文件。

**主要功能：**
- `getCurrentLocale()`: 获取当前语言环境 ('zh-CN' 或 'en-US')
- `getLocalizedFileName(baseName, extension)`: 根据语言生成文件名
- `loadJsonContent(fileName, directory)`: 动态加载 JSON 内容
- `loadMarkdownContent(filePath)`: 动态加载 Markdown 内容
- `loadMultipleJsonContents(fileNames, directory)`: 批量加载多个 JSON 文件
- 快捷加载方法：`loadCategories()`, `loadPosts()`, `loadNotes()` 等

**使用示例：**
```javascript
import { loadPosts } from './contentLoader.js';

const posts = await loadPosts();
// 加载 posts.json 或 posts-en.json 根据当前语言
```

### 2. markdownProcessor.js
Markdown 处理工具，使用 unified 库来解析和渲染 Markdown 内容。

**主要功能：**
- 支持 GitHub Flavored Markdown (GFM)
- 支持数学公式 (KaTeX)
- 支持代码语法高亮 (highlight.js)
- 支持主题切换的语法高亮样式

**使用示例：**
```javascript
import { renderMarkdown } from './markdownProcessor.js';

const html = await renderMarkdown('# Hello World');
// 返回渲染后的 HTML 字符串
```

### 3. runAllGenerators.js
统一生成器运行脚本，用于批量生成项目所需的数据文件。

**执行顺序：**
1. notes, projects, topics (中文)
2. categories (依赖于上述)
3. posts (依赖于 notes)
4. tags (依赖于 posts)
5. 重复上述步骤生成英文版本

**如何使用：**
在项目根目录下运行：
```bash
node src/utils/runAllGenerators.js
```

这个脚本会自动生成所有必要的 JSON 文件到 `src/content/` 目录。

### 4. generators/ 目录
包含各种数据生成器脚本：

- `generateAbout.js`: 生成 about.json
- `generateCategories.js`: 生成 categories.json
- `generateNotes.js`: 生成 notes.json
- `generatePosts.js`: 生成 posts.json
- `generateProjects.js`: 生成 projects.json
- `generateResources.js`: 生成 resources.json
- `generateTags.js`: 生成 tags.json
- `generateTopics.js`: 生成 topics.json

这些脚本通常通过 `runAllGenerators.js` 调用，也可以单独运行：
```bash
node src/utils/generators/generateNotes.js
```

### 5. translator/ 目录
翻译工具，用于处理多语言内容。

- `translator.js`: 翻译核心逻辑
- `translatorCli.js`: 翻译命令行工具
- `translatorConfig.js`: 翻译配置

**如何使用翻译工具：**
```bash
node src/utils/translator/translatorCli.js [options]
```

具体选项请参考 `translatorCli.js` 文件中的帮助信息。

## 注意事项

- 所有工具都设计为不改变现有功能，只统一代码风格（通过 Prettier 格式化）。
- 生成器脚本依赖于 `src/content-src/` 目录中的源文件。
- 内容加载器支持中英文切换，会自动选择对应的文件。
- Markdown 处理器会根据当前主题（明暗）动态加载语法高亮样式。