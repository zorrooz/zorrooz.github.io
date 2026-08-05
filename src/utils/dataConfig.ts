/**
 * 数据目录统一配置（解耦 code/data 的唯一接入点）
 *
 * 数据放在独立的 data 分支（本地通过 git worktree 挂到同级目录 ../blog-data）。
 * 代码分支不包含任何数据文件；所有生成器/翻译器/Vite 别名都从这里取路径。
 *
 * 覆盖方式：环境变量 GBLOG_DATA_DIR（CI 里指向 data 分支的检出目录）。
 */

import path from 'path'

// dataConfig.ts 位于 src/utils/ → 仓库根 = 其上两层
const repoRoot = path.resolve(import.meta.dirname, '../..')

export const dataDir = process.env.GBLOG_DATA_DIR
  ? path.resolve(process.env.GBLOG_DATA_DIR)
  : path.resolve(repoRoot, '../blog-data')

/** 手写源（md/yaml） */
export const contentSrcDir = path.join(dataDir, 'content-src')
/** 网站消费的生成 JSON + html */
export const contentDir = path.join(dataDir, 'content')
/** 机器维护的持久中间状态（tag-mapping、翻译状态、未来一切 LLM 缓存） */
export const cacheDir = path.join(dataDir, 'cache')