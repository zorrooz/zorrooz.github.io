/**
 * 数据目录统一配置（解耦 code/data 的唯一接入点）
 *
 * 数据放在独立的 data 分支（本地通过 git worktree 挂到同级目录 ../blog-data）。
 * 代码分支不包含任何数据文件；所有生成器/翻译器/Vite 别名都从这里取路径。
 *
 * 三层数据模型：
 *   content-src/  第一层 src   — 纯手写中文源（严格不含任何机器产物）
 *   cache/        第二层 cache — 机器维护的持久态（入库）：英文翻译源 en/、
 *                               tag-mapping.json、.translate-state.json
 *   content/      第三层 final — 可再生的最终数据（gitignore，不入库）
 *
 * 覆盖方式：环境变量 GBLOG_DATA_DIR（CI 里指向 data 分支的检出目录）。
 */

import path from 'path'

const repoRoot = path.resolve(import.meta.dirname, '..')

export const dataDir = process.env.GBLOG_DATA_DIR
  ? path.resolve(process.env.GBLOG_DATA_DIR)
  : path.resolve(repoRoot, '../blog-data')

/** 第一层：手写中文源（md/yaml），仅人类编辑 */
export const contentSrcDir = path.join(dataDir, 'content-src')
/** 第三层：网站消费的生成 JSON + html（可再生成，不入库） */
export const contentDir = path.join(dataDir, 'content')
/** 第二层：机器维护的持久中间状态（英文翻译源、tag 映射、翻译状态） */
export const cacheDir = path.join(dataDir, 'cache')

/** 第二层内的英文翻译源：镜像 content-src 结构，文件名沿用 -en 身份后缀 */
export const enSrcDir = path.join(cacheDir, 'en')

/** locale 身份后缀：zh-CN 无后缀，en-US 为 -en（-en 是内容身份而非文件位置） */
export const EN_SUFFIX = '-en'

export const localeSuffix = (locale: string): string => (locale === 'zh-CN' ? '' : EN_SUFFIX)

/** 某一 locale 的源根目录：zh 读手写源，en 读机器翻译缓存 */
export const srcDirFor = (locale: string): string => (locale === 'en-US' ? enSrcDir : contentSrcDir)
