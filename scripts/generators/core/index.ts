/**
 * 内容生成器共享工具（兼容 barrel）。
 * 实现已按职责拆分迁移至 scripts/lib/{fs,text,cli,frontmatter}.ts，
 * 此处仅 re-export，导出名与签名保持不变，消费方 import 路径无需改动。
 */

export {
  ensureDirectoryExistence,
  walk,
  readJson,
  readYaml,
  writeJsonFile,
  requireJsonArray,
  toPosixRelativeNoExt,
  safeArray,
  mdFileFilter,
  findNotesSection,
  walkCategoryArticles,
} from '../../lib/fs.ts'
export { markdownToPlain, countWordsSmart } from '../../lib/text.ts'
export { isDirectRun, runCliScript, logWriteSuccess } from '../../lib/cli.ts'
export { parseFrontMatterAndBody, normalizeTags } from '../../lib/frontmatter.ts'
