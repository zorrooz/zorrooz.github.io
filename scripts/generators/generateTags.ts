import fs from 'fs'
import path from 'path'
import { contentDir, localeSuffix } from '../dataConfig.ts'
import {
  readJson,
  writeJsonFile,
  isDirectRun,
  runCliScript,
  logWriteSuccess,
} from './core/index.ts'
import type { Tag as TagEntry } from '../../src/types/content.ts'

function getFilePaths(locale = 'zh-CN') {
  const suffix = localeSuffix(locale)
  return {
    postsJsonPath: path.join(contentDir, `posts${suffix}.json`),
    tagsJsonPath: path.join(contentDir, `tags${suffix}.json`),
  }
}

function generateTagsFromPosts(posts: Array<{ tags?: unknown }>): TagEntry[] {
  const map = new Map<string, number>()
  for (const p of posts) {
    const tags = Array.isArray(p?.tags) ? p.tags : []
    for (const raw of tags) {
      const name = typeof raw === 'string' ? raw.trim() : ''
      if (!name) continue
      map.set(name, (map.get(name) || 0) + 1)
    }
  }
  return Array.from(map.entries())
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => a.name.localeCompare(b.name, 'zh-Hans-CN'))
}

function generateTagsJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)
  const postsFileName = `posts${localeSuffix(locale)}.json`

  if (!fs.existsSync(paths.postsJsonPath)) {
    console.error(
      `${postsFileName} not found at ${paths.postsJsonPath}. Please generate posts.json first.`,
    )
    process.exitCode = 1
    return
  }

  const parsed = readJson(paths.postsJsonPath)
  if (parsed === null) {
    console.error(`Failed to read/parse ${postsFileName}:`)
    process.exitCode = 1
    return
  }
  if (!Array.isArray(parsed)) {
    console.error(`${postsFileName} is not an array. Abort.`)
    process.exitCode = 1
    return
  }

  const tags = generateTagsFromPosts(parsed as Array<{ tags?: unknown }>)

  try {
    const targetPath = paths.tagsJsonPath
    writeJsonFile(targetPath, tags)
    logWriteSuccess(targetPath, `${tags.length} tags`)
  } catch (e) {
    console.error(`Failed to write ${locale} tags.json:`, e instanceof Error ? e.message : e)
    process.exitCode = 1
  }
}

if (isDirectRun(import.meta)) {
  runCliScript('tags.json', () => {
    generateTagsJson('zh-CN')
    generateTagsJson('en-US')
  })
}

/**
 * 中英标签一致性校验（构建时自动运行）：
 * zh/en tags.json 的数量与名称必须一致（翻译器保留 frontmatter tags 原文，
 * 不一致说明内容数据有问题——如手写 -en.md 或 yaml tags 被改动）。
 */
export function checkTagsConsistency(contentDir: string): void {
  const readNames = (suffix: string): string[] => {
    const arr = readJson(path.join(contentDir, `tags${suffix}.json`))
    return Array.isArray(arr)
      ? arr.map((t: { name?: unknown }) => String(t?.name ?? '')).filter(Boolean)
      : []
  }

  const zh = readNames('')
  const en = readNames('-en')

  // 每篇文章 tags 数量配对检查（posts 按 id 配对）
  const readPosts = (suffix: string): Array<{ id?: unknown; tags?: unknown }> => {
    const arr = readJson(path.join(contentDir, `posts${suffix}.json`))
    return Array.isArray(arr) ? (arr as Array<{ id?: unknown; tags?: unknown }>) : []
  }
  const zhPosts = readPosts('')
  const enPosts = readPosts('-en')
  const mismatches: string[] = []
  for (const p of zhPosts) {
    const e = enPosts.find((x) => x?.id === p?.id)
    const zc = Array.isArray(p?.tags) ? p.tags.length : 0
    const ec = Array.isArray(e?.tags) ? e.tags.length : 0
    if (zc !== ec) mismatches.push(`id=${p?.id}: zh=${zc} en=${ec}`)
  }

  if (zh.length !== en.length || mismatches.length) {
    console.warn(`[Warn] 中英标签不一致: 总数 zh=${zh.length} en=${en.length}`)
    if (mismatches.length) console.warn(`  [文章级差异] ${mismatches.join('; ')}`)
    console.warn(
      '  提示：md frontmatter 的 tags 由 zh→en 映射翻译，数量必须守恒；请检查 -en.md 或映射文件',
    )
  } else {
    console.log(`[OK] 中英标签一致 (${zh.length} = ${en.length})`)
  }
}

export { generateTagsJson }
