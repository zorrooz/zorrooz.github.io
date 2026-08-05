import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
import { contentDir } from '../dataConfig.ts'
import { ensureDirectoryExistence } from './core/index.ts'

function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en'
  return {
    postsJsonPath: path.join(contentDir, `posts${suffix}.json`),
    tagsJsonPath: path.join(contentDir, `tags${suffix}.json`),
  }
}

interface TagEntry {
  name: string
  count: number
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

  if (!fs.existsSync(paths.postsJsonPath)) {
    console.error(
      `posts${locale === 'zh-CN' ? '' : '-en'}.json not found at ${paths.postsJsonPath}. Please generate posts.json first.`,
    )
    process.exitCode = 1
    return
  }

  let posts: unknown = []
  try {
    const raw = fs.readFileSync(paths.postsJsonPath, 'utf-8').trim()
    posts = JSON.parse(raw)
    if (!Array.isArray(posts)) {
      console.error(`posts${locale === 'zh-CN' ? '' : '-en'}.json is not an array. Abort.`)
      process.exitCode = 1
      return
    }
  } catch (e) {
    console.error(
      `Failed to read/parse posts${locale === 'zh-CN' ? '' : '-en'}.json:`,
      e instanceof Error ? e.message : e,
    )
    process.exitCode = 1
    return
  }

  const tags = generateTagsFromPosts(posts as Array<{ tags?: unknown }>)

  try {
    const targetPath = paths.tagsJsonPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(tags, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath} (${tags.length} tags)`)
  } catch (e) {
    console.error(`Failed to write ${locale} tags.json:`, e instanceof Error ? e.message : e)
    process.exitCode = 1
  }
}

function main() {
  console.log('Starting tags.json generation script...')
  generateTagsJson('zh-CN')
  generateTagsJson('en-US')
  console.log('tags.json generation complete.')
}

/**
 * 中英标签一致性校验（构建时自动运行）：
 * zh/en tags.json 的数量与名称必须一致（翻译器保留 frontmatter tags 原文，
 * 不一致说明内容数据有问题——如手写 -en.md 或 yaml tags 被改动）。
 */
export function checkTagsConsistency(contentDir: string): void {
  const readNames = (suffix: string): string[] => {
    const p = path.join(contentDir, `tags${suffix}.json`)
    if (!fs.existsSync(p)) return []
    try {
      const arr = JSON.parse(fs.readFileSync(p, 'utf-8'))
      return Array.isArray(arr)
        ? arr.map((t: { name?: unknown }) => String(t?.name ?? '')).filter(Boolean)
        : []
    } catch {
      return []
    }
  }

  const zh = readNames('')
  const en = readNames('-en')

  // 每篇文章 tags 数量配对检查（posts 按 id 配对）
  const readPosts = (suffix: string): Array<{ id?: unknown; tags?: unknown }> => {
    const p = path.join(contentDir, `posts${suffix}.json`)
    try {
      const arr = JSON.parse(fs.readFileSync(p, 'utf-8'))
      return Array.isArray(arr) ? (arr as Array<{ id?: unknown; tags?: unknown }>) : []
    } catch {
      return []
    }
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
    console.warn('  提示：md frontmatter 的 tags 由 zh→en 映射翻译，数量必须守恒；请检查 -en.md 或映射文件')
  } else {
    console.log(`[OK] 中英标签一致 (${zh.length} = ${en.length})`)
  }
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateTagsJson }
