import path from 'path'
import { contentDir, localeSuffix } from '../dataConfig.ts'
import { findNotesSection, readJson, requireJsonArray, safeArray, writeJsonFile } from '../lib/fs.ts'
import { isDirectRun, logWriteSuccess, runCliScript } from '../lib/cli.ts'
import { logError } from '../lib/log.ts'

function getFilePaths(locale = 'zh-CN') {
  const suffix = localeSuffix(locale)
  return {
    notesJsonPath: path.join(contentDir, `notes${suffix}.json`),
    postsJsonPath: path.join(contentDir, `posts${suffix}.json`),
    categoriesJsonPath: path.join(contentDir, `categories${suffix}.json`),
  }
}

interface CategoryEntry {
  groupTitle: string
  subMap: Record<string, unknown>
}

interface PostItem {
  id: number
  no: number
  title: string
  date: string
  category: string[]
  tags: unknown[]
  preview: string
  wordCount: number
}

function deriveCategory(relativePath: unknown): string[] {
  if (typeof relativePath !== 'string' || !relativePath) return []
  const parts = relativePath.split('/')
  if (parts.length >= 2) return [parts[0], parts[1]]
  if (parts.length === 1) return [parts[0]]
  return []
}

function buildNotesCategoryMap(categoriesArr: unknown): Map<string, CategoryEntry> {
  const map = new Map<string, CategoryEntry>()
  /** 按结构特征定位 notes 段（标题是翻译产物，不可作匹配依据） */
  const notesSection = findNotesSection(categoriesArr)
  if (!notesSection) return map
  const items = notesSection.items as unknown[]
  for (const item of items) {
    if (!item || typeof item !== 'object') continue
    const rec = item as Record<string, unknown>
    const name = rec.name
    if (typeof name !== 'string') continue
    const groupTitle = typeof rec.title === 'string' && rec.title.trim() ? rec.title : name
    const subMap =
      rec.categories && typeof rec.categories === 'object'
        ? (rec.categories as Record<string, unknown>)
        : {}
    map.set(name, { groupTitle, subMap })
  }
  return map
}

function buildPostsFromNotes(
  notesArr: Record<string, unknown>[],
  notesCategoryMap: Map<string, CategoryEntry>,
): PostItem[] {
  const posts: PostItem[] = []
  let id = 1

  for (const it of notesArr) {
    if (!it || typeof it !== 'object') continue
    if (it.draft === true) continue

    const title = typeof it.title === 'string' ? it.title : ''
    const date = typeof it.date === 'string' ? it.date : ''
    const tags = Array.isArray(it.tags) ? it.tags : []
    const preview = typeof it.description === 'string' ? it.description : ''
    const [grp, sub] = deriveCategory(it.relativePath)
    let category: string[] = []
    if (grp) {
      const entry = notesCategoryMap.get(grp)
      if (entry) {
        const first = entry.groupTitle || grp
        const second = sub ? (entry.subMap?.[sub] as string | undefined) || sub : undefined
        category = second ? [first, String(second)] : [first]
      } else {
        category = sub ? [grp, sub] : [grp]
      }
    }
    const wordCount =
      typeof it.wordCount === 'number' && Number.isFinite(it.wordCount) ? it.wordCount : 0

    posts.push({
      id: id,
      no: id,
      title,
      date,
      category,
      tags,
      preview,
      wordCount,
    })

    id += 1
  }

  return posts
}

function generatePostsJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)
  const notesFileName = `notes${localeSuffix(locale)}.json`

  const notesArr = requireJsonArray(
    paths.notesJsonPath,
    notesFileName,
    'Please run generateNotes.ts first.',
  )
  if (notesArr === null) return

  const categoriesArr = safeArray(readJson(paths.categoriesJsonPath))
  const notesCategoryMap = buildNotesCategoryMap(categoriesArr)

  const posts = buildPostsFromNotes(notesArr as Record<string, unknown>[], notesCategoryMap)

  try {
    const targetPath = paths.postsJsonPath
    writeJsonFile(targetPath, posts)
    logWriteSuccess(targetPath, `${posts.length} posts`)
  } catch (e) {
    logError(`写入 ${locale} posts.json 失败:`, e instanceof Error ? e.message : e)
    process.exitCode = 1
  }
}

if (isDirectRun(import.meta)) {
  runCliScript('posts.json', () => {
    generatePostsJson('zh-CN')
    generatePostsJson('en-US')
  })
}

export { generatePostsJson }
