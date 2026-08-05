import fs from 'fs'
import path from 'path'
import { contentDir, localeSuffix } from '../dataConfig.ts'
import {
  readJson,
  safeArray,
  writeJsonFile,
  isDirectRun,
  runCliScript,
  logWriteSuccess,
} from './core/index.ts'

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

export interface PostItem {
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
  if (!Array.isArray(categoriesArr)) return map
  const notesSection = categoriesArr.find((s) => s && s.title === '笔记' && Array.isArray(s.items))
  if (!notesSection) return map
  for (const item of notesSection.items) {
    if (!item || typeof item !== 'object') continue
    const name = (item as Record<string, unknown>).name
    if (typeof name !== 'string') continue
    const groupTitle = typeof item.title === 'string' && item.title.trim() ? item.title : name
    const subMap =
      item.categories && typeof item.categories === 'object'
        ? (item.categories as Record<string, unknown>)
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

  if (!fs.existsSync(paths.notesJsonPath)) {
    console.error(
      `${notesFileName} not found at ${paths.notesJsonPath}. Please run generateNotes.ts first.`,
    )
    process.exitCode = 1
    return
  }

  const parsed = readJson(paths.notesJsonPath)
  if (parsed === null) {
    console.error(`Failed to read/parse ${notesFileName}:`)
    process.exitCode = 1
    return
  }
  if (!Array.isArray(parsed)) {
    console.error(`${notesFileName} is not an array. Abort.`)
    process.exitCode = 1
    return
  }
  const notesArr = parsed as Record<string, unknown>[]

  const categoriesArr = safeArray(readJson(paths.categoriesJsonPath))
  const notesCategoryMap = buildNotesCategoryMap(categoriesArr)

  const posts = buildPostsFromNotes(notesArr, notesCategoryMap)

  try {
    const targetPath = paths.postsJsonPath
    writeJsonFile(targetPath, posts)
    logWriteSuccess(targetPath, `${posts.length} posts`)
  } catch (e) {
    console.error(`Failed to write ${locale} posts.json:`, e instanceof Error ? e.message : e)
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
