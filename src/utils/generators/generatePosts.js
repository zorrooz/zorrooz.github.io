import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const contentDir = path.join(__dirname, '../../content')

function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en'
  return {
    notesJsonPath: path.join(contentDir, `notes${suffix}.json`),
    postsJsonPath: path.join(contentDir, `posts${suffix}.json`),
    categoriesJsonPath: path.join(contentDir, `categories${suffix}.json`),
  }
}

function ensureDirectoryExistence(filePath) {
  const dir = path.dirname(filePath)
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true })
}

function toArray(x) {
  return Array.isArray(x) ? x : x == null ? [] : [x]
}

function deriveCategory(relativePath) {
  if (typeof relativePath !== 'string' || !relativePath) return []
  const parts = relativePath.split('/')
  if (parts.length >= 2) return [parts[0], parts[1]]
  if (parts.length === 1) return [parts[0]]
  return []
}

function buildNotesCategoryMap(categoriesArr) {
  const map = new Map()
  if (!Array.isArray(categoriesArr)) return map
  const notesSection = categoriesArr.find((s) => s && s.title === '笔记' && Array.isArray(s.items))
  if (!notesSection) return map
  for (const item of notesSection.items) {
    if (!item || typeof item !== 'object') continue
    const name = item.name
    if (typeof name !== 'string') continue
    const groupTitle = typeof item.title === 'string' && item.title.trim() ? item.title : name
    const subMap = item.categories && typeof item.categories === 'object' ? item.categories : {}
    map.set(name, { groupTitle, subMap })
  }
  return map
}

function buildPostsFromNotes(notesArr, notesCategoryMap) {
  const posts = []
  let id = 1

  for (const it of notesArr) {
    if (!it || typeof it !== 'object') continue
    if (it.draft === true) continue

    const title = typeof it.title === 'string' ? it.title : ''
    const date = typeof it.date === 'string' ? it.date : ''
    const tags = Array.isArray(it.tags) ? it.tags : []
    const preview = typeof it.description === 'string' ? it.description : ''
    const [grp, sub] = deriveCategory(it.relativePath)
    let category = []
    if (grp) {
      const entry = notesCategoryMap && notesCategoryMap.get(grp)
      if (entry) {
        const first = entry.groupTitle || grp
        const second = sub ? entry.subMap?.[sub] || sub : undefined
        category = second ? [first, second] : [first]
      } else {
        category = sub ? [grp, sub] : [grp]
      }
    }
    const wordCount = Number.isFinite(it.wordCount) ? it.wordCount : 0

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

  if (!fs.existsSync(paths.notesJsonPath)) {
    console.error(
      `notes${locale === 'zh-CN' ? '' : '-en'}.json not found at ${paths.notesJsonPath}. Please run generateNotes.js first.`,
    )
    process.exitCode = 1
    return
  }

  let notesArr = []
  try {
    const raw = fs.readFileSync(paths.notesJsonPath, 'utf-8').trim()
    notesArr = JSON.parse(raw)
    if (!Array.isArray(notesArr)) {
      console.error(`notes${locale === 'zh-CN' ? '' : '-en'}.json is not an array. Abort.`)
      process.exitCode = 1
      return
    }
  } catch (e) {
    console.error(
      `Failed to read/parse notes${locale === 'zh-CN' ? '' : '-en'}.json:`,
      e?.message || e,
    )
    process.exitCode = 1
    return
  }

  let categoriesArr = []
  try {
    if (fs.existsSync(paths.categoriesJsonPath)) {
      const rawCat = fs.readFileSync(paths.categoriesJsonPath, 'utf-8').trim()
      categoriesArr = JSON.parse(rawCat)
    }
  } catch {
    categoriesArr = []
  }
  const notesCategoryMap = buildNotesCategoryMap(categoriesArr)

  const posts = buildPostsFromNotes(notesArr, notesCategoryMap)

  try {
    const targetPath = paths.postsJsonPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(posts, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath} (${posts.length} posts)`)
  } catch (e) {
    console.error(`Failed to write ${locale} posts.json:`, e?.message || e)
    process.exitCode = 1
  }
}

function main() {
  console.log('Starting posts.json generation script...')
  generatePostsJson('zh-CN')
  generatePostsJson('en-US')
  console.log('posts.json generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generatePostsJson }
