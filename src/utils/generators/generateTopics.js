import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const contentSrcDir = path.join(__dirname, '../../content-src')
const topicsSrcDir = path.join(contentSrcDir, 'topics')
const contentOutputDir = path.join(__dirname, '../../content')

function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en'
  return {
    outputPath: path.join(contentOutputDir, `topics${suffix}.json`),
  }
}

function ensureDirectoryExistence(filePath) {
  const dirname = path.dirname(filePath)
  if (!fs.existsSync(dirname)) {
    fs.mkdirSync(dirname, { recursive: true })
  }
}

function walk(dir, predicate = () => true, acc = []) {
  if (!fs.existsSync(dir)) return acc
  const items = fs.readdirSync(dir, { withFileTypes: true })
  for (const it of items) {
    const full = path.join(dir, it.name)
    if (it.isDirectory()) walk(full, predicate, acc)
    else if (predicate(full)) acc.push(full)
  }
  return acc
}

function markdownToPlain(text) {
  let t = text.replace(/```[\s\S]*?```/g, ' ')
  t = t.replace(/`[^`]*`/g, ' ')
  t = t.replace(/!\[[^\]]*]\([^)]+\)/g, ' ')
  t = t.replace(/\[[^\]]*]\([^)]+\)/g, ' ')
  t = t.replace(/<[^>]+>/g, ' ')
  t = t.replace(/^#{1,6}\s+/gm, ' ')
  t = t.replace(/[*_~`>#|-]{1,}/g, ' ')
  t = t.replace(/^\s*\[[^\]]+]:\s+\S+.*$/gm, ' ')
  t = t.replace(/\s+/g, ' ').trim()
  return t
}
function countWordsSmart(text) {
  const cjk = (text.match(/[\u4E00-\u9FFF\u3400-\u4DBF]/g) || []).length
  const words = (text.match(/[A-Za-z0-9]+(?:'[A-Za-z0-9]+)?/g) || []).length
  return cjk + words
}
function toPosixRelativeNoExt(fullPath, baseDir) {
  const rel = path.relative(baseDir, fullPath)
  const noExt = rel.replace(/\.[^/.]+$/, '')
  return noExt.split(path.sep).join('/')
}
function extractTitleFromH1(raw) {
  const lines = raw.split(/\r?\n/)
  for (const line of lines) {
    const m = line.match(/^#\s+(.*)$/)
    if (m) return m[1].trim()
  }
  return ''
}

function buildTopicItem(mdPath) {
  const raw = fs.readFileSync(mdPath, 'utf-8')
  const title = extractTitleFromH1(raw) || path.basename(mdPath).replace(/\.md$/i, '')
  const plain = markdownToPlain(raw)
  const wordCount = countWordsSmart(plain)
  const relativePath = toPosixRelativeNoExt(mdPath, topicsSrcDir)
  return { title, relativePath, wordCount }
}

function generateTopicsJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)

  if (!fs.existsSync(topicsSrcDir)) {
    console.warn(`Warn: topics source directory not found at ${topicsSrcDir}`)
  }

  const mdFiles = walk(topicsSrcDir, (p) => {
    if (locale === 'zh-CN') {
      return /\.md$/i.test(p) && !p.endsWith('-en.md')
    } else {
      return /-en\.md$/i.test(p)
    }
  })

  const items = mdFiles.map(buildTopicItem)
  items.sort((a, b) => a.relativePath.localeCompare(b.relativePath))
  try {
    const targetPath = paths.outputPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(items, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath} (${items.length} topics)`)
  } catch (e) {
    console.error(`Failed to write ${locale} topics.json:`, e)
  }
}

function main() {
  console.log('Starting topics.json generation script...')
  generateTopicsJson('zh-CN')
  generateTopicsJson('en-US')
  console.log('topics.json generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateTopicsJson }
