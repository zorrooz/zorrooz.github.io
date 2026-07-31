import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'

import {
  ensureDirectoryExistence,
  walk,
  normalizeTags,
  parseFrontMatterAndBody,
  markdownToPlain,
  countWordsSmart,
  toPosixRelativeNoExt,
} from './core/index.ts'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const contentSrcDir = path.join(__dirname, '../../content-src')
const notesSrcDir = path.join(contentSrcDir, 'notes')
const contentOutputDir = path.join(__dirname, '../../content')

export interface NoteItem {
  title: string
  date: string
  tags: string[]
  draft: boolean
  description: string
  relativePath: string
  wordCount: number
  tagCount: number
}

function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en'
  return {
    outputPath: path.join(contentOutputDir, `notes${suffix}.json`),
  }
}

function buildNoteItem(mdPath: string): NoteItem {
  const raw = fs.readFileSync(mdPath, 'utf-8')
  const { frontmatter, body } = parseFrontMatterAndBody(raw)

  const title = typeof frontmatter?.title === 'string' ? frontmatter.title : ''
  const date = typeof frontmatter?.date === 'string' ? frontmatter.date : ''
  const description = typeof frontmatter?.description === 'string' ? frontmatter.description : ''
  const draft = typeof frontmatter?.draft === 'boolean' ? frontmatter.draft : false
  const tags = normalizeTags(frontmatter?.tags)

  const plain = markdownToPlain(body)
  const wordCount = countWordsSmart(plain)
  const tagCount = tags.length

  const relativePath = toPosixRelativeNoExt(mdPath, notesSrcDir)

  return {
    title,
    date,
    tags,
    draft,
    description,
    relativePath,
    wordCount,
    tagCount,
  }
}

function generateNotesJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)

  if (!fs.existsSync(notesSrcDir)) {
    console.warn(`Warn: notes source directory not found at ${notesSrcDir}`)
  }

  const mdFiles = walk(notesSrcDir, (p) => {
    if (locale === 'zh-CN') {
      return /\.md$/i.test(p) && !p.endsWith('-en.md')
    } else {
      return /-en\.md$/i.test(p)
    }
  })

  const items = mdFiles.map(buildNoteItem)

  items.sort((a, b) => {
    const ad = Date.parse(a.date || '') || 0
    const bd = Date.parse(b.date || '') || 0
    if (ad !== bd) return bd - ad
    return a.relativePath.localeCompare(b.relativePath)
  })

  try {
    const targetPath = paths.outputPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(items, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath} (${items.length} notes)`)
  } catch (e) {
    console.error(`Failed to write ${locale} notes.json:`, e)
  }
}

function main() {
  console.log('Starting notes.json generation script...')
  generateNotesJson('zh-CN')
  generateNotesJson('en-US')
  console.log('notes.json generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateNotesJson }
