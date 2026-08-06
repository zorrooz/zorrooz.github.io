import fs from 'fs'
import path from 'path'

import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import {
  walk,
  normalizeTags,
  parseFrontMatterAndBody,
  markdownToPlain,
  countWordsSmart,
  toPosixRelativeNoExt,
  mdFileFilter,
  writeJsonFile,
  isDirectRun,
  runCliScript,
  logWriteSuccess,
} from './core/index.ts'
import type { Note as NoteItem } from '../../src/types.ts'

function getFilePaths(locale = 'zh-CN') {
  return {
    // 源目录按 locale 分层：zh 扫手写源 content-src/notes，en 扫机器层 cache/en/notes
    notesSrcDir: path.join(srcDirFor(locale), 'notes'),
    outputPath: path.join(contentDir, `notes${localeSuffix(locale)}.json`),
  }
}

function buildNoteItem(mdPath: string, notesSrcDir: string): NoteItem {
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
  const notesSrcDir = paths.notesSrcDir

  if (!fs.existsSync(notesSrcDir)) {
    console.warn(`Warn: notes source directory not found at ${notesSrcDir}`)
  }

  // 目录即契约：content-src/notes 只含中文源（排除历史遗留 -en），cache/en/notes 全部是英文
  const mdFiles = walk(notesSrcDir, mdFileFilter(locale))

  const items = mdFiles.map((p) => buildNoteItem(p, notesSrcDir))

  items.sort((a, b) => {
    const ad = Date.parse(a.date || '') || 0
    const bd = Date.parse(b.date || '') || 0
    if (ad !== bd) return bd - ad
    return a.relativePath.localeCompare(b.relativePath)
  })

  try {
    const targetPath = paths.outputPath
    writeJsonFile(targetPath, items)
    logWriteSuccess(targetPath, `${items.length} notes`)
  } catch (e) {
    console.error(`Failed to write ${locale} notes.json:`, e)
  }
}

if (isDirectRun(import.meta)) {
  runCliScript('notes.json', () => {
    generateNotesJson('zh-CN')
    generateNotesJson('en-US')
  })
}

export { generateNotesJson }
