// @ts-check
import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'

import {
  ensureDirectoryExistence,
  walk,
  markdownToPlain,
  countWordsSmart,
  toPosixRelativeNoExt,
  extractTitleFromH1,
} from './core/index.js'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const contentSrcDir = path.join(__dirname, '../../content-src')
const projectsSrcDir = path.join(contentSrcDir, 'projects')
const contentOutputDir = path.join(__dirname, '../../content')

function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en'
  return {
    outputPath: path.join(contentOutputDir, `projects${suffix}.json`),
  }
}

/**
 * @param {string} mdPath
 * @returns {{ title: string, relativePath: string, wordCount: number }}
 */
function buildProjectItem(mdPath) {
  const raw = fs.readFileSync(mdPath, 'utf-8')
  const title = extractTitleFromH1(raw) || path.basename(mdPath).replace(/\.md$/i, '')
  const plain = markdownToPlain(raw)
  const wordCount = countWordsSmart(plain)
  const relativePath = toPosixRelativeNoExt(mdPath, projectsSrcDir)

  return { title, relativePath, wordCount }
}

function generateProjectsJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)

  if (!fs.existsSync(projectsSrcDir)) {
    console.warn(`Warn: projects source directory not found at ${projectsSrcDir}`)
  }

  const mdFiles = walk(projectsSrcDir, (p) => {
    if (locale === 'zh-CN') {
      return /\.md$/i.test(p) && !p.endsWith('-en.md')
    } else {
      return /-en\.md$/i.test(p)
    }
  })

  const items = mdFiles.map(buildProjectItem)

  items.sort((a, b) => a.relativePath.localeCompare(b.relativePath))

  try {
    const targetPath = paths.outputPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(items, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath} (${items.length} projects)`)
  } catch (e) {
    console.error(`Failed to write ${locale} projects.json:`, e)
  }
}

function main() {
  console.log('Starting projects.json generation script...')
  generateProjectsJson('zh-CN')
  generateProjectsJson('en-US')
  console.log('projects.json generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateProjectsJson }
