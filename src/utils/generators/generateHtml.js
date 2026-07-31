// @ts-check
/**
 * Generates static HTML for every article at build time.
 * Same file enumeration rules as generateNotes.js:
 *   zh-CN: all .md files except `-en.md`
 *   en-US: only `-en.md` files
 * Output: src/content/html/<relative-path-without-ext>.html
 * (e.g. notes/Omics/genomics/bwa/bwa.html, .../bwa-en.html)
 */

import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'

import { renderMarkdown } from '../markdownProcessor.js'
import { ensureDirectoryExistence, walk } from './core/index.js'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const contentSrcDir = path.join(__dirname, '../../content-src')
const htmlOutputDir = path.join(__dirname, '../../content/html')

const sections = ['notes', 'projects', 'topics']

/** @param {string} p @param {string} locale */
function isTargetFile(p, locale) {
  if (locale === 'zh-CN') return /\.md$/i.test(p) && !p.endsWith('-en.md')
  return /-en\.md$/i.test(p)
}

export async function generateHtml(locale = 'zh-CN') {
  let count = 0
  for (const section of sections) {
    const srcDir = path.join(contentSrcDir, section)
    const mdFiles = walk(srcDir, (p) => isTargetFile(p, locale))

    for (const mdPath of mdFiles) {
      const raw = fs.readFileSync(mdPath, 'utf-8')
      const html = await renderMarkdown(raw)

      const relNoExt = path
        .relative(contentSrcDir, mdPath)
        .replace(/\.[^/.]+$/, '')
        .split(path.sep)
        .join('/')

      const outPath = path.join(htmlOutputDir, `${relNoExt}.html`)
      ensureDirectoryExistence(outPath)
      fs.writeFileSync(outPath, html, 'utf-8')
      count++
    }
  }
  console.log(`Successfully generated: html (${locale}) (${count} files)`)
}
