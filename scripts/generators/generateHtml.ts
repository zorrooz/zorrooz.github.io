/**
 * Generates static HTML for every article at build time.
 * Same file enumeration rules as generateNotes.ts:
 *   zh-CN: all .md files except `-en.md`
 *   en-US: all .md files under cache/en (machine translation layer)
 * Output: content/html/<relative-path-without-ext>.html (dataDir/content)
 * (e.g. notes/Omics/genomics/bwa/bwa.html, .../bwa-en.html)
 *
 * Only `notes` are rendered: projects/topics are metadata-only
 * (their md sources are removed from content-src).
 */

import fs from 'fs'
import path from 'path'

import { contentDir, srcDirFor } from '../dataConfig.ts'
import { renderMarkdown } from '../markdownProcessor.ts'
import {
  ensureDirectoryExistence,
  isDirectRun,
  mdFileFilter,
  runCliScript,
  toPosixRelativeNoExt,
  walk,
} from './core/index.ts'

const htmlOutputDir = path.join(contentDir, 'html')

const sections = ['notes']
const deprecatedSections = ['projects', 'topics']

/** 清理不再生成文章的旧目录残留（projects/topics 已纯元信息化） */
function cleanDeprecatedHtml() {
  for (const section of deprecatedSections) {
    const dir = path.join(htmlOutputDir, section)
    if (fs.existsSync(dir)) {
      fs.rmSync(dir, { recursive: true, force: true })
      console.log(`Cleaned deprecated html output: ${dir}`)
    }
  }
}

export async function generateHtml(locale = 'zh-CN') {
  let count = 0
  cleanDeprecatedHtml()
  for (const section of sections) {
    const srcDir = path.join(srcDirFor(locale), section)
    const mdFiles = walk(srcDir, mdFileFilter(locale))

    for (const mdPath of mdFiles) {
      const raw = fs.readFileSync(mdPath, 'utf-8')
      const html = await renderMarkdown(raw)

      const relNoExt = toPosixRelativeNoExt(mdPath, srcDirFor(locale))

      const outPath = path.join(htmlOutputDir, `${relNoExt}.html`)
      ensureDirectoryExistence(outPath)
      fs.writeFileSync(outPath, html, 'utf-8')
      count++
    }
  }
  console.log(`Successfully generated: html (${locale}) (${count} files)`)
}

if (isDirectRun(import.meta)) {
  runCliScript('html', async () => {
    await generateHtml('zh-CN')
    await generateHtml('en-US')
  })
}
