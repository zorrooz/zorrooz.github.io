/**
 * 构建期把每篇文章渲染为静态 HTML。
 * 文件枚举规则与 generateNotes.ts 一致：
 *   zh-CN：除 `-en.md` 外的全部 .md；en-US：cache/en（机器翻译层）的全部 .md
 * 输出：content/html/<相对路径去扩展名>.html（如 notes/Omics/genomics/bwa/bwa.html、bwa-en.html）
 * 仅渲染 notes：projects/topics 为纯元数据（无 md 文章）。
 */

import fs from 'fs'
import path from 'path'

import { contentDir, srcDirFor } from '../dataConfig.ts'
import { renderMarkdown } from './markdownProcessor.ts'
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
