/**
 * 生成客户端全文检索索引：
 * 读取 categories.json（文章元数据）+ 预渲染文章 HTML（content/html/**），
 * 去标签转纯文本后写入 content/search-index.json（按 locale，-en 后缀）。
 * 在 generateHtml 之后运行，镜像 loadJsonContent 的查找语义。
 */

import fs from 'fs'
import path from 'path'

import { contentDir, localeSuffix } from '../dataConfig.ts'
import { isDirectRun, logWriteSuccess, runCliScript, writeJsonFile } from './core/index.ts'

function stripHtmlToText(html: string): string {
  return html
    .replace(/<style[\s\S]*?<\/style>/gi, ' ')
    .replace(/<script[\s\S]*?<\/script>/gi, ' ')
    .replace(/<[^>]+>/g, ' ')
    .replace(/&amp;/gi, '&')
    .replace(/&lt;/gi, '<')
    .replace(/&gt;/gi, '>')
    .replace(/&quot;/gi, '"')
    .replace(/&#39;|&apos;/gi, "'")
    .replace(/&nbsp;/gi, ' ')
    .replace(/\s+/g, ' ')
    .trim()
}

function loadDescriptions(locale: 'zh-CN' | 'en-US'): Map<string, string> {
  const notesPath = path.join(contentDir, `notes${localeSuffix(locale)}.json`)
  const map = new Map<string, string>()
  if (!fs.existsSync(notesPath)) return map
  const notes = JSON.parse(fs.readFileSync(notesPath, 'utf-8'))
  for (const note of Array.isArray(notes) ? notes : []) {
    if (typeof note?.relativePath !== 'string') continue
    map.set(note.relativePath, typeof note.description === 'string' ? note.description : '')
  }
  return map
}

interface SearchDoc {
  id: string
  title: string
  tags: string[]
  path: string
  description: string
  content: string
}

function collectArticles(locale: 'zh-CN' | 'en-US'): SearchDoc[] {
  const categoriesPath = path.join(contentDir, `categories${localeSuffix(locale)}.json`)
  if (!fs.existsSync(categoriesPath)) {
    console.warn(`Warn: ${categoriesPath} not found, skipping search index.`)
    return []
  }

  const descriptions = loadDescriptions(locale)
  const data = JSON.parse(fs.readFileSync(categoriesPath, 'utf-8'))
  const docs: SearchDoc[] = []

  for (const section of Array.isArray(data) ? data : []) {
    for (const item of Array.isArray(section?.items) ? section.items : []) {
      for (const cat of Array.isArray(item?.categories) ? item.categories : []) {
        for (const art of Array.isArray(cat?.articles) ? cat.articles : []) {
          if (typeof art?.articleUrl !== 'string' || !art.articleUrl) continue
          const rel = art.articleUrl.replace(/^\/article\//, '')
          const htmlPath = path.join(contentDir, 'html', `${rel}.html`)

          let content = ''
          if (fs.existsSync(htmlPath)) {
            content = stripHtmlToText(fs.readFileSync(htmlPath, 'utf-8'))
          }

          const tags: string[] = []
          if (Array.isArray(art.tags)) {
            for (const tag of art.tags) {
              if (typeof tag === 'string') tags.push(tag)
            }
          }

          const key = rel.replace(/^notes\//, '')
          let description = descriptions.get(key) || ''
          if (!description && key.endsWith('-en')) {
            description = descriptions.get(key.slice(0, -3)) || ''
          }

          docs.push({
            id: rel,
            title: typeof art.title === 'string' ? art.title : '',
            tags,
            path: rel,
            description,
            content,
          })
        }
      }
    }
  }
  return docs
}

function generateSearchIndex(locale: 'zh-CN' | 'en-US') {
  const docs = collectArticles(locale)
  const fileName = `search-index${localeSuffix(locale)}.json`
  const outPath = path.join(contentDir, fileName)
  // space=0 有意为之：search-index 由 SearchModal 懒加载，紧凑 JSON 减小首屏 payload
  writeJsonFile(outPath, docs, 0)
  logWriteSuccess(outPath, `${docs.length} docs`)
}

if (isDirectRun(import.meta)) {
  runCliScript('search-index', () => {
    generateSearchIndex('zh-CN')
    generateSearchIndex('en-US')
  })
}

export { generateSearchIndex }
