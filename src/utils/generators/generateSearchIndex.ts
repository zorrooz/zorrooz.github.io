/**
 * Generates the client-side full-text search index.
 * Reads categories.json (article metadata) + pre-rendered article HTML
 * (src/content/html/**), strips tags to plain text, and writes
 * src/content/search-index.json (per locale, -en suffix pattern).
 * Runs after generateHtml, mirrors loadJsonContent lookup semantics.
 */

import fs from 'fs'
import path from 'path'

import { contentDir } from '../dataConfig.ts'

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
  const notesPath = path.join(contentDir, locale === 'en-US' ? 'notes-en.json' : 'notes.json')
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
  const categoriesPath = path.join(
    contentDir,
    locale === 'en-US' ? 'categories-en.json' : 'categories.json',
  )
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
  const fileName = locale === 'en-US' ? 'search-index-en.json' : 'search-index.json'
  const outPath = path.join(contentDir, fileName)
  fs.writeFileSync(outPath, JSON.stringify(docs), 'utf-8')
  console.log(`Successfully generated: ${outPath} (${docs.length} docs)`)
}

export { generateSearchIndex }
