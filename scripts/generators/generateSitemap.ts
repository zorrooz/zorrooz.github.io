/**
 * 由 zh-CN categories 产物生成 sitemap.xml + robots.txt：
 * 读取 content/categories.json（全部文章 URL + 日期），并入四个静态路由，
 * 写入 public/（Vite 原样复制）。
 */

import fs from 'fs'
import path from 'path'

import { contentDir } from '../dataConfig.ts'
import { LOCALE_SEGMENTS, SITE } from '../../src/config.ts'
import { isDirectRun, logWriteSuccess, runCliScript } from './core/index.ts'

const SITE_URL = SITE.url
const publicDir = path.join(import.meta.dirname, '../../public')
const LOCALE_PREFIXES = LOCALE_SEGMENTS.map((segment) => `/${segment}`)

const STATIC_ROUTES = ['/', '/about', '/resource', '/category']

interface ArticleUrlEntry {
  loc: string
  lastmod: string
}

function collectArticleUrls(): ArticleUrlEntry[] {
  const p = path.join(contentDir, 'categories.json')
  if (!fs.existsSync(p)) {
    console.warn('Warn: categories.json not found, skipping sitemap articles.')
    return []
  }
  const data = JSON.parse(fs.readFileSync(p, 'utf-8'))
  const urls: ArticleUrlEntry[] = []
  for (const group of Array.isArray(data) ? data : []) {
    for (const item of Array.isArray(group?.items) ? group.items : []) {
      for (const cat of Array.isArray(item?.categories) ? item.categories : []) {
        for (const art of Array.isArray(cat?.articles) ? cat.articles : []) {
          if (typeof art?.articleUrl === 'string' && art.articleUrl) {
            urls.push({
              loc: art.articleUrl,
              lastmod: typeof art?.date === 'string' ? art.date : '',
            })
          }
        }
      }
    }
  }
  return urls
}

function generateSitemap() {
  const articleUrls = collectArticleUrls()
  const entries = [
    ...LOCALE_PREFIXES.flatMap((prefix) =>
      STATIC_ROUTES.map((loc) => ({ loc: `${prefix}${loc}`, lastmod: '' })),
    ),
    ...LOCALE_PREFIXES.flatMap((prefix) =>
      articleUrls.map((u) => ({ loc: `${prefix}${u.loc}`, lastmod: u.lastmod })),
    ),
  ]

  const urlset = entries
    .map(
      (e) =>
        `  <url>\n    <loc>${SITE_URL}${e.loc}</loc>` +
        (e.lastmod ? `\n    <lastmod>${e.lastmod}</lastmod>` : '') +
        `\n  </url>`,
    )
    .join('\n')

  const sitemap = `<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
${urlset}
</urlset>
`

  const robots = `User-agent: *
Allow: /

Sitemap: ${SITE_URL}/sitemap.xml
`

  const sitemapPath = path.join(publicDir, 'sitemap.xml')
  const robotsPath = path.join(publicDir, 'robots.txt')
  fs.writeFileSync(sitemapPath, sitemap, 'utf-8')
  fs.writeFileSync(robotsPath, robots, 'utf-8')
  logWriteSuccess(sitemapPath, `${entries.length} urls`)
  logWriteSuccess(robotsPath)
}

if (isDirectRun(import.meta)) {
  runCliScript('sitemap', () => generateSitemap())
}

export { generateSitemap }
