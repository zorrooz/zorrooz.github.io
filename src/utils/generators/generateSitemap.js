/**
 * Generates sitemap.xml + robots.txt from the zh-CN categories output.
 * Reads src/content/categories.json (all article URLs + dates), adds the
 * four static routes, and writes to public/ (copied as-is by Vite).
 */

import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const SITE_URL = 'https://zorrooz.github.io'
const contentDir = path.join(__dirname, '../../content')
const publicDir = path.join(__dirname, '../../../public')

const STATIC_ROUTES = ['/', '/about', '/resource', '/category']

function collectArticleUrls() {
  const p = path.join(contentDir, 'categories.json')
  if (!fs.existsSync(p)) {
    console.warn('Warn: categories.json not found, skipping sitemap articles.')
    return []
  }
  const data = JSON.parse(fs.readFileSync(p, 'utf-8'))
  const urls = []
  for (const group of Array.isArray(data) ? data : []) {
    for (const item of Array.isArray(group?.items) ? group.items : []) {
      for (const cat of Array.isArray(item?.categories) ? item.categories : []) {
        for (const art of Array.isArray(cat?.articles) ? cat.articles : []) {
          if (typeof art?.articleUrl === 'string' && art.articleUrl) {
            urls.push({ loc: art.articleUrl, lastmod: typeof art?.date === 'string' ? art.date : '' })
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
    ...STATIC_ROUTES.map((loc) => ({ loc, lastmod: '' })),
    ...articleUrls,
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
  console.log(`Successfully generated: ${sitemapPath} (${entries.length} urls)`)
  console.log(`Successfully generated: ${robotsPath}`)
}

function main() {
  console.log('Starting sitemap generation script...')
  generateSitemap()
  console.log('sitemap generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateSitemap }
