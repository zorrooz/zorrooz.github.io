// @ts-check
/**
 * @typedef {Object} FormattedArticle
 * @property {string} title
 * @property {string} articleUrl
 * @property {number} wordCount
 * @property {string} date
 * @property {string[]} tags
 */

/**
 * @typedef {Object} DetailedSubCategory
 * @property {string} key
 * @property {unknown} title
 * @property {Array<Record<string, unknown>>} articles
 * @property {{ postsCount: number, totalWords: number, latestDate: string }} stats
 */

import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
import yaml from 'js-yaml'
import zhCN from '../../stores/locales/zh-CN.js'
import enUS from '../../stores/locales/en-US.js'
import { ensureDirectoryExistence } from './core/index.js'

const __filename = fileURLToPath(import.meta.url)
const __dirname = path.dirname(__filename)

const contentDir = path.join(__dirname, '../../content')
const contentSrcDir = path.join(__dirname, '../../content-src')

function getFilePaths(locale = 'zh-CN') {
  const suffix = locale === 'zh-CN' ? '' : '-en'
  return {
    notesJsonPath: path.join(contentDir, `notes${suffix}.json`),
    projectsJsonPath: path.join(contentDir, `projects${suffix}.json`),
    topicsJsonPath: path.join(contentDir, `topics${suffix}.json`),
    categoriesYamlPath: path.join(contentSrcDir, `categories${suffix}.yaml`),
    categoriesJsonPath: path.join(contentDir, `categories${suffix}.json`),
  }
}



/** @param {unknown} x @returns {unknown[]} */
function safeArray(x) {
  return Array.isArray(x) ? x : []
}

/**
 * @param {string} p
 * @param {unknown[]} [fallback]
 * @returns {unknown[]}
 */
function readJsonArray(p, fallback = []) {
  try {
    if (!fs.existsSync(p)) return fallback
    const raw = fs.readFileSync(p, 'utf-8').trim()
    const arr = JSON.parse(raw)
    return Array.isArray(arr) ? arr : fallback
  } catch {
    return fallback
  }
}

/**
 * @param {string} p
 * @param {Record<string, unknown>} [fallback]
 * @returns {Record<string, unknown>}
 */
function readYaml(p, fallback = {}) {
  try {
    if (!fs.existsSync(p)) return fallback
    const raw = fs.readFileSync(p, 'utf-8')
    const obj = yaml.load(raw)
    return /** @type {Record<string, unknown>} */ (obj) || fallback
  } catch {
    return fallback
  }
}

/** @param {unknown} url @returns {string} */
function getSubCategoryKeyFromUrl(url) {
  if (typeof url !== 'string') return ''
  const parts = url.replace(/^\/+/, '').split('/')
  const articleIndex = parts[0] === 'article' ? 1 : 0
  return parts[articleIndex + 2] || ''
}

/**
 * @param {Record<string, unknown>} article
 * @param {string} type
 * @returns {FormattedArticle}
 */
function formatArticle(article, type) {
  const relativePath = typeof article?.relativePath === 'string' ? article.relativePath : ''
  return {
    title:
      typeof article.title === 'string' ? article.title : relativePath.split('/').pop() || '未命名',
    articleUrl: `/article/${type}/${relativePath}`,
    wordCount: Number(article?.wordCount) || 0,
    date: typeof article?.date === 'string' ? article.date : '',
    tags: safeArray(article?.tags).filter((t) => typeof t === 'string'),
  }
}

/**
 * @param {Record<string, unknown>} rawDef
 * @returns {{ name: string, title: string, desc: string, date: string, categories: Record<string, unknown> }}
 */
function normalizeNoteConfig(rawDef) {
  const rawCats = rawDef?.categories
  return {
    name: typeof rawDef?.name === 'string' ? rawDef.name : '',
    title: typeof rawDef?.title === 'string' ? rawDef.title : '',
    desc: typeof rawDef?.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef?.date === 'string' ? rawDef.date : '',
    categories:
      typeof rawCats === 'object' && rawCats !== null
        ? /** @type {Record<string, unknown>} */ (rawCats)
        : {},
  }
}

/**
 * @param {Record<string, unknown>} rawDef
 * @returns {{ name: string, title: string, desc: string, date: string, tags: string[], github: string, doi: string, url: string, categories: Record<string, unknown> }}
 */
function normalizeProjectTopicConfig(rawDef) {
  const rawCats = rawDef?.categories
  return {
    name: typeof rawDef?.name === 'string' ? rawDef.name : '',
    title: typeof rawDef?.title === 'string' ? rawDef.title : '',
    desc: typeof rawDef?.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef?.date === 'string' ? rawDef.date : '',
    tags: safeArray(rawDef?.tags).filter((t) => typeof t === 'string'),
    github: typeof rawDef?.github === 'string' ? rawDef.github : '',
    doi: typeof rawDef?.doi === 'string' ? rawDef.doi : '',
    url: typeof rawDef?.url === 'string' ? rawDef.url : '',
    categories:
      typeof rawCats === 'object' && rawCats !== null
        ? /** @type {Record<string, unknown>} */ (rawCats)
        : {},
  }
}

/**
 * @param {Record<string, unknown>[]} noteConfigs
 * @param {Record<string, unknown>[]} noteArticles
 * @param {string} [locale]
 * @returns {Array<{ name: string, title: string, desc: string, tags: string[], stats: DetailedSubCategory['stats'], categories: DetailedSubCategory[], root: string }>}
 */
function buildDetailedNoteCategories(noteConfigs, noteArticles, locale = 'zh-CN') {
  return noteConfigs
    .map((rawConfig) => {
      const config = normalizeNoteConfig(rawConfig)
      if (!config.name) return null

      const currentArticles = noteArticles
        .filter(
          (art) =>
            typeof art?.relativePath === 'string' &&
            art.relativePath.split('/')[0] === config.name,
        )
        .map((art) => formatArticle(art, 'notes'))

      const predefinedSubCats = Object.entries(config.categories)
      const detailedSubCats = predefinedSubCats.map(([key, title]) => {
        const catArticles = currentArticles.filter(
          (art) => getSubCategoryKeyFromUrl(art.articleUrl) === key,
        )
        return {
          key,
          title,
          articles: catArticles.map((art) => {
            /** @type {Record<string, unknown>} */
            const rest = { ...art }
            delete rest.wordCount
            delete rest.date
            return rest
          }),
          stats: {
            postsCount: catArticles.length,
            totalWords: catArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: catArticles.length
              ? catArticles.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())[0].date
              : '',
          },
        }
      })

      const categorizedUrls = detailedSubCats.flatMap((cat) =>
        cat.articles.map((art) => art.articleUrl),
      )
      const uncategorizedArticles = currentArticles.filter(
        (art) => !categorizedUrls.includes(art.articleUrl),
      )
      if (uncategorizedArticles.length > 0) {
        detailedSubCats.push({
          key: 'uncategorized',
          title: getCategoryTitles(locale).uncategorized,
          articles: uncategorizedArticles.map((art) => {
            /** @type {Record<string, unknown>} */
            const rest = { ...art }
            delete rest.wordCount
            delete rest.date
            return rest
          }),
          stats: {
            postsCount: uncategorizedArticles.length,
            totalWords: uncategorizedArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: uncategorizedArticles.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())[0]
              .date,
          },
        })
      }

      const globalStats = {
        postsCount: currentArticles.length,
        totalWords: currentArticles.reduce((sum, art) => sum + art.wordCount, 0),
        latestDate:
          config.date ||
          (currentArticles.length
            ? currentArticles.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())[0].date
            : ''),
      }

      const globalTags = Array.from(
        new Set(
          currentArticles.flatMap((art) =>
            safeArray(art?.tags || []).filter((t) => typeof t === 'string'),
          ),
        ),
      )

      const rootUrl = /** @type {string} */ (
        detailedSubCats.find((cat) => cat.stats.postsCount > 0)?.articles[0]?.articleUrl || ''
      )

      return {
        name: config.name,
        title: config.title || config.name,
        desc: config.desc,
        tags: globalTags,
        stats: globalStats,
        categories: detailedSubCats,
        root: rootUrl,
      }
    })
    .filter((item) => item !== null)
}

/**
 * @param {Record<string, unknown>[]} ptConfigs
 * @param {Record<string, unknown>[]} ptArticles
 * @param {string} type
 * @param {string} [locale]
 * @returns {Array<{ name: string, title: string, desc: string, tags: string[], stats: DetailedSubCategory['stats'], categories: DetailedSubCategory[], root: string, github?: string, doi?: string }>}
 */
function buildDetailedProjectTopicCategories(ptConfigs, ptArticles, type, locale = 'zh-CN') {
  return ptConfigs
    .map((rawConfig) => {
      const config = normalizeProjectTopicConfig(rawConfig)
      if (!config.name) return null

      const nameKey = config.name.toLowerCase().replace(/\s+/g, '')
      const nameAlpha = config.name.toLowerCase().replace(/[^a-z0-9]+/g, '')
      const currentArticles = ptArticles
        .filter((art) => {
          if (typeof art?.relativePath !== 'string') return false
          const artName = art.relativePath.split('/')[0].toLowerCase().replace(/\s+/g, '')
          const artAlpha = art.relativePath
            .split('/')[0]
            .toLowerCase()
            .replace(/[^a-z0-9]+/g, '')
          return [nameKey, nameAlpha].includes(artName) || [nameKey, nameAlpha].includes(artAlpha)
        })
        .map((art) => formatArticle(art, `${type}s`))

      const predefinedSubCats = Object.entries(config.categories)
      const detailedSubCats = predefinedSubCats.map(([key, title]) => {
        const catArticles = currentArticles.filter(
          (art) => getSubCategoryKeyFromUrl(art.articleUrl) === key,
        )
        return {
          key,
          title,
          articles: catArticles.map((art) => {
            /** @type {Record<string, unknown>} */
            const rest = { ...art }
            delete rest.wordCount
            delete rest.date
            return rest
          }),
          stats: {
            postsCount: catArticles.length,
            totalWords: catArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: catArticles.length
              ? catArticles.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())[0].date
              : '',
          },
        }
      })

      const categorizedUrls = detailedSubCats.flatMap((cat) =>
        cat.articles.map((art) => art.articleUrl),
      )
      const uncategorizedArticles = currentArticles.filter(
        (art) => !categorizedUrls.includes(art.articleUrl),
      )
      if (uncategorizedArticles.length > 0) {
        detailedSubCats.push({
          key: 'uncategorized',
          title: getCategoryTitles(locale).uncategorized,
          articles: uncategorizedArticles.map((art) => {
            /** @type {Record<string, unknown>} */
            const rest = { ...art }
            delete rest.wordCount
            delete rest.date
            return rest
          }),
          stats: {
            postsCount: uncategorizedArticles.length,
            totalWords: uncategorizedArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: uncategorizedArticles.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())[0]
              .date,
          },
        })
      }

      const globalStats = {
        postsCount: currentArticles.length,
        totalWords: currentArticles.reduce((sum, art) => sum + art.wordCount, 0),
        latestDate:
          config.date ||
          (currentArticles.length
            ? currentArticles.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())[0].date
            : ''),
      }

      let rootUrl = ''
      if (config.url && config.url.startsWith('/article/')) {
        rootUrl = config.url
      } else {
        rootUrl = /** @type {string} */ (
          detailedSubCats.find((cat) => cat.stats.postsCount > 0)?.articles[0]?.articleUrl || ''
        )
      }

      const extraFields = type === 'project' ? { github: config.github } : { doi: config.doi }

      return {
        name: config.name,
        title: config.title || config.name,
        desc: config.desc,
        tags: config.tags,
        stats: globalStats,
        categories: detailedSubCats,
        root: rootUrl,
        ...extraFields,
      }
    })
    .filter((item) => item !== null)
}

/**
 * @param {string} [locale]
 * @returns {{ notes: string, projects: string, topics: string, uncategorized: string }}
 */
function getCategoryTitles(locale = 'zh-CN') {
  const l = locale === 'en-US' ? enUS : zhCN
  return {
    notes: (l && l.notes) || (locale === 'en-US' ? 'Notes' : '笔记'),
    projects: (l && l.projects) || (locale === 'en-US' ? 'Projects' : '项目'),
    topics: (l && l.topics) || (locale === 'en-US' ? 'Topics' : '课题'),
    uncategorized: (l && l.uncategorized) || (locale === 'en-US' ? 'Uncategorized' : '未分类'),
  }
}

function generateCategoriesJson(locale = 'zh-CN') {
  const paths = getFilePaths(locale)
  try {
    const titles = getCategoryTitles(locale)

    const noteArticles = readJsonArray(paths.notesJsonPath)
    const projectArticles = readJsonArray(paths.projectsJsonPath)
    const topicArticles = readJsonArray(paths.topicsJsonPath)
    const yamlConfig = readYaml(paths.categoriesYamlPath)
    const {
      notes: noteConfigs,
      projects: projectConfigs,
      topics: topicConfigs,
    } = {
      notes: safeArray(yamlConfig?.notes),
      projects: safeArray(yamlConfig?.projects),
      topics: safeArray(yamlConfig?.topics),
    }

    const detailedNotes = buildDetailedNoteCategories(
      /** @type {Record<string, unknown>[]} */ (noteConfigs),
      /** @type {Record<string, unknown>[]} */ (noteArticles),
      locale,
    )
    const detailedProjects = buildDetailedProjectTopicCategories(
      /** @type {Record<string, unknown>[]} */ (projectConfigs),
      /** @type {Record<string, unknown>[]} */ (projectArticles),
      'project',
      locale,
    )
    const detailedTopics = buildDetailedProjectTopicCategories(
      /** @type {Record<string, unknown>[]} */ (topicConfigs),
      /** @type {Record<string, unknown>[]} */ (topicArticles),
      'topic',
      locale,
    )

    const finalStructure = [
      {
        title: titles.notes,
        items: detailedNotes.sort((a, b) =>
          a.name.localeCompare(b.name, locale === 'zh-CN' ? 'zh-Hans-CN' : 'en'),
        ),
      },
      {
        title: titles.projects,
        items: detailedProjects.sort((a, b) =>
          a.name.localeCompare(b.name, locale === 'zh-CN' ? 'zh-Hans-CN' : 'en'),
        ),
      },
      {
        title: titles.topics,
        items: detailedTopics.sort((a, b) =>
          a.name.localeCompare(b.name, locale === 'zh-CN' ? 'zh-Hans-CN' : 'en'),
        ),
      },
    ]

    const targetPath = paths.categoriesJsonPath
    ensureDirectoryExistence(targetPath)
    fs.writeFileSync(targetPath, JSON.stringify(finalStructure, null, 2), 'utf-8')
    console.log(`Successfully generated: ${targetPath}`)
  } catch (error) {
    console.error(`Failed to generate ${paths.categoriesJsonPath}:`, error)
    process.exitCode = 1
  }
}

function main() {
  console.log('Starting categories.json generation script...')
  generateCategoriesJson('zh-CN')
  generateCategoriesJson('en-US')
  console.log('categories.json generation complete.')
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  main()
}

export { generateCategoriesJson }
