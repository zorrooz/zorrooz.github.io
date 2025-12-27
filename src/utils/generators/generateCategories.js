import fs from 'fs'
import path from 'path'
import { fileURLToPath } from 'url'
import yaml from 'js-yaml'

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

function ensureDirectoryExistence(filePath) {
  const dir = path.dirname(filePath)
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true })
}

function safeArray(x) {
  return Array.isArray(x) ? x : []
}

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

function readYaml(p, fallback = {}) {
  try {
    if (!fs.existsSync(p)) return fallback
    const raw = fs.readFileSync(p, 'utf-8')
    const obj = yaml.load(raw)
    return obj || fallback
  } catch {
    return fallback
  }
}


function getSubCategoryKeyFromUrl(url) {
  if (typeof url !== 'string') return ''
  const parts = url.replace(/^\/+/, '').split('/')
  const articleIndex = parts[0] === 'article' ? 1 : 0
  return parts[articleIndex + 2] || ''
}

function formatArticle(article, type) {
  const relativePath = article?.relativePath || ''
  return {
    title:
      typeof article.title === 'string' ? article.title : relativePath.split('/').pop() || '未命名',
    articleUrl: `/article/${type}/${relativePath}`,
    wordCount: Number(article?.wordCount) || 0,
    date: article?.date || '',
    tags: safeArray(article?.tags).filter((t) => typeof t === 'string'),
  }
}

function normalizeNoteConfig(rawDef) {
  return {
    name: typeof rawDef?.name === 'string' ? rawDef.name : '',
    title: typeof rawDef?.title === 'string' ? rawDef.title : '',
    desc: typeof rawDef?.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef?.date === 'string' ? rawDef.date : '',
    categories:
      rawDef?.categories && typeof rawDef.categories === 'object' ? rawDef.categories : {},
  }
}

function normalizeProjectTopicConfig(rawDef) {
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
      rawDef?.categories && typeof rawDef.categories === 'object' ? rawDef.categories : {},
  }
}

function buildDetailedNoteCategories(noteConfigs, noteArticles) {
  return noteConfigs
    .map((rawConfig) => {
      const config = normalizeNoteConfig(rawConfig)
      if (!config.name) return null

      const currentArticles = noteArticles
        .filter(
          (art) =>
            typeof art?.relativePath === 'string' && art.relativePath.split('/')[0] === config.name,
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
          articles: catArticles.map(({ wordCount, date, ...rest }) => rest),
          stats: {
            postsCount: catArticles.length,
            totalWords: catArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: catArticles.length
              ? catArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date
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
          title: '未分类',
          articles: uncategorizedArticles.map(({ wordCount, date, ...rest }) => rest),
          stats: {
            postsCount: uncategorizedArticles.length,
            totalWords: uncategorizedArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: uncategorizedArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0]
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
            ? currentArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date
            : ''),
      }

      const globalTags = Array.from(
        new Set(
          currentArticles.flatMap((art) =>
            safeArray(art?.tags || []).filter((t) => typeof t === 'string'),
          ),
        ),
      )

      const rootUrl =
        detailedSubCats.find((cat) => cat.stats.postsCount > 0)?.articles[0]?.articleUrl || ''

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
    .filter(Boolean)
}

function buildDetailedProjectTopicCategories(ptConfigs, ptArticles, type) {
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
          articles: catArticles.map(({ wordCount, date, ...rest }) => rest),
          stats: {
            postsCount: catArticles.length,
            totalWords: catArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: catArticles.length
              ? catArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date
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
          title: '未分类',
          articles: uncategorizedArticles.map(({ wordCount, date, ...rest }) => rest),
          stats: {
            postsCount: uncategorizedArticles.length,
            totalWords: uncategorizedArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: uncategorizedArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0]
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
            ? currentArticles.sort((a, b) => new Date(b.date) - new Date(a.date))[0].date
            : ''),
      }

      let rootUrl = ''
      if (config.url && config.url.startsWith('/article/')) {
        rootUrl = config.url
      } else {
        rootUrl =
          detailedSubCats.find((cat) => cat.stats.postsCount > 0)?.articles[0]?.articleUrl || ''
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
    .filter(Boolean)
}

function getCategoryTitles(locale = 'zh-CN') {
  if (locale === 'en-US') {
    return {
      notes: 'Notes',
      projects: 'Projects',
      topics: 'Topics',
    }
  }
  return {
    notes: '笔记',
    projects: '项目',
    topics: '课题',
  }
}

function generateCategoriesJson(locale = 'zh-CN') {
  try {
    const paths = getFilePaths(locale)
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

    const detailedNotes = buildDetailedNoteCategories(noteConfigs, noteArticles)
    const detailedProjects = buildDetailedProjectTopicCategories(
      projectConfigs,
      projectArticles,
      'project',
    )
    const detailedTopics = buildDetailedProjectTopicCategories(topicConfigs, topicArticles, 'topic')

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
    console.error(`Failed to generate ${categoriesJsonPath}:`, error)
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
