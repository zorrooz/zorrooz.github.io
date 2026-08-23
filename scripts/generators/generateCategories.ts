import path from 'path'
import zhCN from '../../src/i18n/locales/zh-CN.ts'
import enUS from '../../src/i18n/locales/en-US.ts'
import { ARTICLE_ROUTE_PREFIX } from '../../src/config.ts'
import { contentDir, localeSuffix, srcDirFor } from '../dataConfig.ts'
import { normalizeProjectTopicEntry } from '../lib/yamlEntries.ts'
import { readJson, readYaml, safeArray, writeJsonFile } from '../lib/fs.ts'
import { isDirectRun, logWriteSuccess, runCliScript } from '../lib/cli.ts'

interface FormattedArticle {
  title: string
  articleUrl: string
  wordCount: number
  date: string
  tags: string[]
}

interface DetailedSubCategory {
  key: string
  title: string
  articles: Array<Record<string, unknown>>
  stats: { postsCount: number; totalWords: number; latestDate: string }
}

function getFilePaths(locale = 'zh-CN') {
  const suffix = localeSuffix(locale)
  return {
    notesJsonPath: path.join(contentDir, `notes${suffix}.json`),
    projectsJsonPath: path.join(contentDir, `projects${suffix}.json`),
    topicsJsonPath: path.join(contentDir, `topics${suffix}.json`),
    categoriesYamlPath: path.join(srcDirFor(locale), `categories${suffix}.yaml`),
    categoriesJsonPath: path.join(contentDir, `categories${suffix}.json`),
  }
}

function getSubCategoryKeyFromUrl(url: unknown): string {
  if (typeof url !== 'string') return ''
  const parts = url.replace(/^\/+/, '').split('/')
  const articleIndex = parts[0] === ARTICLE_ROUTE_PREFIX.slice(1) ? 1 : 0
  return parts[articleIndex + 2] || ''
}

function formatArticle(article: Record<string, unknown>, type: string): FormattedArticle {
  const relativePath = typeof article?.relativePath === 'string' ? article.relativePath : ''
  return {
    title:
      typeof article.title === 'string' ? article.title : relativePath.split('/').pop() || '未命名',
    articleUrl: `${ARTICLE_ROUTE_PREFIX}/${type}/${relativePath}`,
    wordCount: Number(article?.wordCount) || 0,
    date: typeof article?.date === 'string' ? article.date : '',
    tags: safeArray(article?.tags).filter((t) => typeof t === 'string'),
  }
}

/** 子分类文章不含 date 字段（与产物契约一致：articles 只保留 title/url/wordCount/tags） */
function withoutDate(art: FormattedArticle): Record<string, unknown> {
  const rest: Record<string, unknown> = { ...art }
  delete rest.date
  return rest
}

/** 最新一篇文章的日期；空列表返回 ''（此处原地排序） */
function latestDateOf(articles: FormattedArticle[]): string {
  return articles.length
    ? articles.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())[0].date
    : ''
}

interface NormalizedNoteConfig {
  name: string
  title: string
  desc: string
  date: string
  categories: Record<string, unknown>
}

function normalizeNoteConfig(rawDef: Record<string, unknown>): NormalizedNoteConfig {
  const rawCats = rawDef?.categories
  return {
    name: typeof rawDef?.name === 'string' ? rawDef.name : '',
    title: typeof rawDef?.title === 'string' ? rawDef.title : '',
    desc: typeof rawDef?.desc === 'string' ? rawDef.desc : '',
    date: typeof rawDef?.date === 'string' ? rawDef.date : '',
    categories:
      typeof rawCats === 'object' && rawCats !== null ? (rawCats as Record<string, unknown>) : {},
  }
}

interface DetailedCategory {
  name: string
  title: string
  desc: string
  tags: string[]
  stats: DetailedSubCategory['stats']
  categories: DetailedSubCategory[]
  root: string
}

function buildDetailedNoteCategories(
  noteConfigs: Record<string, unknown>[],
  noteArticles: Record<string, unknown>[],
  locale = 'zh-CN',
): DetailedCategory[] {
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
      const detailedSubCats: DetailedSubCategory[] = predefinedSubCats.map(([key, title]) => {
        const catArticles = currentArticles.filter(
          (art) => getSubCategoryKeyFromUrl(art.articleUrl) === key,
        )
        return {
          key,
          title: title as string,
          articles: catArticles.map(withoutDate),
          stats: {
            postsCount: catArticles.length,
            totalWords: catArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: latestDateOf(catArticles),
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
          articles: uncategorizedArticles.map(withoutDate),
          stats: {
            postsCount: uncategorizedArticles.length,
            totalWords: uncategorizedArticles.reduce((sum, art) => sum + art.wordCount, 0),
            latestDate: latestDateOf(uncategorizedArticles),
          },
        })
      }

      const globalStats = {
        postsCount: currentArticles.length,
        totalWords: currentArticles.reduce((sum, art) => sum + art.wordCount, 0),
        latestDate: config.date || latestDateOf(currentArticles),
      }

      const globalTags = Array.from(
        new Set(
          currentArticles.flatMap((art) =>
            safeArray(art?.tags || []).filter((t) => typeof t === 'string'),
          ),
        ),
      )

      const rootUrl =
        (detailedSubCats.find((cat) => cat.stats.postsCount > 0)?.articles[0]?.articleUrl as
          | string
          | undefined) || ''

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

interface DetailedProjectTopicCategory extends DetailedCategory {
  github?: string
  doi?: string
}

function buildDetailedProjectTopicCategories(
  ptConfigs: Record<string, unknown>[],
  ptArticles: Record<string, unknown>[],
  type: string,
): DetailedProjectTopicCategory[] {
  return ptConfigs
    .map((rawConfig) => {
      const config = normalizeProjectTopicEntry(rawConfig)
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

      /** projects/topics 为纯元数据（无文章），不生成子分类 */
      const detailedSubCats: DetailedSubCategory[] = []

      const globalStats = {
        postsCount: currentArticles.length,
        totalWords: currentArticles.reduce((sum, art) => sum + art.wordCount, 0),
        latestDate: config.date || latestDateOf(currentArticles),
      }

      let rootUrl = ''
      if (config.url && config.url.startsWith(`${ARTICLE_ROUTE_PREFIX}/`)) {
        rootUrl = config.url
      } else {
        rootUrl =
          (detailedSubCats.find((cat) => cat.stats.postsCount > 0)?.articles[0]?.articleUrl as
            | string
            | undefined) || ''
      }

      const extraFields =
        type === 'project'
          ? {
              github: config.github,
              status: config.status,
              language: config.language,
              stars: config.stars,
              license: config.license,
              version: config.version,
              url: config.url,
            }
          : {
              doi: config.doi,
              status: config.status,
              journal: config.journal,
              year: config.year,
              authors: config.authors,
              url: config.url,
            }

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

interface CategoryTitles {
  notes: string
  projects: string
  topics: string
  uncategorized: string
}

function getCategoryTitles(locale = 'zh-CN'): CategoryTitles {
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

    const noteArticles = readJson(paths.notesJsonPath)
    const projectArticles = readJson(paths.projectsJsonPath)
    const topicArticles = readJson(paths.topicsJsonPath)
    const yamlConfig = (readYaml(paths.categoriesYamlPath) as Record<string, unknown>) || {}
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
      noteConfigs as Record<string, unknown>[],
      (Array.isArray(noteArticles) ? noteArticles : []) as Record<string, unknown>[],
      locale,
    )
    const detailedProjects = buildDetailedProjectTopicCategories(
      projectConfigs as Record<string, unknown>[],
      (Array.isArray(projectArticles) ? projectArticles : []) as Record<string, unknown>[],
      'project',
    )
    const detailedTopics = buildDetailedProjectTopicCategories(
      topicConfigs as Record<string, unknown>[],
      (Array.isArray(topicArticles) ? topicArticles : []) as Record<string, unknown>[],
      'topic',
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
    writeJsonFile(targetPath, finalStructure)
    logWriteSuccess(targetPath)
  } catch (error) {
    console.error(`Failed to generate ${paths.categoriesJsonPath}:`, error)
    process.exitCode = 1
  }
}

if (isDirectRun(import.meta)) {
  runCliScript('categories.json', () => {
    generateCategoriesJson('zh-CN')
    generateCategoriesJson('en-US')
  })
}

export { generateCategoriesJson }
