import type { AboutData, CategoryData, Note, Post, ResourceNode, Tag } from '@/types'
import { currentLocale } from '@/i18n/locale'
import { EN_SUFFIX } from '@/config'

const getLocalizedFileName = (baseName: string, extension = '.json'): string => {
  const locale = currentLocale()
  const isEnglish = locale === 'en-US'

  if (baseName.endsWith(EN_SUFFIX)) {
    return `${baseName}${extension}`
  }

  return isEnglish ? `${baseName}${EN_SUFFIX}${extension}` : `${baseName}${extension}`
}

/**
 * 急切加载核心 JSON（about/categories/notes/posts/projects/resources/tags/topics 含 -en 变体）；
 * search-index* 除外（SearchModal 内懒加载）。@data 别名指向数据分支（GBLOG_DATA_DIR）
 */
const jsonModules = import.meta.glob('@data/content/[abcnprt]*.json', { eager: true })
const htmlModules = import.meta.glob('@data/content/html/**/*.html', {
  query: '?raw',
  import: 'default',
  eager: true,
})
/** 惰性加载 markdown 源（复制文章用，按需拆包）：content-src 为 zh 源，cache/en 为机器层 */
const markdownModules = import.meta.glob('@data/{content-src,cache/en}/**/*.md', {
  query: '?raw',
  import: 'default',
})

interface JsonModule {
  default?: unknown
}

const loadJsonContent = <T>(fileName: string, fallback: T): T => {
  const localizedFileName = getLocalizedFileName(fileName, '.json')

  const matchedKey = Object.keys(jsonModules).find((key) => key.includes(localizedFileName))

  if (matchedKey) {
    return (jsonModules[matchedKey] as JsonModule).default as T
  }

  console.error(`Failed to load JSON content: ${localizedFileName}`)

  if (currentLocale() === 'en-US') {
    const fallbackKey = Object.keys(jsonModules).find((key) => key.includes(`${fileName}.json`))

    if (fallbackKey) {
      return (jsonModules[fallbackKey] as JsonModule).default as T
    }
  }

  return fallback
}

/**
 * 加载预渲染的文章 HTML。
 * `filePath` 为 md 路径（如 `notes/Omics/genomics/bwa/bwa.md`）；优先匹配本地化文件，
 * 回退到中文（en-US: bwa-en.html → bwa.html；zh-CN: bwa.html 或路径自带 -en 时取 -en）。
 */
export const loadHtmlContent = (filePath: string): string => {
  const locale = currentLocale()
  const isEnglish = locale === 'en-US'
  const base = filePath.replace(/\.md$/i, '')

  const candidates = []
  if (base.endsWith('-en')) {
    candidates.push(base)
    if (isEnglish) candidates.push(base.replace(/-en$/, ''))
  } else {
    candidates.push(isEnglish ? `${base}-en` : base)
    candidates.push(base)
  }

  for (const candidate of candidates) {
    const matchedKey = Object.keys(htmlModules).find((key) =>
      key.includes(`/content/html/${candidate}.html`),
    )
    if (matchedKey) {
      return htmlModules[matchedKey] as string
    }
  }

  return '<h1>Content Not Available</h1>\n<p>The requested content could not be loaded.</p>'
}

export const loadCategories = (): CategoryData => loadJsonContent('categories', [])
export const loadPosts = (): Post[] => loadJsonContent('posts', [])
export const loadNotes = (): Note[] => loadJsonContent('notes', [])
export const loadTags = (): Tag[] => loadJsonContent('tags', [])
export const loadResources = (): ResourceNode[] => loadJsonContent('resources', [])

export const EMPTY_ABOUT: AboutData = {
  introduction: '',
  experience: [],
  section: [],
  contacts: [],
}

export const loadAbout = (): AboutData => loadJsonContent('about', EMPTY_ABOUT)

/**
 * 加载文章的完整 markdown 源（含 frontmatter）。
 * `articlePath` 为 md 风格路径（如 `notes/Programming/python/python-basics/python-basics`），
 * 自动按当前 locale 匹配 `-en` 变体。找不到时返回 null。
 */
export const loadMarkdownSource = async (articlePath: string): Promise<string | null> => {
  const suffix = currentLocale() === 'en-US' ? '-en' : ''
  const target = `${articlePath.replace(/\.md$/, '')}${suffix}.md`
  const key = Object.keys(markdownModules).find(
    (k) => k.endsWith(`/${target}`) || k.endsWith(target),
  )
  if (!key) return null
  try {
    return (await markdownModules[key]()) as string
  } catch (err) {
    console.error(`Failed to load markdown source: ${target}`, err)
    return null
  }
}
