import type { AboutData, CategoryData, Note, Post, ResourceNode, Tag } from '@/types/content'

const getCurrentLocale = (): string => {
  // SSR prerender: locale injected by main.ts from the /zh|/en route prefix
  const injected = globalThis.__GBLOG_LOCALE__
  if (injected) return injected
  return (typeof window !== 'undefined' ? localStorage.getItem('locale') : null) || 'zh-CN'
}

const getLocalizedFileName = (baseName: string, extension = '.json'): string => {
  const locale = getCurrentLocale()
  const isEnglish = locale === 'en-US'

  if (baseName.endsWith('-en')) {
    return `${baseName}${extension}`
  }

  return isEnglish ? `${baseName}-en${extension}` : `${baseName}${extension}`
}

// eager for the core JSON set (about/categories/notes/posts/projects/resources/tags/topics,
// incl. -en variants); search-index* is excluded (lazy chunk in SearchModal)
// @data 别名由 vite.config.ts 指向数据分支目录（GBLOG_DATA_DIR）
const jsonModules = import.meta.glob('@data/content/[abcnprt]*.json', { eager: true })
const htmlModules = import.meta.glob('@data/content/html/**/*.html', {
  query: '?raw',
  import: 'default',
  eager: true,
})
// 惰性加载 markdown 源（复制文章用，按需拆包）。
// 覆盖两层源：content-src（zh 手写）与 cache/en（英文机器层），按路径后缀区分语言
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

  if (getCurrentLocale() === 'en-US') {
    const fallbackKey = Object.keys(jsonModules).find((key) => key.includes(`${fileName}.json`))

    if (fallbackKey) {
      return (jsonModules[fallbackKey] as JsonModule).default as T
    }
  }

  return fallback
}

/**
 * Load pre-rendered article HTML.
 * `filePath` uses the markdown-style path (e.g. `notes/Omics/genomics/bwa/bwa.md`);
 * lookup prefers the localized file and falls back to Chinese, mirroring the
 * old `loadMarkdownContent` semantics:
 *   en-US: bwa-en.html -> bwa.html
 *   zh-CN: bwa.html (or bwa-en.html when the path itself carries -en)
 */
export const loadHtmlContent = (filePath: string): string => {
  const locale = getCurrentLocale()
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
  const suffix = getCurrentLocale() === 'en-US' ? '-en' : ''
  // articlePath 可能已带 .md（如 notes/.../python-basics.md），统一去后缀后拼接
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
