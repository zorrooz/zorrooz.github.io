// @ts-check
/** @returns {string} */
const getCurrentLocale = () => {
  return (typeof window !== 'undefined' ? localStorage.getItem('locale') : null) || 'zh-CN'
}

/**
 * @param {string} baseName
 * @param {string} [extension]
 * @returns {string}
 */
const getLocalizedFileName = (baseName, extension = '.json') => {
  const locale = getCurrentLocale()
  const isEnglish = locale === 'en-US'

  if (baseName.endsWith('-en')) {
    return `${baseName}${extension}`
  }

  return isEnglish ? `${baseName}-en${extension}` : `${baseName}${extension}`
}

// eager for the core JSON set (about/categories/notes/posts/projects/resources/tags/topics,
// incl. -en variants); search-index* is excluded (lazy chunk in SearchModal)
const jsonModules = import.meta.glob('../content/[abcnprt]*.json', { eager: true })
const htmlModules = import.meta.glob('../content/html/**/*.html', {
  query: '?raw',
  import: 'default',
  eager: true,
})

/**
 * @param {string} fileName
 * @returns {any}
 */
const loadJsonContent = (fileName) => {
  const localizedFileName = getLocalizedFileName(fileName, '.json')

  const matchedKey = Object.keys(jsonModules).find((key) => key.includes(localizedFileName))

  if (matchedKey) {
    return (/** @type {{ default?: any }} */ (jsonModules[matchedKey])).default || {}
  }

  console.error(`Failed to load JSON content: ${localizedFileName}`)

  if (getCurrentLocale() === 'en-US') {
    const fallbackKey = Object.keys(jsonModules).find((key) => key.includes(`${fileName}.json`))

    if (fallbackKey) {
      return (/** @type {{ default?: any }} */ (jsonModules[fallbackKey])).default || {}
    }
  }

  return {}
}

/**
 * Load pre-rendered article HTML.
 * `filePath` uses the markdown-style path (e.g. `notes/Omics/genomics/bwa/bwa.md`);
 * lookup prefers the localized file and falls back to Chinese, mirroring the
 * old `loadMarkdownContent` semantics:
 *   en-US: bwa-en.html -> bwa.html
 *   zh-CN: bwa.html (or bwa-en.html when the path itself carries -en)
 */
/**
 * @param {string} filePath
 * @returns {string}
 */
export const loadHtmlContent = (filePath) => {
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
      return /** @type {string} */ (htmlModules[matchedKey]) || ''
    }
  }

  return '<h1>Content Not Available</h1>\n<p>The requested content could not be loaded.</p>'
}

export const loadCategories = () => loadJsonContent('categories')
export const loadPosts = () => loadJsonContent('posts')
export const loadNotes = () => loadJsonContent('notes')
export const loadTags = () => loadJsonContent('tags')
export const loadAbout = () => loadJsonContent('about')
export const loadResources = () => loadJsonContent('resources')
