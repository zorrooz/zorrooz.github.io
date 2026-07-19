
const getCurrentLocale = () => {
  return localStorage.getItem('locale') || 'zh-CN'
}

const getLocalizedFileName = (baseName, extension = '.json') => {
  const locale = getCurrentLocale()
  const isEnglish = locale === 'en-US'

  if (baseName.endsWith('-en')) {
    return `${baseName}${extension}`
  }

  return isEnglish ? `${baseName}-en${extension}` : `${baseName}${extension}`
}

const jsonModules = import.meta.glob('../content/**/*.json', { eager: true })
const markdownModules = import.meta.glob('../content-src/**/*.md', {
  query: '?raw',
  import: 'default',
  eager: false,
})

const loadJsonContent = async (fileName) => {
  const localizedFileName = getLocalizedFileName(fileName, '.json')

  try {
    const matchedKey = Object.keys(jsonModules).find((key) => key.includes(localizedFileName))

    if (matchedKey) {
      return jsonModules[matchedKey].default || {}
    }

    throw new Error(`JSON file not found: ${localizedFileName}`)
  } catch (error) {
    console.error(`Failed to load JSON content: ${localizedFileName}`, error)

    if (getCurrentLocale() === 'en-US') {
      const fallbackKey = Object.keys(jsonModules).find((key) => key.includes(`${fileName}.json`))

      if (fallbackKey) {
        return jsonModules[fallbackKey].default || {}
      }
    }

    return {}
  }
}

export const loadMarkdownContent = async (filePath) => {
  const locale = getCurrentLocale()
  const isEnglish = locale === 'en-US'

  const possiblePaths = []

  const localizedPath = isEnglish ? filePath.replace('.md', '-en.md') : filePath
  possiblePaths.push(localizedPath)

  possiblePaths.push(filePath)

  if (isEnglish) {
    const chinesePath = filePath.replace('-en.md', '.md')
    if (chinesePath !== filePath && chinesePath !== localizedPath) {
      possiblePaths.push(chinesePath)
    }
  } else {
    const englishPath = filePath.replace('.md', '-en.md')
    if (englishPath !== filePath && englishPath !== localizedPath) {
      possiblePaths.push(englishPath)
    }
  }

  try {
    for (const path of possiblePaths) {
      const searchPaths = [`../content-src/${path}`, `../content-src/${path}?raw`]

      const matchedKey = Object.keys(markdownModules).find((key) =>
        searchPaths.some((searchPath) => key.endsWith(searchPath)),
      )

      if (matchedKey) {
        const content = await markdownModules[matchedKey]()
        return content
      }
    }

    throw new Error(`Markdown file not found for any of: ${possiblePaths.join(', ')}`)
  } catch (error) {
    console.error(`Failed to load markdown content: ${localizedPath}`, error)

    if (isEnglish) {
      const fallbackPaths = [`../content-src/${filePath}`, `../content-src/${filePath}?raw`]

      const fallbackKey = Object.keys(markdownModules).find((key) =>
        fallbackPaths.some((path) => key.endsWith(path)),
      )

      if (fallbackKey) {
        const content = await markdownModules[fallbackKey]()
        return content
      }
    }

    return '# Content Not Available\n\nThe requested content could not be loaded.'
  }
}

export const loadCategories = () => loadJsonContent('categories')
export const loadPosts = () => loadJsonContent('posts')
export const loadNotes = () => loadJsonContent('notes')
export const loadTags = () => loadJsonContent('tags')
export const loadAbout = () => loadJsonContent('about')
export const loadResources = () => loadJsonContent('resources')
