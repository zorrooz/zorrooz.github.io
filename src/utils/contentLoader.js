
export const getCurrentLocale = () => {
  return localStorage.getItem('locale') || 'zh-CN'
}

export const getLocalizedFileName = (baseName, extension = '.json') => {
  const locale = getCurrentLocale()
  const isEnglish = locale === 'en-US'

  if (baseName.endsWith('-en')) {
    return `${baseName}${extension}`
  }

  return isEnglish ? `${baseName}-en${extension}` : `${baseName}${extension}`
}

export const loadJsonContent = async (fileName, directory = '') => {
  const localizedFileName = getLocalizedFileName(fileName, '.json')

  try {
    const modules = import.meta.glob('../content/**/*.json', { eager: true })

    const fullPath = directory
      ? `../content/${directory}/${localizedFileName}`
      : `../content/${localizedFileName}`

    const matchedKey = Object.keys(modules).find((key) => key.includes(localizedFileName))

    if (matchedKey) {
      return modules[matchedKey].default || {}
    }

    throw new Error(`JSON file not found: ${localizedFileName}`)
  } catch (error) {
    console.error(`Failed to load JSON content: ${localizedFileName}`, error)

    if (getCurrentLocale() === 'en-US') {
      try {
        const modules = import.meta.glob('../content/**/*.json', { eager: true })
        const fallbackKey = Object.keys(modules).find((key) => key.includes(`${fileName}.json`))

        if (fallbackKey) {
          return modules[fallbackKey].default || {}
        }
      } catch (fallbackError) {
        console.error(`Failed to load fallback JSON content: ${fileName}.json`, fallbackError)
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
    const markdownModules = import.meta.glob('../content-src/**/*.md', {
      query: '?raw',
      import: 'default',
      eager: false,
    })

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
      try {
        const markdownModules = import.meta.glob('../content-src/**/*.md', {
          query: '?raw',
          import: 'default',
          eager: false,
        })

        const fallbackPaths = [`../content-src/${filePath}`, `../content-src/${filePath}?raw`]

        const fallbackKey = Object.keys(markdownModules).find((key) =>
          fallbackPaths.some((path) => key.endsWith(path)),
        )

        if (fallbackKey) {
          const content = await markdownModules[fallbackKey]()
          return content
        }
      } catch (fallbackError) {
        console.error(`Failed to load fallback markdown content: ${filePath}`, fallbackError)
      }
    }

    return '# Content Not Available\n\nThe requested content could not be loaded.'
  }
}

export const loadMultipleJsonContents = async (fileNames, directory = '') => {
  const results = {}

  for (const fileName of fileNames) {
    try {
      const data = await loadJsonContent(fileName, directory)
      results[fileName] = data
    } catch (error) {
      console.error(`Failed to load ${fileName}`, error)
      results[fileName] = {}
    }
  }

  return results
}

export const loadCategories = () => loadJsonContent('categories')
export const loadPosts = () => loadJsonContent('posts')
export const loadNotes = () => loadJsonContent('notes')
export const loadProjects = () => loadJsonContent('projects')
export const loadTopics = () => loadJsonContent('topics')
export const loadTags = () => loadJsonContent('tags')
export const loadAbout = () => loadJsonContent('about')
export const loadResources = () => loadJsonContent('resources')

export default {
  getCurrentLocale,
  getLocalizedFileName,
  loadJsonContent,
  loadMarkdownContent,
  loadMultipleJsonContents,
  loadCategories,
  loadPosts,
  loadNotes,
  loadProjects,
  loadTopics,
  loadTags,
  loadAbout,
  loadResources,
}
