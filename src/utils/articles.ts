/**
 * 文章元数据工具：categories.json → 扁平文章列表 / 同组线性文章序列。
 */
import type { CategoryArticle, CategoryData } from '@/types'
import { articlePathFromUrl } from '@/utils/navigation'

/**
 * 遍历 sections → items → categories → articles，扁平化全部文章。
 * `type` 可选：仅保留 articleUrl 含 `/${type}/` 的文章（如 'notes'）。
 */
export function flattenCategoryArticles(
  categoryData: CategoryData,
  type?: string,
): CategoryArticle[] {
  const all: CategoryArticle[] = []
  categoryData.forEach((section) =>
    section.items.forEach((item) => {
      item.articles?.forEach((article) => {
        if (!type || article.articleUrl.includes(`/${type}/`)) all.push(article)
      })
      item.categories.forEach((cat) =>
        cat.articles.forEach((article) => {
          if (!type || article.articleUrl.includes(`/${type}/`)) all.push(article)
        }),
      )
    }),
  )
  return all
}

/** 文章页运行时元数据（categories.json 展开） */
export interface ArticleMeta {
  title: string
  path: string
  date: string
  tags: string[]
  wordCount: number
  description?: string
}

/** 扁平文章索引：date 取所属 item/subcategory 的最新更新日 */
export function buildArticleIndex(categoryData: CategoryData): ArticleMeta[] {
  const latestDateByUrl = new Map<string, string>()
  categoryData.forEach((section) =>
    section.items.forEach((item) => {
      const itemLatest = item.stats.latestDate || ''
      item.articles?.forEach((a) => latestDateByUrl.set(a.articleUrl, itemLatest))
      item.categories.forEach((cat) => {
        const catLatest = cat.stats.latestDate || itemLatest
        cat.articles.forEach((a) => latestDateByUrl.set(a.articleUrl, catLatest))
      })
    }),
  )
  return flattenCategoryArticles(categoryData).map((article) => ({
    title: article.title,
    path: articlePathFromUrl(article.articleUrl),
    date: latestDateByUrl.get(article.articleUrl) ?? '',
    tags: article.tags,
    wordCount: article.wordCount,
  }))
}

export interface LinearArticle {
  title: string
  path: string
}

/** 同一 `type/group`（如 notes/Programming）下的线性文章序列（上一篇/下一篇用） */
export function collectGroupArticles(
  categoryData: CategoryData,
  type: string,
  group: string,
): LinearArticle[] {
  const linear: LinearArticle[] = []
  const pushFromUrl = (title: string, articleUrl: string) => {
    if (!articleUrl.trim()) return
    const path = articlePathFromUrl(articleUrl).replace(/\.md$/, '')
    const [t, g] = path.split('/')
    if (t !== type || g !== group) return
    linear.push({ title, path: `${path}.md` })
  }

  for (const section of categoryData) {
    for (const item of section.items) {
      if (item.name !== group) continue
      item.articles?.forEach((a) => pushFromUrl(a.title, a.articleUrl))
      item.categories.forEach((cat) =>
        cat.articles.forEach((a) => pushFromUrl(a.title, a.articleUrl)),
      )
    }
  }
  return linear
}
