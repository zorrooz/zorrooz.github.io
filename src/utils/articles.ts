/**
 * 文章元数据工具：categories.json → 扁平文章列表。
 */
import type { CategoryArticle, CategoryData } from '@/types'

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
