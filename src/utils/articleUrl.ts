/**
 * 文章路径解析：categories.json 的 articleUrl 与站内路径之间的唯一转换。
 * 覆盖 Article.vue（2 处）与 NavigationTree.vue 的重复拆分逻辑。
 */

/** `/article/notes/Omics/...` → md 风格路径 `notes/Omics/....md` */
export const articlePathFromUrl = (articleUrl: string): string => {
  const parts = articleUrl.replace(/^\/+/, '').split('/')
  const i0 = parts[0] === 'article' ? 1 : 0
  return `${parts[i0]}/${parts[i0 + 1]}/${parts.slice(i0 + 2).join('/')}.md`
}

/** md 风格路径 → 站内文章路由路径（去 .md） */
export const toArticleRoutePath = (path: string): string => `/article/${path.replace(/\.md$/, '')}`

/** 路由 `path` 参数（string | string[]）→ 拼接路径字符串 */
export const joinRoutePathParam = (p: string | string[]): string =>
  Array.isArray(p) ? p.join('/') : typeof p === 'string' && p.length ? p : ''
