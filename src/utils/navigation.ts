/**
 * 站内导航/URL 工具：文章路径转换、locale 前缀路径、标签查询跳转。
 * 合并自原 articleUrl.ts + localePath.ts + tagQuery.ts。
 */
import { nextTick } from 'vue'
import type { Router, RouteLocationNormalizedLoaded } from 'vue-router'

import i18n from '@/i18n/index'
import {
  localeFromPath,
  localeSegmentOf,
  stripLocalePrefix,
  toSupportedLocale,
  type SupportedLocale,
} from '@/config/site'
import { scrollToTop } from '@/utils/scroll'

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

/**
 * Prefix an internal path with the current locale (/zh or /en).
 * Accepts '/category' style paths or bare article paths like
 * 'notes/Omics/genomics/bwa/bwa' (the leading slash is added).
 */
export const toLocalePath = (path: string): string => {
  const segment = localeSegmentOf(toSupportedLocale(i18n.global.locale.value))
  return `/${segment}${path.startsWith('/') ? path : `/${path}`}`
}

/** Flip the /zh|/en prefix of the current path and navigate (used by language switcher). */
export const switchLocale = (router: Router, route: RouteLocationNormalizedLoaded): void => {
  const current = localeFromPath(route.path) ?? toSupportedLocale(i18n.global.locale.value)
  const next = current === 'zh-CN' ? 'en-US' : 'zh-CN'
  router.push({
    path: `/${localeSegmentOf(next)}${stripLocalePrefix(route.path)}`,
    query: route.query,
  })
}

export interface TagNavigationOptions {
  locale: SupportedLocale
  /** 附加 query（如原 route.query / page / from），tag 由本函数统一附加并覆盖 */
  query?: Record<string, unknown>
  scroll?: boolean
}

/**
 * 首页/文章页共用的「跳转到标签筛选页」逻辑。
 * 统一行为：push 到 `${localePrefix}/`，query 合并 tag 与调用方附加项。
 */
export function goToTag(router: Router, tag: string, opts: TagNavigationOptions) {
  if (!tag) return
  router
    .push({ path: `/${localeSegmentOf(opts.locale)}/`, query: { ...(opts.query ?? {}), tag } })
    .catch(() => {})
  if (opts.scroll !== false) {
    nextTick(() => scrollToTop())
  }
}
