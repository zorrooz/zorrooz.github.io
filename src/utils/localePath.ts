import type { Router, RouteLocationNormalizedLoaded } from 'vue-router'
import i18n from '@/i18n/index'

import {
  localeFromPath,
  localeSegmentOf,
  stripLocalePrefix,
  toSupportedLocale,
} from '@/config/site'

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
