/**
 * 首页/文章页共用的「跳转到标签筛选页」逻辑。
 * 统一行为：push 到 `${localePrefix}/`，query 合并 tag 与调用方附加项。
 */
import { nextTick } from 'vue'
import type { Router } from 'vue-router'

import { localeSegmentOf, type SupportedLocale } from '@/config/site'
import { scrollToTop } from '@/utils/scroll'

export interface TagNavigationOptions {
  locale: SupportedLocale
  /** 附加 query（如原 route.query / page / from），tag 由本函数统一附加并覆盖 */
  query?: Record<string, unknown>
  scroll?: boolean
}

export function goToTag(router: Router, tag: string, opts: TagNavigationOptions) {
  if (!tag) return
  router
    .push({ path: `/${localeSegmentOf(opts.locale)}/`, query: { ...(opts.query ?? {}), tag } })
    .catch(() => {})
  if (opts.scroll !== false) {
    nextTick(() => scrollToTop())
  }
}
