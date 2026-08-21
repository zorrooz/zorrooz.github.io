/**
 * 标签筛选跳转统一入口：首页/列表页（保留 query + 页码重置）与文章页（记录来源 + 不滚动）。
 */
import { useI18n } from 'vue-i18n'
import { useRoute, useRouter } from 'vue-router'
import { goToTag } from '@/utils/navigation'
import { toSupportedLocale } from '@/config'

export function useTagNavigation() {
  const { locale } = useI18n()
  const route = useRoute()
  const router = useRouter()

  /** 列表页：跳转标签筛选页，保留现有 query，页码重置为 1 */
  function goTag(tag: string) {
    goToTag(router, tag, {
      locale: toSupportedLocale(locale.value),
      query: { ...route.query, page: '1' },
    })
  }

  /** 文章页：记录来源路径（供「返回文章」），不触发滚动 */
  function goTagFromArticle(tag: string) {
    goToTag(router, tag, {
      locale: toSupportedLocale(locale.value),
      query: { from: route.fullPath },
      scroll: false,
    })
  }

  return { goTag, goTagFromArticle }
}
