/**
 * 路由分页状态：page query 为唯一事实来源（前进/后退/外部跳转自动同步）。
 * `resetOn` 变化（如筛选结果更新）时回到第 1 页；翻页写入 query 并回顶。
 */
import { onMounted, ref, watch, type Ref } from 'vue'
import { useRoute, useRouter } from 'vue-router'
import { scrollToTop } from '@/utils/scroll'

export function useRoutePagination(resetOn: Ref<unknown>, totalPages: () => number) {
  const route = useRoute()
  const router = useRouter()
  const currentPage = ref(1)

  const clampQueryPage = (): number => {
    const p = parseInt(String(route.query.page))
    return Number.isFinite(p) && p >= 1 ? Math.min(p, Math.max(1, totalPages())) : 1
  }

  function pushPage(page: number) {
    currentPage.value = page
    router.push({ path: route.path, query: { ...route.query, page: String(page) } }).catch(() => {})
    scrollToTop()
  }

  function goToPage(page: number) {
    if (page >= 1 && page <= totalPages()) pushPage(page)
  }

  function prevPage() {
    if (currentPage.value > 1) pushPage(currentPage.value - 1)
  }

  function nextPage() {
    if (currentPage.value < totalPages()) pushPage(currentPage.value + 1)
  }

  watch(resetOn, () => {
    currentPage.value = 1
  })

  watch(
    () => route.query.page,
    () => {
      const page = clampQueryPage()
      if (page !== currentPage.value) {
        currentPage.value = page
        scrollToTop()
      }
    },
  )

  onMounted(() => {
    currentPage.value = clampQueryPage()
  })

  return { currentPage, goToPage, prevPage, nextPage }
}
