/**
 * 文章页双 sticky 侧栏（左导航/右 TOC）的视口适配：
 * 按吸顶 header 高度计算可用高度写入 max-height，随窗口 resize 更新。
 */
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'

export function useStickySidebars(getSidebars: () => (HTMLElement | null)[]) {
  const viewportWidth = ref(typeof window !== 'undefined' ? window.innerWidth : 1024)

  const isDesktop = computed(() => viewportWidth.value >= 992)

  /** 内容加载/渲染后调用：重算两栏可用高度 */
  function updateSidebarDimensions() {
    if (typeof window === 'undefined') return
    const headerH = document.querySelector('header')?.offsetHeight || 60
    const availableH = Math.max(200, window.innerHeight - headerH - 24 - 24)
    getSidebars().forEach((el) => {
      if (!el) return
      el.style.maxHeight = `${availableH}px`
      el.style.overflowY = 'auto'
    })
  }

  function onResize() {
    viewportWidth.value = window.innerWidth
    updateSidebarDimensions()
  }

  onMounted(() => {
    viewportWidth.value = window.innerWidth
    window.addEventListener('resize', onResize)
  })

  onBeforeUnmount(() => window.removeEventListener('resize', onResize))

  return { isDesktop, updateSidebarDimensions }
}
