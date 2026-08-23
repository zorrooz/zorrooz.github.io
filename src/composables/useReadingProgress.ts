/**
 * 阅读进度条：监听页面滚动，输出 0-100 的百分比（rAF 节流）。
 */
import { onBeforeUnmount, onMounted, ref } from 'vue'

export function useReadingProgress() {
  const progressPercent = ref(0)
  let ticking = false

  function update() {
    const doc = document.documentElement
    const max = doc.scrollHeight - window.innerHeight
    progressPercent.value = max > 0 ? Math.min(100, (window.scrollY / max) * 100) : 0
  }

  function onScroll() {
    if (ticking) return
    ticking = true
    requestAnimationFrame(() => {
      ticking = false
      update()
    })
  }

  onMounted(() => {
    update()
    window.addEventListener('scroll', onScroll, { passive: true })
  })

  onBeforeUnmount(() => window.removeEventListener('scroll', onScroll))

  return { progressPercent }
}
