/**
 * v-reveal 滚动入场指令：IntersectionObserver 一次性触发，
 * 元素进入视口后加 reveal-visible（样式见 global.scss），reduce-motion 用户直接跳过。
 */
import type { Directive } from 'vue'

interface RevealElement extends HTMLElement {
  __revealIo?: IntersectionObserver | null
}

export const reveal: Directive<RevealElement> = {
  mounted(el) {
    if (typeof window === 'undefined') return
    if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) return
    if (!('IntersectionObserver' in window)) return
    el.classList.add('reveal')
    const io = new IntersectionObserver(
      ([entry]) => {
        if (entry.isIntersecting) {
          el.classList.add('reveal-visible')
          io.disconnect()
        }
      },
      { threshold: 0.12 },
    )
    io.observe(el)
    el.__revealIo = io
  },
  unmounted(el) {
    el.__revealIo?.disconnect()
    el.__revealIo = null
  },
}
