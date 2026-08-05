// reveal.ts — v-reveal 滚动入场指令（IntersectionObserver，一次触发）
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
