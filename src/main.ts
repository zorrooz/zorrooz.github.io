import { createPinia } from 'pinia'
import { VueHeadMixin } from '@unhead/vue'
import { ViteSSG } from 'vite-ssg'
import type { Directive } from 'vue'

import App from './App.vue'
import { routes } from './router'
import i18n from './i18n'
import { useLocaleStore } from './stores/locale'
import { useThemeStore } from './stores/theme'
import { localeFromPath } from './config'

import './assets/styles/global.scss'
import './assets/styles/highlight/github.css'
import './assets/styles/highlight/github-dark-dimmed.css'
import '@fortawesome/fontawesome-free/css/fontawesome.min.css'
import '@fortawesome/fontawesome-free/css/solid.min.css'
import '@fortawesome/fontawesome-free/css/brands.min.css'
import 'katex/dist/katex.min.css'

// bootstrap 在模块顶层触碰 document，SSR 构建不得引入
if (!import.meta.env.SSR) import('bootstrap')

// v-reveal 滚动入场指令（IntersectionObserver，一次触发）
interface RevealElement extends HTMLElement {
  __revealIo?: IntersectionObserver | null
}

const reveal: Directive<RevealElement> = {
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

export const createApp = ViteSSG(
  App,
  { routes },
  ({ app, router, initialState, isClient, routePath }) => {
    const pinia = createPinia()
    app.use(pinia)
    if (import.meta.env.SSR) {
      initialState.pinia = pinia.state.value
    } else if (initialState.pinia) {
      pinia.state.value = initialState.pinia
    }

    app.use(router)
    app.use(i18n)
    app.mixin(VueHeadMixin)
    app.directive('reveal', reveal)

    const themeStore = useThemeStore()
    const localeStore = useLocaleStore()
    themeStore.initTheme()
    // SSR：按 /zh 或 /en 路径前缀逐页预渲染，注入当前 locale
    if (import.meta.env.SSR && typeof routePath === 'string') {
      const locale = localeFromPath(routePath)
      if (locale) {
        localeStore.setLocale(locale)
        // contentLoader（无 router 上下文）在预渲染期读取该全局值
        globalThis.__GBLOG_LOCALE__ = locale
      }
    }
    if (isClient) localeStore.initLocale()
    // PWA：生产环境注册 service worker（index.html 已带 manifest 链接）
    if (isClient && import.meta.env.PROD && 'serviceWorker' in navigator) {
      window.addEventListener('load', () => {
        navigator.serviceWorker.register('/sw.js').catch(() => {})
      })
    }
  },
)
